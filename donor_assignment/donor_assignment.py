import argparse

parser = argparse.ArgumentParser(description="Run donor assignment using google batch")
parser.add_argument("bam_path", type=str)
parser.add_argument("vcf_path", type=str)
parser.add_argument("res_workspace_bucket", type=str)
parser.add_argument("project_name", type=str)
parser.add_argument("project_region", type=str)
args = parser.parse_args()


bam_path = args.bam_path # e.g. 'gs://nnfc-hgrm-output/tes/possorted_genome_bam.bam'
vcf_path = args.vcf_path # e.g 'gs://landerlab-vcf/CIRM_plus_15donors/merged_GTex_iPSCORE_CIRM.vcf.gz'

workspace_bucket = args.res_workspace_bucket # e.g. 'fc-secure-be53f2e7-3be8-4c69-b217-9903a3210729'
temp_path = f'gs://{workspace_bucket}/temp_path'
results_path = f'gs://{workspace_bucket}/results_path/donor_assignment'

project_name = args.project_name # e.g. 'nnfc-hgrm'
project_region = args.project_region # e.g. 'us-central1'


import tempfile
import os
from google.cloud import storage, batch_v1
import uuid
import zlib
import re
import json
import gzip
from collections import Counter
import bamnostic
import pandas as pd
from io import StringIO
import pysam


def split_bucket_blob(path):
    bucket_name = path.split('/')[2]
    blob_name = path[len(bucket_name)+6:]
    return bucket_name, blob_name

def gcs_list_files(path):
    bucket_name, blob_name = split_bucket_blob(path)
    storage_client = storage.Client()
    blobs = storage_client.list_blobs(bucket_name, prefix=blob_name)
    blobs = list(blobs)
    return blobs

def gcs_open(path, mode='r'):
    bucket_name, blob_name = split_bucket_blob(path)
    
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    fid = blob.open(mode)
    return fid

def download_bai_file(bam_path):
    bai_path = bam_path + '.bai'
    local_folder = tempfile.mkdtemp()
    local_bai_path = os.path.join(local_folder, 'donor_assignment.bai')
    bucket_name, blob_name = split_bucket_blob(bai_path)
    
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    blob.download_to_filename(local_bai_path)
    return local_bai_path

def get_contigs(bam_path):
    # reads the first 1E6 bytes of the file to learn the map between contig id and name
    bam_fid = gcs_open(bam_path, 'rb')
    chunk = bam_fid.read(1000000)
    gz = zlib.decompressobj(wbits = zlib.MAX_WBITS | 16) # gzip
    header = gz.decompress(chunk).decode('ascii', errors='ignore')
    contigs = []
    for line in header.split('\n'):
        if 'SN:' in line:
            contig = line.split('\t')[1][3:]
            contigs.append(contig)
    contig_lookup = {idx : contig for idx, contig in enumerate(contigs)}
    return contig_lookup

def include_contig_in_analysis(name):
    try:
        int(name[3:])
        return name[:3] == 'chr'
    except:
        return False
    
def iterate_bai_intervals(bai_file, contig_lookup):
    # a BAI file contains a linear index for each contig with 16384bp intervals
    # for example, chr1 is split into 15195 intervals of 16384bp
    # for each interval, the linear contains the smallest file offset of the
    # alignments that overlap with the interval
    # the value 16384 is available as bai._LINEAR_INDEX_WINDOW
    for ref_id in range(bai_file.n_refs):
        if not include_contig_in_analysis( contig_lookup[ref_id] ):
            continue
        intervals = bai_file.get_ref(ref_id).intervals
        # each array element is int64, with the first 48 bits indicating the
        # position in the compressed file, and the last 16 bits the offset
        # after decompression; we use the compressed file offset
        interval_starts = [interval >> 16 for interval in intervals]
        start = interval_starts[0]
        end = interval_starts[-1]
        size = end - start
        yield ref_id, intervals, interval_starts, start, end, size
        
def get_reference(chr1_size):
    if chr1_size == '249250621':
        reference = 'GRCh37'
    elif chr1_size == '248956422':
        reference = 'GRCh38'
    else:
        reference = f'Unknown (chr1 size is {chr1_size})'
    return reference

def gcs_read_first_bytes(path, num_bytes = 1000000, encoding='ascii'):
    fid = gcs_open(path, 'rb')
    header = fid.read(num_bytes)
    if path.endswith('.gz') or path.endswith('.bgz'):
        gz = zlib.decompressobj(wbits = zlib.MAX_WBITS | 16) # gzip
        header = gz.decompress(header)
    header = header.decode(encoding, errors='ignore')
    return header;

def get_vcf_reference(vcf_path):
    header = gcs_read_first_bytes(vcf_path, 1000000)
    reference = None
    for line in header.split('\n'):
        contig = re.search('ID=((chr)?[0-9]+)[,>]', line)
        contig_size = re.search('length=([0-9]+)', line)
        if line.startswith('##contig=') and contig and contig_size:
            contig = contig.group(1).replace('chr', '')
            contig_size = contig_size.group(1)
            if contig != '1':
                continue
            reference = get_reference(contig_size)
    return reference

def get_bam_reference(bam_path):
    header = gcs_read_first_bytes(bam_path, 1000000)
    reference = None
    for line in header.split("\n"):
        contig = re.search('SN:((chr)?[0-9]+)[\t$]', line)
        contig_size = re.search('LN:([0-9]+)(\t|$)', line)
        if line.startswith('@SQ') and contig and contig_size:
            contig = contig.group(1).replace('chr', '')
            contig_size = contig_size.group(1)
            if contig != '1':
                continue
            reference = get_reference(contig_size)
    return reference



# Step 1: Split up BAM file in equal parts
local_bai_path = download_bai_file(bam_path)


bai = bamnostic.bai.Bai(local_bai_path)
contig_lookup = get_contigs(bam_path)
total_size = sum([size for _, _, _, _, _, size in iterate_bai_intervals(bai, contig_lookup)])
print(f'Total size: {total_size / 1024**3 :.1f} GB')

target_num_jobs = 100
size_per_job = total_size / target_num_jobs
size_done = 0
regions_per_thread = [[]]
assigned_to_thread = [0]
for ref_id, intervals, interval_starts, start, end, size in iterate_bai_intervals(bai, contig_lookup):
    contig_name = contig_lookup[ref_id]
    regions_per_thread[-1].append([contig_name, 1, 0])
    while size_done + size > size_per_job:
        idx = min([idx for idx, interval_start in enumerate(interval_starts) if size_done + interval_start - start > size_per_job])
        locus = idx * bai._LINEAR_INDEX_WINDOW
        size_done += interval_starts[idx] - start
        start = interval_starts[idx]
        size = end - start
        regions_per_thread[-1][-1][2] = locus
        regions_per_thread.append([[contig_name, locus, 0]])
        size_done = 0
    regions_per_thread[-1][-1][2] = len(intervals) * bai._LINEAR_INDEX_WINDOW
    size_done += size
print(len(regions_per_thread))

# Step 2: Use cloud Batch to process BAM file

def create_batch_job(project_id, region, job_name, fstab, script, num_tasks, network_name = None):
    client = batch_v1.BatchServiceClient()

    runnable = batch_v1.Runnable()
    runnable.script = batch_v1.Runnable.Script()
    runnable.script.text = script

    task = batch_v1.TaskSpec()
    task.runnables = [runnable]

    for bucket_name, mount_path in fstab.items():
        gcs_bucket = batch_v1.GCS()
        gcs_bucket.remote_path = bucket_name
        gcs_volume = batch_v1.Volume()
        gcs_volume.gcs = gcs_bucket
        gcs_volume.mount_path = mount_path
        task.volumes.append(gcs_volume)

    resources = batch_v1.ComputeResource()
    resources.cpu_milli = 1000  # in milliseconds per cpu-second
    resources.memory_mib = 512  # in MiB
    task.compute_resource = resources
    
    task.max_retry_count = 5
    task.max_run_duration = "7200s"

    group = batch_v1.TaskGroup()
    group.task_count = num_tasks
    group.task_spec = task

    policy = batch_v1.AllocationPolicy.InstancePolicy()
    policy.machine_type = "e2-highcpu-4"
    policy.provisioning_model = batch_v1.types.AllocationPolicy.ProvisioningModel.SPOT # SPOT or STANDARD
    instances = batch_v1.AllocationPolicy.InstancePolicyOrTemplate()
    instances.policy = policy

    allocation_policy = batch_v1.AllocationPolicy()
    allocation_policy.instances = [instances]
    
    if network_name is not None:
        iface = batch_v1.AllocationPolicy.NetworkInterface(network=network_name)
        network = batch_v1.AllocationPolicy.NetworkPolicy()
        network.network_interfaces = [iface]
        allocation_policy.network = network

    job = batch_v1.Job()
    job.task_groups = [group]
    job.allocation_policy = allocation_policy
    job.logs_policy = batch_v1.LogsPolicy()
    job.logs_policy.destination = batch_v1.LogsPolicy.Destination.CLOUD_LOGGING

    create_request = batch_v1.CreateJobRequest()
    create_request.job = job
    create_request.job_id = job_name
    create_request.parent = f"projects/{project_id}/locations/{region}"

    return client.create_job(create_request)


job_path = f'{temp_path}/{uuid.uuid4()}.json'
python_script_path = f'{temp_path}/{uuid.uuid4()}.py'
job_name = f'donor-assignment-{uuid.uuid4()}'

parameters = {'regions per thread': regions_per_thread, 'bam path': bam_path, 'vcf path': vcf_path, 'results path': results_path, 'python script path': python_script_path, 'job path': job_path}
fstab = {}
for key, value in parameters.items():
    if 'path' not in key:
        continue
    bucket_name, _ = split_bucket_blob(parameters[key])
    fstab[bucket_name] = f'/mnt/{bucket_name}'
    parameters[key] = parameters[key].replace('gs://', '/mnt/')
    
python_script = f"""
import os
import subprocess
import json
import gzip

def read_file(path, mode='r'):
    fid = open(path, mode)
    return fid.read()

class count_variants_on_region:
    def __init__(self, region, bam_path, vcf_path):
        self.region = region
        self.bam_path = bam_path
        self.vcf_path = vcf_path
        self.variants = {{}}
            
    def load_vcf(self):
        import pysam
        vcf = pysam.VariantFile(self.vcf_path)
        region=f'{{self.region[0]}}:{{self.region[1]}}-{{self.region[2]}}'
        try:
            vcf_fetch = vcf.fetch(region=region)
        except:
            vcf_fetch = vcf.fetch(region=region.replace('chr', ''))
        for vcf_line in vcf_fetch:
            pos, ref, alts = vcf_line.pos, vcf_line.ref, vcf_line.alts
            self.variants[pos] = ref

    def count(self, output_file):
        import pysam
        import gzip
        bamfile = pysam.AlignmentFile(self.bam_path, mode="rb")
        fid = gzip.open(output_file, 'at')
        for idx,bam_line in enumerate(bamfile.fetch(region=f'{{self.region[0]}}:{{self.region[1]}}-{{self.region[2]}}')):
            if bam_line.mapping_quality != 255: # unique mapping
                continue
            try: # CB tag isn't always present
                barcode = bam_line.get_tag('CB')
                UMI = bam_line.get_tag('UB')
            except:
                continue

            for read_idx, locus in bam_line.get_aligned_pairs(matches_only=True):
                locus += 1 # AlignmentFile is 0-based while VariantFile is 1-based
                if not locus in self.variants:
                    continue
                read = bam_line.seq[read_idx]                
                fid.write(f'{{self.region[0]}}\\t{{locus}}\\t{{read}}\\t{{barcode}}\\t{{UMI}}\\n')
        fid.close()

task_index = int(os.environ['BATCH_TASK_INDEX'])
parameters = json.loads(read_file('{parameters['job path']}'))
regions_per_thread = parameters['regions per thread']
regions = regions_per_thread[task_index]

output_file = os.path.join(parameters['results path'], f'results_{{task_index}}.tsv.gz')
for region in regions:
    print('processing', region)
    c2d = count_variants_on_region(region, parameters['bam path'], parameters['vcf path'])
    c2d.load_vcf()
    c2d.count(output_file)"""

with gcs_open(job_path, 'w') as f:
    f.write(json.dumps(parameters))
    
with gcs_open(python_script_path, 'w') as f:
    f.write(python_script)

script = f"""
apt -o DPkg::Lock::Timeout=180 -y install python3-pysam
cp {parameters['python script path']} /
python3 /{os.path.basename(python_script_path)}
"""

num_tasks = len(regions_per_thread)
network = "projects/160048904816/global/networks/sc-gcp-network"
batch_job = create_batch_job(project_name, project_region, job_name, fstab, script, num_tasks, network)


# Step 3: Collect and Aggregate Results

blobs = gcs_list_files(results_path)

dfs = []
for blob in blobs:
    match_task_index = re.match('.*_([0-9]+).tsv.gz', blob.name)
    if not match_task_index:
        continue
    task_index = int(match_task_index[1])
    with blob.open('rb') as fid:
        data = fid.read()
    data = gzip.decompress(data)
    data = data.decode('ascii')
    df = pd.read_csv(StringIO(data), sep="\t", header=None, names=('chr', 'pos', 'read', 'barcode', 'UMI'))
    df.sort_values(['chr', 'pos', 'barcode', 'UMI'], inplace=True)
    dfs.append(df)

df = pd.concat(dfs)


# you probably want to move these steps into the parallelized jobs, and access the vcf file through parameters['vcf path']
vcf = pysam.VariantFile('merged_WGS_spikeins.vcf.gz')
variants = {}
regions = regions_per_thread[task_index]
for region in regions:
    vcf_region=f'{region[0]}:{region[1]}-{region[2]}'
    try:
        vcf_fetch = vcf.fetch(region=vcf_region)
    except:
        vcf_fetch = vcf.fetch(region=vcf_region.replace('chr', ''))

    for vcf_line in vcf.fetch():
        pos, ref, alts = vcf_line.pos, vcf_line.ref, vcf_line.alts
        gt_bases = [vcf_line.ref] + list(vcf_line.alts)
        alleles_of_all_samples = [allele for idx in range(len(vcf_line.samples)) for allele in vcf_line.samples[idx]['GT']]
        alleles_of_all_samples = Counter(alleles_of_all_samples)

        likelihoods = []
        for idx, sample in enumerate(vcf_line.samples):
            likelihood_sample = {'C': 0, 'T': 0, 'A': 0, 'G': 0, 'N': 0}
            for allele in vcf_line.samples[idx]['GT']:
                likelihood_sample[ gt_bases[allele] ] += 1 / alleles_of_all_samples[allele]
            likelihoods.append(likelihood_sample)
        variants[f'{region[0]}:{pos}'] = {'ref': ref, 'alt': alts, 'likelihoods' : likelihoods}


n_ref_reads = 0
for row in df.itertuples():
    variant = variants[f'{row.chr}:{row.pos}']
    if row.read == variant['ref']:
        n_ref_reads += 1
print('Fraction of reads that matches the reference:', n_ref_reads / df.shape[0])


from math import log
from collections import defaultdict

samples = list(vcf.header.samples)
N_samples = len(samples)
log_ll = defaultdict(lambda : [0]*len(samples))
for row in df.itertuples():
    variant = variants[f'{row.chr}:{row.pos}']
    p_read_is_noise = 0.05
    for idx, likelihood in enumerate(variant['likelihoods']):
        log_ll[row.barcode][idx] += log((1-p_read_is_noise) * likelihood[row.read] + p_read_is_noise / N_samples)