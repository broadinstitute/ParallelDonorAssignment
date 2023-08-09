from donor_assignment import *

parser = argparse.ArgumentParser(description="Run donor assignment using google batch")


bai = bamnostic.bai.Bai(${BAI_PATH})
contig_lookup = get_contigs(bam_path)
total_size = sum([size for _, _, _, _, _, size in iterate_bai_intervals(bai, contig_lookup)])

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

with open('list_of_regions.txt', 'w') as fh:  
    for region in regions_per_thread:
        region = region[0]
        fh.write(','.join(np.array(region).astype(str)) + '\n')
