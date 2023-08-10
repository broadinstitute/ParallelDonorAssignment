import bamnostic
import argparse
import numpy as np
import pysam

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


def main():
    parser = argparse.ArgumentParser(description="Split the BAM file into equal parts")
    parser.add_argument("contigs", type=argparse.FileType('r'))
    parser.add_argument("BAI_PATH", type=str)
    parser.add_argument("num_splits", type=int)
    args = parser.parse_args()

    BAI_PATH = args.BAI_PATH # 'gs://nnfc-hgrm-output/tes/possorted_genome_bam.bam.bai'
    contig_lookup = dict(enumerate(pysam.AlignmentFile(args.contigs).references)) # arg.contigs == header.sam
    target_num_jobs = args.num_splits 


    bai = bamnostic.bai.Bai(BAI_PATH)
    # contig_lookup = get_contigs(bam_path)
    total_size = sum([size for _, _, _, _, _, size in iterate_bai_intervals(bai, contig_lookup)])

    # target_num_jobs = 100
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
            fh.write(region[0] + ':' + region[1] + '-' + region[2] + '\n')


if __name__ == '__main__':
    main()
