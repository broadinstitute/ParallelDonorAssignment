import bamnostic
import argparse
import pysam
import numpy as np


def include_contig_in_analysis(name):
    return (name.isnumeric() or
            (name[:3] == 'chr' and name[3:].isnumeric()))


def iterate_bai_intervals(bai_file, contig_lookup):
    # a BAI file contains a linear index for each contig with 16384bp intervals
    # for example, chr1 is split into 15195 intervals of 16384bp
    # for each interval, the linear contains the smallest file offset of the
    # alignments that overlap with the interval
    # the value 16384 is available as bai._LINEAR_INDEX_WINDOW
    starts = {}
    for ref_id in range(bai_file.n_refs):
        starts[ref_id] = bai_file.get_ref(ref_id).intervals[0]

    for ref_id in range(bai_file.n_refs):
        if not include_contig_in_analysis(contig_lookup[ref_id]):
            continue
        intervals = list(bai_file.get_ref(ref_id).intervals)
        # tack on the start of next contig if available
        if ref_id + 1 in starts:
            intervals += [starts[ref_id + 1]]
        # each array element is int64, with the first 48 bits indicating the
        # position in the compressed file, and the last 16 bits the offset
        # after decompression; we use the compressed file offset
        interval_starts = [interval >> 16 for interval in intervals]
        yield ref_id, np.array(interval_starts)


def main():
    parser = argparse.ArgumentParser(description="Split the BAM file into equal parts")
    parser.add_argument("contigs", type=argparse.FileType('r'))
    parser.add_argument("BAI_PATH", type=str)
    parser.add_argument("num_splits", type=int)
    args = parser.parse_args()

    BAI_PATH = args.BAI_PATH  # 'gs://nnfc-hgrm-output/test/possorted_genome_bam.bam.bai'
    contig_lookup = dict(enumerate(pysam.AlignmentFile(args.contigs).references))  # arg.contigs == header.sam
    target_num_jobs = args.num_splits

    bai = bamnostic.bai.Bai(BAI_PATH)
    total_size = sum([(v[-1] - v[0]) for _, v in iterate_bai_intervals(bai, contig_lookup)])

    size_per_job = total_size / target_num_jobs

    jobs = []
    for ref_id, intervals in iterate_bai_intervals(bai, contig_lookup):
        contig_name = contig_lookup[ref_id]
        start_interval_idx = 0
        next_interval_idx = 0
        while next_interval_idx < len(intervals):
            chunk_size = intervals[next_interval_idx] - intervals[start_interval_idx]
            if (chunk_size > size_per_job):
                jobs.append([contig_name,
                             start_interval_idx * bai._LINEAR_INDEX_WINDOW,
                             next_interval_idx  * bai._LINEAR_INDEX_WINDOW,
                             intervals[start_interval_idx],
                             intervals[next_interval_idx]])
                start_interval_idx = next_interval_idx
            else:
                next_interval_idx += 1
        if start_interval_idx < len(intervals) - 1
            jobs.append([contig_name,
                         start_interval_idx * bai._LINEAR_INDEX_WINDOW,
                         next_interval_idx  * bai._LINEAR_INDEX_WINDOW,
                         intervals[start_interval_idx],
                         intervals[next_interval_idx - 1]])
    assert len(jobs)

    with open('list_of_regions.txt', 'w') as fh:
        for contig, start, end, byte_start, byte_end in jobs:
            fh.write(f"{contig}:{start}-{end} bytes:{byte_start}-{byte_end - 1}\n")  # inclusive ranges


if __name__ == '__main__':
    main()
