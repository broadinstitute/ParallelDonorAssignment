import gzip
import pysam
import argparse


class count_variants_on_region:
    def __init__(self, bam_path, vcf_path, umi_tag):
        self.bam_path = bam_path
        self.vcf_path = vcf_path
        self.umi_tag = (umi_tag if umi_tag is not None else 'UB')
        self.variants = {}

    def load_vcf(self):
        vcf = pysam.VariantFile(self.vcf_path)
        # this should store the read's position and if it's ref or alt
        # might need to include alt as well in self.variants ??
        for vcf_line in vcf:
            pos, ref, alts = vcf_line.pos, vcf_line.ref, vcf_line.alts  # noqa
            self.variants[pos] = ref

    def count(self, output_file):
        bamfile = pysam.AlignmentFile(self.bam_path, mode="rb", ignore_truncation=True)
        fid = gzip.open(output_file, 'at')
        for bam_line in bamfile:
            if bam_line.mapping_quality < 20:
                continue
            try:  # CB tag isn't always present
                barcode = bam_line.get_tag('CB')
                if self.umi_tag != 'ATAC':
                    UMI = bam_line.get_tag(self.umi_tag)
                else:
                    UMI = bam_line.query_name
            except Exception:
                continue

            for read_idx, genome_pos in bam_line.get_aligned_pairs(matches_only=True):
                # see if a read overlaps a snp site, and if it does print out
                # what you read (location of the snp overlap, base that was
                # read, barcode umi etc)
                genome_pos += 1  # AlignmentFile is 0-based while VariantFile is 1-based
                if genome_pos not in self.variants:
                    continue
                base = bam_line.seq[read_idx]
                fid.write(f'{bam_line.reference_name}\t{genome_pos}\t{base}\t{barcode}\t{UMI}\n')
        fid.close()


def main():
    parser = argparse.ArgumentParser(description="Count the reads on variants")
    parser.add_argument("BAM_region", type=argparse.FileType('r'))
    parser.add_argument("VCF_region", type=argparse.FileType('r'))
    parser.add_argument("--umi-tag", type=str, help='optional argument for UMI tag in BAM, or "ATAC" to use read name')
    args = parser.parse_args()
    bam_file = args.BAM_region
    vcf_file = args.VCF_region

    output_file = 'results.tsv.gz'
    c2d = count_variants_on_region(bam_file, vcf_file, args.umi_tag)
    c2d.load_vcf()
    c2d.count(output_file)


if __name__ == '__main__':
    main()
