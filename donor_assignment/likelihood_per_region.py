import pandas as pd
import numpy as np
import argparse
import subprocess
import gzip


def single_base(read):
    return len(set(read)) == 1

def generate_barcode_lls(barcode_pos_reads, genotypes, donors, num_donors, ref_probs, alt_probs):
    """ Gather loglikelihood a cell came from a donor"""
    barcode_pos_reads.reset_index('pos', inplace=True)

    temp_probs = 0
    for base in 'ACGT':
        # Indicator vector of shape [umi_probs_position_index.shape[0] (barcode-umi-pos)] x 1
        base_is_ref = (genotypes.REF.loc[barcode_pos_reads.pos] == base).values.reshape((-1, 1))
        base_is_alt = (genotypes.ALT.loc[barcode_pos_reads.pos] == base).values.reshape((-1, 1))
        base_is_neither = 1 - (base_is_alt + base_is_ref)

        # temp_probs ([pos] x [donor] with positions repeated in umi_probs_position_index ordering)
        # sum[base] P(donor contributed REF / ALT allele) * Indicator(base is REF / ALT) * P(base | observed bases as SNP)
        temp_probs += ref_probs.loc[barcode_pos_reads.pos] * base_is_ref * barcode_pos_reads[f'prob_{base}'].values.reshape((-1, 1))
        temp_probs += alt_probs.loc[barcode_pos_reads.pos] * base_is_alt * barcode_pos_reads[f'prob_{base}'].values.reshape((-1, 1))
        temp_probs += (1 / num_donors) * base_is_neither * barcode_pos_reads[f'prob_{base}'].values.reshape((-1, 1))
    
    # umi_probs_position_index ([barcode, umi] x [pos, chr, read, prob_A, prob_C, prob_G, prob_T, donor_1 .... donor_k for k donors])
    # has the total probability that a barcode/umi came from a donor per donor
    # for our current subset of CBCs
    assert (barcode_pos_reads.pos == temp_probs.index).all()
    tmp = pd.concat([barcode_pos_reads.reset_index(), temp_probs.reset_index(drop=True)], axis=1)
    umi_probs_position_index = tmp.set_index("barcode UMI".split())

    # Regularize the donor probabilities
    # log transform
    # Sum log likelihoods per barcode
    #
    regularized_log_probs = umi_probs_position_index.copy()
    regularized_log_probs[donors] *= 0.95
    regularized_log_probs[donors] += 0.05 / num_donors
    regularized_log_probs[donors] = np.log(regularized_log_probs[donors])
    
    #
    # Add in umi counts and snp counts
    #
    barcode_log_likelihood = regularized_log_probs.groupby(['barcode'])[donors].sum()
    num_snps = umi_probs_position_index.groupby(['barcode']).size()
    num_umis = umi_probs_position_index.reset_index().groupby('barcode')['UMI'].unique().str.len()
    barcode_log_likelihood['num_snps'] = num_snps
    barcode_log_likelihood['num_umis'] = num_umis

    return barcode_log_likelihood

            

def main():
    #
    # Read in args
    #
    parser = argparse.ArgumentParser(description="Calculate donor likelihood given UMI")
    parser.add_argument("reads_on_variants_results", type=str)
    parser.add_argument("donor_list", type=str)
    parser.add_argument("VCF_region", type=str)
    parser.add_argument("region_name", type=str)
    args = parser.parse_args()
    VCF_str = args.VCF_region

    #
    # Get variants on reads df, and groupby barcode, umi, and position in the genome
    #
    with gzip.open(args.reads_on_variants_results, 'rt') as file:
        df = pd.read_table(file, header=None, names=('chr', 'pos', 'read', 'barcode', 'UMI'))
    unique_read_counts = df.groupby(['barcode', 'UMI', 'pos']).sum()
    unique_read_counts['num_reads'] = unique_read_counts.read.str.len()
    print(f"Read data for {len(unique_read_counts)} UMIs.")

    #
    # Generate probability table for each base
    #
    probs = {"A": list([.9, .1/3, .1/3, .1/3]),
             "C": list([.1/3, .9, .1/3, .1/3]),
             "G": list([.1/3, .1/3, .9, .1/3]),
             "T": list([.1/3, .1/3, .1/3, .9])}
    probs = pd.DataFrame(probs).T
    probs.columns = 'prob_A prob_C prob_G prob_T'.split()

    #
    # Delete all copies that map to different bases (e.g. TA (2 copies, diff bases), TTG (3 copies, diff bases))
    #
    single_base_mask = unique_read_counts.read.map(single_base)
    assert single_base_mask.mean() > 0.95
    single_base_uniq_reads = unique_read_counts[single_base_mask].copy()

    #
    # calculate general probability that each read is a certain base
    #
    for col in probs.columns:
        single_base_uniq_reads[col] = single_base_uniq_reads.read.str[0].map(probs[col])**(single_base_uniq_reads.num_reads)
    single_base_uniq_reads[probs.columns] = single_base_uniq_reads[probs.columns].div(single_base_uniq_reads[probs.columns].sum(axis=1), axis=0)

    #
    # load donors
    #
    donors = [ln.strip() for ln in open(args.donor_list)]
    donors_joined = ",".join(donors)
    #
    # filter to donors
    #
    bcftools_proc = subprocess.Popen(f"bcftools query -s {donors_joined} -f %CHROM\t%POS\t%TYPE\t%REF\t%ALT[\t%GT]\n {VCF_str}".split(' '),
                                     stdout=subprocess.PIPE)

    genotypes = pd.read_table(bcftools_proc.stdout, header=None,
                              names="chrom pos type REF ALT".split() + donors)
    genotypes = genotypes.sort_values("chrom pos".split())
    genotypes = genotypes.set_index(['pos'])
    print(f"Read genotypes for {len(genotypes)} SNPs and {len(donors)} donors")

    #
    # Generate ref and alt counts dfs
    # refs_df and alt_df have the total number of reference and alt alleles for each donor at each SNP position
    #
    refs_df = pd.DataFrame(index=genotypes.index, columns=donors)
    alt_df = pd.DataFrame(index=genotypes.index, columns=donors)
    for col in donors:
        refs_df[col] = genotypes[col].str.count('0')
        alt_df[col] = genotypes[col].str.count('1')
    print("Genotypes to counts done.")

    # ref_probs and alt_probs is the proportion of observed alleles that came
    # from that donor -- normalized to the total counts of REF / ALT at that
    # SNP (it is relative to the total number of reference or alt alleles at
    # that SNP across the donor population)
    ref_probs = refs_df.div(refs_df.sum(axis=1), axis=0)
    alt_probs = alt_df.div(alt_df.sum(axis=1), axis=0)
    assert ref_probs.index.is_unique

    #####
    # Generate barcode loglikelihoods, and write to output file
    # final output is barcode_log_likelihood: [barcode] x [donor] loglikelihood
    #####

    # Split read info into 100 sorted CBC chunks to save memory
    # and get CBC likelihoods  
    barcode_reads = single_base_uniq_reads.reset_index('barcode')
    all_cbcs = barcode_reads.barcode.sort_values().unique()
    num_donors = len(donors)

    simplified_region = args.region_name.replace(":", "_").replace("-", "_")
    with open(f"barcode_log_likelihood_{simplified_region}.txt.gz", "wb") as outf:
        first_block = True
        while len(all_cbcs) > 0:
            cur_cbcs = all_cbcs[:100]
            all_cbcs = all_cbcs[100:]
            # get loglikelihood functions
            barcode_log_likelihood = generate_barcode_lls(barcode_reads[barcode_reads.barcode.isin(cur_cbcs)], genotypes, 
                                                          donors, num_donors, ref_probs, alt_probs)
            # final output is barcode_log_probs: [barcode] x [donor] loglikelihood
            # continuously add chunks of CBC likelihoods to the output file
            barcode_log_likelihood.to_csv(outf, header=first_block, compression='gzip', sep='\t')
            first_block = False

    print("Done.")


if __name__ == '__main__':
    main()
