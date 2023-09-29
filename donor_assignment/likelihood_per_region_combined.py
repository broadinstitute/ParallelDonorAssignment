import pandas as pd
import numpy as np
import argparse
import subprocess
import gzip
import tqdm


def single_base(read):
    return len(set(read)) == 1

@profile
def new_generate_barcode_lls(barcode_reads, donors, fill_donor_alt_cts, error_rate, regularize_factor=0.95):

    # Probability of neither REF or ALT base
    prob_neither = 2 * error_rate / 3 / len(donors)
    
    # Donor population probability adjustment
    # score for donor given genotype = (ref_prob * num_ref_in_donor / cohort_num_ref) + (alt_prob * num_alt_in_donor / cohort_num_alt
    barcode_reads['score_donor_ref'] = barcode_reads.ref_prob * 2 / barcode_reads.adj_num_donor_ref_alleles + prob_neither
    barcode_reads['score_donor_alt'] = barcode_reads.alt_prob * 2 / barcode_reads.adj_num_donor_alt_alleles + prob_neither
    barcode_reads['score_donor_het'] = (
        (barcode_reads.alt_prob / barcode_reads.adj_num_donor_alt_alleles) +
        (barcode_reads.ref_prob / barcode_reads.adj_num_donor_ref_alleles)
    ) + prob_neither
    barcode_reads['score_donor_missing'] = (
        (barcode_reads.alt_prob * barcode_reads.mean_alt_count) / barcode_reads.adj_num_donor_alt_alleles + 
        (barcode_reads.ref_prob * (2 - barcode_reads.mean_alt_count) / barcode_reads.adj_num_donor_ref_alleles)
    ) + prob_neither

    # expand scores to donors. Rows should sum to 1
    score_umi_snp_donor_ref_df = barcode_reads.score_donor_ref.values.reshape(-1, 1) * (fill_donor_alt_cts.loc[barcode_reads.pos, donors] == 0)
    score_umi_snp_donor_alt_df = barcode_reads.score_donor_alt.values.reshape(-1, 1) * (fill_donor_alt_cts.loc[barcode_reads.pos, donors] == 2)
    score_umi_snp_donor_het_df = barcode_reads.score_donor_het.values.reshape(-1, 1) * (fill_donor_alt_cts.loc[barcode_reads.pos, donors] == 1)
    score_umi_snp_donor_missing_df = \
        barcode_reads.score_donor_missing.values.reshape(-1, 1) * \
        (fill_donor_alt_cts.loc[barcode_reads.pos, donors].isna())
    
    score_umi_snp_donors_df = score_umi_snp_donor_ref_df + \
        score_umi_snp_donor_alt_df + \
        score_umi_snp_donor_het_df + \
        score_umi_snp_donor_missing_df

    # regularize
    reg_score_umi_snp_donors_df = score_umi_snp_donors_df.copy()
    reg_score_umi_snp_donors_df *= regularize_factor
    reg_score_umi_snp_donors_df += (1 - regularize_factor) / len(donors)
    reg_score_umi_snp_donors_df = np.log(reg_score_umi_snp_donors_df)
    reg_score_umi_snp_donors_df = pd.concat([barcode_reads, reg_score_umi_snp_donors_df[donors].reset_index(drop=True)], axis=1)

    # group and sum by barcodes
    barcode_lkl_df = reg_score_umi_snp_donors_df.groupby('barcode')[donors].sum()
    barcode_groupby = reg_score_umi_snp_donors_df.groupby('barcode')
    barcode_lkl_df['num_umi_snps'] = barcode_groupby.size()
    barcode_lkl_df['num_umis'] = barcode_groupby['UMI'].nunique()
    barcode_lkl_df['num_snps'] = barcode_groupby['pos'].nunique()

    return barcode_lkl_df, reg_score_umi_snp_donors_df, score_umi_snp_donors_df

def calculate_donor_liks(df, donors):
    """ Calculate donor likelihoods, given a df with columns:
        [barcode] [chr] [pos] [ref_loglikelihood] [alt_loglikelihood] [het_loglikelihood] [donor_i ... donor_k]
            [donors] columns encode "ref" "alt" or "het" for each donor at that position
            [ref_loglikelihood] [alt_loglikelihood] [het_loglikelihood] are the loglikelihoods for each barcode-umi-pos,
                calculated from the base observed in @function dropulation_likelihoods
    """
    # treat donors with no genotype at that location as equally likely for ref, alt, or het
    # this effectively diminishes the overall likelihood that the umi came from that donor
    #   -- I'm not sure that's true.  It would make them more likely than the lowest likelihood set (-Ray)
    donor_ref_liks = df.ref_loglikelihood.values.reshape(-1, 1) * (df[donors] == 0)
    donor_alt_liks = df.alt_loglikelihood.values.reshape(-1, 1) * (df[donors] == 2)
    donor_het_liks = df.het_loglikelihood.values.reshape(-1, 1) * (df[donors] == 1)
    donor_liks = donor_ref_liks + donor_het_liks + donor_alt_liks

    # fill the 0 values (no genotype) with the population average LL, not including those 0s
    donor_liks[df[donors].isna()] = np.nan
    donor_nogt_liks = donor_liks.mean(axis=1).values.reshape(-1, 1) * (df[donors].isna())
    donor_liks = donor_liks.fillna(0) + donor_nogt_liks

    donor_liks['barcode'] = df.barcode
    return donor_liks.groupby(['barcode']).sum()

@profile
def dropulation_likelihoods(barcode_reads, donors, fill_donor_alt_cts, error_rate=0.001, het_rate=0.5):
    """ Calculate the cell-donor loglikelihoods using dropulation methods"""

    # Calculate log likelihoods per UMI and snp
    barcode_reads['ref_loglikelihood'] = np.log10(barcode_reads.ref_prob)
    barcode_reads['alt_loglikelihood'] = np.log10(barcode_reads.alt_prob)

    base_is_either_mask = (
        (barcode_reads.ALT == barcode_reads.base) |
        (barcode_reads.REF == barcode_reads.base)
    )
    barcode_reads.loc[base_is_either_mask, 'het_loglikelihood'] = np.log10(het_rate)

    # Subset to barcode-umi-pos where it's either ref, alt, or het
    barcode_reads = barcode_reads[base_is_either_mask]

    # Sum values over pos - but groupby more than that so we can use .sum()
    per_snp_loglik_grpby = barcode_reads.groupby(['barcode', 'chr', 'pos', 'REF', 'ALT'])  # move  these into index
    # if testing, can sum 'A C T G' as well
    per_snp_loglik = per_snp_loglik_grpby['ref_loglikelihood alt_loglikelihood het_loglikelihood'.split()].sum()
    per_snp_loglik = per_snp_loglik.reset_index()

    # add genotypes, and compute per cbc likelihoods
    donor_likelihoods = calculate_donor_liks(per_snp_loglik.merge(fill_donor_alt_cts.reset_index(), on='pos'),
                                             donors)

    # count UMIs, and num_snps == (UMIS x SNPs)
    umi_snp_counts = barcode_reads.groupby("barcode UMI".split()).pos.nunique()
    umi_snp_counts = umi_snp_counts.reset_index().groupby("barcode")
    donor_likelihoods['num_umis'] = umi_snp_counts.UMI.nunique()
    donor_likelihoods['num_snps'] = umi_snp_counts.pos.sum()
    return donor_likelihoods

@profile
def main():
    #
    # Read in args
    #
    parser = argparse.ArgumentParser(description="Calculate donor likelihood given UMI")
    parser.add_argument("reads_on_variants_results", type=str)
    parser.add_argument("donor_list", type=str)
    parser.add_argument("VCF_region", type=str)
    parser.add_argument("region_name", type=str)
    parser.add_argument("likelihood_method", type=str, choices=['our_method', 'dropulation'])
    parser.add_argument("whitelist_fname", type=str)
    parser.add_argument("error_rate", type=np.float64, help='Error rate of base')
    parser.add_argument("--chunk_size", type=int, default=500)
    parser.add_argument("--error_by_base", help='whether to divide error_rate by 3', action='store_true', default=False)
    args = parser.parse_args()
    VCF_str = args.VCF_region
    simplified_region_name = args.region_name.replace(":", "_").replace("-", "_")
    error_rate = args.error_rate

    #
    # Get variants on reads df, and groupby barcode, umi, and position in the genome
    #
    df = pd.read_table(args.reads_on_variants_results, header=None, names=('chr', 'pos', 'read', 'barcode', 'UMI'))
    # whitelist_barcodes must contain column named 'barcode'
    whitelist_barcodes = pd.read_table(args.whitelist_fname, header=None).squeeze().values
    df = df[df.barcode.isin(whitelist_barcodes)]

    #
    # check for no reads - write empty dataframe
    #
    if len(df) == 0:
        pd.DataFrame(columns="barcode no_donor".split()).to_csv(f"barcode_log_likelihood_{simplified_region_name}.txt.gz",
                                                                sep="\t")
        return


    unique_read_counts = df.groupby(['barcode', 'UMI', 'pos', 'chr']).sum()
    unique_read_counts['num_reads'] = unique_read_counts.read.str.len()
    print(f"Read data for {len(unique_read_counts)} UMIs.")

    #
    # Delete all copies that map to different bases (e.g. TA (2 copies, diff bases), TTG (3 copies, diff bases))
    #
    
    single_base_mask = unique_read_counts.read.map(single_base)
    single_base_uniq_reads = unique_read_counts[single_base_mask].reset_index()
    single_base_uniq_reads['base'] = single_base_uniq_reads['read'].str[0]

    #
    # load donors
    #
    donors = [ln.strip() for ln in open(args.donor_list)]

    #
    # filter to donors
    #
    bcftools_proc = subprocess.Popen(f"bcftools query -s {','.join(donors)} -f %CHROM\t%POS\t%TYPE\t%REF\t%ALT[\t%GT]\n {VCF_str}".split(' '),
                                     stdout=subprocess.PIPE)

    genotypes = pd.read_table(bcftools_proc.stdout, header=None,
                              names="chrom pos type REF ALT".split() + donors)
    genotypes = genotypes.sort_values("chrom pos".split())
    genotypes = genotypes.set_index(['pos'])
    print(f"Read genotypes for {len(genotypes)} SNPs and {len(donors)} donors")

    #
    # Generate alt counts dfs
    # donor_alt_cts have the total number of alt alleles for each donor at each SNP position
    # Missing genotypes uses '.' (ie "./." or ".|.")
    #
    donor_alt_cts = genotypes[donors].applymap(lambda x: x.count('1'))
    donor_missing = genotypes[donors].applymap(lambda x: '.' in x).astype(int).replace({1: np.nan})
    fill_donor_alt_cts = donor_alt_cts + donor_missing # possible 1/.?
    print("Genotypes to counts done.")

    # 
    # Annotate genotype counts per snp
    #
    num_donors = len(donors)
    genotypes['num_donor_alt_alleles'] = fill_donor_alt_cts.sum(axis=1)
    genotypes['num_donor_ref_alleles'] = 2 * fill_donor_alt_cts.notna().sum(axis=1) - genotypes['num_donor_alt_alleles']
    genotypes['num_donor_missing'] = fill_donor_alt_cts.isna().sum(axis=1)
    genotypes['mean_alt_count'] = genotypes['num_donor_alt_alleles'] / (genotypes[['num_donor_alt_alleles', 'num_donor_ref_alleles']].sum(axis=1))
    genotypes['adj_num_donor_alt_alleles'] = genotypes.num_donor_alt_alleles + genotypes.mean_alt_count * genotypes.num_donor_missing
    genotypes['adj_num_donor_ref_alleles'] = genotypes.num_donor_ref_alleles + (2 - genotypes.mean_alt_count) * genotypes.num_donor_missing

    short_genotypes = genotypes.reset_index().rename(
        columns={'chrom': 'chr'}
    )[[
        'chr', 'pos', 'REF', 'ALT', 'num_donor_alt_alleles', 'num_donor_ref_alleles', 
        'num_donor_missing', 'mean_alt_count', 'adj_num_donor_alt_alleles', 'adj_num_donor_ref_alleles'
    ]]

    with gzip.open(f"barcode_log_likelihood_{simplified_region_name}.{args.likelihood_method}.combined.txt.gz", "wb") as outf, gzip.open(f"barcode_log_likelihood_{simplified_region_name}.{args.likelihood_method}.umi.combined.txt.gz", "wb") as umi_outf:
        #
        # Split read info into chunks with all the info for a group of CBCs to
        # save memory, then get per CBC likelihoods and write them out
        #
        all_cbcs = single_base_uniq_reads.barcode.sort_values().unique()
        cbc_chunks = np.array_split(all_cbcs, int(np.ceil(len(all_cbcs) / args.chunk_size)))
        first_block = True
        for cur_cbcs in tqdm.tqdm(cbc_chunks, total=len(cbc_chunks), desc="Unique CBC chunks"):
            
            subset_barcode_reads = single_base_uniq_reads[single_base_uniq_reads.barcode.isin(cur_cbcs)].reset_index().copy()
            subset_barcode_reads = subset_barcode_reads.merge(short_genotypes, on=['chr', 'pos'], how='left')
        
            #
            # Precompute probabilities
            #
            base_is_ref_mask = subset_barcode_reads.REF == subset_barcode_reads.base
            subset_barcode_reads.loc[base_is_ref_mask, 'ref_prob'] = 1 - error_rate
            subset_barcode_reads.loc[base_is_ref_mask, 'alt_prob'] = error_rate / 3 if args.error_by_base else error_rate
            
            base_is_alt_mask = subset_barcode_reads.ALT == subset_barcode_reads.base
            subset_barcode_reads.loc[base_is_alt_mask, 'alt_prob'] = 1 - error_rate
            subset_barcode_reads.loc[base_is_alt_mask, 'ref_prob'] = error_rate / 3 if args.error_by_base else error_rate
            
            if subset_barcode_reads.empty:
                continue
            
            # get loglikelihood functions
            if args.likelihood_method == 'our_method':
                barcode_log_likelihood, reg_donor_scores_df, prereg_donor_scores_df = new_generate_barcode_lls(
                    subset_barcode_reads, donors, fill_donor_alt_cts, error_rate
                )
                prereg_donor_scores_df.to_csv(umi_outf, header=first_block, sep='\t')
            else:
                barcode_log_likelihood = dropulation_likelihoods(
                    subset_barcode_reads, donors, fill_donor_alt_cts
                )
                
            # final output is barcode_log_probs: [barcode] x [donor] loglikelihood
            # continuously add chunks of CBC likelihoods to the output file
            barcode_log_likelihood.to_csv(outf, header=first_block, sep='\t')
            first_block = False

    print("Done.")


if __name__ == '__main__':
    main()
