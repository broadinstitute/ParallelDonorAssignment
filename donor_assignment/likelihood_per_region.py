import pandas as pd
import numpy as np
import argparse
import subprocess
import gzip


def single_base(read):
    return len(set(read)) == 1


def generate_barcode_lls(barcode_pos_reads, genotypes, donors, num_donors,
                         ref_probs, alt_probs, simplified_region_name):
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
        positions = barcode_pos_reads.pos
        temp_probs += ref_probs.loc[positions] * base_is_ref * barcode_pos_reads[f'prob_{base}'].values.reshape((-1, 1))
        temp_probs += alt_probs.loc[positions] * base_is_alt * barcode_pos_reads[f'prob_{base}'].values.reshape((-1, 1))
        temp_probs += (1 / num_donors) * base_is_neither * barcode_pos_reads[f'prob_{base}'].values.reshape((-1, 1))

    # umi_probs_position_index ([barcode, umi] x [pos, chr, read, prob_A, prob_C, prob_G, prob_T,
    #                                             donor_1 .... donor_k for k donors])
    # has the total probability that a barcode/umi came from a donor per donor
    # for our current subset of CBCs
    assert (barcode_pos_reads.pos == temp_probs.index).all()
    tmp = pd.concat([barcode_pos_reads.reset_index(), temp_probs.reset_index(drop=True)], axis=1)
    umi_probs_position_index = tmp.set_index("barcode UMI".split())

    # Regularize the donor probabilities
    # log transform
    # Sum log likelihoods per barcode
    #
    regularized_log_probs = umi_probs_position_index[donors]
    regularized_log_probs *= 0.95
    regularized_log_probs += 0.05 / num_donors
    regularized_log_probs = np.log(regularized_log_probs)

    #
    # Add in umi counts and snp counts
    #
    barcode_log_likelihood = regularized_log_probs.groupby(['barcode']).sum()

    umi_probs_position_index['snp_id'] = umi_probs_position_index['chr'].astype(str) + ':' \
        + umi_probs_position_index['pos'].astype(str) + ':' \
        + umi_probs_position_index['read'].astype(str)

    num_umi_snps = umi_probs_position_index.groupby(['barcode']).size()
    num_snps = umi_probs_position_index.groupby(['barcode'])['snp_id'].nunique()
    num_umis = umi_probs_position_index.reset_index().groupby('barcode')['UMI'].nunique()

    barcode_log_likelihood['num_umi_snps'] = num_umi_snps
    barcode_log_likelihood['num_snps'] = num_snps
    barcode_log_likelihood['num_umis'] = num_umis

    return barcode_log_likelihood


def calculate_donor_liks(df, donors):
    """ Calculate donor likelihoods, given a df with columns:
        [barcode] [chr] [pos] [ref_loglikelihood] [alt_loglikelihood] [het_loglikelihood] [donor_i ... donor_k]
            [donors] columns encode "ref" "alt" or "het" for each donor at that position
            [ref_loglikelihood] [alt_loglikelihood] [het_loglikelihood] are the loglikelihoods for each barcode-umi-pos,
                calculated from the base observed in @function dropulation_likelihoods
    """
    # donors with no genotype get average population likelihood
    donor_ref_liks = df.ref_loglikelihood.values.reshape(-1, 1) * (df[donors] == 'ref')
    donor_alt_liks = df.alt_loglikelihood.values.reshape(-1, 1) * (df[donors] == 'alt')
    donor_het_liks = df.het_loglikelihood.values.reshape(-1, 1) * (df[donors] == 'het')
    donor_liks = donor_ref_liks + donor_het_liks + donor_alt_liks

    # fill the 0 values (no genotype) with the population average LL, not including those 0s
    donor_liks[df[donors] == 'none'] = np.nan
    donor_nogt_liks = donor_liks.mean(axis=1).values.reshape(-1, 1) * (df[donors] == 'none')
    donor_liks = donor_liks.fillna(0) + donor_nogt_liks

    donor_liks['barcode'] = df.barcode
    return donor_liks.groupby(['barcode']).sum()


def dropulation_likelihoods(barcode_reads, genotypes, refs_df, alt_df, donors, error_rate=0.001, het_rate=0.5):
    """ Calculate the cell-donor loglikelihoods using dropulation methods"""
    intermediate_df = barcode_reads.reset_index().copy()

    # Get the single base (assumed mismatching bases from the same UMI are all dropped already)
    intermediate_df['base'] = intermediate_df['read'].str[0]
    # For testing, can get one hot encoding of bases. But not necessary
    # e.g pd.get_dummies(pd.Series(intermediate_df['base']))]

    # Annotate REF and ALT with the genotype data (merge is slow)
    tmp_genotypes = genotypes.reset_index().rename(columns={'chrom': 'chr'})[['chr', 'pos', 'REF', 'ALT']]
    intermediate_df = intermediate_df.merge(tmp_genotypes, on=['chr', 'pos'])

    if intermediate_df.empty:
        return None

    # Calculate log likelihoods per UMI and snp
    base_is_ref_mask = intermediate_df.REF == intermediate_df.base
    intermediate_df.loc[base_is_ref_mask, 'ref_loglikelihood'] = np.log10((1 - error_rate))
    intermediate_df.loc[base_is_ref_mask, 'alt_loglikelihood'] = np.log10((error_rate))

    base_is_alt_mask = intermediate_df.ALT == intermediate_df.base
    intermediate_df.loc[base_is_alt_mask, 'alt_loglikelihood'] = np.log10((1 - error_rate))
    intermediate_df.loc[base_is_alt_mask, 'ref_loglikelihood'] = np.log10((error_rate))

    base_is_either_mask = (
        (intermediate_df.ALT == intermediate_df.base) |
        (intermediate_df.REF == intermediate_df.base)
    )
    intermediate_df.loc[base_is_either_mask, 'het_loglikelihood'] = np.log10(het_rate)

    # Subset to barcode-umi-pos where it's either ref, alt, or het
    intermediate_df = intermediate_df[base_is_either_mask]

    # Sum values over pos - but groupby more than that so we can use .sum()
    per_snp_loglik_grpby = intermediate_df.groupby(['barcode', 'chr', 'pos', 'REF', 'ALT'])  # move  these into index
    # if testing, can sum 'A C T G' as well
    per_snp_loglik = per_snp_loglik_grpby['ref_loglikelihood alt_loglikelihood het_loglikelihood'.split()].sum()
    per_snp_loglik = per_snp_loglik.reset_index()

    assert (refs_df.index == alt_df.index).all()
    # Reencode each donor as ref, alt, or het based on its copies of the reference allele
    replaced_ref = refs_df.replace({1: "het", 2: "ref", 0: "alt"})
    replaced_alt = alt_df.replace({1: "het", 2: "alt", 0: "ref"})
    # Fill values where we don't know the genotype (i.e. 0 copies of ref or alt) with 'none'
    donor_genotypes = replaced_alt[replaced_ref == replaced_alt].fillna('none')

    # add genotypes, and compute per cbc likelihoods
    donor_likelihoods = calculate_donor_liks(per_snp_loglik.merge(donor_genotypes.reset_index(), on='pos'),
                                             donors)

    # count UMIs, and num_snps == (UMIS x SNPs)
    umi_snp_counts = intermediate_df.groupby("barcode UMI".split()).pos.nunique()
    umi_snp_counts = umi_snp_counts.reset_index().groupby("barcode")
    donor_likelihoods['num_umis'] = umi_snp_counts.UMI.nunique()
    donor_likelihoods['num_snps'] = umi_snp_counts.pos.sum()
    return donor_likelihoods


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
    args = parser.parse_args()
    VCF_str = args.VCF_region
    simplified_region_name = args.region_name.replace(":", "_").replace("-", "_")

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
    # Generate probability table for each base
    #
    probs = [[.9, .1/3, .1/3, .1/3],
             [.1/3, .9, .1/3, .1/3],
             [.1/3, .1/3, .9, .1/3],
             [.1/3, .1/3, .1/3, .9]]
    probs = pd.DataFrame(index="A C G T".split(),
                         columns='prob_A prob_C prob_G prob_T'.split(),
                         data=probs)

    #
    # Delete all copies that map to different bases (e.g. TA (2 copies, diff bases), TTG (3 copies, diff bases))
    #
    single_base_mask = unique_read_counts.read.map(single_base)
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

    #
    # filter to donors
    #
    bcftools_proc = subprocess.Popen(f"bcftools query -s {','.join(donors)} -f %CHROM\t%POS\t%TYPE\t%REF\t%ALT[\t%GT]\n {VCF_str}".split(' '),
                                     stdout=subprocess.PIPE)

    genotypes = pd.read_table(bcftools_proc.stdout, header=None,
                              names="chrom pos type REF ALT".split() + donors)
    genotypes = genotypes.sort_values("chrom pos".split())
    genotypes = genotypes.set_index(['pos'])

    # drop INDELs
    genotypes = genotypes[genotypes.type != 'INDEL']
    # check number of repeat positions
    if genotypes.index.shape[0] > 100 and genotypes.index.duplicated(keep=False).sum() / genotypes.index.shape[0] > 0.05:
        raise ValueError("Genotype VCF contains too many (> 5%) repeat positions.")
    # drop duplicates
    genotypes = genotypes[~genotypes.index.duplicated(keep=False)]

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

    num_donors = len(donors)
    with gzip.open(f"barcode_log_likelihood_{simplified_region_name}.txt.gz", "wb") as outf:
        # Split read info into chunks with all the info for a group of CBCs to
        # save memory, then get per CBC likelihoods and write them out
        barcode_reads = single_base_uniq_reads.reset_index('barcode')
        all_cbcs = barcode_reads.barcode.sort_values().unique()
        first_block = True
        while len(all_cbcs) > 0:
            cur_cbcs = all_cbcs[:500]
            all_cbcs = all_cbcs[500:]
            print(len(all_cbcs))
            # get loglikelihood functions
            if args.likelihood_method == 'our_method':
                barcode_log_likelihood = generate_barcode_lls(barcode_reads[barcode_reads.barcode.isin(cur_cbcs)],
                                                              genotypes, donors, num_donors, ref_probs, alt_probs,
                                                              simplified_region_name)
            else:
                barcode_log_likelihood = dropulation_likelihoods(barcode_reads[barcode_reads.barcode.isin(cur_cbcs)],
                                                                 genotypes, refs_df, alt_df, donors)

            # final output is barcode_log_probs: [barcode] x [donor] loglikelihood
            # continuously add chunks of CBC likelihoods to the output file
            if barcode_log_likelihood is not None:
                barcode_log_likelihood.to_csv(outf, header=first_block, sep='\t')
            else:
                continue
            first_block = False

    print("Done.")


if __name__ == '__main__':
    main()
