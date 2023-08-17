import sys
import subprocess
import pandas as pd
from scipy.optimize import nnls

donors_file = sys.argv[1]
genotypes = sys.argv[2]
BAM = sys.argv[3]
min_coverage = int(sys.argv[4])

#
# load donors
#
donors = [ln.strip() for ln in open(donors_file)]

#
# load genotypes
#
genotypes = pd.read_table(genotypes, header=None, names="chrom pos type REF ALT".split() + donors)
genotypes = genotypes.sort_values("chrom pos".split())

#
# extract autosomes, make sure chromosomes start with chr
#
genotypes.chrom = genotypes.chrom.astype(str).str.replace('chr', '')
mask = genotypes.chrom.str.isnumeric()
genotypes = genotypes[mask]
genotypes.chrom = 'chr' + genotypes.chrom
genotypes = genotypes.reset_index(drop=True)
print(f"Computing pileup for {len(genotypes)} sites.")

#
# write positions for pileup
#
genotypes['chrom pos'.split()].to_csv("positions.txt", sep="\t", header=None, index=None)

#
# run pileup
#
mpileup_proc = subprocess.Popen(f"samtools mpileup -l positions.txt --no-output-ins {BAM}".split(' '),
                                stdout=subprocess.PIPE)

pileup = pd.read_table(mpileup_proc.stdout, header=None,
                       names="chrom pos ign cov bases qualities".split())
# recompute pileup coverage
pileup['truecov'] = pileup.bases.str.upper().str.count('[ACGT]')
print(f"Found reads over {len(pileup)} sites.")

# filter to minimum coverage
pileup = pileup[pileup.truecov >= min_coverage]
print(f"Filtered to {len(pileup)} sites with at least {min_coverage} reads.")

#
# filter to overlapped genotypes
#
pileup = pileup.set_index("chrom pos".split())
genotypes = genotypes.set_index("chrom pos".split()).loc[pileup.index]

#
# get dosages from genotypes
#
allele_counts_REF = pd.DataFrame(index=genotypes.index)
allele_counts_ALT = pd.DataFrame(index=genotypes.index)

for donor in genotypes.columns[3:]:
    allele_counts_REF[donor] = genotypes[donor].str.count('0')
    allele_counts_ALT[donor] = genotypes[donor].str.count('1')

#
# filter out lines with missing genotypes
#
mask = (allele_counts_REF + allele_counts_ALT).min(axis=1) < 2
print(f"Removing {mask.sum()} of {len(mask)} sites for missing genotypes.")
pileup = pileup[~ mask]
genotypes = genotypes[~ mask]
allele_counts_REF = allele_counts_REF[~ mask]
allele_counts_ALT = allele_counts_ALT[~ mask]

#
# filter out lines with all donors REF or ALT
#
mask = (allele_counts_REF.min(axis=1) == 2) | (allele_counts_REF.min(axis=1) == 2)
print(f"Removing {mask.sum()} of {len(mask)} sites for idenical genotypes.")
pileup = pileup[~ mask]
genotypes = genotypes[~ mask]
allele_counts_REF = allele_counts_REF[~ mask]
allele_counts_ALT = allele_counts_ALT[~ mask]

#
# count each of ACGT in pileup
#
bases = pileup.bases.str.upper().str
for base in 'ACGT':
    pileup['count_' + base] = bases.count(base)

#
# get reference and alternate read counts in pileup, along with ALT_frac
#
pileup['count_REF'] = 0
pileup['count_ALT'] = 0
for base in 'ACGT':
    pileup['count_REF'] += pileup['count_' + base] * (genotypes.REF == base)
    pileup['count_ALT'] += pileup['count_' + base] * (genotypes.ALT == base)
pileup['ALT_frac'] = (pileup.count_ALT) / (pileup.count_REF + pileup.count_ALT)    
pileup['regularized_ALT_frac'] = (pileup.count_ALT + 1) / (pileup.count_REF + pileup.count_ALT + 2)

#
# set up regression problem
#
good = (pileup.regularized_ALT_frac < 0.9) & (pileup.regularized_ALT_frac > 0.1)  # from cisVar
y = pileup.regularized_ALT_frac[good].values
X = allele_counts_ALT.loc[good].values
weights, residual = nnls(X, y)
print(f"Solved with RMSE = {residual / len(y)}")
print(f"Sum of weights = {weights.sum()}.  Should be near 0.5.")
weights = weights / weights.sum()
with open("donor_weights.txt", "w") as weights_file:
    for k, w in zip(allele_counts_ALT.columns, weights):
        weights_file.write(f"{k}\t{w}\n")


#
# for each base with min_coverage, compute Expected, Observed, and Observed_se
#
X = allele_counts_ALT.values
pileup['pred_regularized_ALT_frac'] = (X @ weights) / 2  # we rescaled weights above
pileup.to_csv("cisvar_estimates.txt.gz", sep="\t")
