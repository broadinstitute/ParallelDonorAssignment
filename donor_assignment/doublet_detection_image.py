import pandas as pd
import matplotlib.pyplot as plt
import pylab
import argparse

params = {'legend.fontsize': '40',
          'figure.figsize': (10, 10),
          'axes.labelsize': '40',
          'axes.titlesize': '50',
          'xtick.labelsize': '40',
          'ytick.labelsize': '40',
          'axes.linewidth': '0.5',
          'pdf.fonttype': '42',
          'font.sans-serif': 'Helvetica'}
pylab.rcParams.update(params)



def main():
    parser = argparse.ArgumentParser(description="Generate a likelihood/umi vs umi count plot for detecting doublets")
    parser.add_argument("cell_donor_likelihoods", type=str)
    parser.add_argument("donor_names", type=str)
    args = parser.parse_args()

    cell_donor_lls = pd.read_table(args.cell_donor_likelihoods)
    full_donor_list = pd.read_csv(args.donor_names, header=None).to_numpy()
    new_donors = cell_donor_lls.columns.intersection(full_donor_list.reshape(-1))

    sample_best_LL = cell_donor_lls.set_index('barcode')[new_donors].agg(['max', 'idxmax'], 
                                                                                    axis='columns').rename(columns = {'max':'bestLikelihood', 'idxmax':'bestSample'})
    sample_best_LL = sample_best_LL.merge(cell_donor_lls.set_index('barcode')[['num_umis', 'num_snps']], left_index=True, right_index=True)
    sample_best_LL['LogLikperUMI'] = sample_best_LL.bestLikelihood / sample_best_LL.num_umis

    fig, ax = plt.subplots(figsize=(9, 7))
    hb = ax.hexbin(sample_best_LL.LogLikperUMI, sample_best_LL.num_umis, yscale='log', cmap='inferno', bins='log')
    fig.colorbar(hb, ax=ax, label='# of cells')
    ax.set_xlabel('LogLik per UMI')
    ax.set_ylabel('num UMIs')

    plt.tight_layout()
    plt.savefig('test_fig.png', dpi=200)

if __name__ == '__main__':
    main()