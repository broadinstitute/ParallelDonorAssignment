import pandas as pd
import matplotlib.pyplot as plt
import pylab
import argparse
import numpy as np

params = {
    "legend.fontsize": "40",
    "figure.figsize": (10, 10),
    "axes.labelsize": "40",
    "axes.titlesize": "50",
    "xtick.labelsize": "40",
    "ytick.labelsize": "40",
    "axes.linewidth": "0.5",
    "pdf.fonttype": "42",
    "font.sans-serif": "Helvetica",
}
pylab.rcParams.update(params)


def generate_loglik_per_umi_fig(sample_best_LL, thresh, prefix=None):
    fig, ax = plt.subplots(figsize=(9, 7))
    hb = ax.hexbin(
        sample_best_LL.LogLikperUMI,
        sample_best_LL.num_umis,
        yscale="log",
        cmap="inferno",
        bins="log",
    )
    fig.colorbar(hb, ax=ax, label="# of cells")
    ax.set_xlabel("LogLik per UMI")
    ax.set_ylabel("num UMIs")
    ax.axvline(thresh, c="r")
    
    plt.xticks(rotation=30)
    plt.tight_layout()
    if prefix is None:
        plt.savefig("loglik_per_umi_plot.png", dpi=200)
    else:
        plt.savefig(prefix + ".loglik_per_umi_plot.png", dpi=200)


def threshold_otsu(x, *args, **kwargs):
    """Find the threshold value for a bimodal histogram using the Otsu method.

    If you have a distribution that is bimodal (AKA with two peaks, with a valley
    between them), then you can use this to find the location of that valley, that
    splits the distribution into two.

    From the SciKit Image threshold_otsu implementation:
    https://github.com/scikit-image/scikit-image/blob/70fa904eee9ef370c824427798302551df57afa1/skimage/filters/thresholding.py#L312
    """
    counts, bin_edges = np.histogram(x, *args, **kwargs)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2

    # class probabilities for all possible thresholds
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(counts * bin_centers) / weight1
    mean2 = (np.cumsum((counts * bin_centers)[::-1]) / weight2[::-1])[::-1]

    # Clip ends to align class 1 and class 2 variables:
    # The last value of ``weight1``/``mean1`` should pair with zero values in
    # ``weight2``/``mean2``, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    idx = np.argmax(variance12)
    threshold = bin_centers[idx]
    return threshold


def get_singlets(sample_best_LL, thresh, prefix=None):

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.hist(sample_best_LL.LogLikperUMI, bins=30)
    ax.axvline(thresh, c="r")
    ax.set_xlabel("LogLikperUMI")
    ax.annotate(
        f"Threshold: {round(thresh,3)}",
        (0.1, 0.9),
        xycoords="axes fraction",
        fontsize=25,
    )
    plt.tight_layout()

    if prefix is None:
        plt.savefig("loglik_per_umi_histogram.png", dpi=200)
    else:
        plt.savefig(prefix + ".loglik_per_umi_histogram.png", dpi=200)

    singlets = sample_best_LL[sample_best_LL.LogLikperUMI > thresh]
    # save singlets df
    if prefix is None:
        singlets.reset_index()["barcode bestSample".split()].to_csv(
            "singlets.txt", index=None, sep="\t")
    else:
        singlets.reset_index()["barcode bestSample".split()].to_csv(
            prefix + ".singlets.txt", index=None, sep="\t")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a likelihood/umi vs umi count plot for detecting doublets"
    )
    parser.add_argument("cell_donor_likelihoods", type=str)
    parser.add_argument("donor_names", type=str)
    parser.add_argument("--threshold", type=float, help='optional argument for pre-set singlet threshold')
    parser.add_argument("--prefix", type=str, help='prefix for output files (implicit . added after prefix)')
    args = parser.parse_args()

    cell_donor_lls = pd.read_table(args.cell_donor_likelihoods)
    full_donor_list = pd.read_csv(args.donor_names, header=None).to_numpy()
    new_donors = cell_donor_lls.columns.intersection(full_donor_list.reshape(-1))

    sample_best_LL = (
        cell_donor_lls.set_index("barcode")[new_donors]
        .agg(["max", "idxmax"], axis="columns")
        .rename(columns={"max": "bestLikelihood", "idxmax": "bestSample"})
    )
    sample_best_LL = sample_best_LL.merge(
        cell_donor_lls.set_index("barcode")[["num_umis", "num_snps"]],
        left_index=True,
        right_index=True,
    )
    sample_best_LL["LogLikperUMI"] = (
        sample_best_LL.bestLikelihood / sample_best_LL.num_umis
    )

    if args.threshold is None:
        thresh = threshold_otsu(sample_best_LL.LogLikperUMI)
    else:
        thresh = args.threshold

    generate_loglik_per_umi_fig(sample_best_LL, thresh, args.prefix)
    get_singlets(sample_best_LL, thresh, args.prefix)


if __name__ == "__main__":
    main()
