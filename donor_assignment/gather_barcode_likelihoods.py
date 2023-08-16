import pandas as pd
import sys
import gzip


def line_iterator(file):
    for chunk in pd.read_table(file, chunksize=1000, index_col=0):
        for barcode, vals in chunk.iterrows():
            yield barcode
            yield vals


def combine_iterators(iters):
    pending_barcodes = {}

    def bump_iter(iter):
        next_bc = next(iter, None)
        if next_bc is None:
            return
        if next_bc in pending_barcodes:
            pending_barcodes[next_bc].append(iter)
        else:
            pending_barcodes[next_bc] = [iter]

    # get the first barcode for each iterator
    for iter in iters:
        bump_iter(iter)

    # find the lowest pending barcode and emit it
    while pending_barcodes:
        lowest_bc = min(pending_barcodes)
        lowest_iterators = pending_barcodes.pop(lowest_bc)
        vals = [next(iter) for iter in lowest_iterators]
        yield lowest_bc, sum(vals)
        # re-activate those iters
        [bump_iter(iter) for iter in lowest_iterators]


if __name__ == '__main__':
    barcode_likelihood_files = sys.argv[1]
    output_file_path = sys.argv[2]
    per_line_iterators = [line_iterator(ln.strip()) for ln in open(barcode_likelihood_files)]
    with gzip.open(output_file_path, "w") as outf:
        for bc, vals in combine_iterators(per_line_iterators):
            outf.write(bc + "\t" + "\t".join(vals.astype(str).values) + "\n")
