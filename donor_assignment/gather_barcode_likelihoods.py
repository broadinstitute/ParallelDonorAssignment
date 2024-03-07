import pandas as pd
import sys
import gzip

if __name__ == '__main__':
    barcode_likelihood_files = sys.argv[1]
    output_file_path = sys.argv[2]
    if len(sys.argv) == 4:
        min_umis = int(sys.argv[3])
    else:
        min_umis = 0

    total_df = 0
    for filepath in [ln.strip() for ln in open(barcode_likelihood_files)]:
        new_table = pd.read_table(filepath, index_col=0)
        total_df = new_table.add(total_df, fill_value=0)
    total_df = total_df[total_df.num_umis >= min_umis]
    total_df.to_csv(output_file_path, sep="\t")
