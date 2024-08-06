
import argparse

import pandas as pd
import numpy as np

from mgnify_pipelines_toolkit.constants.thresholds import MIN_AMPLICON_STRATEGY_CHECK

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="Input")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output")

    args = parser.parse_args()

    INPUT = args.input
    SAMPLE = args.sample
    OUTPUT = args.output

    return INPUT, SAMPLE, OUTPUT


def main():
    
    INPUT, SAMPLE, OUTPUT = parse_args()

    cons_df = pd.read_csv(INPUT, sep='\t')

    cons_values = cons_df.values[0][1:]
    mean_cons = np.mean(cons_values)

    fw = open(f"{OUTPUT}/{SAMPLE}_library_check_out.txt", "w")

    if mean_cons >= MIN_AMPLICON_STRATEGY_CHECK:
        print("This data is likely to be AMPLICON.")
        fw.write("AMPLICON") # File with "AMPLICON" written as a result.

    else:
        print("This data is unlikely to be AMPLICON.")
        # If unlikely to be AMPLICON, the output file will be empty.

    fw.close()

if __name__ == "__main__":
    main()