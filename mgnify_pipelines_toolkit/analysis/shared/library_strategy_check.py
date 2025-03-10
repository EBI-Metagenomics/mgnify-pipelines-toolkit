#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024-2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse

import pandas as pd
import numpy as np

from mgnify_pipelines_toolkit.constants.thresholds import MIN_AMPLICON_STRATEGY_CHECK


def parse_args():

    parser = argparse.ArgumentParser(
        description="Script that checks the output of assess_mcp_proportions.py to guess whether a FASTQ file is AMPLICON or NOT AMPLICON."
    )
    parser.add_argument("-i", "--input", required=True, type=str, help="Input")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output")

    args = parser.parse_args()

    input = args.input
    sample = args.sample
    output = args.output

    return input, sample, output


def main():

    input, sample, output = parse_args()

    cons_df = pd.read_csv(input, sep="\t")

    cons_values = cons_df.values[0][1:]
    mean_cons = np.mean(cons_values)

    fw = open(f"{output}/{sample}_library_check_out.txt", "w")

    if mean_cons >= MIN_AMPLICON_STRATEGY_CHECK:
        print("This data is likely to be AMPLICON.")
        fw.write("AMPLICON")  # File with "AMPLICON" written as a result.

    else:
        print("This data is unlikely to be AMPLICON.")
        # If unlikely to be AMPLICON, the output file will be empty.

    fw.close()


if __name__ == "__main__":
    main()
