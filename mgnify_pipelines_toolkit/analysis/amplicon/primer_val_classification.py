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
from collections import defaultdict
import re

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

from mgnify_pipelines_toolkit.constants.var_region_coordinates import (
    REGIONS_16S_BACTERIA,
    REGIONS_16S_ARCHAEA,
    REGIONS_18S,
)

STRAND_FWD = "fwd"
STRAND_REV = "rev"


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="Path to cmsearch_deoverlap_tblout file",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        type=str,
        help="Path to concatenated primers fasta file",
    )
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument(
        "--se",
        action=argparse.BooleanOptionalAction,
        help="Flag for if run is single-end",
    )
    args = parser.parse_args()

    input = args.input
    fasta = args.fasta
    sample = args.sample
    single_end = args.se

    return input, fasta, sample, single_end


def get_amp_region(beg, end, strand, model):
    prev_region = ""

    margin = -10

    for region, region_coords in model.items():

        region_beg = region_coords[0]
        beg_diff = region_beg - beg
        end_diff = region_beg - end

        if strand == STRAND_FWD:
            if beg_diff >= margin and end_diff >= margin:
                return region
        else:
            if beg_diff >= margin and end_diff >= margin:
                return prev_region

        prev_region = region

    return prev_region


def main():

    input, fasta, sample, single_end = parse_args()
    res_dict = defaultdict(list)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    fwd_primers_fw = open("./fwd_primers.fasta", "w")
    rev_primers_fw = open("./rev_primers.fasta", "w")

    matched_primers_list = []

    with open(input, "r") as fr:
        for line in fr:
            line = line.strip()
            line = re.sub("[ \t]+", "\t", line)
            line_lst = line.split("\t")

            primer_name = line_lst[0]
            rfam = line_lst[3]
            beg = float(line_lst[5])
            end = float(line_lst[6])

            if "variant" not in primer_name:
                continue

            cleaned_primer_name = "_".join(primer_name.split("_")[0:-3])
            if cleaned_primer_name in matched_primers_list:
                continue

            if rfam == "RF00177":
                gene = "16S"
                model = REGIONS_16S_BACTERIA
            elif rfam == "RF01959":
                gene = "16S"
                model = REGIONS_16S_ARCHAEA
            elif rfam == "RF01960":
                gene = "18S"
                model = REGIONS_18S
            else:  # For cases when it's a std primer but for some reason hasn't matched the model
                if cleaned_primer_name == "F_auto" or cleaned_primer_name == "R_auto":
                    continue
                gene = "Unknown"
                amp_region = "Unknown"
                model = ""

            res_dict["Run"].append(sample)
            res_dict["AssertionEvidence"].append("ECO_0000363")
            res_dict["AssertionMethod"].append("automatic assertion")

            strand = ""

            if primer_name[-1] == "F":
                strand = STRAND_FWD
            elif primer_name[-1] == "R":
                strand = STRAND_REV
            else:
                print(f"Not sure what strand this is, exiting: {primer_name}")

            if model:
                amp_region = get_amp_region(beg, end, strand, model)

            primer_seq = str(fasta_dict[cleaned_primer_name].seq)

            res_dict["Gene"].append(gene)
            res_dict["VariableRegion"].append(amp_region)
            res_dict["PrimerName"].append(cleaned_primer_name)
            res_dict["PrimerStrand"].append(strand)
            res_dict["PrimerSeq"].append(primer_seq)

            if strand == STRAND_FWD:
                fwd_primers_fw.write(f">{cleaned_primer_name}\n{primer_seq}\n")
            elif strand == STRAND_REV:
                if single_end:
                    primer_seq = Seq(primer_seq).reverse_complement()
                rev_primers_fw.write(f">{cleaned_primer_name}\n{primer_seq}\n")

            matched_primers_list.append(cleaned_primer_name)

    res_df = pd.DataFrame.from_dict(res_dict)
    res_tsv_name = f"./{sample}_primer_validation.tsv"
    res_df.to_csv(res_tsv_name, sep="\t", index=False) if not res_df.empty else open(res_tsv_name, "w").close()

    fwd_primers_fw.close()
    rev_primers_fw.close()


if __name__ == "__main__":
    main()
