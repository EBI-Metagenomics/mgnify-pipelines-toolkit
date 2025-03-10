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
    args = parser.parse_args()

    input = args.input
    fasta = args.fasta
    sample = args.sample

    return input, fasta, sample


def get_amp_region(beg, end, strand, model):
    prev_region = ""

    for region, region_coords in model.items():

        region_beg = region_coords[0]
        beg_diff = region_beg - beg
        end_diff = region_beg - end

        if strand == STRAND_FWD:
            if beg_diff > 0 and end_diff > 0:
                return region
        else:
            if beg_diff > 0 and end_diff > 0:
                return prev_region

        prev_region = region

    return prev_region


def main():
    input, fasta, sample = parse_args()
    res_dict = defaultdict(list)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    with open(input, "r") as fr:
        for line in fr:
            line = line.strip()
            line = re.sub("[ \t]+", "\t", line)
            line_lst = line.split("\t")

            primer_name = line_lst[0]
            rfam = line_lst[3]
            beg = float(line_lst[5])
            end = float(line_lst[6])

            if rfam == "RF00177":
                gene = "16S"
                model = REGIONS_16S_BACTERIA
            elif rfam == "RF01959":
                gene = "16S"
                model = REGIONS_16S_ARCHAEA
            elif rfam == "RF01960":
                gene = "18S"
                model = REGIONS_18S
            else:
                continue

            res_dict["Run"].append(sample)
            res_dict["AssertionEvidence"].append("ECO_0000363")
            res_dict["AssertionMethod"].append("automatic assertion")

            strand = ""

            if "F" in primer_name:
                strand = STRAND_FWD
            elif "R" in primer_name:
                strand = STRAND_REV

            amp_region = get_amp_region(beg, end, strand, model)
            primer_seq = str(fasta_dict[primer_name].seq)

            res_dict["Gene"].append(gene)
            res_dict["VariableRegion"].append(amp_region)
            res_dict["PrimerName"].append(primer_name)
            res_dict["PrimerStrand"].append(strand)
            res_dict["PrimerSeq"].append(primer_seq)

    res_df = pd.DataFrame.from_dict(res_dict)
    res_df.to_csv(f"./{sample}_primer_validation.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
