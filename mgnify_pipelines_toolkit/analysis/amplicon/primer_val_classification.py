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
import logging
import re
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from mgnify_pipelines_toolkit.constants.var_region_coordinates import (
    REGIONS_16S_ARCHAEA,
    REGIONS_16S_BACTERIA,
    REGIONS_18S,
)

STRAND_FWD = "fwd"
STRAND_REV = "rev"

logging.basicConfig(level=logging.INFO)


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


def get_amp_region(primer_beg: float, primer_end: float, strand: str, model: dict) -> str:
    prev_region = ""

    # some valid primers go inside HV regions a little bit, this margin is to account for that
    margin = -10

    for region, region_coords in model.items():
        # get current region start and end coordinates
        region_beg = region_coords[0]
        region_end = region_coords[1]

        # compute where primer beginning is in relation to current region
        region_beg_primer_beg_diff = region_beg - primer_beg
        region_beg_primer_end_diff = region_beg - primer_end
        primer_beg_near_region_start = region_beg_primer_beg_diff >= margin
        primer_end_near_region_start = region_beg_primer_end_diff >= margin

        # compute where primer end is in relation to current region
        region_end_primer_beg_diff = region_end - primer_beg
        region_end_primer_end_diff = region_end - primer_end
        primer_beg_before_region_end = region_end_primer_beg_diff >= margin
        primer_end_before_region_end = region_end_primer_end_diff >= margin

        if primer_beg_near_region_start and primer_end_near_region_start:
            # if both these statements are true then primer is before a HV region
            # i.e. validation = true
            if strand == STRAND_FWD:
                return region
            else:
                # if primer strand is REV then we return the previous region
                return prev_region
        elif primer_beg_before_region_end and primer_end_before_region_end:
            # if the previous if statement is FALSE
            # AND if both these statements are true then primer is within a HV region
            # i.e. validation = false
            logging.warning(f"This primer is within HV region {region}: {str(int(primer_beg))}-{str(int(primer_end))} vs {region_beg}-{region_end}")
            return ""
        # keep iterating through HV regions otherwise

        prev_region = region

    return prev_region


def main():
    input, fasta, sample, single_end = parse_args()
    res_dict = defaultdict(list)

    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    logging.info(f"Total primers read (including permutations): {len(fasta_dict)}")

    fwd_primers_fw = open("./fwd_primers.fasta", "w")
    rev_primers_fw = open("./rev_primers.fasta", "w")

    matched_primers_list = []

    with open(input, "r") as fr:
        logging.info(f"Reading deoverlap file: {input}")
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

            strand = ""

            if primer_name[-1] == "F":
                strand = STRAND_FWD
            elif primer_name[-1] == "R":
                strand = STRAND_REV
            else:
                logging.warning(f"Not sure what strand this is, skipping: {primer_name}")
                continue

            if model:
                logging.info(f"Checking match coordinates for primer {primer_name}")
                amp_region = get_amp_region(beg, end, strand, model)

            if not amp_region:
                logging.warning(f"Primer validation failed for {primer_name}, skipping")
                continue

            primer_seq = str(fasta_dict[cleaned_primer_name].seq)

            res_dict["Run"].append(sample)
            res_dict["AssertionEvidence"].append("ECO_0000363")
            res_dict["AssertionMethod"].append("automatic assertion")
            res_dict["Gene"].append(gene)
            res_dict["VariableRegion"].append(amp_region)
            res_dict["PrimerName"].append(cleaned_primer_name)
            res_dict["PrimerStrand"].append(strand)

            if strand == STRAND_FWD:
                fwd_primers_fw.write(f">{cleaned_primer_name}\n{primer_seq}\n")
            elif strand == STRAND_REV:
                if single_end:
                    primer_seq = Seq(primer_seq).reverse_complement()
                rev_primers_fw.write(f">{cleaned_primer_name}\n{primer_seq}\n")

            res_dict["PrimerSeq"].append(primer_seq)

            matched_primers_list.append(cleaned_primer_name)
            logging.info(f"Added {cleaned_primer_name} to list of matched primers")

    res_tsv_name = f"./{sample}_primer_validation.tsv"
    if res_dict:
        res_df = pd.DataFrame.from_dict(res_dict)
        res_df.to_csv(res_tsv_name, sep="\t", index=False) if not res_df.empty else open(res_tsv_name, "w").close()
        logging.info(f"{len(res_df)} primers validated, generating output")

    else:
        logging.warning("No primers were successfully validated, generating empty outputs")
        primer_val_fw = open(res_tsv_name, "w")
        primer_val_fw.close()

    fwd_primers_fw.close()
    rev_primers_fw.close()


if __name__ == "__main__":
    main()
