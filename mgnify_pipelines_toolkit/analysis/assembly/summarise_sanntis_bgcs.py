#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
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
import fileinput
import logging
import pandas as pd

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s"
)

SANNTIS_VERSION = "0.9.4.1"

f"""
Script parses SanntiS GFF output and adds descriptions of annotated MIBiGs classes.
Descriptions were pre-parsed for version {SANNTIS_VERSION} and stored as a dictionary.
"""

DESCRIPTIONS = {
    "Polyketide": "Built from iterative condensation of acetate units derived from acetyl-CoA",
    "Terpene": "Composed of isoprene (C5) units derived from isopentenyl pyrophosphate",
    "Alkaloid": "Nitrogen-containing compounds derived from amino acids (e.g., ornithine, lysine, tyrosine, tryptophan)",
    "RiPP": "Ribosomally synthesised and Post-translationally modified Peptide",
    "NRP": "Nonribosomal Peptide",
    "Saccharide": "Carbohydrate-based natural products (e.g., aminoglycoside antibiotics)",
    "Other": "Catch-all class for clusters encoding metabolites outside main classes (e.g., cyclitols, indolocarbazoles, and phosphonates)",
}


def parse_args():
    description = (
        "Sanntis output summary generator. "
        "Script takes SanntiS GFF and counts pairs of (nearest_MiBIG, nearest_MiBIG_class)."
        "It also adds pre-parsed descriptions of classes stored in that script as a dictionary. "
        f"Descriptions were taken from SanntiS docs v{SANNTIS_VERSION}."
    )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--sanntis-gff", help="SanntiS GFF", required=True)
    parser.add_argument(
        "-o", "--output", help="SanntiS summary TSV output file.", required=True
    )
    args = parser.parse_args()
    return args.sanntis_gff, args.output


def main():
    input_gff, output_filename = parse_args()
    dict_list = []
    with fileinput.hook_compressed(input_gff, "r") as file_in:
        # TODO: to be merged with the GFF toolkit
        for line in file_in:
            if line.startswith("#"):
                continue
            info = line.strip().split("\t")[8].split(";")
            entry_dict = {}
            # TODO: merge this with the GFF toolkit GFF reader
            for pair in info:
                key, value = pair.split(
                    "=", 1
                )  # Ensure split only occurs at the first '=' occurrence
                entry_dict[key] = value
            dict_list.append(entry_dict)

        # Convert to DataFrame
        df = pd.DataFrame(dict_list)
        df = df.rename(
            columns={
                "nearest_MiBIG": "nearest_mibig",
                "nearest_MiBIG_class": "nearest_mibig_class",
            }
        )
        df_grouped = (
            df.groupby(["nearest_mibig", "nearest_mibig_class"])
            .size()
            .reset_index(name="count")
        )
        df_grouped = df_grouped.sort_values(by="count", ascending=False)

        df_desc = pd.DataFrame(
            list(DESCRIPTIONS.items()), columns=["mibig_class", "description"]
        )
        df_desc = df_desc.set_index("mibig_class")
        df_merged = df_grouped.merge(
            df_desc, left_on="nearest_mibig_class", right_index=True, how="left"
        )
        df_merged["description"] = df_merged.apply(
            lambda row: (
                row["nearest_mibig_class"].replace(
                    "NRP", df_desc.loc["NRP"]["description"]
                )
                if pd.isna(row["description"]) and "NRP" in row["nearest_mibig_class"]
                else row["description"]
            ),
            axis=1,
        )
        df_merged = df_merged[
            ["nearest_mibig", "nearest_mibig_class", "description", "count"]
        ]
        df_merged = df_merged.rename(columns={
            "Description": "description",
            "Count": "count"
        })
        df_merged.to_csv(output_filename, sep="\t", index=False)


if __name__ == "__main__":
    main()
