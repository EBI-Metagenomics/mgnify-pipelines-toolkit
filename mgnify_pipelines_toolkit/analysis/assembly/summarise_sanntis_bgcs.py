#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024-2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import fileinput
import logging
import os.path

import pandas as pd

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s"
)

"""
Script parses sanntis GFF output and adds descriptions of annotated MIBiGs classes.

Descriptions were provided by Santiago Fragoso in doc with links to papers
https://docs.google.com/document/d/1IIGokcEpX_yLjg1H-w4_VAbekGRl6LgjSvafym3KMz0/edit?tab=t.0
Those were copied to TSV file constants/bgc_descriptions/sanntis.txt
"""


def parse_args():
    description = "Sanntis output summary."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--sanntis-gff", help="Sanntis gff", required=True)
    parser.add_argument(
        "-o", "--output", help="Sanntis summary output file.", required=True
    )
    parser.add_argument(
        "-d",
        "--descriptions",
        help="Sanntis descriptions file.",
        required=False,
        default="mgnify_pipelines_toolkit/constants/bgc_descriptions/sanntis.txt",
    )
    args = parser.parse_args()
    return args.sanntis_gff, args.output, args.descriptions


def main():
    input_gff, output_filename, descriptions = parse_args()
    dict_list = []
    with fileinput.hook_compressed(input_gff, "r") as file_in:
        for line in file_in:
            if line.startswith("#"):
                continue
            info = line.strip().split("\t")[8].split(";")
            entry_dict = {}
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
                "nearest_MiBIG": "nearest_MIBiG",
                "nearest_MiBIG_class": "nearest_MIBiG_class",
            }
        )
        df_grouped = (
            df.groupby(["nearest_MIBiG", "nearest_MIBiG_class"])
            .size()
            .reset_index(name="Count")
        )
        df_grouped = df_grouped.sort_values(by="Count", ascending=False)

        if descriptions and os.path.exists(descriptions):
            df_desc = pd.read_csv(descriptions, sep="\t", header=0)
            df_desc = df_desc.set_index("MIBiG_class")
            df_merged = df_grouped.merge(
                df_desc, left_on="nearest_MIBiG_class", right_index=True, how="left"
            )
            df_merged["Description"] = df_merged.apply(
                lambda row: row["nearest_MIBiG_class"].replace(
                    "NRP", df_desc.loc["NRP"]["Description"]
                )
                if pd.isna(row["Description"]) and "NRP" in row["nearest_MIBiG_class"]
                else row["Description"],
                axis=1,
            )
            df_merged = df_merged[
                ["nearest_MIBiG", "nearest_MIBiG_class", "Description", "Count"]
            ]
            df_merged.to_csv(output_filename, sep="\t", index=False)
        else:
            df_grouped.to_csv(output_filename, sep="\t", index=False)


if __name__ == "__main__":
    main()
