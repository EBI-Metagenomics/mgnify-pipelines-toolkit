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
import pandas as pd

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s"
)


def parse_args():
    description = "AntiSMASH output summary."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--antismash-gff", help="Antismash gff", required=True)
    parser.add_argument(
        "-o", "--output", help="Antisamsh summary output file.", required=True
    )
    args = parser.parse_args()
    return args.antismash_gff, args.output


def main():
    input_gff, output_filename = parse_args()
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
        df = df[df["product"].notna()]
        df_grouped = (
            df.groupby(["product"]).size().reset_index(name="Count")
        ).sort_values(by="Count", ascending=False)

        df_grouped = df_grouped.rename(
            columns={
                "product": "ClassID",
            }
        )

        df_grouped.to_csv(output_filename, sep="\t", index=False)


if __name__ == "__main__":
    main()
