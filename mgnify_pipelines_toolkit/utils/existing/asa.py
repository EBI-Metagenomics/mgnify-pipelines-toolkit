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

import pyfastx
import argparse
import csv


def rename_fasta_sequences(input_fasta, output_fasta, prefix, mapping_file):
    fasta = pyfastx.Fasta(input_fasta, build_index=False)

    with (
        open(output_fasta, "w") as out_fasta,
        open(mapping_file, "w", newline="") as map_file,
    ):
        writer = csv.writer(map_file, delimiter="\t")
        writer.writerow(["original", "renamed"])

        for index, (name, seq) in enumerate(fasta, start=1):
            new_name = f"{prefix}{index}"
            out_fasta.write(f">{new_name}\n{seq}\n")
            writer.writerow([name, new_name])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rename FASTA sequences with a prefix and auto-incrementing index."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument(
        "-o", "--output", required=True, help="Output renamed FASTA file"
    )
    parser.add_argument(
        "-p", "--prefix", required=True, help="Prefix for new sequence names"
    )
    parser.add_argument(
        "-m",
        "--mapping",
        required=True,
        help="Output mapping file for old and new names",
    )

    args = parser.parse_args()

    rename_fasta_sequences(args.input, args.output, args.prefix, args.mapping)