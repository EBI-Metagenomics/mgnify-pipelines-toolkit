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
import os
import sys

from mgnify_pipelines_toolkit.constants.ncrna import TRNA


def main():
    parser = argparse.ArgumentParser(description="Script that extarcts tRNA from tRNAscan-SE")
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="The trnas_stats.out output file.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Name of the file to save the output to",
        required=False,
    )

    args = parser.parse_args()

    with open(args.input, "r") as input_handle:
        trnas = 0
        flag = 0
        for line in input_handle:
            if "Isotype / Anticodon" in line:
                flag = 1
            elif flag == 1:
                cols = line.split()
                if len(cols) > 1:
                    aa_pred = line.split(":")[0].split()[0]
                    counts = int(line.split(":")[1].split()[0])
                    if (aa_pred in TRNA or "Met" in aa_pred) and counts > 0:
                        trnas += 1

    new_line = "{name}\t{trnas}".format(name=os.path.basename(args.input).split("_stats")[0], trnas=trnas)
    with open(args.output or sys.stdout, "w") as file_out:
        file_out.write(new_line)


if __name__ == "__main__":
    main()
