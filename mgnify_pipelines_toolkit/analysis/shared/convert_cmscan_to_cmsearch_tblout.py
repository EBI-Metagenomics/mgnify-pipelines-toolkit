#!/usr/bin/env python3
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

"""
    Script to convert cmscan-table to cmsearch-table (swap columns 1 and 2 with 3 and 4)

    input example:
    #target name         accession query name           accession mdl mdl from   mdl to seq from   seq to strand ..
    #------------------- --------- -------------------- --------- --- -------- -------- -------- -------- ------ ..
    SSU_rRNA_eukarya     RF01960   SRR17062740.1        -          cm      582     1025        1      452      + ..

    expected output:
    #------------------- --------- -------------------- --------- --- -------- -------- -------- -------- ------ ..
    #target name         accession query name           accession mdl mdl from   mdl to seq from   seq to strand ..
    SRR17062740.1        -         SSU_rRNA_eukarya     RF01960    cm      582     1025        1      452      + ..

"""

import sys
import argparse
import fileinput
from itertools import accumulate


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Convert cmscan table to cmsearch table"
    )
    parser.add_argument(
        "-i", "--input", dest="input", help="Input cmscan file", required=True
    )
    parser.add_argument(
        "-o", "--output", dest="output", help="Output filename", required=True
    )
    return parser.parse_args(argv)


class TableModifier:
    def __init__(
        self,
        input_file: str,
        output_file: str,
    ):
        """
        Output of cmsearch-table has columns separated with different number of spaces (to keep humanreadable format)
        :param input_file: output of cmscan-table
        :param output_file: name of cmsearch table
        """
        self.input_file = input_file
        self.output_file = output_file

    def modify_table(self):
        with (
            fileinput.hook_compressed(
                self.input_file, "r", encoding="utf-8"
            ) as file_in,
            open(self.output_file, "w") as file_out,
        ):
            header_written = False
            separator_line, header = "", ""
            for line in file_in:
                if line.startswith("#"):
                    if "--" in line:
                        separator_line = line.split(" ")
                        separator_line[0] = separator_line[0].replace("#", "-")
                        lengths = [0] + list(
                            accumulate(len(s) + 1 for s in separator_line)
                        )
                    else:
                        header = line
                else:
                    coord_to_keep = len(" ".join(separator_line[0:4]))
                    if not header_written:
                        file_out.write(header)
                        file_out.write(
                            " ".join(
                                [
                                    "#" + separator_line[2][1:],
                                    separator_line[3],
                                    separator_line[0].replace("#", ""),
                                    separator_line[1],
                                ]
                                + separator_line[4:]
                            )
                        )
                        header_written = True
                    new_line = (
                        line[lengths[2] : lengths[3]]
                        + line[lengths[3] : lengths[4]]
                        + line[lengths[0] : lengths[1]]
                        + line[lengths[1] : lengths[2]]
                        + line[coord_to_keep + 1 :]
                    )
                    file_out.write(new_line)


def main():
    args = parse_args(sys.argv[1:])
    table_modifier = TableModifier(
        input_file=args.input,
        output_file=args.output,
    )
    table_modifier.modify_table()


if __name__ == "__main__":
    main()
