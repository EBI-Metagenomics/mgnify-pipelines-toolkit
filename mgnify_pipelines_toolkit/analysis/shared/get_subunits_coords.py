#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024 EMBL - European Bioinformatics Institute
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
import sys
import os
import gzip
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="Extract lsu, ssu and 5s")
    parser.add_argument("-i", "--input", dest="input", help="Input fasta file", required=True)
    parser.add_argument("-l", "--lsu", dest="lsu", help="LSU pattern", required=True)
    parser.add_argument("-s", "--ssu", dest="ssu", help="SSU pattern", required=True)

    SSU_coords = "SSU_coords"
    LSU_coords = "LSU_coords"
    SSU_count = 0
    LSU_count = 0

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()

        with open(SSU_coords, 'w') as out_ssu, open(LSU_coords, 'w') as out_lsu, open(args.input, 'r') as input:
            for line in input:
                if args.lsu in line:
                    out_lsu.write(line)
                    LSU_count += 1
                elif args.ssu in line:
                    out_ssu.write(line)
                    SSU_count += 1
        with open("RNA-counts", 'w') as count:
            count.write("LSU count\t" + str(LSU_count) + "\nSSU count\t" + str(SSU_count))

    out_ssu.close()
    out_lsu.close()
    count.close()

if __name__ == "__main__":
    main()
