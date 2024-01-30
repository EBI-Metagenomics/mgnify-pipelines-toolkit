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

from Bio import Seq, SeqIO

def parse_args(argv=None):

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="Path to finalised primer list fasta file")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path")
    args = parser.parse_args(argv)
  
    _INPUT = args.input
    _SAMPLE = args.sample
    _OUTPUT = args.output

    return _INPUT, _SAMPLE, _OUTPUT

def main(argv=None):
    
    _INPUT, _SAMPLE, _OUTPUT = parse_args(argv)

    primers_dict = SeqIO.to_dict(SeqIO.parse(_INPUT, "fasta"))
    
    for primer_key in primers_dict.keys():

        primer = primers_dict[primer_key]
        primer_name = primer.name

        if "R" in primer_name:
            primers_dict[primer_key].seq = primer.seq.reverse_complement()

    SeqIO.write(primers_dict.values(), f"{_OUTPUT}/{_SAMPLE}_rev_comp_se_primers.fasta", "fasta")


if __name__ == "__main__":
    main()