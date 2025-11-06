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
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compare predicted proteins (transcripts) between BRAKER and MetaEuk GFF files using exon/CDS overlap."
    )
    parser.add_argument('--fasta_aa', type=str, required=True, help="BRAKER protein fasta")
    parser.add_argument('--fasta_fn', type=str, required=True, help="BRAKER nucleotide fasta")
    parser.add_argument('--output', type=str, required=True, help="Output prefix (sample name)")
    return parser.parse_args()


def rename_fasta_ids(fasta_file, output, suffix):
    '''
    Rename BRAKER protein fasta IDs to append '|braker_only' suffix
    '''
    records = SeqIO.parse(fasta_file, "fasta")

    for rec in records:
        new_id = f"{rec.id}|braker_only"
        rec.id = new_id

    SeqIO.write(records, output + suffix, "fasta" )

def main():
    args = parse_args()
    print(f"Renaming {args.output}...")

    rename_fasta_ids(args.fast_aa, args.output, "_braker_unique.faa")
    rename_fasta_ids(args.fasta_fn, args.output, "_braker_unique.ffn")

    print("Renaming done")


if __name__ == "__main__":
    main()

