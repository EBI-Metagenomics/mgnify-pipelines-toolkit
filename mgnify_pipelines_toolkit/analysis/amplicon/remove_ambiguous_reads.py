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

# Script removes any reads with ambiguous bases (Ns) for the purpose of DADA2

import argparse
import gzip

from Bio import SeqIO, bgzf

def parse_args(argv=None):

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fwd", required=True, type=str, help="Path to forward (or single-end) fastq file")
    parser.add_argument("-r", "--rev", required=False, type=str, help="Path to reverse fastq file")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    args = parser.parse_args(argv)
    
    _FWD = args.fwd
    _REV = args.rev
    _SAMPLE = args.sample

    return _FWD, _REV, _SAMPLE


def main(argv=None):

    _FWD, _REV, _SAMPLE = parse_args(argv)

    fwd_handle = gzip.open(_FWD, "rt")
    fwd_reads = SeqIO.to_dict(SeqIO.parse(fwd_handle, "fastq"))
    fwd_handle.close()

    paired_end = True

    if _REV == None:
        paired_end = False
    else:
        rev_handle = gzip.open(_REV, "rt")
        rev_reads = SeqIO.to_dict(SeqIO.parse(rev_handle, "fastq"))
        rev_handle.close()
    
    remove_lst = []

    for read_id in fwd_reads.keys():

        if "N" in str(fwd_reads[read_id].seq):
            print(read_id)
            remove_lst.append(read_id)
            continue
        elif paired_end and "N" in str(rev_reads[read_id].seq):
            print(read_id)
            remove_lst.append(read_id)
            continue

    [ fwd_reads.pop(read_id) for read_id in remove_lst ]
    if paired_end:
        [ rev_reads.pop(read_id) for read_id in remove_lst ]

    if paired_end:
        
        fwd_handle = bgzf.BgzfWriter(f"./{_SAMPLE}_noambig_1.fastq.gz", "wb")
        rev_handle = bgzf.BgzfWriter(f"./{_SAMPLE}_noambig_2.fastq.gz", "wb")
        
        SeqIO.write(sequences=fwd_reads.values(), handle=fwd_handle, format="fastq")
        SeqIO.write(sequences=rev_reads.values(), handle=rev_handle, format="fastq")

        fwd_handle.close()
        rev_handle.close()
    else:
        fwd_handle = bgzf.BgzfWriter(f"./{_SAMPLE}_noambig.fastq.gz", "wb")
        SeqIO.write(sequences=fwd_reads.values(), handle=fwd_handle, format="fastq")
        fwd_handle.close()

if __name__ == "__main__":
    main()