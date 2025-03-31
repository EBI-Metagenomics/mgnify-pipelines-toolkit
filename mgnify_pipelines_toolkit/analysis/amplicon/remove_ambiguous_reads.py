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

# Script removes any reads with ambiguous bases (Ns) for the purpose of DADA2

import argparse
import gzip

from Bio import SeqIO, bgzf


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f",
        "--fwd",
        required=True,
        type=str,
        help="Path to forward (or single-end) fastq file",
    )
    parser.add_argument(
        "-r", "--rev", required=False, type=str, help="Path to reverse fastq file"
    )
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    args = parser.parse_args()

    fwd = args.fwd
    rev = args.rev
    sample = args.sample

    return fwd, rev, sample


def main():

    fwd, rev, sample = parse_args()

    fwd_handle = gzip.open(fwd, "rt")
    fwd_reads = SeqIO.to_dict(SeqIO.parse(fwd_handle, "fastq"))
    fwd_handle.close()

    paired_end = True

    if rev is None:
        paired_end = False
    else:
        rev_handle = gzip.open(rev, "rt")
        rev_reads = SeqIO.to_dict(SeqIO.parse(rev_handle, "fastq"))
        rev_handle.close()

    remove_set = set()

    for read_id in fwd_reads.keys():

        fwd_read_seq = str(fwd_reads[read_id].seq)
        if len(fwd_read_seq) < 100:
            remove_set.add(read_id)
            continue
        elif "N" in fwd_read_seq:
            print(read_id)
            remove_set.add(read_id)
            continue

        if paired_end:
            rev_read_seq = str(rev_reads[read_id].seq)
            if len(rev_read_seq) < 100:
                print(read_id)
                remove_set.add(read_id)
                continue
            elif "N" in rev_read_seq:
                print(read_id)
                remove_set.add(read_id)
                continue

    [fwd_reads.pop(read_id) for read_id in remove_set]
    if paired_end:
        [rev_reads.pop(read_id) for read_id in remove_set]

    if paired_end:
        fwd_handle = bgzf.BgzfWriter(f"./{sample}_noambig_1.fastq.gz", "wb")
        rev_handle = bgzf.BgzfWriter(f"./{sample}_noambig_2.fastq.gz", "wb")

        SeqIO.write(sequences=fwd_reads.values(), handle=fwd_handle, format="fastq")
        SeqIO.write(sequences=rev_reads.values(), handle=rev_handle, format="fastq")

        fwd_handle.close()
        rev_handle.close()
    else:
        fwd_handle = bgzf.BgzfWriter(f"./{sample}_noambig.fastq.gz", "wb")
        SeqIO.write(sequences=fwd_reads.values(), handle=fwd_handle, format="fastq")
        fwd_handle.close()


if __name__ == "__main__":
    main()
