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

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def get_tblout_column_indices(tool):
    cmsearch_indices = {
        "contig": 0,
        "gene": 2,
        "model": 3,
        "strand": 9,
        "start": 7,
        "end": 8,
    }
    cmscan_indices = {
        "contig": 3,
        "gene": 1,
        "model": 2,
        "strand": 11,
        "start": 9,
        "end": 10,
    }
    if tool == "cmsearch":
        return cmsearch_indices
    else:
        return cmscan_indices


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract rRNAs from the Infernal cmsearch/cmscan dereplicated tblout output.\n"
            "Models:\n"
            "RF00001 -> 5S ribosomal RNA (5S rRNA)\n"
            "RF00177 -> This family is a member of clan (CL00111), which contains the following 5 members: "
            "SSU_rRNA_archaea, SSU_rRNA_bacteria, SSU_rRNA_eukarya, "
            "SSU_rRNA_microsporidia, SSU_trypano_mito\n"
            "RF02541 -> This family is a member of clan (CL00112), which contains the following 5 members: "
            "5_8S_rRNA, LSU_rRNA_archaea, LSU_rRNA_bacteria, LSU_rRNA_eukarya, LSU_trypano_mito"
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="Input fasta, the name of the file will be used to prefix the entries in the extracted fasta file",
        required=True,
    )
    parser.add_argument(
        "-d", "--deoverlapped", dest="deoverlapped", help="tblout.deoverlapped", required=True
    )
    parser.add_argument(
        "-s",
        "--source",
        dest="source",
        help="Software used to generate the tblout file",
        choices=["cmsearch", "cmscan"],
        required=False,
        default="cmsearch",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        help="Path to file where the output FASTA will be saved to",
        required=False,
    )
    args = parser.parse_args()
    hits = {}
    added = {}
    idx = get_tblout_column_indices(args.source)
    with open(args.deoverlapped, "r") as f:
        for line in f:
            line = line.strip("\n")
            if line.startswith("#"):
                continue
            cols = line.split()
            model = cols[idx["model"]]

            # RF00001 -> 5S ribosomal RNA (5S rRNA)
            # RF00177 -> This family is a member of clan (CL00111),
            #   It contains the following 5 members:
            #   - SSU_rRNA_archaea
            #   - SSU_rRNA_bacteria
            #   - SSU_rRNA_eukarya
            #   - SSU_rRNA_microsporidia
            #   - SSU_trypano_mito
            # RF02541 -> This family is a member of clan (CL00112).
            #   It contains the following 5 members:
            #   - 5_8S_rRNA LSU_rRNA_archaea
            #   - LSU_rRNA_bacteria
            #   - LSU_rRNA_eukarya
            #   - LSU_trypano_mito
            if model not in ["RF00001", "RF00177", "RF02541"]:
                continue

            contig = cols[idx["contig"]]
            gene = cols[idx["gene"]]
            strand = cols[idx["strand"]]
            if strand == "+":
                start = int(cols[idx["start"]])
                end = int(cols[idx["end"]])
            else:
                start = int(cols[idx["end"]])
                end = int(cols[idx["start"]])
            if contig not in added.keys():
                added[contig] = 1
            else:
                added[contig] += 1
            contig = f"{contig}__{gene}_hit-{added[contig]}__{start}-{end}_len={end - start + 1}"
            hits[contig] = [start, end]

        records_to_write = []  # List to hold SeqRecord objects
        with open(args.input, "r") as input_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                for contig in hits.keys():
                    if contig.split("__")[0] == record.id:
                        start = hits[contig][0] - 1
                        end = hits[contig][1]
                        seq = record.seq[start:end]
                        name = (
                            os.path.basename(args.input).split(".")[0] + "__" + contig
                        )
                        records_to_write.append(SeqRecord(seq, id=name, description=""))

        with open(args.outfile or sys.stdout, "w") as file_out:
            SeqIO.write(records_to_write, file_out, "fasta")


if __name__ == "__main__":
    main()
