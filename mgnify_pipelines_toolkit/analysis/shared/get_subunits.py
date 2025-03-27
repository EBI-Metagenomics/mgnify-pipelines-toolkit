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
from Bio import SeqIO
from mgnify_pipelines_toolkit.constants.ncrna import (
    DIRECTORY_SEQ_CAT,
    SSU,
    LSU,
    Seq5S,
    Seq5_8S,
    SSU_rRNA_archaea,
    SSU_rRNA_bacteria,
    SSU_rRNA_eukarya,
    SSU_rRNA_microsporidia,
    LSU_rRNA_archaea,
    LSU_rRNA_bacteria,
    LSU_rRNA_eukarya,
    NON_CODING_RNA,
    SSU_MODELS,
    LSU_MODELS,
    RFAM_MODELS,
)


def set_model_names(prefix, name, directory, separate_subunits):
    pattern_dict = {}
    pattern_dict[SSU] = os.path.join(directory, f"{name}_SSU.fasta")
    pattern_dict[LSU] = os.path.join(directory, f"{name}_LSU.fasta")
    pattern_dict[Seq5S] = os.path.join(directory, f"{name}_5S.fasta")
    pattern_dict[Seq5_8S] = os.path.join(directory, f"{name}_5_8S.fasta")
    if separate_subunits:
        pattern_dict[SSU_rRNA_archaea] = os.path.join(
            directory,
            f"{prefix}{name}_{SSU_rRNA_archaea}.{RFAM_MODELS[SSU_rRNA_archaea]}.fasta",
        )
        pattern_dict[SSU_rRNA_bacteria] = os.path.join(
            directory,
            f"{prefix}{name}_{SSU_rRNA_bacteria}.{RFAM_MODELS[SSU_rRNA_bacteria]}.fasta",
        )
        pattern_dict[SSU_rRNA_eukarya] = os.path.join(
            directory,
            f"{prefix}{name}_{SSU_rRNA_eukarya}.{RFAM_MODELS[SSU_rRNA_eukarya]}.fasta",
        )
        pattern_dict[SSU_rRNA_microsporidia] = os.path.join(
            directory,
            f"{prefix}{name}_{SSU_rRNA_microsporidia}.{RFAM_MODELS[SSU_rRNA_microsporidia]}.fasta",
        )
        pattern_dict[LSU_rRNA_archaea] = os.path.join(
            directory,
            f"{prefix}{name}_{LSU_rRNA_archaea}.{RFAM_MODELS[LSU_rRNA_archaea]}.fasta",
        )
        pattern_dict[LSU_rRNA_bacteria] = os.path.join(
            directory,
            f"{prefix}{name}_{LSU_rRNA_bacteria}.{RFAM_MODELS[LSU_rRNA_bacteria]}.fasta",
        )
        pattern_dict[LSU_rRNA_eukarya] = os.path.join(
            directory,
            f"{prefix}{name}_{LSU_rRNA_eukarya}.{RFAM_MODELS[LSU_rRNA_eukarya]}.fasta",
        )
    return pattern_dict


def main():
    parser = argparse.ArgumentParser(
        description="Extract lsu, ssu and 5s and other models"
    )
    parser.add_argument(
        "-i", "--input", dest="input", help="Input fasta file", required=True
    )
    parser.add_argument(
        "-p", "--prefix", dest="prefix", help="prefix for models", required=False
    )
    parser.add_argument("-n", "--name", dest="name", help="Accession", required=True)
    parser.add_argument(
        "--separate-subunits-by-models",
        action="store_true",
        help="Create separate files for each kingdon example: sample_SSU_rRNA_eukarya.RF01960.fasta",
    )

    args = parser.parse_args()
    prefix = args.prefix if args.prefix else ""
    name = args.name if args.name else "accession"

    directory = DIRECTORY_SEQ_CAT
    if not os.path.exists(directory):
        os.makedirs(directory)

    print("Start fasta mode")
    pattern_dict = set_model_names(
        prefix, name, directory, args.separate_subunits_by_models
    )

    open_files = {}
    for record in SeqIO.parse(args.input, "fasta"):
        model = "-".join("/".join(record.id.split("/")[:-1]).split("-")[-1:])
        if model in SSU_MODELS:
            if SSU not in open_files:
                file_out = open(pattern_dict[SSU], "w")
                open_files[SSU] = file_out
            SeqIO.write(record, open_files[SSU], "fasta")
        elif model in LSU_MODELS:
            if LSU not in open_files:
                file_out = open(pattern_dict[LSU], "w")
                open_files[LSU] = file_out
            SeqIO.write(record, open_files[LSU], "fasta")

        if model in NON_CODING_RNA:
            if model in pattern_dict:
                filename = pattern_dict[model]
            else:
                filename = None
        else:
            filename = os.path.join(directory, f"{name}_other_ncRNA.fasta")
        if filename:
            if model not in open_files:
                file_out = open(filename, "w")
                open_files[model] = file_out
            SeqIO.write(record, open_files[model], "fasta")

    for item in open_files:
        open_files[item].close()


if __name__ == "__main__":
    main()
