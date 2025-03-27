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
import hashlib
import logging
import sys
from pathlib import Path

from Bio import SeqIO
import pandas as pd


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)


def process_lines(lines, output_handler, rhea2reaction_dict, protein_hashes):
    current_protein = None
    for line in lines:
        parts = line.strip().split("\t")
        protein_id = parts[0]
        if protein_id != current_protein:
            current_protein = protein_id
            protein_rheas = set()
        rhea_list = parts[-1].split("RheaID=")[1].split()
        top_hit = "top hit" if rhea_list and not protein_rheas else ""

        for rhea in rhea_list:
            if rhea not in protein_rheas:
                chebi_reaction, reaction = rhea2reaction_dict[rhea]
                contig_id = protein_id.split("_")[0]
                protein_hash = protein_hashes[protein_id]

                print(
                    contig_id,
                    protein_id,
                    protein_hash,
                    rhea,
                    chebi_reaction,
                    reaction,
                    top_hit,
                    sep="\t",
                    file=output_handler,
                )
                protein_rheas.add(rhea)


def main():
    parser = argparse.ArgumentParser(
        "Use diamond output file to create a table with Rhea and CHEBI reaction annotation for every protein."
    )
    parser.add_argument(
        "-d",
        "--diamond_hits",
        required=True,
        type=str,
        help="DIAMOND results file, use '-' for stdin",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help=(
            "Output TSV file with columns: contig_id, protein_id, protein hash, "
            "Rhea IDs, CHEBI reaction, reaction definition, 'top hit' if it is "
            "the first hit for the protein"
        ),
    )
    parser.add_argument(
        "-p",
        "--proteins",
        required=True,
        type=Path,
        help="Protein fasta file used as DIAMOND input",
    )
    parser.add_argument(
        "--rhea2chebi",
        required=True,
        type=Path,
        help="File that maps rhea_ids to CHEBI",
    )

    args = parser.parse_args()

    diamond_hits = args.diamond_hits
    output = args.output
    proteins = args.proteins
    rhea2chebi = args.rhea2chebi

    logging.info(
        f"Step 1/3: Parse protein fasta and calculating SHA256 hash from {proteins.resolve()}"
    )
    protein_hashes = {}
    with open(proteins, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_hash = hashlib.sha256(str(record.seq).encode("utf-8")).hexdigest()
            protein_hashes[record.id] = protein_hash

    logging.info(f"Step 2/3: Load reactions from provided file {rhea2chebi.resolve()}")
    df = pd.read_csv(rhea2chebi, delimiter="\t")
    rhea2reaction_dict = dict(zip(df["ENTRY"], zip(df["EQUATION"], df["DEFINITION"])))

    logging.info(
        f"Step 3/3: Read DIAMOND results from {'STDIN' if diamond_hits == '-' else Path(diamond_hits).resolve()} and write output"
    )
    with open(output, "w") as output_handler:
        if diamond_hits == "-":
            process_lines(sys.stdin, output_handler, rhea2reaction_dict, protein_hashes)
        else:
            with open(diamond_hits, "r") as input_file:
                process_lines(
                    input_file, output_handler, rhea2reaction_dict, protein_hashes
                )

    logging.info("Processed successfully. Exiting.")


if __name__ == "__main__":
    main()
