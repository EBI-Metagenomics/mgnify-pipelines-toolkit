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
import gzip
import hashlib
import logging
import re
import sys
from pathlib import Path

from Bio import SeqIO
import pandas as pd
import requests


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)


def main(proteins: Path, output: Path, rhea2chebi: Path):
    logging.info("Step 0/3: Check Rhea-CHEBI mapping file...")
    if not rhea2chebi:
        logging.info("Rhea-CHEBI mapping is not provided. Starting download...")
        download_path = Path(".")
        rhea2chebi = download_and_convert_to_tsv(download_path)
        logging.info(
            f"File downloaded successfully and saved to {rhea2chebi.resolve()}"
        )

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

    logging.info("Step 3/3: Read DIAMOND results from STDIN and write output")
    with open(output, "w") as fh:
        current_protein = None
        for line in sys.stdin:
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
                        file=fh,
                    )
                    protein_rheas.add(rhea)

    logging.info("Processed successfully. Exiting.")


def download_and_convert_to_tsv(download_path: Path) -> Path:
    url = "https://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz"

    # Ensure the download directory exists
    download_path.mkdir(parents=True, exist_ok=True)

    # Define paths for the created files
    gz_file_path = download_path / "rhea-reactions.txt.gz"
    tsv_file_path = download_path / "rhea2chebi.tsv"

    # Download the .gz file
    response = requests.get(url)
    response.raise_for_status()
    with gz_file_path.open("wb") as file:
        file.write(response.content)

    # Process the .gz file and convert it to TSV
    with gzip.open(gz_file_path, "rt") as file_in, open(tsv_file_path, "w") as file_out:
        input_text = file_in.read()
        pattern = re.compile(
            r"ENTRY\s+(.*?)\nDEFINITION\s+(.*?)\nEQUATION\s+(.*?)\n(?:ENZYME\s+(.*?))?\n*",
            re.DOTALL,
        )
        matches = pattern.findall(input_text)
        file_out.write("ENTRY\tDEFINITION\tEQUATION\tENZYME\n")
        for match in matches:
            entry, definition, equation, enzyme = match
            enzyme = enzyme.strip() if enzyme else ""
            file_out.write(f"{entry}\t{definition}\t{equation}\t{enzyme}\n")

    # Clean up the downloaded .gz file
    gz_file_path.unlink()
    return tsv_file_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        """
                                     Use diamond output file to create a table with Rhea and CHEBI
                                     reaction annotation for every protein.
                                     """
    )
    parser.add_argument(
        "-p",
        "--proteins",
        required=True,
        type=Path,
        help="Protein fasta file used as DIAMOND input",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output TSV file with columns: contig_id, protein_id, UniRef90 cluster, rhea_ids, CHEBI reaction participants",
    )
    parser.add_argument(
        "--rhea2chebi",
        default=None,
        type=Path,
        help="File that maps rhea_ids to CHEBI. Will be downloaded if not provided",
    )

    args = parser.parse_args()
    main(args.proteins, args.output, args.rhea2chebi)
