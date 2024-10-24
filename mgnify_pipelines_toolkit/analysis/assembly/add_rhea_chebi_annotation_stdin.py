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
from collections import defaultdict
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


def main(proteins: Path, output: Path, rhea2chebi: Path, up2rhea: Path):
    logging.info("Step 0/5: Checking Rhea-CHEBI mapping file...")
    if not rhea2chebi:
        logging.info("Rhea-CHEBI mapping is not provided. Starting download...")
        download_path = Path(".")
        rhea2chebi = download_and_convert_to_tsv(download_path)
        logging.info(
            f"File downloaded successfully and saved to {rhea2chebi.resolve()}"
        )

    logging.info(f"Step 1/5: Loading RHEA ids from provided file {up2rhea.resolve()}")
    df = pd.read_csv(up2rhea, delimiter="\t", header=None)
    up2rhea_dict = dict(zip(df.iloc[:, 0], df.iloc[:, 1].str.split()))

    logging.info(
        f"Step 2/5: Loading reactions from provided file {rhea2chebi.resolve()}"
    )
    df = pd.read_csv(rhea2chebi, delimiter="\t")
    rhea2reaction_dict = dict(zip(df["ENTRY"], zip(df["EQUATION"], df["DEFINITION"])))

    logging.info("Step 3/5: Reading DIAMOND results from STDIN")
    query2rhea = defaultdict(dict)
    for line in sys.stdin:
        parts = line.strip().split("\t")
        protein_id, uniref90_ID = parts[0], parts[1].split("_")[1]
        rhea_list = up2rhea_dict.get(uniref90_ID, [])
        top_hit = "top hit" if rhea_list and protein_id not in query2rhea else ""

        for rhea in rhea_list:
            chebi_reaction, reaction = rhea2reaction_dict[rhea]
            query2rhea[protein_id][rhea] = (chebi_reaction, reaction, top_hit)

    logging.info(
        f"Step 4/5: Parsing protein fasta and calculating SHA256 hash from {proteins.resolve()}"
    )
    protein_hashes = {}
    with open(proteins, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_hash = hashlib.sha256(str(record.seq).encode("utf-8")).hexdigest()
            protein_hashes[record.id] = protein_hash

    logging.info(f"Step 5/5: Finalising and saving output table to {output.resolve()}")
    with open(output, "w") as fh:
        for protein_id in query2rhea:
            contig_id = protein_id.split("_")[0]
            protein_hash = protein_hashes.get(protein_id, "N/A")
            for rhea in query2rhea[protein_id]:
                print(
                    contig_id,
                    protein_id,
                    protein_hash,
                    rhea,
                    "\t".join(query2rhea[protein_id][rhea]),
                    sep="\t",
                    file=fh,
                )

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
    parser.add_argument(
        "--up2rhea",
        required=True,
        type=Path,
        help="File that maps UniProt IDs to Rhea. Must contain at least 2 columns 'Entry' and 'Rhea ID'",
    )
    args = parser.parse_args()
    main(args.proteins, args.output, args.rhea2chebi, args.up2rhea)
