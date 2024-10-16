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
from pathlib import Path

from Bio import SeqIO
import pandas as pd
import requests


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)


def main(input: Path, proteins: Path, output: Path, rhea2chebi: Path, up2rhea: Path):
    logging.info("Step 0/5: Checking Rhea-CHEBI mapping file...")
    if not rhea2chebi:
        logging.info("Rhea-CHEBI mapping is not provided. Starting download...")
        download_path = Path(".")
        rhea2chebi = download_and_convert_to_tsv(download_path)
        logging.info(
            f"File downloaded successfully and saved to {rhea2chebi.resolve()}"
        )

    logging.info(f"Step 1/5: Reading input file {input.resolve()}")
    # TODO File is quite big, inadequate memory load
    diamond_df = pd.read_csv(
        input,
        sep="\t",
        usecols=[0, 1, 10],
        names=["protein_id", "uniref90_ID", "e_value"],
        header=None,
    )
    # UniRef90 cluster name contains a representative protein name
    diamond_df["uniref90_rep"] = diamond_df["uniref90_ID"].str.split("_").str[1]

    logging.info(
        f"Step 2/5: Adding RHEA IDs based on provided file {up2rhea.resolve()}"
    )
    # TODO File is quite big, inadequate memory load
    up2rhea_df = pd.read_csv(up2rhea, sep="\t", header=None)
    up2rhea_df.columns = ["unirefKB_id", "rhea_id"]
    diamond_df = diamond_df.merge(
        up2rhea_df, left_on="uniref90_rep", right_on="unirefKB_id", how="left"
    )
    diamond_df = diamond_df[diamond_df["rhea_id"].notna()]

    split_protein_id = diamond_df["protein_id"].str.split("_", n=1, expand=True)
    diamond_df["contig_id"] = split_protein_id[0]

    diamond_df["rank"] = diamond_df.groupby("protein_id")["e_value"].rank(
        method="first"
    )

    diamond_df["rhea_id"] = diamond_df["rhea_id"].str.split()
    diamond_df = diamond_df.explode("rhea_id")

    logging.info(
        f"Step 3/5: Adding CHEBI reactions based on provided file {rhea2chebi.resolve()}"
    )
    rhea2chebi_df = pd.read_csv(rhea2chebi, sep="\t")
    rhea2chebi_df.columns = [
        "rhea_id",
        "reaction_definition",
        "chebi_reaction",
        "enzyme_id",
    ]
    diamond_df = diamond_df.merge(
        rhea2chebi_df[["rhea_id", "chebi_reaction", "reaction_definition"]],
        on="rhea_id",
        how="left",
    )

    grouped_df = (
        diamond_df.groupby(
            ["protein_id", "rhea_id", "reaction_definition", "chebi_reaction"]
        )
        .agg(
            {
                "uniref90_ID": lambda x: ",".join(x),  # Combine all UniRef90 IDs
                "rank": "min",
            }
        )
        .reset_index()
    )

    # Assign top hit based on first match with RHEA IDs
    grouped_df["top_hit"] = grouped_df["rank"].apply(
        lambda x: "top hit" if x == 1 else ""
    )

    # Merge contig ID back into the grouped dataframe
    grouped_df = grouped_df.merge(
        diamond_df[["protein_id", "contig_id"]].drop_duplicates(),
        on="protein_id",
        how="left",
    )

    logging.info(
        f"Step 4/5: Parsing protein fasta and calculating SHA256 hash from {proteins.resolve()}"
    )
    protein_hashes = {}
    with open(proteins, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_hash = hashlib.sha256(str(record.seq).encode("utf-8")).hexdigest()
            protein_hashes[record.id] = protein_hash

    grouped_df["checksum"] = grouped_df["protein_id"].map(protein_hashes)

    logging.info(f"Step 5/5: Saving output table to {output.resolve()}")
    grouped_df = grouped_df[
        [
            "contig_id",
            "protein_id",
            "checksum",
            # "uniref90_ID",
            "rhea_id",
            "reaction_definition",
            "chebi_reaction",
            "top_hit",
        ]
    ]
    grouped_df.to_csv(output, sep="\t", index=False)

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
        "-i", "--input", required=True, type=Path, help="DIAMOND results file"
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
    main(args.input, args.proteins, args.output, args.rhea2chebi, args.up2rhea)
