#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2023-2025 EMBL - European Bioinformatics Institute
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
#

import argparse
import fileinput
import logging
from pathlib import Path
import re

logging.basicConfig(level=logging.INFO)


def main():

    args = parse_args()
    hmm_file, overview_file, genome_gff, outfile, dbcan_ver = (
        args.hmm_file,
        args.overview_file,
        args.genome_gff,
        args.outfile,
        args.dbcan_ver,
    )

    hmm_path = Path(hmm_file)
    overview_path = Path(overview_file)

    if not hmm_path.is_file():
        raise FileNotFoundError(f"Input hmm path does not exist: {hmm_file}")

    if not overview_path.is_file():
        raise FileNotFoundError(f"Input overview path does not exist: {overview_file}")

    substrates = load_substrates(hmm_path)
    genome_gff_lines = load_gff(genome_gff)

    print_gff(overview_file, outfile, dbcan_ver, substrates, genome_gff_lines)


def load_gff(gff):
    genome_gff_lines = dict()
    with fileinput.hook_compressed(gff, "r", encoding="utf-8") as gff:
        for line in gff:
            if line.startswith("##FASTA"):
                return genome_gff_lines

            fields = line.strip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue

            if "Parent=" in line:
                # Get transcript name from the 9th column for mettannotator
                match = re.search(r"Parent=([^;]+)", fields[8])
            elif "ID=" in line:
                # Get transcript name from the 9th column for ASA
                match = re.search(r"ID=([^;]+)", fields[8])
            else:
                logging.error(
                    "Not sure what gff annotation delimiter is in use. Exiting"
                )
                exit(1)

            transcript_name = match.group(1)
            genome_gff_lines.setdefault(transcript_name, []).append(line)
    return genome_gff_lines


def print_gff(overview_file, outfile, dbcan_version, substrates, genome_gff_lines):
    with open(outfile, "w") as file_out:
        file_out.write("##gff-version 3\n")
        with fileinput.hook_compressed(overview_file, "r", encoding="utf-8") as file_in:
            for line in file_in:
                if line.startswith("MGYG") or line.startswith("ERZ"):
                    (
                        transcript,
                        ec_number_raw,
                        dbcan_hmmer,
                        dbcan_sub_ecami,
                        diamond,
                        num_of_tools,
                    ) = line.strip().split("\t")
                    # EC is reported as 2.4.99.-:5 with :5 meaning 5 proteins in the subfamily have EC 2.4.99.-

                    ec_number = ""
                    ec_list = ec_number_raw.split("|")
                    for ec in ec_list:
                        if ec != "-":
                            ec_number += ec.split(":")[0] + "|"

                    ec_number = ec_number.strip("|")

                    # Dbcan recommends to use subfamily preference as dbcan_hmmer > dbcan_sub_ecami > diamond
                    # diamond is messier, so we don't report it here
                    if dbcan_hmmer != "-":
                        # the field dbcan_hmmer reports match positions in parentheses, clear them out first:
                        subfamily = dbcan_hmmer.split("(")[0]
                    elif dbcan_sub_ecami != "-":
                        subfamily = dbcan_sub_ecami
                    else:
                        continue
                    cleaned_substrates = ",".join(
                        sorted(
                            {
                                subsrate.strip()
                                for subsrate in substrates.get(transcript, "N/A").split(
                                    ","
                                )
                            }
                        )
                    )
                    # Assemble information to add to the 9th column
                    col9_parts = [
                        f"protein_family={subfamily}",
                        f"substrate_dbcan-sub={cleaned_substrates}",
                    ]

                    if ec_number:
                        col9_parts.append(f"eC_number={ec_number}")

                    col9_parts.append(f"num_tools={num_of_tools}")
                    col9_text = ";".join(col9_parts)

                    for gff_line in genome_gff_lines[transcript]:
                        fields = gff_line.strip().split("\t")
                        # Replace the tool
                        fields[1] = f"dbCAN:{dbcan_version}"
                        # Replace the feature
                        fields[2] = "CAZyme"
                        # Replace the confidence value
                        fields[5] = "."
                        # Keep only the ID in the 9th column
                        attributes = fields[8].split(";")[0]
                        # Add dbcan information to the 9th column
                        attributes = f"{attributes};{col9_text};"
                        fields[8] = attributes
                        file_out.write("\t".join(fields) + "\n")


def load_substrates(hmm_path):
    substrates = dict()
    with fileinput.hook_compressed(hmm_path, "r", encoding="utf-8") as file_in:
        header = next(file_in)
        header_fields = header.strip().split("\t")
        substrate_idx = header_fields.index("Substrate")
        gene_idx = header_fields.index("Gene ID")
        evalue_idx = header_fields.index("E Value")
        for line in file_in:
            fields = line.strip().split("\t")
            if float(fields[evalue_idx]) < 1e-15:  # evalue is the default from dbcan
                substrate = fields[substrate_idx]
                if not substrate == "-":
                    gene_id = fields[gene_idx]
                    substrates.setdefault(gene_id, []).append(substrate)
    # resolve cases with multiple substrates
    for gene_id, substrate_list in substrates.items():
        substrate_list = list(set(substrate_list))
        if len(substrate_list) == 1:
            substrates[gene_id] = substrate_list[0]
        else:
            substrates[gene_id] = ",".join(substrate_list)
    return substrates


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script takes dbCAN output for a eukaryotic genome and parses it to create a standalone GFF."
        )
    )
    parser.add_argument(
        "-hmm",
        dest="hmm_file",
        required=True,
        help="Path to the hmm file.",
    )
    parser.add_argument(
        "-ov",
        dest="overview_file",
        required=True,
        help="Path to the overview file.",
    )
    parser.add_argument(
        "-g",
        dest="genome_gff",
        required=True,
        help="Path to the genome GFF.",
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=True,
        help="Path to the output file.",
    )
    parser.add_argument(
        "-v",
        dest="dbcan_ver",
        required=True,
        help="dbCAN version used.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
