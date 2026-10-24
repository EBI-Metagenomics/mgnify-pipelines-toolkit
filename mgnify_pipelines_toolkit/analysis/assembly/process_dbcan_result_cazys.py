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
import csv
import fileinput
import logging
import re
from pathlib import Path

logging.basicConfig(level=logging.INFO)


HMM_HIT_EVALUE: float = 1e-15


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
    """Group CDS lines by Parent (mettannotator) or ID (assembly analysis).

    Ignore other feature types and stop before an embedded FASTA section.
    Keep a list because one transcript can have several CDS lines.
    """
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
                logging.error("Not sure what gff annotation delimiter is in use. Exiting")
                exit(1)

            transcript_name = match.group(1)
            genome_gff_lines.setdefault(transcript_name, []).append(line)
    return genome_gff_lines


def print_gff(overview_file: Path, outfile: str, dbcan_version: str, substrates: dict, genome_gff_lines: dict):
    """
    Write recommended CAZyme annotations using the original CDS coordinates.

    Overview columns are matched by name, so added or reordered columns are safe.
    EC evidence counts are removed; empty values, '-' and 'N/A' are omitted.
    Substrates come from the HMM-derived mapping, not the overview's Substrate column.
    Tool counts are retained for downstream filtering by gff_toolkit.

    :param overview_file: Path to the input overview file containing annotations and dbCAN tool information.
    :param outfile: Path to the output GFF3 file to be generated.
    :param dbcan_version: Version of the dbCAN tool used to annotate.
    :param substrates: A dictionary mapping transcript identifiers to sets of substrate names.
    :param genome_gff_lines: A dictionary where keys are transcript identifiers and values are lists of GFF3 lines corresponding to those transcripts.
    """
    with open(outfile, "w") as file_out:
        file_out.write("##gff-version 3\n")
        with fileinput.hook_compressed(overview_file, "r", encoding="utf-8") as file_in:
            # Only these named columns are needed; extra columns are ignored.
            reader = csv.DictReader(file_in, delimiter="\t")
            required_columns = {"Gene ID", "EC#", "#ofTools", "Recommend Results"}
            missing_columns = required_columns - set(reader.fieldnames or [])
            if missing_columns:
                raise ValueError(f"Missing required dbCAN overview columns: {', '.join(sorted(missing_columns))}")
            for row in reader:
                # DictReader uses None for cells absent from a truncated row, unlike an empty cell ('').
                missing_values = [column for column in sorted(required_columns) if row[column] is None]
                if missing_values:
                    raise ValueError(f"Incomplete dbCAN overview row {reader.line_num}: {', '.join(missing_values)}")
                transcript = row["Gene ID"].strip()
                # TODO: why do we have such strict rule?
                if not transcript.startswith("MGYG") and not transcript.startswith("ERZ"):
                    continue

                ec_number_raw = row["EC#"]
                num_of_tools = row["#ofTools"].strip()
                recommended_family = row["Recommend Results"].strip()

                # No recommended family means no CAZyme feature should be written.
                if recommended_family in {"", "-", "N/A"}:
                    continue

                # dbCAN 5.1.2 uses '|'; dbCAN 5.2.9 may mix ';' within groups and '|' between groups.
                # "3.2.1.24:1;3.2.1.113:2" -> "3.2.1.24|3.2.1.113"
                # The suffix ':2' counts proteins with that EC in the subfamily; it is not part of the EC.
                ec_list = []
                for entry in re.split(r"[|;]", ec_number_raw):
                    ec = entry.split(":", 1)[0].strip()
                    # A standalone '-' is missing; a partial EC such as '3.2.1.-' must be kept.
                    if ec not in {"", "-", "N/A"}:
                        ec_list.append(ec)
                # Preserve the original order and repeated ECs for compatibility with dbCAN 5.1.2 output.
                ec_number = "|".join(ec_list)

                # Sorted, deduplicated HMM substrates become a comma-separated GFF value.
                # Genes without a named substrate get N/A, rather than an empty attribute.
                substrates_str = ",".join(sorted(substrates.get(transcript, ["N/A"])))

                col9_parts = [
                    f"protein_family={recommended_family}",
                    f"substrate_dbcan-sub={substrates_str}",
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


def load_substrates(hmm_path: Path) -> dict:
    """
    Collect named substrates from subfamily HMM hits with i-Evalue < 1e-15.

    Read columns by header name. Split substrate lists on ';' or ',', discard
    missing-value tokens, and combine all qualifying hits for each gene into a set.
    Genes with no named substrates are absent, so print_gff can supply N/A.

    :param hmm_path: The path to the HMM output file. The file must be in a tab-delimited format and may be compressed.
    :return: A dictionary mapping gene IDs to sets of substrate names.
    """
    substrates = dict()
    with fileinput.hook_compressed(hmm_path, "r", encoding="utf-8") as file_in:
        header = next(file_in)
        header_fields = header.strip().split("\t")
        substrate_idx = header_fields.index("Substrate")
        gene_idx = header_fields.index("Target Name")
        evalue_idx = header_fields.index("i-Evalue")
        for line in file_in:
            # Preserve empty TSV cells; stripping tabs would shift the remaining columns.
            fields = line.rstrip("\r\n").split("\t")
            if float(fields[evalue_idx]) < HMM_HIT_EVALUE:
                gene_id = fields[gene_idx]
                # "host glycan;-;N/A,host glycan" -> {"host glycan"}
                # Filter each token, since '-' can appear alongside a named substrate.
                for entry in re.split(r"[;,]", fields[substrate_idx]):
                    substrate = entry.strip()
                    if substrate not in {"", "-", "N/A"}:
                        # Sets also remove duplicates across multiple hits for the same gene.
                        substrates.setdefault(gene_id, set()).add(substrate)
    return substrates


def parse_args():
    parser = argparse.ArgumentParser(description="The script takes dbCAN outputs and parses them to create a standalone GFF.")
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
