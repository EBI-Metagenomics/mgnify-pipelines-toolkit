#!/usr/bin/env python3

# Copyright 2023-2024 EMBL - European Bioinformatics Institute
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

logging.basicConfig(level=logging.INFO)


def main():

    args = parse_args()
    standard_file, substrate_file, outfile, dbcan_ver = (
        args.standard_file,
        args.substrate_file,
        args.outfile,
        args.dbcan_ver,
    )
    standard_path = Path(standard_file)
    substrate_path = Path(substrate_file)

    if not standard_path.exists():
        raise FileNotFoundError(f"Input standards path does not exist: {standard_file}")

    if not substrate_path.exists():
        raise FileNotFoundError(
            f"Input substrate path does not exist: {substrate_file}"
        )

    substrates = load_substrates(substrate_path)
    cgc_locations = load_cgcs(standard_path)
    print_gff(standard_path, outfile, dbcan_ver, substrates, cgc_locations)


def load_cgcs(standard_path):
    cgc_locations = dict()
    with fileinput.hook_compressed(standard_path, "r", encoding="utf-8") as file_in:
        for line in file_in:
            if not line.startswith("CGC#"):
                cgc, _, contig, _, start, end, _, _ = line.strip().split("\t")
                cgc_id = f"{contig}_{cgc}"
                if cgc_id in cgc_locations:
                    if cgc_locations[cgc_id]["start"] > int(start):
                        cgc_locations[cgc_id]["start"] = int(start)
                    if cgc_locations[cgc_id]["end"] < int(end):
                        cgc_locations[cgc_id]["end"] = int(end)
                else:
                    cgc_locations[cgc_id] = {
                        "start": int(start),
                        "end": int(end),
                        "contig": contig,
                    }
    return cgc_locations


def print_gff(standard_path, outfile, dbcan_version, substrates, cgc_locations):
    with open(outfile, "w") as file_out:
        file_out.write("##gff-version 3\n")
        cgcs_printed = list()
        with fileinput.hook_compressed(standard_path, "r", encoding="utf-8") as file_in:
            for line in file_in:
                if not line.startswith("CGC#"):
                    cgc, gene_type, contig, prot_id, start, end, strand, protein_fam = (
                        line.strip().split("\t")
                    )
                    cgc_id = f"{contig}_{cgc}"
                    protein_fam = protein_fam.replace(" ", "")
                    if cgc_id not in cgcs_printed:
                        substrate = (
                            substrates[cgc_id]
                            if cgc_id in substrates
                            else "substrate_dbcan-pul=N/A;substrate_dbcan-sub=N/A"
                        )
                        file_out.write(
                            "{}\tdbCAN:{}\tpredicted PUL\t{}\t{}\t.\t.\t.\tID={};{}\n".format(
                                contig,
                                dbcan_version,
                                cgc_locations[cgc_id]["start"],
                                cgc_locations[cgc_id]["end"],
                                cgc_id,
                                substrate,
                            )
                        )
                        cgcs_printed.append(cgc_id)
                    file_out.write(
                        (
                            f"{contig}\tdbCAN:{dbcan_version}\t{gene_type}\t{start}"
                            + f"\t{end}\t.\t{strand}\t.\tID={prot_id};Parent={cgc_id};protein_family={protein_fam}\n"
                        )
                    )


def load_substrates(substrate_path):
    substrates = dict()
    with fileinput.hook_compressed(substrate_path, "r", encoding="utf-8") as file_in:
        for line in file_in:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                cgc_parts = parts[0].rsplit("|", 1)
                cgc = "_".join(cgc_parts)
                try:
                    substrate_pul = parts[2]
                except IndexError:
                    substrate_pul = "N/A"
                try:
                    substrate_ecami = parts[5]
                except IndexError:
                    substrate_ecami = "N/A"
                if not substrate_pul:
                    substrate_pul = "N/A"
                if not substrate_ecami:
                    substrate_ecami = "N/A"
                substrates[cgc] = (
                    f"substrate_dbcan-pul={substrate_pul};substrate_dbcan-sub={substrate_ecami}"
                )

    return substrates


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script takes dbCAN output and parses it to create a standalone GFF."
        )
    )
    parser.add_argument(
        "-st",
        dest="standard_file",
        required=True,
        help="Path to the standard file (*cgc_standard.out)",
    )
    parser.add_argument(
        "-sb",
        dest="substrate_file",
        required=True,
        help="Path to the substrate file (*substrate.out)",
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
