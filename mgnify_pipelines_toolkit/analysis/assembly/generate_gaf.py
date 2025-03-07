#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024-2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import os
import logging

from mgnify_pipelines_toolkit.analysis.assembly.go_utils import parse_interproscan_tsv


def parse_args():

    description = "Go slim pipeline for processing InterProScan results"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-i", "--ips_input", help="InterProScan result file", required=True
    )
    parser.add_argument("-o", "--output", help="GO summary output file", required=True)
    args = parser.parse_args()

    ips_input = args.ips_input
    output = args.output

    return ips_input, output


# Constants
PROJECT_NAME = "EBI Metagenomics"
PROJECT_URL = "http://www.ebi.ac.uk/metagenomics"
PROJECT_CONTACT = "metagenomics-help@ebi.ac.uk"
FIXED_TIMESTAMP = "20160528"  # What is this timestamp?

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s"
)


def write_gaf_file(gaf_input_file_path: str, go_id_set: set[str]) -> None:
    """
    Create a GO Annotation File (GAF) from a set of GO IDs.

    :param gaf_input_file_path: Path to output GAF file
    :param go_id_set: Set of GO IDs to include in the file
    """
    with open(gaf_input_file_path, "w") as fw:
        # Write GAF header
        fw.write("!gaf-version: 2.1\n")
        fw.write(f"!Project_name: {PROJECT_NAME}\n")
        fw.write(f"!URL: {PROJECT_URL}\n")
        fw.write(f"!Contact Email: {PROJECT_CONTACT}\n")

        # Write GO entries
        for go_id in go_id_set:
            gaf_entry = "\t".join(
                [
                    "EMG",
                    go_id,
                    "GO",
                    "",
                    go_id,
                    "PMID:12069591",
                    "IEA",
                    "",
                    "P",
                    "",
                    "",
                    "protein",
                    "taxon:1310605",
                    FIXED_TIMESTAMP,
                    "InterPro",
                    "",
                ]
            )
            fw.write(gaf_entry + "\n")

    logging.info(f"GAF file created successfully: {gaf_input_file_path}")


def main():
    """
    Process the InterProScan TSV output and generate a GO annotation file (GAF)).
    """

    ips_input, output = parse_args()

    # Validate input file
    if not os.path.exists(ips_input):
        raise FileNotFoundError(f"Input file not found: {ips_input}")

    if os.path.getsize(ips_input) == 0:
        logging.warning("Input file is empty. Skipping processing.")
        return

    # Parse InterProScan result file
    logging.info(f"Parsing InterProScan input: {ips_input}")
    go2protein_count_dict = parse_interproscan_tsv(ips_input)
    logging.info("Finished parsing InterProScan file")

    logging.info("Writing the GAF file")
    write_gaf_file(f"{output}_ips_annotations.gaf", go2protein_count_dict.keys())


if __name__ == "__main__":
    main()
