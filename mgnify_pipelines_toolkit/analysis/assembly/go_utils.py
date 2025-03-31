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

from collections import defaultdict
import logging
import os
from pathlib import Path
import re

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s"
)


def count_and_assign_go_annotations(
    go2protein_count: defaultdict[int],
    go_annotations: set[str],
    num_of_proteins: int,
    mapped_go_terms: defaultdict[set] = None,
) -> defaultdict[int]:
    """Increments counts dictionary for GO terms found on a protein.
        If used for GO-slim terms, then a mapped_go_terms dictionary is required
        (with default value of None).
    :param go2protein_count: Current state of the count dictionary
    :type go2protein_count: defaultdict[int]
    :param go_annotations: GO-terms to be incremented
    :type go_annotations: set[str]
    :param num_of_proteins: Number of proteins to be incremented (not sure if we need this, see TODO below)
    :type num_of_proteins: int
    :param mapped_go_terms: Dictionary containin the GO-slim conversion
    :type mapped_go_terms: defaultdict(set)
    :return: _description_
    :rtype: _type_
    """

    if not mapped_go_terms:
        for go_id in go_annotations:
            go2protein_count[go_id] += num_of_proteins
    else:
        slim_go_ids_set = set()
        for go_annotation in go_annotations:
            mapped_go_ids = mapped_go_terms.get(go_annotation)
            if mapped_go_ids:
                slim_go_ids_set.update(mapped_go_ids)
        for slim_go_id in slim_go_ids_set:
            go2protein_count[slim_go_id] += num_of_proteins

    return go2protein_count


def parse_interproscan_tsv(ips_file: Path, mapped_go_terms: dict = None) -> dict:
    """Parses an InterProScan output line by line and return a dictionary of counts for the different GO terms.
        The structure of the IPS file is one annotation per line, some of which will be GO terms. If a protein
        has multiple annotations, then those annotations will follow one by one in order. This function therefore
        parses the file by keeping some flags to track which proteins it's currently on, and which GO terms were found
        for said protein. It then finally increments the count of said protein's GO terms when it's done being parsed.
    :param ips_file: InterProScan .tsv file
    :type ips_file: Path
    :return: Dictionary containing GO term counts in the input InterProScan file
    :rtype: dict
    """

    go2protein_count = defaultdict(int)
    if not os.path.exists(ips_file):
        logging.error(f"The InterProScan file {ips_file} could not be found. Exiting.")
        exit(1)

    num_of_proteins_with_go = 0
    total_num_of_proteins = 0
    line_counter = 0
    previous_protein_acc = None
    go_annotations_single_protein = set()

    go_pattern = re.compile("GO:\\d+")

    with open(ips_file, "r") as fr:

        for line in fr:
            # IPS files are parsed line by line - the same protein accession will appear multiple lines in a row with different annotation
            line_counter += 1
            line = line.strip()
            chunks = line.split("\t")
            # Get protein accession
            current_protein_acc = chunks[0]

            # TODO: not sure if this line is needed - do we ever have more than one protein in a single line of IPS?
            # Will keep just in case
            num_of_proteins = len(current_protein_acc.split("|"))

            # If we're at a new protein accession in the IPS file then we finally increment
            # the go2protein_count dictionary for each term that was found in that protein
            if current_protein_acc != previous_protein_acc:
                total_num_of_proteins += 1
                if len(go_annotations_single_protein) > 0:
                    num_of_proteins_with_go += 1
                    go2protein_count = count_and_assign_go_annotations(
                        go2protein_count,
                        go_annotations_single_protein,
                        num_of_proteins,
                        mapped_go_terms,
                    )
                # reset GO id set because we hit a new protein accession
                go_annotations_single_protein = set()
                previous_protein_acc = current_protein_acc

            # Parse out GO annotations
            # GO annotations are associated to InterPro entries (InterPro entries start with 'IPR')
            # Than use the regex to extract the GO Ids (e.g. GO:0009842)
            if len(chunks) >= 13 and chunks[11].startswith("IPR"):
                for go_annotation in go_pattern.findall(line):
                    go_annotations_single_protein.add(go_annotation)

        # Do final counting for the last protein
        go2protein_count = count_and_assign_go_annotations(
            go2protein_count,
            go_annotations_single_protein,
            num_of_proteins,
            mapped_go_terms,
        )

    return go2protein_count
