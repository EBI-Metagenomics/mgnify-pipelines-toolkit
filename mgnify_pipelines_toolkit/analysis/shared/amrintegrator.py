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
import fileinput
import logging
import os
from collections import defaultdict
from typing import Dict, List, Optional, Union, DefaultDict

# Constant to simplify long strings
DRUG_REPLACEMENT_STRINGS: Dict[str, str] = {
    "macrolide antibiotic; lincosamide antibiotic; streptogramin antibiotic": "MLS",
    "lincosamide/macrolide/streptogramin": "MLS",
}

# Set up logger
logger = logging.getLogger(__name__)


def validate_inputs(optional_inputs: Dict[str, Optional[str]]) -> List[str]:
    """
    Validate that input files exist and return a list of valid input names.

    Args:
        optional_inputs: Dictionary mapping input names to file paths (can be None)

    Returns:
        List of valid input names whose files exist
    """
    valid_inputs: List[str] = []
    for name, path in optional_inputs.items():
        if path and os.path.exists(path):
            valid_inputs.append(name)
        elif path:
            logger.warning(f"File not found for '{name}' → {path}")
    return valid_inputs


def normalize_drug_class(drug_class: str) -> str:
    """
    Normalize drug class names by applying replacements and removing suffixes.

    Args:
        drug_class: Raw drug class string

    Returns:
        Normalized drug class string
    """
    # Apply known full-string replacements
    for original, replacement in DRUG_REPLACEMENT_STRINGS.items():
        if original in drug_class:
            drug_class = drug_class.replace(original, replacement)

    # Remove generic suffixes after replacement
    drug_class = drug_class.replace(" antibiotic", "")

    return drug_class


def parse_hamronized(
    amr_annotation: DefaultDict[str, Dict[str, Dict[str, Union[List[str], str]]]], input_file: str
) -> DefaultDict[str, Dict[str, Dict[str, Union[List[str], str]]]]:
    """
    Parse hamronized AMR annotation files (deeparg, rgi).

    Args:
        amr_annotation: Nested defaultdict to store AMR annotations
        input_file: Path to the hamronized input file

    Returns:
        Updated amr_annotation defaultdict
    """
    with fileinput.hook_compressed(input_file, "r", encoding="utf-8") as input_table:
        # Skipping the header line
        next(input_table)
        for line in input_table:
            line_l: List[str] = line.rstrip().split("\t")
            analysis_software_name: str = line_l[6]
            drug_class: str = normalize_drug_class(line_l[13])
            protein_id: str = line_l[20].split(" ")[0]
            seq_identity: str = line_l[35]
            drug_class_list: List[str] = [cls.strip() for cls in drug_class.split(";")]

            amr_annotation[protein_id][analysis_software_name] = {"drug_class": drug_class_list, "seq_identity": seq_identity}

    return amr_annotation


def parse_amrfinderplus(
    amr_annotation: DefaultDict[str, Dict[str, Dict[str, Union[List[str], str]]]], input_file: str
) -> DefaultDict[str, Dict[str, Dict[str, Union[List[str], str]]]]:
    """
    Parse AMRFinderPlus output files.

    Args:
        amr_annotation: Nested defaultdict to store AMR annotations
        input_file: Path to the AMRFinderPlus output file

    Returns:
        Updated amr_annotation defaultdict
    """
    with fileinput.hook_compressed(input_file, "r", encoding="utf-8") as input_table:
        # Skipping the header line
        next(input_table)
        for line in input_table:
            line_l: List[str] = line.rstrip().split("\t")
            analysis_software_name: str = "amrfinderplus"
            drug_class: str = normalize_drug_class(line_l[6].lower())
            drug_class_list: List[str] = [cls.strip() for cls in drug_class.split(";")]
            protein_id: str = line_l[0]
            seq_identity: str = line_l[12]

            amr_annotation[protein_id][analysis_software_name] = {"drug_class": drug_class_list, "seq_identity": seq_identity}

    return amr_annotation


def parse_amr_dict(amr_annotation: DefaultDict[str, Dict[str, Dict[str, Union[List[str], str]]]]) -> Dict[str, List[str]]:
    """
    Parse AMR annotation dictionary and format protein attributes.

    Args:
        amr_annotation: Nested defaultdict containing AMR annotations

    Returns:
        Dictionary mapping protein IDs to lists of formatted attribute strings
    """
    protein_attributes: Dict[str, List[str]] = {}

    for protein, tools in amr_annotation.items():
        all_drugs: List[str] = []
        tool_names: List[str] = []
        tool_identities: List[str] = []

        for tool, info in tools.items():
            # Collect drug classes
            drug_classes: Union[List[str], str] = info.get("drug_class", [])
            if isinstance(drug_classes, list):
                all_drugs.extend(drug_classes)
            else:
                all_drugs.append(drug_classes)

            # Collect tool name and identity value
            tool_names.append(tool)
            seq_identity: Union[List[str], str] = info.get("seq_identity", "NA")
            if isinstance(seq_identity, str):
                tool_identities.append(seq_identity)
            else:
                tool_identities.append("NA")

        # Clean up and deduplicate drug classes
        unique_drugs: List[str] = list(dict.fromkeys(all_drugs))
        unique_drugs = [item.replace(" ", "_") for item in unique_drugs]

        # Build attribute strings
        drugs_string: str = "drug_class=" + ",".join(unique_drugs)
        tools_string: str = "amr_tool=" + ",".join(tool_names)
        idents_string: str = "amr_tool_ident=" + ",".join(tool_identities)

        # Append to dictionary
        protein_attributes[protein] = [drugs_string, tools_string, idents_string]

    return protein_attributes


def parse_gff(cds_gff: str, output_file: str, protein_attributes: Dict[str, List[str]]) -> None:
    """
    Parse GFF file and integrate AMR annotations into output GFF.

    Args:
        cds_gff: Path to input GFF file containing CDS coordinates
        output_file: Path to output GFF file
        protein_attributes: Dictionary mapping protein IDs to attribute lists
    """
    with (
        fileinput.hook_compressed(cds_gff, "r", encoding="utf-8") as input_table,
        open(output_file, "w") as output_gff,
    ):
        output_gff.write("##gff-version 3\n")
        for line in input_table:
            line = line.rstrip()
            line_l: List[str] = line.split("\t")
            # Annotation lines have exactly 9 columns
            if len(line_l) == 9:
                (
                    contig,
                    seq_source,
                    seq_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attr,
                ) = line.rstrip().split("\t")
                if seq_type == "CDS":
                    features_list: List[str] = attr.split(";")
                    feature_id: str = features_list[0].split("=")[1]
                    if feature_id in protein_attributes:
                        new_attribute: str = attr + ";" + ";".join(protein_attributes[feature_id])
                        line_l.pop(-1)
                        line_l.append(new_attribute)
                        to_print: str = "\t".join(line_l)
                        output_gff.write(to_print + "\n")


def main() -> None:
    """
    Main function to orchestrate AMR annotation integration.
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="Integration of antimicrobial resistance genes annotation with deeparg, rgi and amrfinderplus into a single gff"
    )
    parser.add_argument("-d", "--deeparg_hamr", dest="deeparg_hamr", help="Result of deeparg tool after hamronization", required=False, default=None)
    parser.add_argument("-r", "--rgi_hamr", dest="rgi_hamr", help="Result of rgi tool after hamronization", required=False, default=None)
    parser.add_argument("-a", "--amrfp_out", dest="amrfp_out", help="Result of amrfinderplus tool", required=False, default=None)
    parser.add_argument(
        "-c", "--cds_gff", dest="cds_gff", help="GFF file containing the coordinates of the CDSs used for amr prediction", required=True, default=None
    )
    parser.add_argument("-o", "--output", dest="output", help="Name of the output file", required=False, default="integrated_result.gff")
    parser.add_argument("-v", "--verbose", dest="verbose", help="Enable verbose logging (DEBUG level)", action="store_true", default=False)
    args: argparse.Namespace = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(),  # Console output
        ],
    )

    optional_inputs: Dict[str, Optional[str]] = {
        "deeparg": args.deeparg_hamr,
        "rgi": args.rgi_hamr,
        "amrfinderplus": args.amrfp_out,
    }

    valid_inputs: List[str] = validate_inputs(optional_inputs)
    if valid_inputs:
        amr_annotation: DefaultDict[str, Dict[str, Dict[str, Union[List[str], str]]]] = defaultdict(dict)
        logger.info(f"The valid inputs provided for integration are: {', '.join(valid_inputs)}")
        logger.info("Parsing valid inputs")

        for name in valid_inputs:
            logger.debug(f"Processing {name} input file: {optional_inputs[name]}")
            if name == "amrfinderplus":
                amr_annotation = parse_amrfinderplus(amr_annotation, optional_inputs[name])
            else:
                amr_annotation = parse_hamronized(amr_annotation, optional_inputs[name])

        logger.info("Parsing the annotation information per protein")
        protein_attributes: Dict[str, List[str]] = parse_amr_dict(amr_annotation)
        logger.debug(f"Found AMR annotations for {len(protein_attributes)} proteins")

        logger.info("Parsing gff file and writing output file")
        parse_gff(args.cds_gff, args.output, protein_attributes)
        logger.info(f"Output written to: {args.output}")

    else:
        logger.info("No inputs to parse, the output file is not generated.")


if __name__ == "__main__":
    main()
