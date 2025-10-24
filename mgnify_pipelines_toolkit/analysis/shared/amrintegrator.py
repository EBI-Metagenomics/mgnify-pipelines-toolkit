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
import gzip
import os


def validate_inputs(optional_inputs):
    valid_inputs = []
    for name, path in optional_inputs.items():
        if path and os.path.exists(path):
            valid_inputs.append(name)
        elif path:
            print(f"Warning: file not found for '{name}' → {path}")
    return valid_inputs


def open_file(filename):
    """
    Open a file, handling both compressed (.gz) and uncompressed files.
    Returns a file handle that can be used for reading.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")  # 'rt' for text mode
    else:
        return open(filename, "r")


def parse_hamronized(amr_annotation, input_file):
    with open_file(input_file) as input_table:
        next(input_table)
        for line in input_table:
            line_l = line.rstrip().split("\t")
            analysis_software_name = line_l[6]
            drug_class = (
                line_l[13].replace("macrolide antibiotic; lincosamide antibiotic; streptogramin antibiotic", "MLS").replace(" antibiotic", "")
            )
            protein_id = line_l[20].split(" ")[0]
            seq_identity = line_l[35]

            drug_class_list = [cls.strip() for cls in drug_class.split(";")]

            if protein_id not in amr_annotation:
                amr_annotation[protein_id] = {}

            amr_annotation[protein_id][analysis_software_name] = {"drug_class": drug_class_list, "seq_identity": seq_identity}

    return amr_annotation


def parse_amrfinderplus(amr_annotation, input_file):
    with open_file(input_file) as input_table:
        next(input_table)
        for line in input_table:
            drug_class_list = []
            line_l = line.rstrip().split("\t")
            analysis_software_name = "amrfinderplus"
            drug_class_list.append(line_l[6].lower().replace("lincosamide/macrolide/streptogramin", "MLS"))
            protein_id = line_l[0]
            seq_identity = line_l[12]

            if protein_id not in amr_annotation:
                amr_annotation[protein_id] = {}

            amr_annotation[protein_id][analysis_software_name] = {"drug_class": drug_class_list, "seq_identity": seq_identity}

    return amr_annotation


def parse_amr_dict(amr_annotation):
    protein_attributes = {}

    for protein, tools in amr_annotation.items():
        all_drugs = []
        tool_names = []
        tool_identities = []

        for tool, info in tools.items():
            # Collect drug classes
            all_drugs.extend(info.get("drug_class", []))
            # Collect tool name and identity value
            tool_names.append(tool)
            tool_identities.append(info.get("seq_identity", "NA"))

        # Clean up and deduplicate drug classes
        unique_drugs = list(dict.fromkeys(all_drugs))
        unique_drugs = [item.replace(" ", "_") for item in unique_drugs]

        # Build attribute strings
        drugs_string = "drug_class=" + ",".join(unique_drugs)
        tools_string = "amr_tool=" + ",".join(tool_names)
        idents_string = "amr_tool_ident=" + ",".join(tool_identities)

        # Append to dictionary
        protein_attributes[protein] = [drugs_string, tools_string, idents_string]

    return protein_attributes


def parse_gff(cds_gff, output_file, protein_attributes):
    with (
        open_file(cds_gff) as input_table,
        open(output_file, "w") as output_gff,
    ):
        output_gff.write("##gff-version 3\n")
        for line in input_table:
            line = line.rstrip()
            line_l = line.split("\t")
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
                    features_list = attr.split(";")
                    feature_id = features_list[0].split("=")[1]
                    if feature_id in protein_attributes:
                        new_attribute = attr + ";".join(protein_attributes[feature_id])
                        line_l.pop(-1)
                        line_l.append(new_attribute)
                        to_print = "\t".join(line_l)
                        output_gff.write(to_print + "\n")

    # Example of the output:
    ##gff-version 3
    # b20_05_1.circ	Prodigal_v2.6.3	CDS	2925164	2927089	152.2	-	0	ID=b20_05_1.circ_2009;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=5bp;gc_cont=0.398;conf=100.00;score=152.15;cscore=148.05;sscore=4.10;rscore=-0.03;uscore=0.76;tscore=4.02;drug_class=tetracycline;amr_tool=rgi,amrfinderplus;amr_tool_ident=96.41,99.53
    # b20_05_1.circ	Prodigal_v2.6.3	CDS	3085110	3088241	420.2	+	0	ID=b20_05_1.circ_2130;partial=00;start_type=ATG;rbs_motif=TAAAA;rbs_spacer=9bp;gc_cont=0.462;conf=99.99;score=420.21;cscore=405.19;sscore=15.02;rscore=8.54;uscore=1.45;tscore=4.02;drug_class=fluoroquinolone,tetracycline;amr_tool=rgi;amr_tool_ident=42.24


def main():
    parser = argparse.ArgumentParser(
        description="Integration of antimicrobial resistance genes annotation with deeparg, rgi and amrfinderplus into a single gff"
    )
    parser.add_argument("-d", "--deeparg_hamr", dest="deeparg_hamr", help="Result of deeparg tool after hamronization", required=False, default=None)
    parser.add_argument("-r", "--rgi_hamr", dest="rgi_hamr", help="Result of rgi tool after hamronization", required=False, default=None)
    parser.add_argument("-a", "--amrfp_out", dest="amrfp_out", help="Result of amrfinderplus tool", required=False, default=None)
    parser.add_argument(
        "-c", "--cds_gff", dest="cds_gff", help="GFF file containing the coordinates of the CDSs used for amr prediction", required=True, default=None
    )
    parser.add_argument("-o", "--output", dest="output", help="Name of the output file", required=False, default="integrated_result.gff")
    args = parser.parse_args()

    optional_inputs = {
        "deeparg": args.deeparg_hamr,
        "rgi": args.rgi_hamr,
        "amrfinderplus": args.amrfp_out,
    }

    valid_inputs = validate_inputs(optional_inputs)
    if valid_inputs:
        amr_annotation = {}
        print(f"The valid inputs provided for integration are: {', '.join(valid_inputs)}")
        print("Parsing valid inputs")

        for name in valid_inputs:
            if name == "amrfinderplus":
                amr_annotation = parse_amrfinderplus(amr_annotation, optional_inputs[name])
            else:
                amr_annotation = parse_hamronized(amr_annotation, optional_inputs[name])

        print("Parsing the annottion information per protein")
        protein_attributes = parse_amr_dict(amr_annotation)

        print("Parsing gff file and writing output file")
        parse_gff(args.cds_gff, args.output, protein_attributes)

    else:
        print("No inputs to parse, generating empty output.")
        with open(args.output, "w") as out:
            out.write("##gff-version 3\n")


if __name__ == "__main__":
    main()
