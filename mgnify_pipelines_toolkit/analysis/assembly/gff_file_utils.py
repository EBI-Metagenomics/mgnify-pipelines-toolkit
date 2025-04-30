#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
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


import csv


def write_results_to_file(
    outfile, header, main_gff_extended, fasta, ncrnas, trnas, crispr_annotations
):
    with open(outfile, "w") as file_out:
        file_out.write("\n".join(header) + "\n")
        contig_list = list(main_gff_extended.keys())
        # check if there are any contigs that don't have CDS; if so add them in
        contig_list = check_for_additional_keys(
            ncrnas, trnas, crispr_annotations, contig_list
        )
        # sort contigs by digit at the end of contig/genome accession
        if contig_list[0].startswith(
            "MGYG"
        ):  # e.g. 'MGYG000500002_1', 'MGYG000500002_2', 'MGYG000500002_3'
            contig_list = sorted(list(contig_list), key=lambda x: int(x.split("_")[-1]))
        elif contig_list[0].startswith(
            "ERZ"
        ):  # e.g. 'ERZ1049444', 'ERZ1049445', 'ERZ1049446'
            contig_list = sorted(
                list(contig_list), key=lambda x: int(x.split("ERZ")[-1])
            )
        for contig in contig_list:
            sorted_pos_list = sort_positions(
                contig, main_gff_extended, ncrnas, trnas, crispr_annotations
            )
            for pos in sorted_pos_list:
                for my_dict in (ncrnas, trnas, crispr_annotations, main_gff_extended):
                    if contig in my_dict and pos in my_dict[contig]:
                        for line in my_dict[contig][pos]:
                            if type(line) is str:
                                file_out.write(f"{line}\n")
                            else:
                                for element in line:
                                    file_out.write(element)
        for line in fasta:
            file_out.write(f"{line}\n")


def sort_positions(contig, main_gff_extended, ncrnas, trnas, crispr_annotations):
    sorted_pos_list = list()
    for my_dict in (main_gff_extended, ncrnas, trnas, crispr_annotations):
        if contig in my_dict:
            sorted_pos_list += list(my_dict[contig].keys())
    return sorted(list(set(sorted_pos_list)))


def check_for_additional_keys(ncrnas, trnas, crispr_annotations, contig_list):
    for my_dict in (ncrnas, trnas, crispr_annotations):
        dict_keys = set(my_dict.keys())
        absent_keys = dict_keys - set(contig_list)
        if absent_keys:
            contig_list = contig_list + list(absent_keys)
    return contig_list


def print_pseudogene_report(pseudogene_report_dict, pseudogene_report_file):
    with open(pseudogene_report_file, "w") as file_out:
        writer = csv.writer(file_out, delimiter="\t", lineterminator="\n")
        # Print header
        writer.writerow(
            [
                "ID",
                "Pseudogene according to Bakta/Prokka",
                "Pseudogene according to Pseudofinder",
                "AntiFam hit",
            ]
        )

        all_keys = ["gene_caller", "pseudofinder", "antifams"]
        for protein, attributes in pseudogene_report_dict.items():
            # Fill in missing attributes with False
            line = [protein] + [str(attributes.get(key, False)) for key in all_keys]
            writer.writerow(line)
