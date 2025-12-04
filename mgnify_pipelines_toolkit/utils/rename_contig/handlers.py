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

from typing import Dict

from mgnify_pipelines_toolkit.utils.io import open_file


class GFFHandler:
    """Handler for renaming contig IDs in GFF files."""

    @staticmethod
    def rename(input_file: str, output_file: str, mapping: Dict[str, str]) -> None:
        """
        Rename contig IDs in a GFF file.

        :param input_file: Input GFF file path
        :type input_file: str
        :param output_file: Output GFF file path
        :type output_file: str
        :param mapping: Dictionary mapping old_name -> new_name
        :type mapping: dict
        """
        with open_file(input_file, "rt") as in_f, open(output_file, "w") as out_f:
            in_fasta = False

            for line in in_f:
                # Handle FASTA section in GFF
                if line.startswith("##FASTA"):
                    in_fasta = True
                    out_f.write(line)
                    continue

                if in_fasta:
                    if line.startswith(">"):
                        seqid = line[1:].strip().split()[0]
                        new_id = mapping.get(seqid, seqid)
                        out_f.write(f">{new_id}\n")
                    else:
                        out_f.write(line)
                    continue

                # Handle sequence-region directive
                if line.startswith("##sequence-region"):
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        seqid = parts[1]
                        new_id = mapping.get(seqid, seqid)
                        parts[1] = new_id
                        out_f.write(" ".join(parts) + "\n")
                    else:
                        out_f.write(line)
                    continue

                # Handle comment lines
                if line.startswith("#"):
                    out_f.write(line)
                    continue

                # Handle feature lines
                columns = line.rstrip("\n").split("\t")
                if len(columns) >= 9:
                    seqid = columns[0]
                    new_id = mapping.get(seqid, seqid)
                    columns[0] = new_id
                    out_f.write("\t".join(columns) + "\n")
                else:
                    out_f.write(line)


class GenBankHandler:
    """Handler for renaming contig IDs in GenBank files."""

    @staticmethod
    def rename(input_file: str, output_file: str, mapping: Dict[str, str]) -> None:
        """
        Rename contig IDs in a GenBank file.

        :param input_file: Input GenBank file path
        :type input_file: str
        :param output_file: Output GenBank file path
        :type output_file: str
        :param mapping: Dictionary mapping old_name -> new_name
        :type mapping: dict
        """
        with open(input_file, "r") as in_f, open(output_file, "w") as out_f:
            for line in in_f:
                if line.startswith("LOCUS"):
                    parts = line.split()
                    if len(parts) >= 2:
                        old_id = parts[1]
                        new_id = mapping.get(old_id, old_id)
                        parts[1] = new_id

                        # Preserve Prokka LOCUS line format
                        if len(parts) == 7:
                            loc_field = f"{parts[0]:<12}"
                            id_field = f"{parts[1]:<21} "
                            length_field = f"{parts[2]} {parts[3]}    "
                            type_field = f"{parts[4]}     "
                            topology_field = f"{parts[5]}       "
                            date_field = parts[6]
                            out_f.write(f"{loc_field}{id_field}{length_field}{type_field}{topology_field}{date_field}\n")
                        else:
                            out_f.write(" ".join(parts) + "\n")
                    else:
                        out_f.write(line)
                else:
                    out_f.write(line)
