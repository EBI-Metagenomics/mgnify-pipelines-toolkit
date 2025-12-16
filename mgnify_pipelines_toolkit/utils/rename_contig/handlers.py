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

"""
File format handlers for renaming sequences in different formats.

Provides consistent class-based interface for handling:
- FASTA files (FASTAHandler)
- GFF files (GFFHandler)
- GenBank files (GenBankHandler)
"""

import sys
from typing import Callable, Dict, List, Optional, Tuple

from Bio import SeqIO

from mgnify_pipelines_toolkit.utils.io import open_file

from .parsers import ParsedHeader, parse_header
from .writers import format_header


class FASTAHandler:
    """Handler for renaming sequences in FASTA files."""

    @staticmethod
    def read(file_path: str) -> List:
        """
        Read all sequences from a FASTA file.

        Uses BioPython's SeqIO for efficient parsing.

        :param file_path: Path to FASTA file (may be gzipped)
        :type file_path: str
        :returns: List of Bio.SeqRecord.SeqRecord objects
        :rtype: list
        """
        records = []
        with open_file(file_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                records.append(record)
        return records

    @staticmethod
    def rename(
        input_file: str,
        output_file: str,
        prefix: str = "contig",
        parser_fn: Callable[[str], ParsedHeader] = parse_header,
        formatter_mode: str = "simple",
        mapping: Optional[Dict[str, str]] = None,
        separate_by_size: bool = False,
    ) -> Dict[str, Tuple[str, str]]:
        """
        Rename sequences in FASTA file.

        Creates output file with renamed sequences. Can optionally separate
        output by size thresholds or use pre-defined mapping.

        :param input_file: Input FASTA file path
        :type input_file: str
        :param output_file: Output FASTA file path (or prefix if separate_by_size=True)
        :type output_file: str
        :param prefix: Prefix for new names (default: 'contig')
        :type prefix: str
        :param parser_fn: Header parser function (default: parse_header)
        :type parser_fn: callable
        :param formatter_mode: Metadata formatter mode ('simple', 'generic', or 'virify')
        :type formatter_mode: str
        :param mapping: Pre-defined mapping of original_name -> new_name (optional)
        :type mapping: dict or None
        :param separate_by_size: Separate output by size thresholds (1kb, 5kb, 100kb)
        :type separate_by_size: bool
        :returns: Mapping of new_name -> (old_name, short_name)
        :rtype: dict
        """
        records = FASTAHandler.read(input_file)
        name_mapping = {}
        counter = 1

        if separate_by_size:
            # Create three output files with size-based separation
            output_1kb = f"{output_file}_1kb_contigs.fasta"
            output_5kb = f"{output_file}_5kb_contigs.fasta"
            output_100kb = f"{output_file}_100kb_contigs.fasta"

            with open(output_1kb, "w") as f1kb, open(output_5kb, "w") as f5kb, open(output_100kb, "w") as f100kb:
                for record in records:
                    # Parse header
                    parsed = parser_fn(record.id)
                    seq_length = len(record.seq)

                    # Generate new name (use mapping if available)
                    if mapping and parsed.seq_name in mapping:
                        new_name = mapping[parsed.seq_name]
                    elif mapping and parsed.short_name in mapping:
                        new_name = mapping[parsed.short_name]
                    else:
                        new_name = f"{prefix}{counter}"
                        counter += 1

                    # Store mapping
                    name_mapping[new_name] = (parsed.seq_name, parsed.short_name)

                    # Format header
                    header = format_header(new_name, parsed, formatter_mode)

                    # Write to appropriate files based on length
                    if seq_length > 1000:
                        f1kb.write(header + "\n")
                        f1kb.write(str(record.seq).upper() + "\n")
                    if seq_length > 5000:
                        f5kb.write(header + "\n")
                        f5kb.write(str(record.seq).upper() + "\n")
                    if seq_length >= 100000:
                        f100kb.write(header + "\n")
                        f100kb.write(str(record.seq).upper() + "\n")

        else:
            # Single output file
            with open(output_file, "w") as out_f:
                for record in records:
                    # Parse header
                    parsed = parser_fn(record.id)

                    # Generate new name (use mapping if available)
                    if mapping and parsed.seq_name in mapping:
                        new_name = mapping[parsed.seq_name]
                    elif mapping and parsed.short_name in mapping:
                        new_name = mapping[parsed.short_name]
                    else:
                        new_name = f"{prefix}{counter}"
                        counter += 1

                    # Store mapping
                    name_mapping[new_name] = (parsed.seq_name, parsed.short_name)

                    # Format header
                    header = format_header(new_name, parsed, formatter_mode)

                    # Write sequence
                    out_f.write(header + "\n")
                    out_f.write(str(record.seq) + "\n")

        return name_mapping

    @staticmethod
    def restore(
        input_file: str,
        output_file: str,
        mapping: Dict[str, str],
        parser_fn: Callable[[str], ParsedHeader] = parse_header,
        formatter_mode: str = "simple",
    ) -> None:
        """
        Restore original sequence names using a mapping.

        :param input_file: Input FASTA file with renamed sequences
        :type input_file: str
        :param output_file: Output FASTA file with original names
        :type output_file: str
        :param mapping: Dictionary mapping current_name -> original_name
        :type mapping: dict
        :param parser_fn: Header parser function (default: parse_header)
        :type parser_fn: callable
        :param formatter_mode: Metadata formatter mode ('simple', 'generic', or 'virify')
        :type formatter_mode: str
        """
        from .parsers import _extract_short_name

        records = FASTAHandler.read(input_file)

        with open(output_file, "w") as out_f:
            for record in records:
                # Parse header
                parsed = parser_fn(record.id)

                # Try to find mapping
                original = mapping.get(parsed.seq_name, mapping.get(parsed.short_name, None))
                if not original:
                    print(
                        f"Warning: No mapping found for {parsed.seq_name}. Using current name.",
                        file=sys.stderr,
                    )
                    original = parsed.seq_name

                # Create a new ParsedHeader with original name for formatting
                restored_parsed = ParsedHeader(
                    seq_name=original,
                    short_name=_extract_short_name(original),
                    metadata=parsed.metadata,
                    viral_id=parsed.viral_id,
                )

                # Format and write header
                header = format_header(original, restored_parsed, formatter_mode)
                out_f.write(header + "\n")

                # Write sequence
                out_f.write(str(record.seq) + "\n")


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
