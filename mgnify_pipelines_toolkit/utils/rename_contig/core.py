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


import csv
import sys
from typing import Callable, Dict, List, Optional, Tuple

from Bio import SeqIO

from mgnify_pipelines_toolkit.utils.io import open_file

from .parsers import ParsedHeader, parse_header
from .writers import format_header


def read_fasta(file_path: str) -> List:
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


def rename_fasta(
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
    records = read_fasta(input_file)
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


def restore_fasta(
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

    records = read_fasta(input_file)

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


def read_mapping_file(map_file: str, from_col: str = "original", to_col: str = "renamed") -> Dict[str, str]:
    """
    Read a mapping file and return a dictionary.

    :param map_file: Path to TSV mapping file
    :type map_file: str
    :param from_col: Source column name
    :type from_col: str
    :param to_col: Target column name
    :type to_col: str
    :returns: Dictionary mapping from_col values to to_col values
    :rtype: dict
    :raises ValueError: If required columns are missing
    """
    mapping = {}
    with open(map_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if from_col not in row or to_col not in row:
                raise ValueError(f"Mapping file must contain '{from_col}' and '{to_col}' columns")
            mapping[row[from_col]] = row[to_col]
    return mapping


def write_mapping_file(
    mapping: Dict[str, Tuple[str, str]],
    output_file: str,
) -> None:
    """
    Write a mapping file with three columns: original, renamed, short.

    :param mapping: Dictionary mapping new_name -> (old_name, short_name)
    :type mapping: dict
    :param output_file: Output TSV file path
    :type output_file: str
    """
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["original", "renamed", "short"])
        for new_name, (old_name, short_name) in mapping.items():
            writer.writerow([old_name, new_name, short_name])
