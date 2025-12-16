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
from typing import Callable, Dict, List, Optional, Tuple

from .parsers import ParsedHeader, parse_header
from .handlers import FASTAHandler, GFFHandler, GenBankHandler


# Wrapper functions for backward compatibility
# All FASTA operations now delegate to FASTAHandler


def read_fasta(file_path: str) -> List:
    """
    Read all sequences from a FASTA file.

    Wrapper for FASTAHandler.read() for backward compatibility.

    :param file_path: Path to FASTA file (may be gzipped)
    :type file_path: str
    :returns: List of Bio.SeqRecord.SeqRecord objects
    :rtype: list
    """
    return FASTAHandler.read(file_path)


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

    Wrapper for FASTAHandler.rename() for backward compatibility.

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
    return FASTAHandler.rename(
        input_file=input_file,
        output_file=output_file,
        prefix=prefix,
        parser_fn=parser_fn,
        formatter_mode=formatter_mode,
        mapping=mapping,
        separate_by_size=separate_by_size,
    )


def restore_fasta(
    input_file: str,
    output_file: str,
    mapping: Dict[str, str],
    parser_fn: Callable[[str], ParsedHeader] = parse_header,
    formatter_mode: str = "simple",
) -> None:
    """
    Restore original sequence names using a mapping.

    Wrapper for FASTAHandler.restore() for backward compatibility.

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
    return FASTAHandler.restore(
        input_file=input_file,
        output_file=output_file,
        mapping=mapping,
        parser_fn=parser_fn,
        formatter_mode=formatter_mode,
    )


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


# API Wrapper Functions for backward compatibility


def _parse_fasta_header(header: str) -> Tuple[str, str, List[str], Optional[str]]:
    """
    Parse FASTA header and return components as a tuple.

    Wrapper function for backward compatibility with test API.

    :param header: Raw FASTA header
    :type header: str
    :returns: Tuple of (seq_name, short_name, metadata, viral_id)
    :rtype: tuple
    """
    parsed = parse_header(header)
    return parsed.seq_name, parsed.short_name, parsed.metadata, parsed.viral_id


def rename_gff(input_file: str, output_file: str, mapping: Dict[str, str]) -> None:
    """
    Rename contig IDs in a GFF file.

    Wrapper for GFFHandler.rename for API compatibility.

    :param input_file: Input GFF file path
    :type input_file: str
    :param output_file: Output GFF file path
    :type output_file: str
    :param mapping: Dictionary mapping old_name -> new_name
    :type mapping: dict
    """
    return GFFHandler.rename(input_file, output_file, mapping)


def rename_genbank(input_file: str, output_file: str, mapping: Dict[str, str]) -> None:
    """
    Rename contig IDs in a GenBank file.

    Wrapper for GenBankHandler.rename for API compatibility.

    :param input_file: Input GenBank file path
    :type input_file: str
    :param output_file: Output GenBank file path
    :type output_file: str
    :param mapping: Dictionary mapping old_name -> new_name
    :type mapping: dict
    """
    return GenBankHandler.rename(input_file, output_file, mapping)


def rename_fasta_by_size(
    input_file: str,
    output_prefix: str,
    prefix: str = "contig_",
    parser_fn: Callable[[str], ParsedHeader] = parse_header,
    formatter_mode: str = "simple",
    mapping: Optional[Dict[str, str]] = None,
) -> Dict[str, Tuple[str, str]]:
    """
    Rename sequences in FASTA file and separate by size thresholds.

    Creates three output files:
    - {output_prefix}_1kb_contigs.fasta (sequences > 1000 bp)
    - {output_prefix}_5kb_contigs.fasta (sequences > 5000 bp)
    - {output_prefix}_100kb_contigs.fasta (sequences >= 100000 bp)

    This matches the functionality of the map.py script for the Mobilome pipeline.

    :param input_file: Input FASTA file path
    :type input_file: str
    :param output_prefix: Output file prefix (without extension)
    :type output_prefix: str
    :param prefix: Prefix for new names (default: 'contig_')
    :type prefix: str
    :param parser_fn: Header parser function (default: parse_header)
    :type parser_fn: callable
    :param formatter_mode: Metadata formatter mode ('simple', 'generic', or 'virify')
    :type formatter_mode: str
    :param mapping: Pre-defined mapping of original_name -> new_name (optional)
    :type mapping: dict or None
    :returns: Mapping of new_name -> (old_name, short_name)
    :rtype: dict
    """
    return rename_fasta(
        input_file=input_file,
        output_file=output_prefix,
        prefix=prefix,
        parser_fn=parser_fn,
        formatter_mode=formatter_mode,
        mapping=mapping,
        separate_by_size=True,
    )
