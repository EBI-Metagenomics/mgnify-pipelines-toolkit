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
Sequence/Contig Renaming Toolkit
=================================

A unified, modular toolkit for renaming sequences in FASTA, GFF, and GenBank files.
This module replaces and consolidates multiple pipeline-specific renaming scripts:
- metaviraverse.py (viral identifier preservation)
- virify.py (VirSorter metadata handling)
- mett.py (multi-format support)
- asa.py (simple renaming)
- map.py (size-based separation)

Features
--------
- Rename sequences in FASTA, GFF, and GenBank files
- Preserve viral metadata (viral IDs, phage-circular, prophage coordinates)
- Size-based file separation (1kb, 5kb, 100kb thresholds)
- Restore original names using mapping files
- Flexible parser and formatter modes for different pipelines

Supported Pipelines
-------------------

**Metaviraverse Pipeline**
    Preserves viral identifiers in format: contig|viral_id
    Usage: --parser virify --formatter virify

**Virify Pipeline**
    Handles VirSorter metadata (phage-circular, prophage coordinates)
    Usage: --parser virify --formatter virify

**METT Pipeline**
    Multi-format renaming (FASTA, GFF, GenBank)
    Usage: --parser base --formatter simple

**ASA Pipeline**
    Simple sequence renaming with mapping file generation
    Usage: --parser base --formatter simple

**MAP/Mobilome Pipeline**
    Size-based file separation (>1kb, >5kb, >=100kb)
    Usage: --separate-by-size

Quick Start - CLI
-----------------

Basic renaming::

    mpt rename -i input.fasta -o output_dir --prefix contig_

With viral metadata (Metaviraverse/Virify)::

    mpt rename -i input.fasta -o output_dir \\
        --parser virify --formatter virify --prefix contig_

Multi-format renaming (FASTA + GFF + GenBank)::

    mpt rename -i assembly.fasta prokka.gff prokka.gbk -o output_dir

Size-based separation (MAP pipeline)::

    mpt rename -i assembly.fasta -o output_dir --separate-by-size

Restore original names::

    mpt restore -i renamed.fasta -o output_dir -m mapping.tsv

Quick Start - Python API
-------------------------

Basic renaming::

    from mgnify_pipelines_toolkit.utils.rename_contig import rename_fasta

    mapping = rename_fasta(
        "input.fasta",
        "output.fasta",
        prefix="contig_"
    )

With viral metadata::

    from mgnify_pipelines_toolkit.utils.rename_contig import (
        rename_fasta,
        parse_virify_header
    )

    mapping = rename_fasta(
        "input.fasta",
        "output.fasta",
        prefix="contig_",
        parser_fn=parse_virify_header,
        formatter_mode="virify"
    )

Size-based separation::

    from mgnify_pipelines_toolkit.utils.rename_contig import rename_fasta_by_size

    mapping = rename_fasta_by_size(
        "assembly.fasta",
        "output_prefix",
        prefix="contig_"
    )

Multi-format renaming::

    from mgnify_pipelines_toolkit.utils.rename_contig import (
        rename_fasta,
        rename_gff,
        rename_genbank
    )

    # Rename FASTA first
    mapping = rename_fasta("input.fasta", "output.fasta", prefix="seq")

    # Convert mapping and rename other formats
    old_to_new = {old: new for new, (old, _) in mapping.items()}
    rename_gff("input.gff", "output.gff", old_to_new)
    rename_genbank("input.gbk", "output.gbk", old_to_new)

Module Structure
----------------
- cli.py: Command-line interface
- core.py: Core renaming functions
- parsers.py: Header parsing strategies
- writers.py: Header formatting strategies
- handlers.py: File format handlers (GFF, GenBank)
"""

# CLI entry point
from .cli import main

# Core renaming functions
from .core import (
    _parse_fasta_header,
    read_fasta,
    read_mapping_file,
    rename_fasta,
    rename_fasta_by_size,
    rename_genbank,
    rename_gff,
    restore_fasta,
    write_mapping_file,
)

# File format handlers
from .handlers import GenBankHandler, GFFHandler

# Parsers
from .parsers import FormatterMode, ParsedHeader, parse_header, parse_virify_header

# Writers
from .writers import format_header

__all__ = [
    # Core functions
    "rename_fasta",
    "rename_fasta_by_size",
    "restore_fasta",
    "read_fasta",
    "read_mapping_file",
    "write_mapping_file",
    # API compatibility wrappers
    "_parse_fasta_header",
    "rename_gff",
    "rename_genbank",
    # Parsers
    "parse_header",
    "parse_virify_header",
    "ParsedHeader",
    "FormatterMode",
    # Writers
    "format_header",
    # File format handlers
    "GFFHandler",
    "GenBankHandler",
    # CLI
    "main",
]
