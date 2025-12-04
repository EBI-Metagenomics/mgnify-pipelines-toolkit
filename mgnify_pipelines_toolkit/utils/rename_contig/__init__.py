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
Sequence/Contig renaming toolkit.

A simple, modular library for renaming sequences in FASTA, GFF, and GenBank files.

Example usage:

    from mgnify_pipelines_toolkit.utils.rename_contig import (
        rename_fasta,
        parse_viral_header,
        read_mapping_file,
    )

    # Rename sequences with viral metadata parsing
    mapping = rename_fasta(
        "input.fasta",
        "output.fasta",
        prefix="contig_",
        parser_fn=parse_viral_header,
        formatter_mode="viral"
    )
"""

# CLI entry point
from .cli import main

# Core renaming functions
from .core import (
    read_fasta,
    read_mapping_file,
    rename_fasta,
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
    "restore_fasta",
    "read_fasta",
    "read_mapping_file",
    "write_mapping_file",
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
