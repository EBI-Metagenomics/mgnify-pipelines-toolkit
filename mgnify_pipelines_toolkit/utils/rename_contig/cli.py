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
CLI entry point for the contig renaming toolkit.

This module serves as the main entry point and imports commands from:
- cli_rename.py: Rename command implementation
- cli_restore.py: Restore command implementation
"""

import click


@click.group()
def main():
    """
    Universal contig renaming tool for FASTA, GFF, and GenBank files.

    A unified toolkit that replaces pipeline-specific renaming scripts:
    - metaviraverse (viral ID preservation)
    - virify (VirSorter metadata)
    - mett (multi-format support)
    - asa (simple renaming)
    - map (size-based separation)

    \b
    PIPELINE USAGE EXAMPLES:

    Metaviraverse (preserves viral IDs):
      mpt rename -i input.fasta -o outdir --parser virify --formatter virify

    Virify (preserves VirSorter metadata):
      mpt rename -i input.fasta -o outdir --parser virify --formatter virify

    METT (multi-format renaming):
      mpt rename -i assembly.fasta prokka.gff prokka.gbk -o outdir

    ASA (simple renaming):
      mpt rename -i input.fasta -o outdir --prefix seq

    MAP/Mobilome (size-based separation):
      mpt rename -i assembly.fasta -o outdir --separate-by-size

    \b
    For more info: mpt rename --help
    """
    pass


# Import and register commands from separate modules
from .cli_rename import rename
from .cli_restore import restore

main.add_command(rename)
main.add_command(restore)


if __name__ == "__main__":
    main()
