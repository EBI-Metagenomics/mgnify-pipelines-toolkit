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

from pathlib import Path

import click

from .constants import SUPPORTED_FASTA_EXTS
from .core import read_mapping_file, restore_fasta
from .parsers import parse_header, parse_virify_header


@click.command()
@click.argument("input", nargs=-1, required=True, type=click.Path(exists=True))
@click.option("-o", "--outdir", required=True, type=click.Path(), help="Output directory")
@click.option("-m", "--map", required=True, type=click.Path(exists=True), help="Mapping file (TSV)")
@click.option("--from-col", default="renamed", help="Source column in mapping file (default: renamed)")
@click.option("--to-col", default="original", help="Target column in mapping file (default: original)")
@click.option(
    "--parser",
    type=click.Choice(["base", "virify"]),
    default="base",
    help="Header parser strategy (default: base). Use 'virify' for Virify metadata",
)
@click.option(
    "--formatter",
    type=click.Choice(["simple", "generic", "virify"]),
    default="simple",
    help="Metadata formatter strategy (default: simple)",
)
def restore(input, outdir, map, from_col, to_col, parser, formatter):
    """
    Restore original contig names using a mapping file.

    Reverses the renaming process by looking up original names in the mapping
    file. Preserves metadata when using appropriate parser/formatter modes.

    \b
    USAGE EXAMPLES:

    Basic restore (renamed -> original):
      mpt restore -i renamed.fasta -o outdir -m mapping.tsv

    Restore with viral metadata:
      mpt restore -i renamed.fasta -o outdir -m mapping.tsv \\
        --parser virify --formatter virify

    Restore from short names:
      mpt restore -i renamed.fasta -o outdir -m mapping.tsv \\
        --from-restore temporary --to-restore short

    \b
    MAPPING FILE COLUMNS:
      original:  Full original sequence name
      renamed:   New sequence name (e.g., contig_1)
      short:     Short name (first token before space)

    \b
    TYPICAL WORKFLOW:
      1. Rename sequences: mpt rename -i input.fasta -o outdir -m mapping.tsv
      2. Process files...
      3. Restore names: mpt restore -i processed.fasta -o restored -m mapping.tsv

    :param input: Input FASTA file(s) with renamed sequences
    :param outdir: Output directory
    :param map: Mapping file (TSV)
    :param from_col: Source column in mapping file
    :param to_col: Target column in mapping file
    :param parser: Header parser to use ('base' or 'virify')
    :param formatter: Metadata formatter to use ('simple', 'generic', or 'virify')
    """
    # Create output directory if it doesn't exist
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)

    # Select parser function
    parser_fn = parse_virify_header if parser == "virify" else parse_header

    # Read mapping file
    mapping = read_mapping_file(map, from_col=from_col, to_col=to_col)

    click.echo(f"Restoring {len(input)} file(s) using mapping from {map}")
    click.echo(f"Mapping: {from_col} -> {to_col}\n")

    # Process each input file
    for input_file in input:
        input_path = Path(input_file)
        file_ext = "".join(input_path.suffixes).lower()

        is_fasta = any(file_ext.endswith(ext) for ext in SUPPORTED_FASTA_EXTS)

        if not is_fasta:
            click.secho(
                f"Warning: Restore mode currently only supports FASTA files, skipping {input_file}",
                fg="yellow",
            )
            continue

        output_file = outdir_path / input_path.name

        click.echo(f"  Restoring: {input_file}")

        restore_fasta(
            input_file,
            str(output_file),
            mapping=mapping,
            parser_fn=parser_fn,
            formatter_mode=formatter,
        )
        click.echo(f"  → {output_file}")

    click.echo(f"\nAll files restored successfully. Output directory: {outdir_path}")
