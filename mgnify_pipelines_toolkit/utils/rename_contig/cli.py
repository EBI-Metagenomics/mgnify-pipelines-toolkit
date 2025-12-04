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

from .core import (
    read_mapping_file,
    rename_fasta,
    restore_fasta,
    write_mapping_file,
)
from .handlers import GenBankHandler, GFFHandler
from .parsers import parse_header, parse_virify_header

SUPPORTED_FASTA_EXTS = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")
SUPPORTED_GBK_EXTS = (".gbk", ".gb", ".genbank")
SUPPORTED_GFF_EXTS = (".gff", ".gff3", ".gff.gz", ".gff3.gz")


@click.group()
def main():
    """
    Universal contig renaming tool for FASTA, GFF, and GenBank files.

    Rename sequences in multiple file formats using composable strategies.
    """
    pass


@main.command()
@click.argument("input", nargs=-1, required=True, type=click.Path(exists=True))
@click.option("-o", "--outdir", required=True, type=click.Path(), help="Output directory")
@click.option("-p", "--prefix", default="contig", help="Prefix for new names (default: contig)")
@click.option("-m", "--map", type=click.Path(), help="Output mapping file (TSV)")
@click.option(
    "--use-mapping",
    type=click.Path(exists=True),
    help="Use existing mapping file instead of generating new names",
)
@click.option("--from-col", default="original", help="Source column in mapping file (default: original)")
@click.option("--to-col", default="renamed", help="Target column in mapping file (default: renamed)")
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
    help="Metadata formatter strategy (default: simple). 'generic' preserves all metadata, 'virify' handles Virify IDs",
)
@click.option("--separate-by-size", is_flag=True, help="Separate output by contig size (_1kb, _5kb, _100kb)")
def rename(
    input,
    outdir,
    prefix,
    map,
    use_mapping,
    from_col,
    to_col,
    parser,
    formatter,
    separate_by_size,
):
    """
    Rename contigs in FASTA, GFF, and/or GenBank files.

    Supports multiple input files of different formats.

    Example with viral metadata:

        mpt rename -i input.fasta -o output_dir --parser viral --formatter viral

    :param input: Input file(s) (FASTA, GFF, and/or GenBank)
    :param outdir: Output directory
    :param prefix: Prefix for new sequence names
    :param map: Output mapping file (TSV)
    :param use_mapping: Use existing mapping file instead of generating new names
    :param from_col: Source column in mapping file
    :param to_col: Target column in mapping file
    :param parser: Header parser to use ('base' or 'viral')
    :param formatter: Metadata formatter to use ('simple', 'generic', or 'viral')
    :param separate_by_size: Separate output by contig size
    """
    # Create output directory if it doesn't exist
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)

    # Select parser function
    parser_fn = parse_virify_header if parser == "virify" else parse_header

    # Categorize input files by type
    fasta_files = []
    gff_files = []
    gbk_files = []

    for input_file in input:
        input_path = Path(input_file)
        file_ext = "".join(input_path.suffixes).lower()

        is_fasta = any(file_ext.endswith(ext) for ext in SUPPORTED_FASTA_EXTS)
        is_gff = any(file_ext.endswith(ext) for ext in SUPPORTED_GFF_EXTS)
        is_gbk = any(file_ext.endswith(ext) for ext in SUPPORTED_GBK_EXTS)

        if is_fasta:
            fasta_files.append(input_file)
        elif is_gff:
            gff_files.append(input_file)
        elif is_gbk:
            gbk_files.append(input_file)
        else:
            raise click.BadParameter(
                f"Unsupported file type: {input_file}. Supported: {SUPPORTED_FASTA_EXTS + SUPPORTED_GFF_EXTS + SUPPORTED_GBK_EXTS}"
            )

    # Read existing mapping if provided
    old_to_new_mapping = None
    if use_mapping:
        old_to_new_mapping = read_mapping_file(use_mapping, from_col=from_col, to_col=to_col)

    # Combined mapping from all FASTA files
    combined_name_mapping = {}

    # First pass: Rename all FASTA files
    if fasta_files:
        click.echo(f"Processing {len(fasta_files)} FASTA file(s)...")

        for fasta_file in fasta_files:
            input_path = Path(fasta_file)

            if separate_by_size:
                if len(fasta_files) > 1:
                    raise click.BadParameter("--separate-by-size only supports a single FASTA input file")

                output_prefix = str(outdir_path / input_path.stem)

                click.echo(f"  Renaming FASTA with size separation: {fasta_file}")

                name_mapping = rename_fasta(
                    fasta_file,
                    output_prefix,
                    prefix=prefix,
                    parser_fn=parser_fn,
                    formatter_mode=formatter,
                    mapping=old_to_new_mapping,
                    separate_by_size=True,
                )

                combined_name_mapping.update(name_mapping)
                click.echo(f"  → {output_prefix}_1kb_contigs.fasta")
                click.echo(f"  → {output_prefix}_5kb_contigs.fasta")
                click.echo(f"  → {output_prefix}_100kb_contigs.fasta")
                click.echo(f"  Total sequences renamed: {len(name_mapping)}")

            else:
                output_file = outdir_path / input_path.name

                click.echo(f"  Renaming FASTA: {fasta_file}")

                name_mapping = rename_fasta(
                    fasta_file,
                    str(output_file),
                    prefix=prefix,
                    parser_fn=parser_fn,
                    formatter_mode=formatter,
                    mapping=old_to_new_mapping,
                    separate_by_size=False,
                )

                combined_name_mapping.update(name_mapping)
                click.echo(f"  → {output_file} ({len(name_mapping)} sequences)")

        # Convert combined mapping to old -> new format for GFF/GenBank
        if not old_to_new_mapping:
            old_to_new_mapping = {old: new for new, (old, short) in combined_name_mapping.items()}
            # Also add short names as keys
            for new, (old, short) in combined_name_mapping.items():
                if short != old:
                    old_to_new_mapping[short] = new

        # Write combined mapping file
        if map and combined_name_mapping:
            write_mapping_file(combined_name_mapping, map)
            click.echo(f"Wrote mapping to: {map}")

    # Second pass: Rename GFF files
    if gff_files:
        if not old_to_new_mapping:
            raise click.UsageError("GFF renaming requires either FASTA file(s) in input or --use-mapping")

        click.echo(f"\nProcessing {len(gff_files)} GFF file(s)...")
        for gff_file in gff_files:
            input_path = Path(gff_file)
            output_file = outdir_path / input_path.name

            click.echo(f"  Renaming GFF: {gff_file}")
            GFFHandler.rename(gff_file, str(output_file), old_to_new_mapping)
            click.echo(f"  → {output_file}")

    # Third pass: Rename GenBank files
    if gbk_files:
        if not old_to_new_mapping:
            raise click.UsageError("GenBank renaming requires either FASTA file(s) in input or --use-mapping")

        click.echo(f"\nProcessing {len(gbk_files)} GenBank file(s)...")
        for gbk_file in gbk_files:
            input_path = Path(gbk_file)
            output_file = outdir_path / input_path.name

            click.echo(f"  Renaming GenBank: {gbk_file}")
            GenBankHandler.rename(gbk_file, str(output_file), old_to_new_mapping)
            click.echo(f"  → {output_file}")

    click.echo(f"\nAll files processed successfully. Output directory: {outdir_path}")


@main.command()
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

    Example with Virify metadata:

        mpt restore -i renamed.fasta -o output_dir -m mapping.tsv --parser virify --formatter virify

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


if __name__ == "__main__":
    main()
