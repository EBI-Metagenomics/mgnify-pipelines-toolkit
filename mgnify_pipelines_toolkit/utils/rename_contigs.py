#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
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
Universal contig renaming script for FASTA, GFF, and GenBank files.
Combines functionality from multiple pipeline-specific scripts.
"""

import argparse
import csv
import gzip
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

try:
    from Bio import SeqIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

try:
    import pyfastx
    HAS_PYFASTX = True
except ImportError:
    HAS_PYFASTX = False


SUPPORTED_FASTA_EXTS = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")
SUPPORTED_GBK_EXTS = (".gbk", ".gb", ".genbank")
SUPPORTED_GFF_EXTS = (".gff", ".gff3", ".gff.gz", ".gff3.gz")


def _open_file(file_path, mode="r"):
    """Open a file, handling both compressed and uncompressed formats."""
    if file_path.endswith(".gz"):
        if "b" not in mode and "t" not in mode:
            mode = mode + "t"
        return gzip.open(file_path, mode)
    return open(file_path, mode)


def _parse_fasta_header(header: str) -> Tuple[str, List[str]]:
    """
    Parse a FASTA header to extract sequence name and metadata.

    Handles:
    - Viral identifiers with "|" separator
    - VirSorter metadata (phage-circular, prophage-<start>:<end>)

    Args:
        header: FASTA header line (with or without ">")

    Returns:
        Tuple of (sequence_name, metadata_list)
    """
    clean = header.replace(">", "").replace("\n", "").strip()

    # Check for viral identifier with "|"
    parts = clean.split("|")
    seq_name = parts[0]

    metadata = []

    # Check for phage-circular in header
    if "phage-circular" in clean:
        metadata.append("phage-circular")

    # Check for viral identifier
    viral_identifier = None
    if len(parts) > 1:
        # Check if second part is a viral identifier (not prophage metadata)
        if not parts[1].startswith("prophage-") and not parts[1] == "phage-circular":
            viral_identifier = parts[1]
            if viral_identifier not in metadata:
                metadata.append(viral_identifier)

    # Parse prophage metadata
    for part in parts:
        match = re.search(r"prophage-\d+:\d+", part)
        if match and match[0] not in metadata:
            metadata.append(match[0])

    # Get short name (first token before space)
    short_name = seq_name.split()[0]

    return seq_name, short_name, metadata, viral_identifier


def read_mapping_file(map_file: str, from_col: str = "original", to_col: str = "renamed") -> Dict[str, str]:
    """
    Read a mapping file and return a dictionary.

    Args:
        map_file: Path to TSV mapping file
        from_col: Source column name
        to_col: Target column name

    Returns:
        Dictionary mapping from_col values to to_col values
    """
    mapping = {}
    with open(map_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if from_col not in row or to_col not in row:
                raise ValueError(f"Mapping file must contain '{from_col}' and '{to_col}' columns")
            mapping[row[from_col]] = row[to_col]
    return mapping


def rename_fasta(
    input_file: str,
    output_file: str,
    prefix: str = "contig",
    mapping: Optional[Dict[str, str]] = None,
    preserve_metadata: bool = False,
    use_pyfastx: bool = False
) -> Dict[str, Tuple[str, str]]:
    """
    Rename sequences in a FASTA file.

    Args:
        input_file: Input FASTA file path
        output_file: Output FASTA file path
        prefix: Prefix for new sequence names (default: "contig")
        mapping: Optional pre-defined mapping dict (old_name -> new_name)
        preserve_metadata: Whether to preserve metadata in headers
        use_pyfastx: Use pyfastx library if available (faster for large files)

    Returns:
        Dictionary mapping new_name -> (old_name, short_name)
    """
    name_mapping = {}  # new_name -> (old_name, short_name)
    counter = 1

    # Use pyfastx if requested and available
    if use_pyfastx and HAS_PYFASTX and not input_file.endswith(".gz"):
        fasta = pyfastx.Fasta(input_file, build_index=False)
        with open(output_file, "w") as out_f:
            for name, seq in fasta:
                if mapping and name in mapping:
                    new_name = mapping[name]
                else:
                    new_name = f"{prefix}{counter}" if not prefix.endswith("_") else f"{prefix}{counter}"
                    counter += 1

                short_name = name.split()[0]
                name_mapping[new_name] = (name, short_name)
                out_f.write(f">{new_name}\n{seq}\n")

    # Use BioPython if available
    elif HAS_BIOPYTHON and not preserve_metadata:
        records = []
        with _open_file(input_file, "r") as in_f:
            for i, rec in enumerate(SeqIO.parse(in_f, "fasta"), 1):
                old_name = rec.id
                if mapping and old_name in mapping:
                    new_name = mapping[old_name]
                else:
                    new_name = f"{prefix}{counter}" if not prefix.endswith("_") else f"{prefix}{counter}"
                    counter += 1

                short_name = old_name.split()[0]
                name_mapping[new_name] = (old_name, short_name)

                rec.id = new_name
                rec.name = new_name
                rec.description = ""
                records.append(rec)

        with open(output_file, "w") as out_f:
            SeqIO.write(records, out_f, "fasta")

    # Fallback: manual parsing (handles all cases including metadata preservation)
    else:
        with _open_file(input_file, "r") as in_f, open(output_file, "w") as out_f:
            for line in in_f:
                if line.startswith(">"):
                    original_header = line[1:].strip()
                    seq_name, short_name, metadata, viral_id = _parse_fasta_header(original_header)

                    if mapping and seq_name in mapping:
                        new_name = mapping[seq_name]
                    elif mapping and short_name in mapping:
                        new_name = mapping[short_name]
                    else:
                        new_name = f"{prefix}{counter}" if not prefix.endswith("_") else f"{prefix}{counter}"
                        counter += 1

                    name_mapping[new_name] = (seq_name, short_name)

                    # Write header with or without metadata
                    if preserve_metadata and metadata:
                        # Add viral identifier first if it exists
                        if viral_id:
                            out_f.write(f">{new_name}|{viral_id}")
                            # Add other metadata
                            other_metadata = [m for m in metadata if m != viral_id]
                            if other_metadata:
                                out_f.write("|" + "|".join(other_metadata))
                            out_f.write("\n")
                        else:
                            out_f.write(f">{new_name}|{'|'.join(metadata)}\n")
                    else:
                        out_f.write(f">{new_name}\n")
                else:
                    out_f.write(line)

    return name_mapping


def rename_fasta_by_size(
    input_file: str,
    output_prefix: str,
    prefix: str = "contig",
    mapping: Optional[Dict[str, str]] = None,
    preserve_metadata: bool = False
) -> Dict[str, Tuple[str, str]]:
    """
    Rename sequences in a FASTA file and separate by size thresholds.
    Creates three output files: _1kb_contigs.fasta, _5kb_contigs.fasta, _100kb_contigs.fasta

    Args:
        input_file: Input FASTA file path
        output_prefix: Prefix for output files (will create prefix_1kb_contigs.fasta, etc.)
        prefix: Prefix for new sequence names (default: "contig")
        mapping: Optional pre-defined mapping dict (old_name -> new_name)
        preserve_metadata: Whether to preserve metadata in headers

    Returns:
        Dictionary mapping new_name -> (old_name, short_name)
    """
    output_1kb = f"{output_prefix}_1kb_contigs.fasta"
    output_5kb = f"{output_prefix}_5kb_contigs.fasta"
    output_100kb = f"{output_prefix}_100kb_contigs.fasta"

    name_mapping = {}
    counter = 1

    # Use BioPython for length filtering
    if HAS_BIOPYTHON:
        with _open_file(input_file, "r") as in_f:
            with open(output_1kb, "w") as out_1kb, \
                 open(output_5kb, "w") as out_5kb, \
                 open(output_100kb, "w") as out_100kb:

                for rec in SeqIO.parse(in_f, "fasta"):
                    old_name = rec.id
                    seq_length = len(rec.seq)
                    my_chain = str(rec.seq).upper()

                    # Get new name
                    if mapping and old_name in mapping:
                        new_name = mapping[old_name]
                    else:
                        new_name = f"{prefix}{counter}" if not prefix.endswith("_") else f"{prefix}{counter}"
                        counter += 1

                    short_name = old_name.split()[0]
                    name_mapping[new_name] = (old_name, short_name)

                    # Handle metadata if needed
                    if preserve_metadata:
                        _, _, metadata, viral_id = _parse_fasta_header(rec.description)
                        if metadata:
                            if viral_id:
                                header = f">{new_name}|{viral_id}"
                                other_metadata = [m for m in metadata if m != viral_id]
                                if other_metadata:
                                    header += "|" + "|".join(other_metadata)
                                header += "\n"
                            else:
                                header = f">{new_name}|{'|'.join(metadata)}\n"
                        else:
                            header = f">{new_name}\n"
                    else:
                        header = f">{new_name}\n"

                    # Write to appropriate files based on length
                    if seq_length > 1000:
                        out_1kb.write(header)
                        out_1kb.write(my_chain + "\n")
                    if seq_length > 5000:
                        out_5kb.write(header)
                        out_5kb.write(my_chain + "\n")
                    if seq_length >= 100000:
                        out_100kb.write(header)
                        out_100kb.write(my_chain + "\n")

    # Fallback: manual parsing
    else:
        with _open_file(input_file, "r") as in_f:
            with open(output_1kb, "w") as out_1kb, \
                 open(output_5kb, "w") as out_5kb, \
                 open(output_100kb, "w") as out_100kb:

                current_seq = []
                current_header = None

                def write_sequence():
                    """Helper to write accumulated sequence to appropriate files."""
                    if current_header and current_seq:
                        seq = "".join(current_seq).upper()
                        seq_length = len(seq)

                        if seq_length > 1000:
                            out_1kb.write(current_header)
                            out_1kb.write(seq + "\n")
                        if seq_length > 5000:
                            out_5kb.write(current_header)
                            out_5kb.write(seq + "\n")
                        if seq_length >= 100000:
                            out_100kb.write(current_header)
                            out_100kb.write(seq + "\n")

                for line in in_f:
                    if line.startswith(">"):
                        # Write previous sequence
                        write_sequence()

                        # Process new header
                        original_header = line[1:].strip()
                        seq_name, short_name, metadata, viral_id = _parse_fasta_header(original_header)

                        if mapping and seq_name in mapping:
                            new_name = mapping[seq_name]
                        elif mapping and short_name in mapping:
                            new_name = mapping[short_name]
                        else:
                            new_name = f"{prefix}{counter}" if not prefix.endswith("_") else f"{prefix}{counter}"
                            counter += 1

                        name_mapping[new_name] = (seq_name, short_name)

                        # Build header
                        if preserve_metadata and metadata:
                            if viral_id:
                                current_header = f">{new_name}|{viral_id}"
                                other_metadata = [m for m in metadata if m != viral_id]
                                if other_metadata:
                                    current_header += "|" + "|".join(other_metadata)
                                current_header += "\n"
                            else:
                                current_header = f">{new_name}|{'|'.join(metadata)}\n"
                        else:
                            current_header = f">{new_name}\n"

                        current_seq = []
                    else:
                        current_seq.append(line.strip())

                # Write last sequence
                write_sequence()

    return name_mapping


def restore_fasta(
    input_file: str,
    output_file: str,
    mapping: Dict[str, str],
    preserve_metadata: bool = False
) -> None:
    """
    Restore original sequence names in a FASTA file using a mapping.

    Args:
        input_file: Input FASTA file with renamed sequences
        output_file: Output FASTA file with original names
        mapping: Dictionary mapping current_name -> original_name
        preserve_metadata: Whether to preserve VirSorter metadata
    """
    with _open_file(input_file, "r") as in_f, open(output_file, "w") as out_f:
        for line in in_f:
            if line.startswith(">"):
                header = line[1:].strip()
                seq_name, short_name, metadata, _ = _parse_fasta_header(header)

                # Try to find mapping
                original = mapping.get(seq_name, mapping.get(short_name, None))
                if not original:
                    print(f"Warning: No mapping found for {seq_name}. Using current name.", file=sys.stderr)
                    original = seq_name

                # Write header with or without metadata
                if preserve_metadata and metadata:
                    out_f.write(f">{original}|{'|'.join(metadata)}\n")
                else:
                    out_f.write(f">{original}\n")
            else:
                out_f.write(line)


def rename_gff(
    input_file: str,
    output_file: str,
    mapping: Dict[str, str]
) -> None:
    """
    Rename contig IDs in a GFF file.

    Args:
        input_file: Input GFF file path
        output_file: Output GFF file path
        mapping: Dictionary mapping old_name -> new_name
    """
    with _open_file(input_file, "rt") as in_f, open(output_file, "w") as out_f:
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


def rename_genbank(
    input_file: str,
    output_file: str,
    mapping: Dict[str, str]
) -> None:
    """
    Rename contig IDs in a GenBank file.

    Args:
        input_file: Input GenBank file path
        output_file: Output GenBank file path
        mapping: Dictionary mapping old_name -> new_name
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


def write_mapping_file(
    mapping: Dict[str, Tuple[str, str]],
    output_file: str
) -> None:
    """
    Write a mapping file with three columns: original, renamed, short.

    Args:
        mapping: Dictionary mapping new_name -> (old_name, short_name)
        output_file: Output TSV file path
    """
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["original", "renamed", "short"])
        for new_name, (old_name, short_name) in mapping.items():
            writer.writerow([old_name, new_name, short_name])


def rename_mode(args):
    """Execute rename mode with multiple input files."""
    # Create output directory if it doesn't exist
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Categorize input files by type
    fasta_files = []
    gff_files = []
    gbk_files = []

    for input_file in args.input:
        if not os.path.exists(input_file):
            print(f"Error: File {input_file} does not exist", file=sys.stderr)
            sys.exit(1)

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
            print(f"Error: Unsupported file type: {input_file}", file=sys.stderr)
            print(f"Supported: {SUPPORTED_FASTA_EXTS + SUPPORTED_GFF_EXTS + SUPPORTED_GBK_EXTS}", file=sys.stderr)
            sys.exit(1)

    # Read existing mapping if provided
    old_to_new_mapping = None
    if args.use_mapping:
        if args.from_col and args.to_col:
            old_to_new_mapping = read_mapping_file(args.use_mapping, from_col=args.from_col, to_col=args.to_col)
        else:
            print('Mapping file was provided without --from-col/--to-col, default values would be used')

    # Combined mapping from all FASTA files
    combined_name_mapping = {}

    # First pass: Rename all FASTA files
    if fasta_files:
        print(f"Processing {len(fasta_files)} FASTA file(s)...")

        # Handle size-based separation
        if args.separate_by_size:
            if len(fasta_files) > 1:
                print("Error: --separate-by-size only supports a single FASTA input file", file=sys.stderr)
                sys.exit(1)

            fasta_file = fasta_files[0]
            input_path = Path(fasta_file)
            # Use the stem (filename without extension) as the output prefix
            output_prefix = str(outdir / input_path.stem)

            print(f"  Renaming FASTA with size separation: {fasta_file}")
            name_mapping = rename_fasta_by_size(
                fasta_file,
                output_prefix,
                prefix=args.prefix,
                mapping=old_to_new_mapping,
                preserve_metadata=args.preserve_metadata
            )

            combined_name_mapping.update(name_mapping)
            print(f"  → {output_prefix}_1kb_contigs.fasta")
            print(f"  → {output_prefix}_5kb_contigs.fasta")
            print(f"  → {output_prefix}_100kb_contigs.fasta")
            print(f"  Total sequences renamed: {len(name_mapping)}")

        # Normal renaming (single output file)
        else:
            for fasta_file in fasta_files:
                input_path = Path(fasta_file)
                output_file = outdir / input_path.name

                print(f"  Renaming FASTA: {fasta_file}")
                name_mapping = rename_fasta(
                    fasta_file,
                    str(output_file),
                    prefix=args.prefix,
                    mapping=old_to_new_mapping,
                    preserve_metadata=args.preserve_metadata,
                    use_pyfastx=args.use_pyfastx
                )

                # Merge mappings
                combined_name_mapping.update(name_mapping)
                print(f"  → {output_file} ({len(name_mapping)} sequences)")

        # Convert combined mapping to old -> new format for GFF/GenBank
        if not old_to_new_mapping:
            old_to_new_mapping = {old: new for new, (old, short) in combined_name_mapping.items()}
            # Also add short names as keys
            for new, (old, short) in combined_name_mapping.items():
                if short != old:
                    old_to_new_mapping[short] = new

        # Write combined mapping file
        if args.map and combined_name_mapping:
            write_mapping_file(combined_name_mapping, args.map)
            print(f"Wrote mapping to: {args.map}")

    # Second pass: Rename GFF files
    if gff_files:
        if not old_to_new_mapping:
            print("Error: GFF renaming requires either FASTA file(s) in input or --use-mapping", file=sys.stderr)
            sys.exit(1)

        print(f"\nProcessing {len(gff_files)} GFF file(s)...")
        for gff_file in gff_files:
            input_path = Path(gff_file)
            output_file = outdir / input_path.name

            print(f"  Renaming GFF: {gff_file}")
            rename_gff(gff_file, str(output_file), old_to_new_mapping)
            print(f"  → {output_file}")

    # Third pass: Rename GenBank files
    if gbk_files:
        if not old_to_new_mapping:
            print("Error: GenBank renaming requires either FASTA file(s) in input or --use-mapping", file=sys.stderr)
            sys.exit(1)

        print(f"\nProcessing {len(gbk_files)} GenBank file(s)...")
        for gbk_file in gbk_files:
            input_path = Path(gbk_file)
            output_file = outdir / input_path.name

            print(f"  Renaming GenBank: {gbk_file}")
            rename_genbank(gbk_file, str(output_file), old_to_new_mapping)
            print(f"  → {output_file}")

    print(f"\nAll files processed successfully. Output directory: {outdir}")


def restore_mode(args):
    """Execute restore mode with multiple input files."""
    if not args.map:
        print("Error: Restore mode requires --map", file=sys.stderr)
        sys.exit(1)

    # Create output directory if it doesn't exist
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Read mapping file
    mapping = read_mapping_file(args.map, from_col=args.from_col, to_col=args.to_col)

    print(f"Restoring {len(args.input)} file(s) using mapping from {args.map}")
    print(f"Mapping: {args.from_col} -> {args.to_col}\n")

    # Process each input file
    for input_file in args.input:
        # Check if input file exists
        if not os.path.exists(input_file):
            print(f"Error: File {input_file} does not exist", file=sys.stderr)
            sys.exit(1)

        # Currently only supports FASTA restoration
        input_path = Path(input_file)
        file_ext = "".join(input_path.suffixes).lower()

        is_fasta = any(file_ext.endswith(ext) for ext in SUPPORTED_FASTA_EXTS)

        if not is_fasta:
            print(f"Warning: Restore mode currently only supports FASTA files, skipping {input_file}", file=sys.stderr)
            continue

        output_file = outdir / input_path.name

        print(f"  Restoring: {input_file}")
        restore_fasta(input_file, str(output_file), mapping, preserve_metadata=args.preserve_metadata)
        print(f"  → {output_file}")

    print(f"\nAll files restored successfully. Output directory: {outdir}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Universal contig renaming tool for FASTA, GFF, and GenBank files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Rename FASTA with custom prefix
  %(prog)s rename -i input.fasta -o output_dir -p "contig_" -m mapping.tsv

  # Rename FASTAs with custom prefix
  %(prog)s rename -i input1.fasta input2.fasta -o output_dir -p "contig_" -m mapping.tsv

  # Rename FASTA and GFF together
  %(prog)s rename -i input.fasta input.gff -o output_dir -p "seq" -m mapping.tsv

  # Rename multiple files (FASTA, GFF, and GenBank)
  %(prog)s rename -i input.fasta input.gff input.gbk -o output_dir -p "contig_"

  # Rename with size-based separation (creates _1kb, _5kb, _100kb output files)
  %(prog)s rename -i input.fasta -o output_dir -p "contig_" --separate-by-size -m mapping.tsv

  # Restore original names from temporary names
  %(prog)s restore -i renamed.fasta -o output_dir -m mapping.tsv --from-col renamed --to-col original

  # Restore multiple files
  %(prog)s restore -i renamed1.fasta renamed2.fasta -o output_dir -m mapping.tsv
        """
    )

    subparsers = parser.add_subparsers(dest="mode", help="Operation mode")

    # Rename mode
    rename_parser = subparsers.add_parser("rename", help="Rename contigs in files")
    rename_parser.add_argument("-i", "--input", nargs="+", required=True, help="Input file(s) (FASTA, GFF, and/or GenBank)")
    rename_parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    rename_parser.add_argument("-p", "--prefix", default="contig", help="Prefix for new names (default: contig)")
    rename_parser.add_argument("-m", "--map", help="Output mapping file (TSV)")
    rename_parser.add_argument("--use-mapping", help="Use existing mapping file instead of generating new names")
    rename_parser.add_argument("--from-col", default="original", help="Source column in mapping file (use if existing mapping was provided)")
    rename_parser.add_argument("--to-col", default="renamed", help="Target column in mapping file (use if existing mapping was provided)")
    rename_parser.add_argument("--preserve-metadata", action="store_true", help="Preserve metadata in headers (VirSorter, viral IDs)")
    rename_parser.add_argument("--use-pyfastx", action="store_true", help="Use pyfastx library for faster FASTA processing")
    rename_parser.add_argument("--separate-by-size", action="store_true", help="Separate output by contig size (creates _1kb_contigs.fasta, _5kb_contigs.fasta, _100kb_contigs.fasta)")
    rename_parser.set_defaults(func=rename_mode)

    # Restore mode
    restore_parser = subparsers.add_parser("restore", help="Restore original contig names")
    restore_parser.add_argument("-i", "--input", nargs="+", required=True, help="Input FASTA file(s) with renamed sequences")
    restore_parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    restore_parser.add_argument("-m", "--map", required=True, help="Mapping file (TSV)")
    restore_parser.add_argument("--from-col", default="renamed", help="Source column in mapping file")
    restore_parser.add_argument("--to-col", default="original", help="Target column in mapping file")
    restore_parser.add_argument("--preserve-metadata", action="store_true", help="Preserve metadata in headers")
    restore_parser.set_defaults(func=restore_mode)

    args = parser.parse_args()

    if not args.mode:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
