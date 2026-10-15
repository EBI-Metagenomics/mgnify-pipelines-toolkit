#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path


def generate_reads_manifest(
    output_path: Path,
    study: str,
    sample: str,
    name: str,
    platform: str,
    instrument: str,
    library_source: str,
    library_selection: str,
    library_strategy: str,
    fastq_files: list[str] | None = None,
    bam_file: str | None = None,
    insert_size: int | None = None,
    library_name: str | None = None,
    description: str | None = None,
) -> None:
    """Write a tab-delimited ENA reads submission manifest file.

    Reads are declared either as FASTQ file(s) or as a single BAM file (e.g. long reads), not both.

    Appends a timestamp to `name` so re-submitting the same run data
    produces a distinct NAME field and avoids ENA duplicate-rejection errors.
    """
    if bool(fastq_files) == bool(bam_file):
        raise ValueError("Exactly one of fastq_files or bam_file must be provided")

    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    run_name = f"{name}_{timestamp}"

    fields: list[tuple[str, str]] = [
        ("STUDY", study),
        ("SAMPLE", sample),
        ("NAME", run_name),
        ("PLATFORM", platform),
        ("INSTRUMENT", instrument),
        ("LIBRARY_SOURCE", library_source),
        ("LIBRARY_SELECTION", library_selection),
        ("LIBRARY_STRATEGY", library_strategy),
    ]

    if insert_size is not None:
        fields.append(("INSERT_SIZE", str(insert_size)))
    if library_name:
        fields.append(("LIBRARY_NAME", library_name))
    if description:
        fields.append(("DESCRIPTION", description))

    if bam_file:
        fields.append(("BAM", bam_file))
    else:
        for fq in fastq_files:
            fields.append(("FASTQ", fq))

    with output_path.open("w") as fh:
        for key, value in fields:
            fh.write(f"{key}\t{value}\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate an ENA reads submission manifest file.")
    parser.add_argument("--study", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--name", required=True, help="Run name base; a YYYYmmddHHMMSS timestamp is appended")
    parser.add_argument("--platform", required=True)
    parser.add_argument("--instrument", required=True)
    parser.add_argument("--library-source", required=True)
    parser.add_argument("--library-selection", required=True)
    parser.add_argument("--library-strategy", required=True)
    reads_group = parser.add_mutually_exclusive_group(required=True)
    reads_group.add_argument("--fastq", action="append", dest="fastq_files", help="FASTQ file (repeat for paired-end)")
    reads_group.add_argument("--bam", dest="bam_file", help="BAM file (e.g. unaligned or basecalled long reads)")
    parser.add_argument("--output", required=True, type=Path, help="Path for the output manifest file")
    parser.add_argument("--insert-size", type=int, default=None)
    parser.add_argument("--library-name", default=None)
    parser.add_argument("--description", default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    generate_reads_manifest(
        output_path=args.output,
        study=args.study,
        sample=args.sample,
        name=args.name,
        platform=args.platform,
        instrument=args.instrument,
        library_source=args.library_source,
        library_selection=args.library_selection,
        library_strategy=args.library_strategy,
        fastq_files=args.fastq_files,
        bam_file=args.bam_file,
        insert_size=args.insert_size,
        library_name=args.library_name,
        description=args.description,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
