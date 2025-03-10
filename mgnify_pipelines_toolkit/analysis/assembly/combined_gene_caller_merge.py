#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024-2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import json
import logging
from collections import defaultdict
import csv
import re

from intervaltree import Interval, IntervalTree
from Bio import SeqIO


MASK_OVERLAP_THRESHOLD = 5


def parse_gff(gff_file):
    """
    Parse a GFF file and extract CDS features as Interval objects.

    Args:
        gff_file (str): Path to the GFF file.

    Returns:
        dict: A nested dictionary with sequence IDs as keys, and within each,
            strand (+/-) as keys, containing a list of Intervals for CDS regions.
            Each Interval object apart from the start and end positions of the CDS region
            also stores the protein ID.
    """
    predictions = defaultdict(lambda: defaultdict(list))
    with open(gff_file, "r") as gff_in:
        for line in gff_in:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            seq_id, _, feature_type, start, end, _, strand, _, attributes = fields
            if feature_type == "CDS":
                # Parse attributes to get the ID value
                attr_dict = dict(
                    attr.split("=") for attr in attributes.split(";") if "=" in attr
                )
                protein_id = attr_dict["ID"]
                predictions[seq_id][strand].append(
                    Interval(int(start), int(end), data={"protein_id": protein_id})
                )
    if not predictions:
        raise ValueError("Zero gene predictions was read from the GFF file")
    return predictions


def parse_pyrodigal_output(file):
    """
    Parse Pyrodigal *.out file to extract gene predictions as Interval objects.
    Example of *.out file:
    # Sequence Data: seqnum=1;seqlen=25479;seqhdr="Bifidobacterium-longum-subsp-infantis-MC2-contig1"
    # Model Data: version=Pyrodigal.v2.6.3;run_type=Single;model="Ab initio";gc_cont=59.94;transl_table=11;uses_sd=1
    >1_1_279_+

    Args:
        file (str): Path to the Pyrodigal *.out file.

    Returns:
        dict: A nested dictionary with sequence IDs as keys, and within each,
            strand (+/-) as keys, containing a list of Intervals for CDS regions.
            Each Interval object apart from the start and end positions of the CDS region
            also stores the protein ID.
    """
    predictions = defaultdict(lambda: defaultdict(list))
    with open(file) as file_in:
        for line in file_in:
            if line.startswith("# Model Data"):
                continue
            if line.startswith("# Sequence Data"):
                matches = re.search(r'seqhdr="(\S+)"', line)
                if matches:
                    seq_id = matches.group(1)
            else:
                fields = line[1:].strip().split("_")
                # Fragment_id is an index of the fragment
                # Pyrodigal uses these (rather than coordinates) to identify sequences in the fasta output
                fragment_id, start, end, strand = fields
                protein_id = f"{seq_id}_{fragment_id}"
                predictions[seq_id][strand].append(
                    Interval(int(start), int(end), data={"protein_id": protein_id})
                )
    if not predictions:
        raise ValueError("Zero gene predictions was read from the *.out file")
    return predictions


def parse_fgsrs_output(file):
    """
    Parse FragGeneScanRS *.out file to extract gene predictions as Interval objects.
    Example of *.out file:
    >Bifidobacterium-longum-subsp-infantis-MC2-contig1
    256	2133	-	1	1.263995	I:	D:

    Args:
        file (str): Path to the FragGeneScanRS *.out file.

    Returns:
        dict: A nested dictionary with sequence IDs as keys, and within each,
            strand (+/-) as keys, containing a list of Intervals for CDS regions.
            Each Interval object apart from the start and end positions of the CDS region
            also stores the protein ID.
    """
    predictions = defaultdict(lambda: defaultdict(list))
    with open(file) as file_in:
        for line in file_in:
            if line.startswith(">"):
                seq_id = line.split()[0][1:]
            else:
                fields = line.strip().split("\t")
                start, end, strand, *_ = fields
                protein_id = f"{seq_id}_{start}_{end}_{strand}"
                predictions[seq_id][strand].append(
                    Interval(int(start), int(end), data={"protein_id": protein_id})
                )
    if not predictions:
        raise ValueError("Zero gene predictions was read from the *.out file")
    return predictions


def parse_cmsearch_output(mask_file):
    """
    Parse masking regions from a cmsearch output file and store them as Intervals.

    Args:
        mask_file (str): Path to the masking file (possibly BED or GFF-like format).

    Returns:
        dict: A dictionary with sequence IDs as keys, and a list of Intervals representing masked regions.
    """
    regions = defaultdict(list)
    with open(mask_file) as file_in:
        for line in file_in:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split()
            seq_id = fields[0]
            start = int(fields[7])
            end = int(fields[8])
            if start > end:
                start, end = end, start
            regions[seq_id].append(Interval(start, end))
    if not regions:
        raise ValueError("Zero intervals was read from the input masking file")
    return regions


def mask_regions(predictions, mask):
    """
    Apply masking to predictions by removing regions that overlap significantly
    (more than MASK_OVERLAP_THRESHOLD)
    with masked regions.

    Args:
        predictions (dict): A nested dictionary with sequence IDs as keys, and within each,
            strand (+/-) as keys, containing a list of Intervals as values.
        mask (dict): A dictionary with sequence IDs as keys, and a list of Intervals as values.

    Returns:
        dict: Updated predictions with masked regions removed.
    """
    masked = defaultdict(lambda: defaultdict(list))

    for seq_id, strand_dict in predictions.items():
        if seq_id in mask:
            mask_tree = create_interval_tree(mask[seq_id])
            for strand, regions in strand_dict.items():
                tree = create_interval_tree(regions)
                masked_intervals = []
                for region in tree:
                    # Check for overlaps greater than 5 base pairs
                    overlapping_intervals = mask_tree.overlap(region.begin, region.end)
                    overlap = False
                    for mask_region in overlapping_intervals:
                        # If overlap is more than 5 base pairs, mark for masking
                        # Add 1 to make boundaries inclusive
                        overlap_len = 1 + abs(
                            min(region.end, mask_region.end)
                            - max(region.begin, mask_region.begin)
                        )
                        if overlap_len > MASK_OVERLAP_THRESHOLD:
                            overlap = True
                            break
                    if not overlap:
                        masked_intervals.append(region)
                masked[seq_id][strand] = sorted(masked_intervals)
        else:
            # If no mask information exists, add the predictions directly
            masked[seq_id] = strand_dict
    return masked


def merge_predictions(predictions, priority):
    """
    Merge gene predictions from two sources, applying a priority order.

    Args:
        predictions (dict): Nested dictionary containing gene predictions from both sources.
        priority (list): List specifying the order of priority for merging the predictions.

    Returns:
        dict: Nested dictionary with all predictions of the first priority source merged with non-overlapping predictions
            the secondary source.
    """
    merged = defaultdict(lambda: defaultdict((lambda: defaultdict(list))))
    primary, secondary = priority

    # Primary merge
    merged[primary] = predictions[primary]

    # Secondary merge: add non-overlapping regions from the secondary gene caller
    for seq_id in predictions[secondary]:
        for strand in ["+", "-"]:
            secondary_regions = predictions[secondary][seq_id][strand]
            if seq_id in predictions[primary]:
                primary_regions = merged[primary][seq_id][strand]
                merged[secondary][seq_id][strand].extend(
                    check_against_gaps(primary_regions, secondary_regions)
                )
            else:
                merged[secondary][seq_id][strand] = secondary_regions
    return merged


def check_against_gaps(regions, candidates):
    """
    Check candidate regions against existing regions and select those
    that do not overlap with any existing ones.

    Args:
        regions (list): Interval objects for existing regions.
        candidates (list): Interval objects for candidate regions.

    Returns:
        list: Selected candidate Intervals that do not overlap with existing ones.
    """
    regions_tree = create_interval_tree(regions)
    selected_candidates = []
    for candidate in candidates:
        # Check if the candidate overlaps with any existing region
        if not regions_tree.overlap(candidate.begin, candidate.end):
            selected_candidates.append(candidate)
    return selected_candidates


def output_fasta_files(predictions, files_dict, output_faa, output_ffn):
    """
    Write FASTA output files containing protein and transcript sequences for
    the predicted genes after merging.

    Args:
        predictions (dict): Nested dictionary with merged gene predictions as Interval objects.
            Each Interval object stores a protein ID in the data attribute.
        files_dict (dict): Dictionary containing input FASTA files for both Pyrodigal and FragGeneScanRS.
        output_faa (str): Path to output protein FASTA file.
        output_ffn (str): Path to output transcript FASTA file.
    """
    with (
        open(output_faa, "w") as output_faa_fh,
        open(output_ffn, "w") as output_ffn_fh,
    ):
        for caller, seq_data in predictions.items():
            proteins = set()
            for seq_id, strand_dict in seq_data.items():
                for strand, regions in strand_dict.items():
                    for region in regions:
                        protein_id = region.data["protein_id"]
                        proteins.add(protein_id)

            for input_file, output_file in [
                (files_dict[caller]["proteins"], output_faa_fh),
                (files_dict[caller]["transcripts"], output_ffn_fh),
            ]:
                sequences = []
                for record in SeqIO.parse(input_file, "fasta"):
                    if record.id in proteins:
                        # Prodigal appends * to the end of a truncated sequence
                        # FGS uses * to mark an ambiguous amino acid
                        # Replace ending * and replace any other "*" with "X"
                        record.seq = record.seq.rstrip("*").replace("*", "X")
                        sequences.append(record)
                SeqIO.write(sequences, output_file, "fasta")


def output_gff(predictions, output_gff):
    """
    Write merged gene predictions to a GFF output file.

    Args:
        predictions (dict): Nested dictionary with merged gene predictions as Interval objects.
            Each Interval object stores a protein ID in the data attribute.
        output_gff (str): Path to the output GFF file.
    """
    with open(output_gff, "w") as gff_out:
        writer = csv.writer(gff_out, delimiter="\t")
        gff_out.write("##gff-version 3\n")
        for caller, seq_data in predictions.items():
            for seq_id, strand_dict in seq_data.items():
                for strand, regions in strand_dict.items():
                    for region in regions:
                        writer.writerow(
                            [
                                seq_id,  # Sequence ID
                                caller,  # Source
                                "CDS",  # Feature type
                                region.begin,  # Start position
                                region.end,  # End position
                                ".",  # Score (not used, hence '.')
                                strand,  # Strand (+/-)
                                ".",  # Phase (not used, hence '.')
                                f"ID={region.data['protein_id']}",  # Attributes
                            ]
                        )


def output_summary(summary, output_file):
    """
    Write a summary of gene counts to a text file in JSON format.

    Args:
        summary (dict): Summary of gene counts.
        output_file (str): Path to the summary output file.
    """
    with open(output_file, "w") as sf:
        sf.write(json.dumps(summary, sort_keys=True, indent=4) + "\n")


def get_counts(predictions):
    """
    Count the number of gene predictions for each caller.

    Args:
        predictions (dict): Nested dictionary with gene predictions for each caller.

    Returns:
        dict: Total count of genes for each caller.
    """
    total = {}
    for caller, seq_data in predictions.items():
        count = sum(
            len(seq_data[seq_id]["+"] + seq_data[seq_id]["-"]) for seq_id in seq_data
        )
        total[caller] = count
    return total


def create_interval_tree(regions):
    """
    Create an IntervalTree from a list of regions.

    Args:
        regions (list): List of Interval objects.

    Returns:
        IntervalTree: An interval tree for efficient overlap checking.
    """
    tree = IntervalTree()
    for region in regions:
        tree.add(region)
    return tree


def main():
    parser = argparse.ArgumentParser(
        """
        MGnify gene caller combiner.
        This script merges gene predictions made by Pyrodigal and FragGeneScanRS (FGS)
        and outputs FASTA and GFF files.
        For each gene caller, the script expects a set of files:
        - GFF file with gene predictions OR *.out file
        - FASTA file with protein sequences
        - FASTA file with transcript sequences
        """
    )
    parser.add_argument(
        "--name", "-n", required=True, help="Base name for output files"
    )
    parser.add_argument(
        "--priority",
        "-P",
        choices=["Pyrodigal_FragGeneScanRS", "FragGeneScanRS_Pyrodigal"],
        default="Pyrodigal_FragGeneScanRS",
        help="Merge priority",
    )
    parser.add_argument(
        "--mask",
        "-m",
        help="Regions for masking (Infernal cmsearch output file)",
    )
    parser.add_argument("--pyrodigal-gff", "-pg", help="Pyrodigal *.gff file")
    parser.add_argument("--pyrodigal-out", "-po", help="Pyrodigal *.out file")
    parser.add_argument(
        "--pyrodigal-ffn",
        "-pt",
        required=True,
        help="Pyrodigal *.ffn file with transcripts",
    )
    parser.add_argument(
        "--pyrodigal-faa",
        "-pp",
        required=True,
        help="Pyrodigal *.faa file with proteins",
    )
    parser.add_argument("--fgsrs-gff", "-fg", help="FragGeneScanRS *.gff file")
    parser.add_argument("--fgsrs-out", "-fo", help="FragGeneScanRS *.out file")
    parser.add_argument(
        "--fgsrs-ffn",
        "-ft",
        required=True,
        help="FragGeneScanRS *.ffn file with transcripts",
    )
    parser.add_argument(
        "--fgsrs-faa",
        "-fp",
        required=True,
        help="FragGeneScanRS *.faa file with proteins",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Increase verbosity level to debug"
    )
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(levelname)s %(asctime)s - %(message)s",
        datefmt="%Y/%m/%d %H:%M:%S",
    )

    if not args.pyrodigal_out and not args.pyrodigal_gff:
        parser.error(
            "For Pyrodigal, you must provide either --pyrodigal-out or --pyrodigal-gff"
        )

    if not args.fgsrs_out and not args.fgsrs_gff:
        parser.error(
            "For FragGeneScanRS, you must provide either --fgsrs-out or --fgsrs-gff"
        )

    summary = {}
    all_predictions = {}

    caller_priority = args.priority.split("_")
    logging.info(f"Caller priority: 1. {caller_priority[0]}, 2. {caller_priority[1]}")

    logging.info("Parsing Pyrodigal annotations...")
    if args.pyrodigal_out:
        all_predictions["Pyrodigal"] = parse_pyrodigal_output(args.pyrodigal_out)
    elif args.pyrodigal_gff:
        all_predictions["Pyrodigal"] = parse_gff(args.pyrodigal_gff)

    logging.info("Parsing FragGeneScanRS annotations...")
    if args.fgsrs_out:
        all_predictions["FragGeneScanRS"] = parse_fgsrs_output(args.fgsrs_out)
    elif args.fgsrs_gff:
        all_predictions["FragGeneScanRS"] = parse_gff(args.fgsrs_gff)

    summary["all"] = get_counts(all_predictions)

    if args.mask:
        logging.info("Masking of non-coding RNA regions was enabled")
        logging.info(f"Parsing masking intervals from file {args.mask}")
        mask_regions_file = parse_cmsearch_output(args.mask)
        for caller in all_predictions:
            logging.info(f"Masking {caller} outputs...")
            all_predictions[caller] = mask_regions(
                all_predictions[caller], mask_regions_file
            )
        summary["after_masking"] = get_counts(all_predictions)

    logging.info("Merging combined gene caller results")
    merged_predictions = merge_predictions(all_predictions, caller_priority)
    summary["merged"] = get_counts(merged_predictions)

    logging.info("Writing output files...")
    output_summary(summary, f"{args.name}.summary.txt")
    output_gff(merged_predictions, f"{args.name}.gff")
    files = {
        "Pyrodigal": {
            "proteins": args.pyrodigal_faa,
            "transcripts": args.pyrodigal_ffn,
        },
        "FragGeneScanRS": {"proteins": args.fgsrs_faa, "transcripts": args.fgsrs_ffn},
    }
    output_fasta_files(
        merged_predictions,
        files,
        f"{args.name}.faa",
        f"{args.name}.ffn",
    )


if __name__ == "__main__":
    main()
