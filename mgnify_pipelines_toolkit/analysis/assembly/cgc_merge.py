#!/usr/bin/env python3

import argparse
import logging
from collections import defaultdict

from intervaltree import Interval, IntervalTree


def create_interval_tree(regions):
    tree = IntervalTree()
    for region in regions:
        tree.add(region)  # 'region' is already an Interval
    return tree


# Flatten regions (merge overlapping intervals)
def flatten_regions(regions):
    """Take a list of intervals, merge overlapping intervals and return modified list"""
    tree = create_interval_tree(regions)
    tree.merge_overlaps()
    return sorted(tree)


def check_against_gaps(regions, candidates):
    # TODO there is no check if regions is empty, is it a problem?
    regions_tree = create_interval_tree(regions)
    selected_candidates = []
    for candidate in candidates:
        # Check if the candidate overlaps with any existing region
        if not regions_tree.overlap(candidate.begin, candidate.end):
            selected_candidates.append(candidate)
    return selected_candidates


def output_prodigal(predictions, files, outputs):
    """From the combined predictions output the prodigal data"""
    pass


def output_fgs(predictions, files, outputs):
    """From the combined predictions output the FGS data"""
    pass


def output_files(predictions, summary, files):
    """Output all files"""
    pass


# Parsing GFF files into a structure containing IntervalTrees
def parse_gff(gff_file):
    predictions = defaultdict(lambda: defaultdict(list))
    with open(gff_file, "r") as gff_in:
        for line in gff_in:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            seq_id, _, feature_type, start, end, _, strand, _, _ = fields
            if feature_type == "gene":
                predictions[seq_id][strand].append(Interval(int(start), int(end)))
    return predictions


def parse_cmsearch_output(mask_file):
    regions = defaultdict(list)
    with open(mask_file) as file_in:
        for line in file_in:
            if line.startswith("#"):
                continue
            # TODO maybe it's TSV?
            fields = line.rstrip().split()
            seq_id = fields[0]
            start = int(fields[7])
            end = int(fields[8])
            if start > end:
                start, end = end, start
            regions[seq_id].append(Interval(start, end))
    return regions


def mask_regions(predictions, mask):
    masked = defaultdict(lambda: defaultdict(list))

    # Convert mask file into an IntervalTree for each sequence
    mask_tree = create_interval_tree(mask)

    for seq_id, strand_dict in predictions.items():
        for strand, regions in strand_dict.items():
            tree = create_interval_tree(regions)
            # Subtract mask intervals from predicted regions
            masked_intervals = tree - mask_tree
            masked[seq_id][strand] = sorted(masked_intervals)

    return masked


def merge_predictions(predictions, priority):
    merged = {}
    primary, secondary = priority

    # Primary merge
    for seq_id in predictions[primary]:
        merged[seq_id] = defaultdict(list)
        for strand in ["+", "-"]:
            merged[seq_id][strand] = flatten_regions(
                predictions[primary][seq_id][strand]
            )
    # Secondary merge: add non-overlapping regions from the secondary gene caller
    for seq_id in predictions[secondary]:
        if seq_id not in merged:
            merged[seq_id] = defaultdict(list)
        for strand in ["+", "-"]:
            if seq_id in predictions[primary]:
                primary_regions = merged[seq_id][strand]
                secondary_regions = predictions[secondary][seq_id][strand]
                merged[seq_id][strand].extend(
                    check_against_gaps(primary_regions, secondary_regions)
                )
            else:
                merged[seq_id][strand] = predictions[secondary][seq_id][strand]
    return merged


def get_counts(predictions):
    total = {}
    for caller, seq_data in predictions.items():
        count = sum(
            len(seq_data[seq_id]["+"] + seq_data[seq_id]["-"]) for seq_id in seq_data
        )
        total[caller] = count
    return total


def combine_main():
    parser = argparse.ArgumentParser(
        "MGnify gene caller combiner. This script will merge the gene called by prodigal and fraggenescan (in any order)"
    )
    parser.add_argument(
        "-n", "--name", action="store", dest="name", required=True, help="basename"
    )
    parser.add_argument(
        "-k",
        "--mask",
        action="store",
        dest="mask",
        required=False,
        help="Sequence mask file",
    )

    parser.add_argument(
        "-a",
        "--prodigal-out",
        action="store",
        dest="prodigal_out",
        required=False,
        help="Stats out prodigal",
    )
    parser.add_argument(
        "-b",
        "--prodigal-ffn",
        action="store",
        dest="prodigal_ffn",
        required=False,
        help="Stats ffn prodigal",
    )
    parser.add_argument(
        "-c",
        "--prodigal-faa",
        action="store",
        dest="prodigal_faa",
        required=False,
        help="Stats faa prodigal",
    )

    parser.add_argument(
        "-d",
        "--fgs-out",
        action="store",
        dest="fgs_out",
        required=False,
        help="Stats out FGS",
    )
    parser.add_argument(
        "-e",
        "--fgs-ffn",
        action="store",
        dest="fgs_ffn",
        required=False,
        help="Stats ffn FGS",
    )
    parser.add_argument(
        "-f",
        "--fgs-faa",
        action="store",
        dest="fgs_faa",
        required=False,
        help="Stats faa FGS",
    )

    parser.add_argument(
        "-p",
        "--caller-priority",
        action="store",
        dest="caller_priority",
        required=False,
        choices=["prodigal_fgs", "fgs_prodigal"],
        default="prodigal_fgs",
        help="Caller priority.",
    )

    parser.add_argument(
        "-v", "--verbose", action="count", help="Increase verbosity level"
    )

    args = parser.parse_args()

    # Set up logging system
    verbose_mode = args.verbose or 0

    log_level = logging.WARNING
    if verbose_mode:
        log_level = logging.DEBUG if verbose_mode > 1 else logging.INFO

    logging.basicConfig(
        level=log_level,
        format="%(levelname)s %(asctime)s - %(message)s",
        datefmt="%Y/%m/%d %I:%M:%S %p",
    )

    summary = {}
    all_predictions = {}
    files = {}
    caller_priority = []
    if args.caller_priority:
        caller_priority = args.caller_priority.split("_")
    else:
        caller_priority = ["prodigal", "fgs"]

    logging.info(f"Caller priority: 1. {caller_priority[0]}, 2. {caller_priority[1]}")

    if args.prodigal_out:
        logging.info("Prodigal presented")
        logging.info("Getting Prodigal regions...")
        all_predictions["prodigal"] = parse_gff(args.prodigal_out)

        files["prodigal"] = [args.prodigal_out, args.prodigal_ffn, args.prodigal_faa]

    if args.fgs_out:
        logging.info("FGS presented")
        logging.info("Getting FragGeneScan regions ...")
        all_predictions["fgs"] = parse_gff(args.fgs_out)

        files["fgs"] = [args.fgs_out, args.fgs_ffn, args.fgs_faa]

    summary["all"] = get_counts(all_predictions)

    # Apply mask of ncRNA search
    logging.info("Masking non coding RNA regions...")
    if args.mask:
        logging.info("Reading regions for masking...")
        mask = parse_cmsearch_output(args.mask)
        if "prodigal" in all_predictions:
            logging.info("Masking Prodigal outputs...")
            all_predictions["prodigal"] = mask_regions(
                all_predictions["prodigal"], mask
            )
        if "fgs" in all_predictions:
            logging.info("Masking FragGeneScan outputs...")
            all_predictions["fgs"] = mask_regions(all_predictions["fgs"], mask)
        summary["masked"] = get_counts(all_predictions)

    # Run the merging step
    if len(all_predictions) > 1:
        logging.info("Merging combined gene caller results...")
        merged_predictions = merge_predictions(all_predictions, caller_priority)
    else:
        logging.info("Skipping merging step...")
        merged_predictions = all_predictions
    summary["merged"] = get_counts(merged_predictions)

    # Output fasta files and summary (json)
    logging.info("Writing output files...")

    files["merged"] = [args.name + ext for ext in [".out", ".ffn", ".faa"]]

    output_files(merged_predictions, summary, files)


if __name__ == "__main__":
    combine_main()
