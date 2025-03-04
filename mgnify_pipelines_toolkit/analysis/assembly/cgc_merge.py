#!/usr/bin/env python3

import argparse
import json
import logging
from collections import defaultdict
import re

from intervaltree import Interval, IntervalTree
from Bio import SeqIO


MASK_OVERLAP_THRESHOLD = 5


def parse_gff(gff_file):
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


def parse_prodigal_output(file):
    """
    Parse Prodigal *.out file.
    Example:
    # Sequence Data: seqnum=1;seqlen=25479;seqhdr="Bifidobacterium-longum-subsp-infantis-MC2-contig1"
    # Model Data: version=Prodigal.v2.6.3;run_type=Single;model="Ab initio";gc_cont=59.94;transl_table=11;uses_sd=1
    >1_1_279_+
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
                # Prodigal uses these (rather than coordinates) to identify sequences in the fasta output
                fragment_id, start, end, strand = fields
                protein_id = f"{seq_id}_{fragment_id}"
                predictions[seq_id][strand].append(
                    Interval(int(start), int(end), data={"protein_id": protein_id})
                )
    if not predictions:
        raise ValueError("Zero gene predictions was read from the *.out file")
    return predictions


def parse_fgs_output(file):
    """
    Parse FGS *.out file.
    Example:
    >Bifidobacterium-longum-subsp-infantis-MC2-contig1
    256	2133	-	1	1.263995	I:	D:
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
    if not regions:
        raise ValueError("Zero intervals was read from the input masking file")
    return regions


def mask_regions(predictions, mask):
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
                        # Add 1 to make bondaries inclusive
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
    return masked


def merge_predictions(predictions, priority):
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
    # TODO there is no check if region is empty, is it a problem?
    regions_tree = create_interval_tree(regions)
    selected_candidates = []
    for candidate in candidates:
        # Check if the candidate overlaps with any existing region
        if not regions_tree.overlap(candidate.begin, candidate.end):
            selected_candidates.append(candidate)
    return selected_candidates


def output_fasta_files(predictions, files_dict, output_faa, output_ffn):
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
                        record.seq = record.seq.rstrip("*")
                        sequences.append(record)
                SeqIO.write(sequences, output_file, "fasta")


def output_gff(predictions, output_gff):
    with open(output_gff, "w") as gff_out:
        for caller, seq_data in predictions.items():
            for seq_id, strand_dict in seq_data.items():
                for strand, regions in strand_dict.items():
                    for region in regions:
                        gff_out.write(
                            f"{seq_id}\tMGnify\tgene\t{region.begin}\t{region.end}\t.\t{strand}\t.\tID=gene_{seq_id}_{region.begin}_{region.end}\n"
                        )


def output_summary(summary, output_file):
    with open(output_file, "w") as sf:
        sf.write(json.dumps(summary, sort_keys=True, indent=4) + "\n")


def get_counts(predictions):
    total = {}
    for caller, seq_data in predictions.items():
        count = sum(
            len(seq_data[seq_id]["+"] + seq_data[seq_id]["-"]) for seq_id in seq_data
        )
        total[caller] = count
    return total


def create_interval_tree(regions):
    tree = IntervalTree()
    for region in regions:
        tree.add(region)
    return tree


def main():
    parser = argparse.ArgumentParser(
        """
        MGnify gene caller combiner.
        This script merges gene predictions made by Prodigal and FragGeneScan
        and outputs FASTA and GFF files.
        """
    )
    parser.add_argument(
        "--name", "-n", required=True, help="Base name for output files"
    )
    parser.add_argument(
        "--priority",
        "-P",
        choices=["prodigal_fgs", "fgs_prodigal"],
        default="prodigal_fgs",
        help="Merge priority",
    )
    parser.add_argument(
        "--mask",
        "-m",
        help="Masked regions (in GFF or BED format)",  # TODO why GFF or BED?
    )
    parser.add_argument("--prodigal-gff", "-pg", help="Prodigal *.gff file")
    parser.add_argument("--prodigal-out", "-po", help="Prodigal *.out file")
    parser.add_argument(
        "--prodigal-ffn",
        "-pt",
        required=True,
        help="Prodigal *.ffn file with transcripts",
    )
    parser.add_argument(
        "--prodigal-faa", "-pp", required=True, help="Prodigal *.faa file with proteins"
    )
    parser.add_argument("--fgs-gff", "-fg", help="FragGeneScan *.gff file")
    parser.add_argument("--fgs-out", "-fo", help="FragGeneScan *.out file")
    parser.add_argument(
        "--fgs-ffn",
        "-ft",
        required=True,
        help="FragGeneScan *.ffn file with transcripts",
    )
    parser.add_argument(
        "--fgs-faa", "-fp", required=True, help="FragGeneScan *.faa file with proteins"
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

    if not args.prodigal_out and not args.prodigal_gff:
        parser.error(
            "For Prodigal, you must provide either --prodigal-out or --prodigal-gff"
        )

    if not args.fgs_out and not args.fgs_gff:
        parser.error("For FragGeneScan, you must provide either --fgs-out or --fgs-gff")

    summary = {}
    all_predictions = {}

    caller_priority = args.priority.split("_")
    logging.info(f"Caller priority: 1. {caller_priority[0]}, 2. {caller_priority[1]}")

    logging.info("Parsing Prodigal annotations...")
    if args.prodigal_out:
        all_predictions["prodigal"] = parse_prodigal_output(args.prodigal_out)
    elif args.prodigal_gff:
        all_predictions["prodigal"] = parse_gff(args.prodigal_gff)

    logging.info("Parsing FragGeneScan annotations...")
    if args.fgs_out:
        all_predictions["fgs"] = parse_fgs_output(args.fgs_out)
    elif args.fgs_gff:
        all_predictions["fgs"] = parse_gff(args.fgs_gff)

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
        "prodigal": {"proteins": args.prodigal_faa, "transcripts": args.prodigal_ffn},
        "fgs": {"proteins": args.fgs_faa, "transcripts": args.fgs_ffn},
    }
    output_fasta_files(
        merged_predictions,
        files,
        f"{args.name}.faa",
        f"{args.name}.ffn",
    )


if __name__ == "__main__":
    main()
