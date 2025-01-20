#!/usr/bin/env python3

import argparse
import json
import logging
import os
import re

from Bio import SeqIO

__version__ = "1.0.4"


class Region:
    def __init__(self, start, end):
        # if end < start: # assuming that for +/- start always lower
        #    start, end = end, start
        self.start = int(start)
        self.end = int(end)

    def __str__(self):
        return "[" + str(self.start) + "," + str(self.end) + "]"

    def __ge__(self, other):
        return self.start >= other.end

    def __gt__(self, other):
        return self.start > other.end

    def __le__(self, other):
        return self.end <= other.start

    def __lt__(self, other):
        return self.end < other.start

    def length(self):
        return self.end - self.start + 1

    # If 'other' overlaps and has a greater end position
    def extends_right(self, other):
        if self.overlaps(other) and self.end > other.end:
            return True
        return False

    # For overlapping fragments extend start and end to match other
    def extend(self, other):
        if self.overlaps(other):
            if other.end > self.end:
                self.end = other.end
            if other.start < self.start:
                self.start = other.start

    def within(self, other):
        if self.start >= other.start and self.end <= other.end:
            return True
        return False

    # Return length of overlap between regions
    def overlaps(self, other):
        if self > other or other > self:
            return False
        # overlap = sum of the individual lengths ...
        ltot = self.length() + other.length()
        # ... minus length of the combined region (i.e. min start to max end)
        lmax = max(self.end, other.end) - min(self.start, other.start) + 1
        return ltot - lmax


# FGS has seq_id/start/end in the fasta files - use those to extract the sequences we want to keep;
# for prodigal it uses a seq_id/index_number, so need to add an extra field
class NumberedRegion(Region):
    def __init__(self, start, end, nid):
        super().__init__(start, end)
        self.nid = nid


def flatten_regions(regions):
    """Take a list of regions (possibly overlapping) and return the non-overlapping set"""
    if len(regions) < 2:
        return regions

    flattened = []
    regions = sorted(regions, key=lambda x: x.start)  # sort by start
    flattened = [regions[0]]
    regions = regions[1:]  # store the first
    for region in regions:
        if not region.overlaps(flattened[-1]):  # doesn't overlap: store new region
            flattened.append(region)
        elif region.extends_right(flattened[-1]):  # overlaps to the right: extend previous region
            flattened[-1].extend(region)
            # else end < prev end => new region within old: do nothing
    return flattened


def check_against_gaps(regions, candidates):
    """Given a set of non-overlapping gaps and a list of candidate regions, return the candidates that do not overlap"""
    regions = sorted(regions, key=lambda line: line.start)
    candidates = sorted(candidates, key=lambda line: line.start)
    selected = []
    r = 0
    if not len(regions):
        return candidates  # no existing predictions - all candidates accepted

    for c in candidates:
        if c < regions[0] or c > regions[-1]:  # outside any of the regions: just append
            selected.append(c)
        else:
            while r < len(regions) - 1 and c >= regions[r]:
                r += 1
            if c < regions[r]:  # found a gap
                selected.append(c)

    return selected


def output_prodigal(predictions, files, outputs):
    """From the combined predictions output the prodigal data"""

    sequence_set = set()
    for seq in predictions:
        for strand in ["-", "+"]:
            for region in predictions[seq][strand]:
                sequence_set.add("_".join([seq, str(region.nid)]))

    # files contains the .faa and .ffn fasta files
    for index in [1, 2]:
        sequences = []
        for record in SeqIO.parse(files[index], "fasta"):
            # remove anything after the first space
            seq_name = record.id.split(" ")[0]
            # Replace ending * #
            record.seq = record.seq.rstrip("*")
            if seq_name in sequence_set:
                sequences.append(record)

        with open(outputs[index], "a") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")


def output_fgs(predictions, files, outputs):
    """From the combined predictions output the FGS data"""
    sequence_set = set()
    for seq in predictions:
        for strand in ["-", "+"]:
            for region in predictions[seq][strand]:
                sequence_set.add("_".join([seq, str(region.start), str(region.end), strand]))

    # files contains the .faa and .ffn fasta files
    for index in [1, 2]:
        sequences = []
        for record in SeqIO.parse(files[index], "fasta"):
            # remove anything after the first space
            seq_name = record.id.split(" ")[0]
            # Replace "*" with "X"
            record.seq = record.seq.replace("*", "X")
            if seq_name in sequence_set:
                sequences.append(record)

        with open(outputs[index], "a") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")


def output_files(predictions, summary, files):
    """Output all files"""
    # To avoid that sequences get appended to the merged output files after restart,
    # make sure the files get deleted if they exist
    logging.info("Removing output files if they exist.")
    for file_ in files["merged"]:
        if os.path.exists(file_):
            logging.info(f"Removing {file_}")
            os.remove(file_)

    for caller in predictions:
        if caller == "fgs":
            output_fgs(predictions["fgs"], files["fgs"], files["merged"])
        if caller == "prodigal":
            output_prodigal(predictions["prodigal"], files["prodigal"], files["merged"])

    with open(files["merged"][0], "w") as sf:
        sf.write(json.dumps(summary, sort_keys=True, indent=4) + "\n")


def get_regions_fgs(fn):
    """Parse FGS output.
    Example:
    # >Bifidobacterium-longum-subsp-infantis-MC2-contig1
    # 256	2133	-	1	1.263995	I:	D:
    """
    regions = {}
    with open(fn) as f:
        for line in f:
            if line[0] == ">":
                id_ = line.split()[0][1:]
                regions[id_] = {}
                regions[id_]["+"] = []
                regions[id_]["-"] = []
            else:
                r = line.split()  # start end strand
                s = int(r[0])
                e = int(r[1])
                regions[id_][r[2]].append(Region(s, e))
    return regions


"""
# noqa: E501
This is from cmsearch
ERR855786.1000054-HWI-M02024:111:000000000-A8H14:1:1115:23473:14586-1 -         LSU_rRNA_bacteria    RF02541   hmm     1224     1446        5      227      +     -    6 0.61   0.8  135.2   2.8e-38 !   -
"""


def get_regions_mask(mask_file):
    """Parse masked region file (i.e. ncRNA)"""
    regions = {}
    with open(mask_file) as f:
        for line in f:
            if line[:1] == "#":
                continue
            r = line.rstrip().split()
            id_ = r[0]
            start = int(r[7])
            end = int(r[8])
            if id_ not in regions:
                regions[id_] = []
            if start > end:
                start, end = end, start
            regions[id_].append(Region(start, end))
    return regions


# # Sequence Data: seqnum=1;seqlen=25479;seqhdr="Bifidobacterium-longum-subsp-infantis-MC2-contig1"
# # Model Data: version=Prodigal.v2.6.3;run_type=Single;model="Ab initio";gc_cont=59.94;transl_table=11;uses_sd=1
# >1_1_279_+
def get_regions_prodigal(fn):
    """Parse prodigal output"""
    regions = {}
    with open(fn) as f:
        for line in f:
            if line[:12] == "# Model Data":
                continue
            if line[:15] == "# Sequence Data":
                m = re.search(r'seqhdr="(\S+)"', line)
                if m:
                    id_ = m.group(1)
                regions[id_] = {}
                regions[id_]["+"] = []
                regions[id_]["-"] = []
            else:
                r = line[1:].rstrip().split("_")
                n = int(
                    r[0]
                )  # also store the index of the fragment - prodigal uses these (rather than coords) to identify sequences in the fasta output
                s = int(r[1])
                e = int(r[2])
                regions[id_][r[3]].append(NumberedRegion(s, e, n))
    return regions


def mask_regions(regions, mask):
    """Look for overlaps of more than 5 base pairs of the supplied regions against a set of masks
    This is probably O(N^2) but, in theory, there shouldn't be many mask regions
    """
    new_regions = {}
    for seq in regions:
        new_regions[seq] = {}
        for strand in ["-", "+"]:
            new_regions[seq][strand] = []
            for r in regions[seq][strand]:
                if seq in mask:
                    overlap = 0
                    for r2 in mask[seq]:
                        if r.overlaps(r2) > 5:
                            overlap = 1
                    if not overlap:
                        new_regions[seq][strand].append(r)
                else:
                    new_regions[seq][strand].append(r)

    return new_regions


# FIXME - This won't work if we have only a single set of predictions, but then
# there's no point in trying to merge
def merge_predictions(predictions, callers):
    """Check that we have priorities set of for all callers we have data for"""
    p = set(callers)
    new_predictions = {}
    for type_ in predictions:
        if type_ not in p:
            return None
            # throw here? - if we've used a caller that we don't have a priority for

    # first set of predictions takes priority - just transfer them
    new_predictions[callers[0]] = predictions[callers[0]]

    # for now assume only two callers, but can be extended
    new_predictions[callers[1]] = {}  # empty set for second priority caller
    for seq in predictions[callers[1]]:
        new_predictions[callers[1]][seq] = {}
        for strand in ["-", "+"]:
            new_predictions[callers[1]][seq][strand] = []
            if seq in predictions[callers[0]]:  # if this sequence already has predictions
                prev_predictions = flatten_regions(
                    predictions[callers[0]][seq][strand]
                )  # non-overlapping set of existing predictions/regions
                new_predictions[callers[1]][seq][strand] = check_against_gaps(
                    prev_predictions, predictions[callers[1]][seq][strand]
                )  # plug new predictions/regions into gaps
            else:  # no existing predictions: just add them
                new_predictions[callers[1]][seq][strand] = predictions[callers[1]][seq][strand]

    return new_predictions


def get_counts(predictions):
    total = {}
    for caller in predictions:
        total[caller] = 0
        for sample in predictions[caller]:
            for strand in ["-", "+"]:
                total[caller] += len(predictions[caller][sample][strand])
    return total


def combine_main():
    parser = argparse.ArgumentParser(
        "MGnify gene caller combiner. This script will merge the gene called by prodigal and fraggenescan (in any order)"
    )
    parser.add_argument("-n", "--name", action="store", dest="name", required=True, help="basename")
    parser.add_argument("-k", "--mask", action="store", dest="mask", required=False, help="Sequence mask file")

    parser.add_argument("-a", "--prodigal-out", action="store", dest="prodigal_out", required=False, help="Stats out prodigal")
    parser.add_argument("-b", "--prodigal-ffn", action="store", dest="prodigal_ffn", required=False, help="Stats ffn prodigal")
    parser.add_argument("-c", "--prodigal-faa", action="store", dest="prodigal_faa", required=False, help="Stats faa prodigal")

    parser.add_argument("-d", "--fgs-out", action="store", dest="fgs_out", required=False, help="Stats out FGS")
    parser.add_argument("-e", "--fgs-ffn", action="store", dest="fgs_ffn", required=False, help="Stats ffn FGS")
    parser.add_argument("-f", "--fgs-faa", action="store", dest="fgs_faa", required=False, help="Stats faa FGS")

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

    parser.add_argument("-v", "--verbose", help="verbose output", dest="verbose", action="count", required=False)

    parser.add_argument("--version", action="version", version=f"{__version__}")

    args = parser.parse_args()

    # Set up logging system
    verbose_mode = args.verbose or 0

    log_level = logging.WARNING
    if verbose_mode:
        log_level = logging.DEBUG if verbose_mode > 1 else logging.INFO

    logging.basicConfig(level=log_level, format="%(levelname)s %(asctime)s - %(message)s", datefmt="%Y/%m/%d %I:%M:%S %p")

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
        all_predictions["prodigal"] = get_regions_prodigal(args.prodigal_out)

        files["prodigal"] = [args.prodigal_out, args.prodigal_ffn, args.prodigal_faa]

    if args.fgs_out:
        logging.info("FGS presented")
        logging.info("Getting FragGeneScan regions ...")
        all_predictions["fgs"] = get_regions_fgs(args.fgs_out)

        files["fgs"] = [args.fgs_out, args.fgs_ffn, args.fgs_faa]

    summary["all"] = get_counts(all_predictions)

    # Apply mask of ncRNA search
    logging.info("Masking non coding RNA regions...")
    if args.mask:
        logging.info("Reading regions for masking...")
        mask = get_regions_mask(args.mask)
        if "prodigal" in all_predictions:
            logging.info("Masking Prodigal outputs...")
            all_predictions["prodigal"] = mask_regions(all_predictions["prodigal"], mask)
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
