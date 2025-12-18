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


import argparse
import os

rRNAs_exp = {  # noqa: N816
    "5S_rRNA": 119,
    "SSU_rRNA_bacteria": 1533,
    "LSU_rRNA_bacteria": 2925,
}
rRNAs_obs = {  # noqa: N816
    "5S_rRNA": [],
    "SSU_rRNA_bacteria": [],
    "LSU_rRNA_bacteria": [],
}
rRNAs_merged = {}  # noqa: N816


def get_tblout_column_indices(tool):
    cmsearch_indices = {
        "gene": 2,
        "start": 5,
        "end": 6,
    }
    cmscan_indices = {
        "gene": 1,
        "start": 7,
        "end": 8,
    }
    if tool == "cmsearch":
        return cmsearch_indices
    else:
        return cmscan_indices


def main():
    parser = argparse.ArgumentParser(description="Script detects bacterial 5S/SSU and LSU rRNA")
    parser.add_argument("-i", "--input", dest="input", help="rrna.tblout.deoverlapped", required=True)
    parser.add_argument(
        "-s",
        "--source",
        dest="source",
        help="Program that generated the tblout file",
        choices=["cmsearch", "cmscan"],
        required=False,
        default="cmsearch",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        help="Path to file where the output will be saved to",
        required=False,
    )

    args = parser.parse_args()

    idx = get_tblout_column_indices(args.source)

    # store start and end position of each hit
    with open(args.input, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            rfam = cols[idx["gene"]]
            if rfam not in rRNAs_exp:
                continue
            rfam_start = int(cols[idx["start"]])
            rfam_end = int(cols[idx["end"]])
            rRNAs_obs[rfam].append([rfam_start, rfam_end])

    # sort intervals by start position and merge
    for ele in rRNAs_obs.keys():
        rRNAs_merged[ele] = []
        try:
            saved = sorted(rRNAs_obs[ele])[0]
            for i in sorted(rRNAs_obs[ele]):
                if i[0] <= saved[-1]:
                    saved[-1] = max(saved[-1], i[-1])
                else:
                    rRNAs_merged[ele].append(tuple(saved))
                    saved[0] = i[0]
                    saved[1] = i[1]
            rRNAs_merged[ele].append(tuple(saved))
        except (KeyError, IndexError, TypeError):
            rRNAs_merged[ele] = [0, 0]

    # calculate total length based on merged intervals
    # The name of the file is: MGYGXX.tblout.deoverlapped#
    # and the name we need is the MGYGXX (or a different style genome accession)
    run_name = os.path.basename(args.input).split(".")[0]

    if args.outfile:
        file_out = open(args.outfile, "w")

    for rna in rRNAs_merged.keys():
        total_length = 0
        for interval in rRNAs_merged[rna]:
            try:
                total_length += interval[1] - interval[0]
            except (KeyError, IndexError, TypeError):
                total_length = 0
        new_line = "{}\t{}\t{:.2f}\n".format(run_name, rna, float(total_length) / rRNAs_exp[rna] * 100)
        if args.outfile:
            file_out.write(new_line)
        else:
            print(new_line)

    if args.outfile:
        file_out.close()


if __name__ == "__main__":
    main()
