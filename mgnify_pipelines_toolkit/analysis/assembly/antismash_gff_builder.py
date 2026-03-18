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
import re
from collections import defaultdict

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="Input JSON from antiSMASH")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output GFF3 file name")
    parser.add_argument(
        "--cds_tag",
        default="ID",
        type=str,
        help="Type of CDS ID tag to use in the GFF3 (default: locus_tag)",
    )  # The CDS' identifier changes from tool to tool.

    args = parser.parse_args()

    return args.input, args.output, args.cds_tag


def _strip_fuzzy(coord: str) -> int:
    """Remove GenBank-style fuzzy markers such as < or > and cast to int."""
    return int(coord.lstrip("<>"))


def parse_antismash_location(location: str, record_length: int | None = None) -> tuple[int, int, str]:
    """
    Parse antiSMASH JSON location strings and return:
        (start_0_based, end_1_based_like_antismash, strand)

    Supported formats:
        [81883:82231](+)
        [<81883:>82231](-)
        join{[0:345](-), [4646840:4647689](-)}
        join{[4646840:4647689](-), [0:345](-)}

    Notes
    -----
    antiSMASH JSON uses:
      - 0-based starts
      - end as right boundary
      - strand in (+)/(-)

    For origin-crossing circular joins, this function linearizes the feature:
      join{[0:C](strand), [A:L](strand)}  ->  start=A, end=L+C
    where L is the contig length.

    The caller should then convert to GFF with:
      gff_start = start + 1
      gff_end = end
    """
    location = location.strip()

    # Simple location: [start:end](strand)
    m = re.fullmatch(r"\[(<??\d+):(>??\d+)\]\(([+-])\)", location)
    if m:
        start = _strip_fuzzy(m.group(1))
        end = _strip_fuzzy(m.group(2))
        strand = m.group(3)
        return start, end, strand

    # Joined location: join{[s1:e1](strand), [s2:e2](strand)}
    m = re.fullmatch(
        r"join\{\[(<??\d+):(>??\d+)\]\(([+-])\),\s*\[(<??\d+):(>??\d+)\]\(([+-])\)\}",
        location,
    )
    if m:
        s1 = _strip_fuzzy(m.group(1))
        e1 = _strip_fuzzy(m.group(2))
        strand1 = m.group(3)
        s2 = _strip_fuzzy(m.group(4))
        e2 = _strip_fuzzy(m.group(5))
        strand2 = m.group(6)

        if strand1 != strand2:
            raise ValueError(f"Inconsistent strand in joined location: {location}")
        strand = strand1

        seg1 = (s1, e1)
        seg2 = (s2, e2)

        # Detect circular origin-crossing join regardless of segment order
        # antiSMASH JSON example:
        #   join{[0:345](-), [4646840:4647689](-)}
        # or
        #   join{[4646840:4647689](-), [0:345](-)}
        if seg1[0] == 0 and seg2[1] == record_length:
            origin_seg = seg1
            end_seg = seg2
            start = end_seg[0]
            end = record_length + origin_seg[1]
            return start, end, strand

        if seg2[0] == 0 and seg1[1] == record_length:
            origin_seg = seg2
            end_seg = seg1
            start = end_seg[0]
            end = record_length + origin_seg[1]
            return start, end, strand

        # Non-origin join fallback: collapse to min/max
        start = min(s1, s2)
        end = max(e1, e2)
        return start, end, strand

    raise ValueError(f"Unsupported antiSMASH location string: {location}")


def get_record_length(record: dict) -> int:
    """
    Infer record length from antiSMASH JSON.

    Preferred sources:
      1. 'source' feature spanning the full contig
      2. sequence length from record['seq']['data']
    """
    for feature in record.get("features", []):
        if feature.get("type") == "source":
            start, end, _ = parse_antismash_location(feature["location"])
            if start != 0:
                raise ValueError(f"Unexpected source feature start for record {record.get('id', 'unknown')}: {feature['location']}")
            return end

    seq_data = record.get("seq", {}).get("data")
    if isinstance(seq_data, str) and seq_data:
        return len(seq_data)

    raise ValueError(f"Could not infer record length for record {record.get('id', 'unknown')}: no source feature and no sequence data found")


def main():
    """Transform an antiSMASH JSON into a GFF3 with 'regions' and CDS within those regions"""

    json_input, output_file, cds_tag = parse_args()

    with open(json_input, "r") as json_data:
        antismash_analysis = json.load(json_data)

    res_dict = defaultdict(list)
    attributes_dict = defaultdict(dict)

    antismash_ver = antismash_analysis["version"]

    for record in antismash_analysis["records"]:
        record_id = record["id"]

        iter_cds = "antismash.detection.genefunctions" in record["modules"].keys()  # Flag to iterate CDS
        region_name = None
        record_length = get_record_length(record)

        for feature in record["features"]:
            if feature["type"] == "region":
                # Annotate region features
                region_name = f"{record_id}_region{feature['qualifiers']['region_number'][0]}"
                region_start, region_end, region_strand = parse_antismash_location(
                    feature["location"],
                    record_length,
                )

                res_dict["contig"].append(record_id)
                res_dict["version"].append(f"antiSMASH:{antismash_ver}")
                res_dict["type"].append("region")
                res_dict["start"].append(region_start + 1)
                res_dict["end"].append(region_end)
                res_dict["score"].append(".")
                res_dict["strand"].append(".")
                res_dict["phase"].append(".")

                product = ",".join(feature["qualifiers"].get("product", []))
                attributes_dict[region_name].update({"ID": region_name, "product": product})

            if iter_cds and feature["type"] == "CDS":
                # Annotate CDS features
                # The > and < are removed to work with pseudogene outputs in Bakta
                # A feature["location"] example that can be seen in Bakta outputs: "[81883:>82231](+)"
                start, end, strand = parse_antismash_location(
                    feature["location"],
                    record_length,
                )

                if not region_name or not (region_start <= end and start <= region_end):
                    continue

                res_dict["contig"].append(record_id)
                res_dict["version"].append(f"antiSMASH:{antismash_ver}")
                res_dict["type"].append("gene")
                res_dict["start"].append(start + 1)  # Correct for 1-based indexing
                res_dict["end"].append(end)
                res_dict["score"].append(".")
                res_dict["strand"].append(strand)
                res_dict["phase"].append(".")

                locus_tag = feature["qualifiers"][cds_tag][0]
                attributes_dict[locus_tag].update(
                    {
                        "ID": locus_tag,
                        "as_type": ",".join(feature["qualifiers"].get("gene_kind", ["other"])),
                        "gene_functions": ",".join(feature["qualifiers"].get("gene_functions", []))
                        .replace(" ", "_")
                        .replace(":_", ":")
                        .replace(";_", "%3B"),
                        "Parent": region_name,
                    }
                )

        # Extended CDS attributes
        if "antismash.detection.hmm_detection" in record["modules"].keys():
            cds_by_protocluster = record["modules"]["antismash.detection.hmm_detection"]["rule_results"]["cds_by_protocluster"]

            if not cds_by_protocluster:
                continue

            for feature in cds_by_protocluster[0][1]:
                if locus_tag := feature.get("cds_name"):
                    as_clusters = ",".join(list(feature["definition_domains"].keys()))
                    if locus_tag in attributes_dict:
                        attributes_dict[locus_tag].update({"as_gene_clusters": as_clusters})

        if "antismash.detection.genefunctions" in record["modules"].keys():
            gene_function_tools = record["modules"]["antismash.detection.genefunctions"]["tools"]
            if tool_data := gene_function_tools.get("smcogs"):
                for locus_tag in tool_data["best_hits"]:
                    smcog_id = tool_data["best_hits"][locus_tag]["reference_id"]
                    smcog_description = tool_data["best_hits"][locus_tag]["description"]

                    score = tool_data["best_hits"][locus_tag]["bitscore"]
                    e_value = tool_data["best_hits"][locus_tag]["evalue"]

                    smcog_note = f"smCOG:{smcog_id}:{smcog_description.replace(' ', '_')}(Score:{score}%3BE-value:{e_value})"
                    if locus_tag in attributes_dict.keys():
                        attributes_dict[locus_tag].update({"as_notes": smcog_note})

    attributes = [";".join(f"{k}={v}" for k, v in attrib_data.items() if v) for attrib_data in attributes_dict.values()]
    res_dict["attributes"] = attributes

    res_df = pd.DataFrame.from_dict(res_dict)

    with open(output_file, "w") as f_out:
        f_out.write("##gff-version 3\n")  # Save data to the GFF3 file with the proper header
        res_df.to_csv(f_out, header=False, index=False, sep="\t")


if __name__ == "__main__":
    main()
