#!/usr/bin/env python3
# Copyright 2024-2026 EMBL - European Bioinformatics Institute
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

from __future__ import annotations

import logging
from collections.abc import Sequence

from mgnify_pipelines_toolkit.analysis.shared.bgc.models import (
    META_KEYS_ORDER,
    BGCRegion,
    CDSRec,
    MergedRegion,
)
from mgnify_pipelines_toolkit.analysis.shared.gff.io import (
    extract_id_from_attr,
    join_unique,
    replace_or_append_attr,
)

log = logging.getLogger(__name__)


def merge_overlaps(regions: Sequence[BGCRegion]) -> list[MergedRegion]:
    """Merge any overlapping intervals per contig (data-portal rule).

    :param regions: Sequence of BGCRegion objects to merge.
    :type regions: Sequence[BGCRegion]
    :returns: List of MergedRegion objects with overlapping members combined.
    :rtype: list[MergedRegion]
    """
    by_contig: dict[str, list[BGCRegion]] = {}
    for r in regions:
        by_contig.setdefault(r.contig, []).append(r)

    merged: list[MergedRegion] = []
    for contig, regs in by_contig.items():
        if not regs:
            continue
        regs_sorted = sorted(regs, key=lambda r: r.start)
        cur_start, cur_end = regs_sorted[0].start, regs_sorted[0].end
        cur_members = [regs_sorted[0]]
        for r in regs_sorted[1:]:
            if r.start <= cur_end:
                cur_end = max(cur_end, r.end)
                cur_members.append(r)
            else:
                merged.append(MergedRegion(contig=contig, start=cur_start, end=cur_end, members=cur_members))
                cur_start, cur_end, cur_members = r.start, r.end, [r]
        merged.append(MergedRegion(contig=contig, start=cur_start, end=cur_end, members=cur_members))
    return merged


def cds_within_region(start: int, end: int, rstart: int, rend: int) -> bool:
    """Return True if the CDS interval [start, end] is fully inside [rstart, rend].

    :param start: CDS start coordinate.
    :param end: CDS end coordinate.
    :param rstart: Region start coordinate.
    :param rend: Region end coordinate.
    :returns: True if CDS is fully covered by the region.
    :rtype: bool
    """
    return start >= rstart and end <= rend


def _tools_covering_cds(members: list[BGCRegion], cds_start: int, cds_end: int) -> list[str]:
    """Return sorted unique tool names for member predictions that fully cover this CDS.

    :param members: List of BGCRegion members in a merged region.
    :param cds_start: CDS start coordinate.
    :param cds_end: CDS end coordinate.
    :returns: Sorted list of tool names that cover the CDS.
    :rtype: list[str]
    """
    return sorted({m.tool for m in members if cds_within_region(cds_start, cds_end, m.start, m.end)})


def _collect_member_meta_for_cds(members: list[BGCRegion], cds_start: int, cds_end: int) -> dict[str, str]:
    """Collect metadata keys for a CDS from member predictions that cover the CDS.

    Values are unioned (unique values joined by commas).

    NOTE: antiSMASH BGC members only carry antismash_product at region-level (by design).
          antiSMASH gene-level keys are applied separately via CDS ID matching.

    :param members: List of BGCRegion members in a merged region.
    :param cds_start: CDS start coordinate.
    :param cds_end: CDS end coordinate.
    :returns: Merged metadata dict for the CDS.
    :rtype: dict[str, str]
    """
    collected: dict[str, list[str]] = {k: [] for k in META_KEYS_ORDER}
    for m in members:
        if not cds_within_region(cds_start, cds_end, m.start, m.end):
            continue
        for k in META_KEYS_ORDER:
            if k in m.attrs and m.attrs[k]:
                collected[k].append(m.attrs[k])

    out: dict[str, str] = {}
    for k, vals in collected.items():
        if vals:
            out[k] = join_unique(vals)
    return out


def support_and_filter_cds(
    contig_to_cds: dict[str, list[CDSRec]],
    merged_regions: Sequence[MergedRegion],
    antismash_gene_ann_by_id: dict[str, dict[str, str]],
    n_tools: int,
) -> dict[str, list[str]]:
    """Filter CDS to those inside BGC regions and annotate with support scores.

    For each CDS fully inside a merged region:
      - compute bgc_support = (# DISTINCT tools covering CDS) / (n_tools provided as input)
      - add bgc_tools = comma-separated list of tools covering the CDS
      - add tool metadata keys (from covering member regions)
      - apply antiSMASH gene-level keys ONLY for the matching CDS ID
      - output only those CDS (filtering out non-BGC CDS)

    :param contig_to_cds: Base CDS records grouped by contig.
    :type contig_to_cds: dict[str, list[CDSRec]]
    :param merged_regions: Merged BGC regions to check CDS coverage against.
    :type merged_regions: Sequence[MergedRegion]
    :param antismash_gene_ann_by_id: Gene-level antiSMASH annotations keyed by gene ID.
    :type antismash_gene_ann_by_id: dict[str, dict[str, str]]
    :param n_tools: Total number of tools provided (denominator for support score).
    :type n_tools: int
    :returns: Annotated CDS lines grouped by contig.
    :rtype: dict[str, list[str]]
    :raises ValueError: If n_tools < 1.
    """
    if n_tools < 1:
        raise ValueError("n_tools must be >= 1")

    out: dict[str, list[str]] = {}
    for mr in merged_regions:
        cds_list = contig_to_cds.get(mr.contig, [])
        if not cds_list or not mr.members:
            continue

        for cds in cds_list:
            # cds_list is sorted by start (guaranteed by load_base_cds); early exit is safe.
            if cds.start > mr.end:
                break
            if not cds_within_region(cds.start, cds.end, mr.start, mr.end):
                continue

            tools = _tools_covering_cds(mr.members, cds.start, cds.end)
            if not tools:
                continue

            bgc_support = len(tools) / n_tools

            cols = cds.line.split("\t")
            attr = cols[8]

            attr = replace_or_append_attr(attr, "bgc_tools", ",".join(tools))

            meta = _collect_member_meta_for_cds(mr.members, cds.start, cds.end)

            # antiSMASH gene-level metadata ONLY for this CDS (ID match)
            cds_id = extract_id_from_attr(attr)
            if cds_id and cds_id in antismash_gene_ann_by_id:
                gene_meta = antismash_gene_ann_by_id[cds_id]
                for k in (
                    "antismash_gene_function",
                    "antismash_as_type",
                    "antismash_as_gene_clusters",
                    "antismash_product",
                ):
                    if k in gene_meta and gene_meta[k]:
                        meta[k] = gene_meta[k]

            for k in META_KEYS_ORDER:
                if k in meta:
                    attr = replace_or_append_attr(attr, k, meta[k])

            attr = replace_or_append_attr(attr, "bgc_support", f"{bgc_support:.2f}")

            cols[8] = attr
            out.setdefault(mr.contig, []).append("\t".join(cols))

    for contig in out:
        out[contig].sort(key=lambda ln: int(ln.split("\t")[3]))

    return out
