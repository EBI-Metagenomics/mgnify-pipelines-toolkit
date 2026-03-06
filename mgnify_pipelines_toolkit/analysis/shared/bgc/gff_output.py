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

from collections.abc import Sequence

from mgnify_pipelines_toolkit.analysis.shared.bgc.models import (
    META_KEYS_ORDER,
    BGCRegion,
    MergedRegion,
)
from mgnify_pipelines_toolkit.analysis.shared.gff.io import (
    attrs_to_str_with_id_first,
    join_unique,
)


def _collect_member_meta_for_region(members: list[BGCRegion]) -> dict[str, str]:
    """Union metadata keys across member regions (unique values joined by commas).

    antiSMASH contributes only antismash_product at region-level.

    :param members: List of BGCRegion members.
    :type members: list[BGCRegion]
    :returns: Merged metadata dict for the region.
    :rtype: dict[str, str]
    """
    collected: dict[str, list[str]] = {k: [] for k in META_KEYS_ORDER}
    for m in members:
        for k in META_KEYS_ORDER:
            if m.attrs.get(k):
                collected[k].append(m.attrs[k])

    out: dict[str, str] = {}
    for k, vals in collected.items():
        if vals:
            out[k] = join_unique(vals)
    return out


def build_region_lines(
    merged_regions: Sequence[MergedRegion],
) -> list[tuple[str, int, str]]:
    """Build GFF3 bgc_region feature lines for all merged regions.

    - members > 1  => source 'bgc_merged'
    - members == 1 => source is the EXACT original source string from the predictor GFF
    - type 'bgc_region'
    - ID first: contig|bgc:start-end
    - bgc_tools= added for all cases

    :param merged_regions: Sequence of MergedRegion objects to emit lines for.
    :type merged_regions: Sequence[MergedRegion]
    :returns: List of (contig, start, line) tuples for sorting.
    :rtype: list[tuple[str, int, str]]
    """
    lines: list[tuple[str, int, str]] = []
    for mr in merged_regions:
        if not mr.members:
            continue
        tools = sorted({m.tool for m in mr.members})
        attrs: dict[str, str] = {
            "ID": f"{mr.contig}|bgc:{mr.start}-{mr.end}",
            "bgc_tools": ",".join(tools),
        }

        if len(mr.members) > 1:
            source = "bgc_merged"
            attrs["member_bgcs"] = str(len(mr.members))
            attrs.update(_collect_member_meta_for_region(mr.members))
        else:
            source = mr.members[0].source
            attrs.update(mr.members[0].attrs)

        line = "\t".join(
            [
                mr.contig,
                source,
                "bgc_region",
                str(mr.start),
                str(mr.end),
                ".",
                ".",
                ".",
                attrs_to_str_with_id_first(attrs),
            ]
        )
        lines.append((mr.contig, mr.start, line))
    return lines
