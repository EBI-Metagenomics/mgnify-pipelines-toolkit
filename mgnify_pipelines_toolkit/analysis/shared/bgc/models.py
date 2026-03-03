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

from pydantic import BaseModel, Field

# Keys that may be added to CDS/BGC attributes.
# Order is preserved when appending to GFF attribute strings.
META_KEYS_ORDER = [
    "antismash_gene_function",
    "antismash_product",
    "antismash_as_type",
    "antismash_as_gene_clusters",
    "gecco_bgc_type",
    "nearest_MiBIG",
    "nearest_MiBIG_class",
]


class BGCRegion(BaseModel):
    """A single BGC region prediction from one tool.

    :param contig: Contig/sequence identifier.
    :param start: Start coordinate (1-based, inclusive).
    :param end: End coordinate (1-based, inclusive).
    :param tool: Normalised tool key (gecco|sanntis|antismash).
    :param source: Exact column 2 source string from the predictor GFF.
    :param attrs: Tool-specific metadata attributes for this region.
    :param label: Short human-readable label for this region (set by the parser).
    """

    contig: str
    start: int
    end: int
    tool: str
    source: str
    attrs: dict[str, str] = Field(default_factory=dict)
    label: str = ""


class CDSRec(BaseModel):
    """A CDS record from the base GFF, preserving the original line.

    :param contig: Contig/sequence identifier.
    :param start: Start coordinate (1-based, inclusive).
    :param end: End coordinate (1-based, inclusive).
    :param line: Original base GFF line (no trailing newline).
    """

    contig: str
    start: int
    end: int
    line: str


class MergedRegion(BaseModel):
    """A merged BGC region spanning one or more overlapping BGCRegion predictions.

    :param contig: Contig/sequence identifier.
    :param start: Merged start coordinate.
    :param end: Merged end coordinate.
    :param members: List of constituent BGCRegion predictions.
    """

    contig: str
    start: int
    end: int
    members: list[BGCRegion] = Field(default_factory=list)
