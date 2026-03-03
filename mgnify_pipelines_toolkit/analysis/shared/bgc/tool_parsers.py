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
from abc import ABC, abstractmethod
from pathlib import Path

from mgnify_pipelines_toolkit.analysis.shared.bgc.models import BGCRegion
from mgnify_pipelines_toolkit.analysis.shared.gff.io import (
    iter_gff_rows,
    parse_attr_str,
)

logger = logging.getLogger(__name__)

# Gene-level annotation dict keyed by gene ID. Only antiSMASH produces these;
# GECCO and SanntiS return an empty dict as the second element of parse_regions.
GeneAnnotations = dict[str, dict[str, str]]


class BGCToolParser(ABC):
    """Abstract base class for BGC tool output parsers."""

    @property
    @abstractmethod
    def tool_name(self) -> str:
        """Normalised tool key (gecco|sanntis|antismash).

        :returns: Tool identifier string.
        :rtype: str
        """
        ...

    @abstractmethod
    def parse_regions(self, path: Path) -> tuple[list[BGCRegion], GeneAnnotations]:
        """Parse tool output GFF and return regions plus any gene-level annotations.

        :param path: Path to the tool's GFF output file.
        :type path: Path
        :returns: Tuple of (BGCRegion list, gene annotation dict). The gene
            annotation dict is empty for tools that don't produce gene-level data.
        :rtype: tuple[list[BGCRegion], GeneAnnotations]
        """
        ...


class GECCOParser(BGCToolParser):
    """Parser for GECCO BGC prediction GFF output."""

    @property
    def tool_name(self) -> str:
        """Normalised tool key.

        :returns: "gecco"
        :rtype: str
        """
        return "gecco"

    def parse_regions(self, path: Path) -> tuple[list[BGCRegion], GeneAnnotations]:
        """Parse GECCO output GFF. Only rows with a Type attribute are collected.

        :param path: Path to GECCO output GFF file (may be compressed).
        :type path: Path
        :returns: (regions, {}) — GECCO has no gene-level annotations.
        :rtype: tuple[list[BGCRegion], GeneAnnotations]
        """
        regs: list[BGCRegion] = []
        for _, cols in iter_gff_rows(path):
            contig, source, _ftype, start_s, end_s, *_rest, attr_s = cols
            attrs = parse_attr_str(attr_s)
            if "Type" in attrs:
                regs.append(
                    BGCRegion(
                        contig=contig,
                        start=int(start_s),
                        end=int(end_s),
                        tool="gecco",
                        source=source,
                        attrs={"gecco_bgc_type": attrs["Type"]},
                        label=attrs["Type"],
                    )
                )
        logger.info("Parsed GECCO regions: %d", len(regs))
        return regs, {}


class SanntiSParser(BGCToolParser):
    """Parser for SanntiS BGC prediction GFF output."""

    @property
    def tool_name(self) -> str:
        """Normalised tool key.

        :returns: "sanntis"
        :rtype: str
        """
        return "sanntis"

    def parse_regions(self, path: Path) -> tuple[list[BGCRegion], GeneAnnotations]:
        """Parse SanntiS output GFF. Only rows with MiBIG attributes are collected.

        :param path: Path to SanntiS output GFF file (may be compressed).
        :type path: Path
        :returns: (regions, {}) — SanntiS has no gene-level annotations.
        :rtype: tuple[list[BGCRegion], GeneAnnotations]
        """
        regs: list[BGCRegion] = []
        for _, cols in iter_gff_rows(path):
            contig, source, _ftype, start_s, end_s, *_rest, attr_s = cols
            attrs = parse_attr_str(attr_s)
            if "nearest_MiBIG" in attrs or "nearest_MiBIG_class" in attrs:
                ra = {k: attrs[k] for k in ("nearest_MiBIG", "nearest_MiBIG_class") if k in attrs}
                regs.append(
                    BGCRegion(
                        contig=contig,
                        start=int(start_s),
                        end=int(end_s),
                        tool="sanntis",
                        source=source,
                        attrs=ra,
                        label=attrs.get("nearest_MiBIG_class", ""),
                    )
                )
        logger.info("Parsed SanntiS regions: %d", len(regs))
        return regs, {}


class AntiSMASHParser(BGCToolParser):
    """Parser for antiSMASH BGC prediction GFF output."""

    @property
    def tool_name(self) -> str:
        """Normalised tool key.

        :returns: "antismash"
        :rtype: str
        """
        return "antismash"

    def parse_regions(self, path: Path) -> tuple[list[BGCRegion], GeneAnnotations]:
        """Parse antiSMASH output GFF, returning regions and gene-level annotations.

        Gene annotations (gene_functions, as_type, as_gene_clusters, product) are
        keyed by gene ID and returned directly rather than stored as instance state.

        :param path: Path to antiSMASH output GFF file (may be compressed).
        :type path: Path
        :returns: (regions, gene_annotations) where gene_annotations maps gene ID
            to a dict of annotation keys.
        :rtype: tuple[list[BGCRegion], GeneAnnotations]
        """
        regions_by_id: dict[str, tuple[str, int, int, str | None, str]] = {}
        gene_parent_by_id: dict[str, str] = {}
        gene_ann_by_id: GeneAnnotations = {}

        for _, cols in iter_gff_rows(path):
            contig, source, ftype, start_s, end_s, *_rest, attr_s = cols
            attrs = parse_attr_str(attr_s)

            if ftype in ("region", "biosynthetic-gene-cluster"):
                rid = attrs.get("ID")
                if not rid:
                    continue
                regions_by_id[rid] = (contig, int(start_s), int(end_s), attrs.get("product"), source)
                continue

            if ftype == "gene":
                gid = attrs.get("ID")
                parent = attrs.get("Parent")
                if not gid or not parent:
                    continue
                gene_parent_by_id[gid] = parent
                g: dict[str, str] = {}
                if "gene_functions" in attrs:
                    g["antismash_gene_function"] = attrs["gene_functions"]
                if "as_type" in attrs:
                    g["antismash_as_type"] = attrs["as_type"]
                if "as_gene_clusters" in attrs:
                    g["antismash_as_gene_clusters"] = attrs["as_gene_clusters"]
                if g:
                    gene_ann_by_id[gid] = g

        regions: list[BGCRegion] = []
        for rid, (contig, start, end, product, source) in regions_by_id.items():
            ra: dict[str, str] = {"antismash_region_id": rid}
            if product:
                ra["antismash_product"] = product
            regions.append(
                BGCRegion(
                    contig=contig,
                    start=start,
                    end=end,
                    tool="antismash",
                    source=source,
                    attrs=ra,
                    label=product or "",
                )
            )

        for gid, parent in gene_parent_by_id.items():
            if parent in regions_by_id:
                _c, _s, _e, product, _src = regions_by_id[parent]
                if product:
                    gene_ann_by_id.setdefault(gid, {})
                    gene_ann_by_id[gid]["antismash_product"] = product

        logger.info("Parsed antiSMASH regions: %d", len(regions))
        return regions, gene_ann_by_id
