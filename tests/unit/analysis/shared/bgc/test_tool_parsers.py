#!/usr/bin/env python3
"""
test_tool_parsers.py

Tests for mgnify_pipelines_toolkit.analysis.shared.bgc.tool_parsers –
BGCToolParser ABC and concrete parser implementations.
BGCRegion model tests are in test_bgc_mapper.py (models section).
"""

import pytest

from mgnify_pipelines_toolkit.analysis.shared.bgc.tool_parsers import (
    AntiSMASHParser,
    GECCOParser,
    SanntiSParser,
)

# ──────────────────────────────────────────────────────────────────────────────
# GECCOParser
# ──────────────────────────────────────────────────────────────────────────────


def test_gecco_parser_tool_name():
    """GECCOParser.tool_name returns 'gecco'."""
    assert GECCOParser().tool_name == "gecco"


def test_gecco_parser_parse_regions(tmp_path):
    """GECCOParser.parse_regions returns one BGCRegion per row with a Type attribute."""
    gff = tmp_path / "gecco.gff"
    gff.write_text(
        "##gff-version 3\n"
        "ctg1\tGECCO v0.9.8\tBGC\t100\t500\t.\t.\t.\tID=cluster1;Type=NRP,Polyketide\n"
        "ctg1\tGECCO v0.9.8\tBGC\t600\t900\t.\t.\t.\tID=cluster2;Type=NRP\n"
        "ctg1\tGECCO v0.9.8\tCDS\t100\t200\t.\t+\t0\tID=cds1\n"  # no Type => skip
    )
    regions, gene_ann = GECCOParser().parse_regions(gff)
    assert len(regions) == 2
    assert gene_ann == {}
    assert all(r.tool == "gecco" for r in regions)
    assert all(r.source == "GECCO v0.9.8" for r in regions)
    assert regions[0].attrs["gecco_bgc_type"] == "NRP,Polyketide"
    assert regions[0].label == "NRP,Polyketide"
    assert regions[1].attrs["gecco_bgc_type"] == "NRP"
    assert regions[1].label == "NRP"
    assert regions[0].start == 100
    assert regions[0].end == 500


def test_gecco_parser_empty_file(tmp_path):
    """GECCOParser.parse_regions returns empty list for a header-only file."""
    gff = tmp_path / "gecco.gff"
    gff.write_text("##gff-version 3\n")
    regions, gene_ann = GECCOParser().parse_regions(gff)
    assert regions == []
    assert gene_ann == {}


# ──────────────────────────────────────────────────────────────────────────────
# SanntiSParser
# ──────────────────────────────────────────────────────────────────────────────


def test_sanntis_parser_tool_name():
    """SanntiSParser.tool_name returns 'sanntis'."""
    assert SanntiSParser().tool_name == "sanntis"


def test_sanntis_parser_parse_regions(tmp_path):
    """SanntiSParser.parse_regions returns one region per row with MiBIG attributes."""
    gff = tmp_path / "sanntis.gff"
    gff.write_text(
        "##gff-version 3\n"
        "ctg1\tSanntiSv0.9.3.3\tCLUSTER\t100\t500\t.\t.\t.\tnearest_MiBIG=BGC001;nearest_MiBIG_class=Other\n"
        "ctg1\tSanntiSv0.9.3.3\tCLUSTER\t600\t900\t.\t.\t.\tnearest_MiBIG=BGC002\n"
        "ctg1\tSanntiSv0.9.3.3\tCLUSTER\t1000\t1200\t.\t.\t.\tsome_other_attr=val\n"  # no MiBIG => skip
    )
    regions, gene_ann = SanntiSParser().parse_regions(gff)
    assert len(regions) == 2
    assert gene_ann == {}
    assert all(r.tool == "sanntis" for r in regions)
    assert all(r.source == "SanntiSv0.9.3.3" for r in regions)
    assert regions[0].attrs["nearest_MiBIG"] == "BGC001"
    assert regions[0].attrs["nearest_MiBIG_class"] == "Other"
    assert regions[0].label == "Other"
    assert "nearest_MiBIG_class" not in regions[1].attrs
    assert regions[1].label == ""  # no class attribute => empty label


# ──────────────────────────────────────────────────────────────────────────────
# AntiSMASHParser
# ──────────────────────────────────────────────────────────────────────────────


def test_antismash_parser_tool_name():
    """AntiSMASHParser.tool_name returns 'antismash'."""
    assert AntiSMASHParser().tool_name == "antismash"


def _write_antismash_gff(path, regions_data, genes_data):
    """Helper: write a minimal antiSMASH-style GFF for tests."""
    lines = ["##gff-version 3\n"]
    for r in regions_data:
        attrs = ";".join(f"{k}={v}" for k, v in r["attrs"].items())
        lines.append(
            "\t".join(
                [
                    r["contig"],
                    r.get("source", "antiSMASH:8.0.1"),
                    r.get("ftype", "region"),
                    str(r["start"]),
                    str(r["end"]),
                    ".",
                    ".",
                    ".",
                    attrs,
                ]
            )
            + "\n"
        )
    for g in genes_data:
        attrs_parts = [f"ID={g['id']}", f"Parent={g['parent']}"]
        attrs_parts.extend(f"{k}={v}" for k, v in g.get("attrs", {}).items())
        lines.append(
            "\t".join(
                [
                    g.get("contig", "ctg1"),
                    "antiSMASH:8.0.1",
                    "gene",
                    "100",
                    "200",
                    ".",
                    "+",
                    ".",
                    ";".join(attrs_parts),
                ]
            )
            + "\n"
        )
    path.write_text("".join(lines))


def test_antismash_parser_parse_regions(tmp_path):
    """AntiSMASHParser.parse_regions extracts regions from 'region' feature rows."""
    gff = tmp_path / "antismash.gff"
    _write_antismash_gff(
        gff,
        regions_data=[
            {
                "contig": "ctg1",
                "source": "antiSMASH:8.0.1",
                "start": 100,
                "end": 500,
                "attrs": {"ID": "ctg1_region1", "product": "sactipeptide"},
            },
            {
                "contig": "ctg1",
                "source": "antiSMASH:8.0.1",
                "start": 600,
                "end": 900,
                "attrs": {"ID": "ctg1_region2", "product": "NRPS,T1PKS"},
            },
        ],
        genes_data=[],
    )
    regions, gene_ann = AntiSMASHParser().parse_regions(gff)
    assert len(regions) == 2
    assert gene_ann == {}
    assert all(r.tool == "antismash" for r in regions)
    assert all(r.source == "antiSMASH:8.0.1" for r in regions)
    products = {r.attrs.get("antismash_product") for r in regions}
    assert "sactipeptide" in products
    assert "NRPS,T1PKS" in products
    labels = {r.label for r in regions}
    assert "sactipeptide" in labels
    assert "NRPS,T1PKS" in labels


def test_antismash_parser_gene_annotations_returned_in_tuple(tmp_path):
    """AntiSMASHParser.parse_regions returns gene-level annotations in the tuple."""
    gff = tmp_path / "antismash.gff"
    _write_antismash_gff(
        gff,
        regions_data=[
            {
                "contig": "ctg1",
                "source": "antiSMASH:8.0.1",
                "start": 100,
                "end": 500,
                "attrs": {"ID": "ctg1_region1", "product": "sactipeptide"},
            },
        ],
        genes_data=[
            {
                "id": "gene001",
                "parent": "ctg1_region1",
                "attrs": {
                    "gene_functions": "biosynthetic_(rule-based-clusters)_sactipeptide:subtilosin",
                    "as_type": "biosynthetic",
                    "as_gene_clusters": "sactipeptide",
                },
            },
        ],
    )
    regions, gene_ann = AntiSMASHParser().parse_regions(gff)
    assert len(regions) == 1
    assert "gene001" in gene_ann
    assert gene_ann["gene001"]["antismash_gene_function"] == "biosynthetic_(rule-based-clusters)_sactipeptide:subtilosin"
    assert gene_ann["gene001"]["antismash_as_type"] == "biosynthetic"
    assert gene_ann["gene001"]["antismash_as_gene_clusters"] == "sactipeptide"
    # product should be inherited from parent region
    assert gene_ann["gene001"]["antismash_product"] == "sactipeptide"


def test_antismash_parser_empty_file(tmp_path):
    """AntiSMASHParser.parse_regions returns empty list and dict for header-only file."""
    gff = tmp_path / "antismash.gff"
    gff.write_text("##gff-version 3\n")
    regions, gene_ann = AntiSMASHParser().parse_regions(gff)
    assert regions == []
    assert gene_ann == {}


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
