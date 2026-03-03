#!/usr/bin/env python3
"""
test_bgc_mapper.py

Pytest suite for bgc_mapper.py script.
Tests BGC integration from GECCO/antiSMASH/SanntiS into a base GFF.
"""

from pathlib import Path

import pytest

from mgnify_pipelines_toolkit.analysis.shared.bgc_mapper import (
    BGCRegion,
    CDSRec,
    MergedRegion,
    build_region_lines,
    build_sideload_json_payload,
    load_base_cds,
    merge_overlaps,
    parse_antismash_regions_and_genes,
    parse_args,
    parse_gecco_regions,
    parse_sanntis_regions,
    support_and_filter_cds,
    validate_inputs,
    write_gff,
)

# ──────────────────────────────────────────────────────────────────────────────
# Test data fixtures - using dictionaries to build inputs
# ──────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def base_gff_data():
    """Base GFF test data as dictionary."""
    return {
        "header": "##gff-version 3\n##sequence-region b1_04_2.circ 1 412031\n",
        "entries": [
            {
                "contig": "b1_04_2.circ",
                "source": "Prodigal:002006",
                "type": "CDS",
                "start": 142110,
                "end": 142544,
                "score": ".",
                "strand": "+",
                "phase": "0",
                "attributes": "ID=b1_04_05961;inference=ab initio prediction:Prodigal:002006;locus_tag=b1_04_05961;product=hypothetical protein",
            },
            {
                "contig": "b1_04_2.circ",
                "source": "Prodigal:002006",
                "type": "CDS",
                "start": 144244,
                "end": 144501,
                "score": ".",
                "strand": "+",
                "phase": "0",
                "attributes": "ID=b1_04_05963;inference=ab initio prediction:Prodigal:002006;locus_tag=b1_04_05963;product=hypothetical protein",
            },
            {
                "contig": "b1_04_2.circ",
                "source": "Prodigal:002006",
                "type": "CDS",
                "start": 217314,
                "end": 218099,
                "score": ".",
                "strand": "+",
                "phase": "0",
                "attributes": "ID=b1_04_06033;eC_number=1.1.1.47;db_xref=COG:COG1028;inference=ab initio prediction:Prodigal:002006;locus_tag=b1_04_06033;product=Glucose 1-dehydrogenase",
            },
            {
                "contig": "b1_04_2.circ",
                "source": "Prodigal:002006",
                "type": "CDS",
                "start": 220661,
                "end": 222040,
                "score": ".",
                "strand": "-",
                "phase": "0",
                "attributes": "ID=b1_04_06036;Name=albA;gene=albA;inference=ab initio prediction:Prodigal:002006;locus_tag=b1_04_06036;product=Antilisterial bacteriocin subtilosin biosynthesis protein AlbA",
            },
            {
                "contig": "b1_04_2.circ",
                "source": "Prodigal:002006",
                "type": "CDS",
                "start": 241259,
                "end": 242098,
                "score": ".",
                "strand": "-",
                "phase": "0",
                "attributes": "ID=b1_04_06058;inference=ab initio prediction:Prodigal:002006;locus_tag=b1_04_06058;product=IS3 family transposase ISBce14",
            },
        ],
    }


@pytest.fixture
def gecco_regions_data():
    """GECCO regions test data as dictionary."""
    return {
        "regions": [
            {
                "contig": "b1_04_2.circ",
                "source": "GECCO v0.9.8",
                "start": 241259,
                "end": 308608,
                "attributes": {
                    "Type": "NRP,Polyketide",
                    "ID": "b1_04_2.circ_cluster_1",
                },
            },
            {
                "contig": "b1_04_2.circ",
                "source": "GECCO v0.9.8",
                "start": 315331,
                "end": 326252,
                "attributes": {"Type": "NRP", "ID": "b1_04_2.circ_cluster_2"},
            },
        ]
    }


@pytest.fixture
def sanntis_regions_data():
    """SanntiS regions test data as dictionary."""
    return {
        "regions": [
            {
                "contig": "b1_04_2.circ",
                "source": "SanntiSv0.9.3.3",
                "start": 142110,
                "end": 148437,
                "attributes": {
                    "nearest_MiBIG": "BGC0000849",
                    "nearest_MiBIG_class": "Other",
                },
            },
            {
                "contig": "b1_04_2.circ",
                "source": "SanntiSv0.9.3.3",
                "start": 217314,
                "end": 226605,
                "attributes": {
                    "nearest_MiBIG": "BGC0000600",
                    "nearest_MiBIG_class": "RiPP",
                },
            },
            {
                "contig": "b1_04_2.circ",
                "source": "SanntiSv0.9.3.3",
                "start": 234258,
                "end": 330050,
                "attributes": {
                    "nearest_MiBIG": "BGC0001059",
                    "nearest_MiBIG_class": "NRP Polyketide",
                },
            },
        ]
    }


@pytest.fixture
def antismash_regions_data():
    """antiSMASH regions test data as dictionary."""
    return {
        "regions": [
            {
                "contig": "b1_04_2.circ",
                "source": "antiSMASH:8.0.1",
                "start": 210661,
                "end": 232040,
                "attributes": {"ID": "b1_04_2.circ_region1", "product": "sactipeptide"},
            },
            {
                "contig": "b1_04_2.circ",
                "source": "antiSMASH:8.0.1",
                "start": 237180,
                "end": 342929,
                "attributes": {
                    "ID": "b1_04_2.circ_region2",
                    "product": "NRPS,lanthipeptide-class-ii,T1PKS",
                },
            },
        ],
        "genes": [
            {
                "id": "b1_04_06033",
                "parent": "b1_04_2.circ_region1",
                "attributes": {
                    "gene_functions": "biosynthetic-additional_(rule-based-clusters)_adh_short",
                    "as_type": "biosynthetic-additional",
                    "as_gene_clusters": "sactipeptide",
                },
            },
            {
                "id": "b1_04_06036",
                "parent": "b1_04_2.circ_region1",
                "attributes": {
                    "gene_functions": "biosynthetic_(rule-based-clusters)_sactipeptide:subtilosin",
                    "as_type": "biosynthetic",
                    "as_gene_clusters": "sactipeptide",
                },
            },
            {
                "id": "b1_04_06058",
                "parent": "b1_04_2.circ_region2",
                "attributes": {
                    "gene_functions": "other_(smcogs)_SMCOG1026:transposase",
                    "as_type": "other",
                },
            },
        ],
    }


@pytest.fixture
def temp_gff_files(tmp_path, base_gff_data):
    """Create temporary GFF files for testing."""
    # Base GFF
    base_gff = tmp_path / "test_base.gff"
    with base_gff.open("w") as f:
        f.write(base_gff_data["header"])
        for entry in base_gff_data["entries"]:
            line = "\t".join(
                [
                    entry["contig"],
                    entry["source"],
                    entry["type"],
                    str(entry["start"]),
                    str(entry["end"]),
                    entry["score"],
                    entry["strand"],
                    entry["phase"],
                    entry["attributes"],
                ]
            )
            f.write(line + "\n")

    return {"base": base_gff}


# ──────────────────────────────────────────────────────────────────────────────
# Tests for data structures
# ──────────────────────────────────────────────────────────────────────────────


def test_cds_rec_creation():
    """Test CDSRec dataclass creation."""
    cds = CDSRec(
        contig="test_contig",
        start=100,
        end=200,
        line="test_contig\tProdigal\tCDS\t100\t200\t.\t+\t0\tID=test001",
    )
    assert cds.contig == "test_contig"
    assert cds.start == 100
    assert cds.end == 200
    assert "ID=test001" in cds.line


def test_bgc_region_creation():
    """Test BGCRegion dataclass creation."""
    attrs = {"gecco_bgc_type": "NRP", "nearest_MiBIG": "BGC0000001"}
    region = BGCRegion(
        contig="test_contig",
        start=1000,
        end=5000,
        tool="gecco",
        source="GECCO v0.9.8",
        attrs=attrs,
    )
    assert region.contig == "test_contig"
    assert region.start == 1000
    assert region.end == 5000
    assert region.tool == "gecco"
    assert region.source == "GECCO v0.9.8"
    assert region.attrs["gecco_bgc_type"] == "NRP"


def test_merged_region_creation():
    """Test MergedRegion dataclass creation."""
    region1 = BGCRegion("contig1", 100, 500, "gecco", "GECCO v0.9.8", {"Type": "NRP"})
    region2 = BGCRegion("contig1", 400, 800, "sanntis", "SanntiSv0.9.3.3", {"nearest_MiBIG": "BGC001"})

    merged = MergedRegion(contig="contig1", start=100, end=800, members=[region1, region2])
    assert merged.contig == "contig1"
    assert merged.start == 100
    assert merged.end == 800
    assert len(merged.members) == 2


# ──────────────────────────────────────────────────────────────────────────────
# Tests for parsing functions
# ──────────────────────────────────────────────────────────────────────────────


def test_load_base_cds(temp_gff_files):
    """Test loading base CDS from GFF."""
    contig_to_cds = load_base_cds(temp_gff_files["base"])

    assert "b1_04_2.circ" in contig_to_cds
    assert len(contig_to_cds["b1_04_2.circ"]) == 5

    # Check first CDS
    first_cds = contig_to_cds["b1_04_2.circ"][0]
    assert first_cds.contig == "b1_04_2.circ"
    assert first_cds.start == 142110
    assert first_cds.end == 142544
    assert "ID=b1_04_05961" in first_cds.line


def test_parse_gecco_regions(tmp_path, gecco_regions_data):
    """Test parsing GECCO regions."""
    # Create temporary GECCO GFF
    gecco_gff = tmp_path / "test_gecco.gff"
    with gecco_gff.open("w") as f:
        f.write("##gff-version 3\n")
        for region in gecco_regions_data["regions"]:
            attrs = ";".join([f"{k}={v}" for k, v in region["attributes"].items()])
            line = "\t".join(
                [
                    region["contig"],
                    region["source"],
                    "BGC",
                    str(region["start"]),
                    str(region["end"]),
                    ".",
                    ".",
                    ".",
                    attrs,
                ]
            )
            f.write(line + "\n")

    regions = parse_gecco_regions(gecco_gff)
    assert len(regions) == 2
    assert regions[0].tool == "gecco"
    assert regions[0].source == "GECCO v0.9.8"
    assert regions[0].attrs["gecco_bgc_type"] == "NRP,Polyketide"
    assert regions[1].attrs["gecco_bgc_type"] == "NRP"


def test_parse_sanntis_regions(tmp_path, sanntis_regions_data):
    """Test parsing SanntiS regions."""
    # Create temporary SanntiS GFF
    sanntis_gff = tmp_path / "test_sanntis.gff"
    with sanntis_gff.open("w") as f:
        f.write("##gff-version 3\n")
        for region in sanntis_regions_data["regions"]:
            attrs = ";".join([f"{k}={v}" for k, v in region["attributes"].items()])
            line = "\t".join(
                [
                    region["contig"],
                    region["source"],
                    "CLUSTER",
                    str(region["start"]),
                    str(region["end"]),
                    ".",
                    ".",
                    ".",
                    attrs,
                ]
            )
            f.write(line + "\n")

    regions = parse_sanntis_regions(sanntis_gff)
    assert len(regions) == 3
    assert regions[0].tool == "sanntis"
    assert regions[0].source == "SanntiSv0.9.3.3"
    assert regions[0].attrs["nearest_MiBIG"] == "BGC0000849"
    assert regions[0].attrs["nearest_MiBIG_class"] == "Other"


def test_parse_antismash_regions_and_genes(tmp_path, antismash_regions_data):
    """Test parsing antiSMASH regions and gene annotations."""
    # Create temporary antiSMASH GFF
    antismash_gff = tmp_path / "test_antismash.gff"
    with antismash_gff.open("w") as f:
        f.write("##gff-version 3\n")
        # Write regions
        for region in antismash_regions_data["regions"]:
            attrs = ";".join([f"{k}={v}" for k, v in region["attributes"].items()])
            line = "\t".join(
                [
                    region["contig"],
                    region["source"],
                    "region",
                    str(region["start"]),
                    str(region["end"]),
                    ".",
                    ".",
                    ".",
                    attrs,
                ]
            )
            f.write(line + "\n")
        # Write genes
        for gene in antismash_regions_data["genes"]:
            attrs_list = [f"ID={gene['id']}", f"Parent={gene['parent']}"]
            attrs_list.extend([f"{k}={v}" for k, v in gene["attributes"].items()])
            attrs = ";".join(attrs_list)
            line = "\t".join(
                [
                    "b1_04_2.circ",
                    "antiSMASH:8.0.1",
                    "gene",
                    "100",
                    "200",
                    ".",
                    "+",
                    ".",
                    attrs,
                ]
            )
            f.write(line + "\n")

    regions, gene_ann = parse_antismash_regions_and_genes(antismash_gff)
    assert len(regions) == 2
    assert regions[0].tool == "antismash"
    assert regions[0].source == "antiSMASH:8.0.1"
    assert regions[0].attrs["antismash_product"] == "sactipeptide"

    # Check gene annotations
    assert "b1_04_06033" in gene_ann
    assert "antismash_gene_function" in gene_ann["b1_04_06033"]
    assert "antismash_product" in gene_ann["b1_04_06033"]


# ──────────────────────────────────────────────────────────────────────────────
# Tests for merging and filtering
# ──────────────────────────────────────────────────────────────────────────────


def test_merge_overlaps_no_overlap():
    """Test merging when regions don't overlap."""
    regions = [
        BGCRegion("contig1", 100, 200, "gecco", "GECCO v0.9.8", {"Type": "NRP"}),
        BGCRegion("contig1", 300, 400, "sanntis", "SanntiSv0.9.3.3", {"nearest_MiBIG": "BGC001"}),
    ]
    merged = merge_overlaps(regions)
    assert len(merged) == 2
    assert merged[0].start == 100
    assert merged[0].end == 200
    assert len(merged[0].members) == 1


def test_merge_overlaps_with_overlap():
    """Test merging when regions overlap."""
    regions = [
        BGCRegion("contig1", 100, 300, "gecco", "GECCO v0.9.8", {"Type": "NRP"}),
        BGCRegion("contig1", 250, 400, "sanntis", "SanntiSv0.9.3.3", {"nearest_MiBIG": "BGC001"}),
    ]
    merged = merge_overlaps(regions)
    assert len(merged) == 1
    assert merged[0].start == 100
    assert merged[0].end == 400
    assert len(merged[0].members) == 2


def test_merge_overlaps_adjacent():
    """Test merging when regions are adjacent (touching)."""
    regions = [
        BGCRegion("contig1", 100, 200, "gecco", "GECCO v0.9.8", {"Type": "NRP"}),
        BGCRegion("contig1", 200, 300, "sanntis", "SanntiSv0.9.3.3", {"nearest_MiBIG": "BGC001"}),
    ]
    merged = merge_overlaps(regions)
    # Adjacent regions should merge
    assert len(merged) == 1
    assert merged[0].start == 100
    assert merged[0].end == 300


def test_support_and_filter_cds_single_tool():
    """Test CDS filtering and support calculation with one tool."""
    contig_to_cds = {
        "contig1": [
            CDSRec("contig1", 150, 250, "contig1\tProdigal\tCDS\t150\t250\t.\t+\t0\tID=cds001"),
            CDSRec("contig1", 500, 600, "contig1\tProdigal\tCDS\t500\t600\t.\t+\t0\tID=cds002"),
        ]
    }

    merged_regions = [
        MergedRegion(
            "contig1",
            100,
            300,
            [BGCRegion("contig1", 100, 300, "gecco", "GECCO v0.9.8", {"Type": "NRP"})],
        )
    ]

    result = support_and_filter_cds(contig_to_cds, merged_regions, {}, n_tools=1)

    assert "contig1" in result
    assert len(result["contig1"]) == 1  # Only one CDS inside region
    assert "bgc_support=1.00" in result["contig1"][0]
    assert "bgc_tools=gecco" in result["contig1"][0]


def test_support_and_filter_cds_multiple_tools():
    """Test CDS filtering and support calculation with multiple tools."""
    contig_to_cds = {
        "contig1": [
            CDSRec("contig1", 150, 250, "contig1\tProdigal\tCDS\t150\t250\t.\t+\t0\tID=cds001"),
        ]
    }

    merged_regions = [
        MergedRegion(
            "contig1",
            100,
            300,
            [
                BGCRegion("contig1", 100, 300, "gecco", "GECCO v0.9.8", {"Type": "NRP"}),
                BGCRegion("contig1", 100, 300, "sanntis", "SanntiSv0.9.3.3", {"nearest_MiBIG": "BGC001"}),
            ],
        )
    ]

    result = support_and_filter_cds(contig_to_cds, merged_regions, {}, n_tools=3)

    assert "contig1" in result
    assert len(result["contig1"]) == 1
    # 2 tools out of 3 = 0.67
    assert "bgc_support=0.67" in result["contig1"][0]
    assert "bgc_tools=gecco,sanntis" in result["contig1"][0]


# ──────────────────────────────────────────────────────────────────────────────
# Tests for region line building
# ──────────────────────────────────────────────────────────────────────────────


def test_build_region_lines_single_member():
    """Test building region lines for non-merged regions."""
    merged_regions = [
        MergedRegion(
            "contig1",
            100,
            300,
            [BGCRegion("contig1", 100, 300, "gecco", "GECCO v0.9.8", {"gecco_bgc_type": "NRP"})],
        )
    ]

    lines = build_region_lines(merged_regions)
    assert len(lines) == 1

    contig, start, line = lines[0]
    assert contig == "contig1"
    assert start == 100
    assert "bgc_region" in line
    assert "GECCO v0.9.8" in line  # Original source preserved
    assert "ID=contig1|bgc:100-300" in line
    assert "bgc_tools=gecco" in line


def test_build_region_lines_merged():
    """Test building region lines for merged regions."""
    merged_regions = [
        MergedRegion(
            "contig1",
            100,
            400,
            [
                BGCRegion("contig1", 100, 300, "gecco", "GECCO v0.9.8", {"gecco_bgc_type": "NRP"}),
                BGCRegion("contig1", 250, 400, "sanntis", "SanntiSv0.9.3.3", {"nearest_MiBIG": "BGC001"}),
            ],
        )
    ]

    lines = build_region_lines(merged_regions)
    assert len(lines) == 1

    contig, start, line = lines[0]
    assert "bgc_merged" in line  # Merged source
    assert "ID=contig1|bgc:100-400" in line
    assert "bgc_tools=gecco,sanntis" in line
    assert "member_bgcs=2" in line


# ──────────────────────────────────────────────────────────────────────────────
# Tests for JSON sideloader
# ──────────────────────────────────────────────────────────────────────────────


def test_build_sideload_json_payload():
    """Test building antiSMASH sideloader JSON payload."""
    merged_regions = [
        MergedRegion(
            "contig1",
            100,
            300,
            [
                BGCRegion("contig1", 100, 300, "gecco", "GECCO v0.9.8", {"gecco_bgc_type": "NRP"}),
            ],
        ),
        MergedRegion(
            "contig1",
            500,
            800,
            [
                BGCRegion("contig1", 500, 700, "gecco", "GECCO v0.9.8", {"gecco_bgc_type": "Polyketide"}),
                BGCRegion("contig1", 650, 800, "sanntis", "SanntiSv0.9.3.3", {"nearest_MiBIG": "BGC002"}),
            ],
        ),
    ]

    payload = build_sideload_json_payload(
        merged_regions,
        tool_name="test_tool",
        tool_version="1.0.0",
        tool_description="Test description",
    )

    assert payload["tool"]["name"] == "test_tool"
    assert payload["tool"]["version"] == "1.0.0"
    assert len(payload["records"]) == 1
    assert payload["records"][0]["name"] == "contig1"
    assert len(payload["records"][0]["subregions"]) == 2

    # Check first subregion (0-based, end-exclusive)
    subregion1 = payload["records"][0]["subregions"][0]
    assert subregion1["start"] == 99
    assert subregion1["end"] == 300
    assert "gecco" in subregion1["details"]["bgc_tools"]

    # Check second subregion (merged)
    subregion2 = payload["records"][0]["subregions"][1]
    assert subregion2["start"] == 499
    assert subregion2["end"] == 800
    assert "member_bgcs" in subregion2["details"]
    assert subregion2["details"]["member_bgcs"] == "2"


# ──────────────────────────────────────────────────────────────────────────────
# Tests for command-line argument parsing
# ──────────────────────────────────────────────────────────────────────────────


def test_parse_args_minimal():
    """Test minimal argument parsing."""
    args = parse_args(
        [
            "--base_gff",
            "base.gff",
            "--gecco_gff",
            "gecco.gff",
            "--output_gff",
            "output.gff",
        ]
    )
    assert args.base_gff == Path("base.gff")
    assert args.gecco_gff == Path("gecco.gff")
    assert args.output_gff == Path("output.gff")
    assert args.antismash_gff is None
    assert args.sanntis_gff is None


def test_parse_args_all_tools():
    """Test argument parsing with all tools."""
    args = parse_args(
        [
            "--base_gff",
            "base.gff",
            "--gecco_gff",
            "gecco.gff",
            "--antismash_gff",
            "antismash.gff",
            "--sanntis_gff",
            "sanntis.gff",
            "--output_gff",
            "output.gff",
        ]
    )
    assert args.gecco_gff == Path("gecco.gff")
    assert args.antismash_gff == Path("antismash.gff")
    assert args.sanntis_gff == Path("sanntis.gff")


def test_validate_inputs_missing_base(tmp_path):
    """Test validation with missing base GFF."""
    args_dict = {
        "base_gff": tmp_path / "nonexistent.gff",
        "gecco_gff": None,
        "antismash_gff": None,
        "sanntis_gff": None,
        "output_gff": tmp_path / "output.gff",
    }

    class Args:
        pass

    args = Args()
    for k, v in args_dict.items():
        setattr(args, k, v)

    with pytest.raises(FileNotFoundError, match="Base GFF not found"):
        validate_inputs(args)


def test_validate_inputs_no_predictors(tmp_path):
    """Test validation with no predictor GFFs."""
    base_gff = tmp_path / "base.gff"
    base_gff.write_text("##gff-version 3\n")

    args_dict = {
        "base_gff": base_gff,
        "gecco_gff": None,
        "antismash_gff": None,
        "sanntis_gff": None,
        "output_gff": tmp_path / "output.gff",
    }

    class Args:
        pass

    args = Args()
    for k, v in args_dict.items():
        setattr(args, k, v)

    with pytest.raises(ValueError, match="At least one optional predictor"):
        validate_inputs(args)


# ──────────────────────────────────────────────────────────────────────────────
# Integration tests
# ──────────────────────────────────────────────────────────────────────────────


def test_full_workflow(tmp_path):
    """Test full workflow from input files to output."""
    # Create test files from dictionaries
    base_gff_dict = {
        "entries": [
            {
                "contig": "contig1",
                "source": "Prodigal",
                "type": "CDS",
                "start": 150,
                "end": 250,
                "score": ".",
                "strand": "+",
                "phase": "0",
                "attributes": "ID=cds001;product=test protein",
            },
            {
                "contig": "contig1",
                "source": "Prodigal",
                "type": "CDS",
                "start": 350,
                "end": 450,
                "score": ".",
                "strand": "+",
                "phase": "0",
                "attributes": "ID=cds002;product=another protein",
            },
        ]
    }

    gecco_dict = {
        "regions": [
            {
                "contig": "contig1",
                "source": "GECCO v0.9.8",
                "start": 100,
                "end": 500,
                "attributes": {"Type": "NRP"},
            }
        ]
    }

    # Write base GFF
    base_gff = tmp_path / "base.gff"
    with base_gff.open("w") as f:
        f.write("##gff-version 3\n")
        for entry in base_gff_dict["entries"]:
            line = "\t".join(
                [
                    entry["contig"],
                    entry["source"],
                    entry["type"],
                    str(entry["start"]),
                    str(entry["end"]),
                    entry["score"],
                    entry["strand"],
                    entry["phase"],
                    entry["attributes"],
                ]
            )
            f.write(line + "\n")

    # Write GECCO GFF
    gecco_gff = tmp_path / "gecco.gff"
    with gecco_gff.open("w") as f:
        f.write("##gff-version 3\n")
        for region in gecco_dict["regions"]:
            attrs = ";".join([f"{k}={v}" for k, v in region["attributes"].items()])
            line = "\t".join(
                [
                    region["contig"],
                    region["source"],
                    "BGC",
                    str(region["start"]),
                    str(region["end"]),
                    ".",
                    ".",
                    ".",
                    attrs,
                ]
            )
            f.write(line + "\n")

    # Load and process
    contig_to_cds = load_base_cds(base_gff)
    gecco_regions = parse_gecco_regions(gecco_gff)
    merged = merge_overlaps(gecco_regions)
    filtered_cds = support_and_filter_cds(contig_to_cds, merged, {}, n_tools=1)
    region_lines = build_region_lines(merged)

    # Write output
    output_gff = tmp_path / "output.gff"
    write_gff(output_gff, region_lines, filtered_cds)

    # Verify output
    assert output_gff.exists()
    content = output_gff.read_text()
    assert "##gff-version 3" in content
    assert "bgc_region" in content
    assert "bgc_support=1.00" in content
    assert "bgc_tools=gecco" in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
