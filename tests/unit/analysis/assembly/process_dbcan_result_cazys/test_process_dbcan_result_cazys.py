#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2023-2025 EMBL - European Bioinformatics Institute
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
#

import hashlib
import tempfile
from pathlib import Path

import pytest

from mgnify_pipelines_toolkit.analysis.assembly.process_dbcan_result_cazys import (
    load_gff,
    load_substrates,
    print_gff,
)


@pytest.fixture
def fixtures_path() -> Path:
    """Get the path to test fixtures.

    :return: Path to the fixtures directory
    :rtype: Path
    """
    # Go up from this file to tests/ directory, then to fixtures/process_dbcan
    tests_dir = Path(__file__).parents[4]
    return tests_dir / "fixtures" / "process_dbcan"


def test_load_substrates_with_semicolons(fixtures_path: Path) -> None:
    """Test that substrates separated by semicolons are parsed correctly.

    :param fixtures_path: Path to the fixtures directory
    :type fixtures_path: Path
    """
    hmm_file = fixtures_path / "5.1.2" / "dbcan_sub_hmm.tsv"
    substrates = load_substrates(hmm_file)

    # Test gene with semicolon-separated substrates
    assert "ERZ23872251_126_1" in substrates
    assert substrates["ERZ23872251_126_1"] == {"cellulose", "chitin", "xylan"}

    # Test gene with semicolon-separated substrates
    assert "ERZ23872251_132_1" in substrates
    assert substrates["ERZ23872251_132_1"] == {"beta-glucan", "plant glycoside"}

    # Test gene with a single substrate
    assert "ERZ23872251_29_1" in substrates
    assert substrates["ERZ23872251_29_1"] == {"lignin"}

    # Test gene with multiple HMM hits (same substrate)
    assert "ERZ23872251_451_1" in substrates
    assert substrates["ERZ23872251_451_1"] == {"chitin", "peptidoglycan"}


def test_print_gff_substrate_parsing(fixtures_path: Path) -> None:
    """Test that substrates are correctly parsed and written to GFF output.

    This test verifies that substrates with various separators (semicolons, commas)
    are correctly parsed and appear in the GFF output with proper formatting.

    :param fixtures_path: Path to the fixtures directory
    :type fixtures_path: Path
    """
    hmm_file = fixtures_path / "5.1.2" / "dbcan_sub_hmm.tsv"
    overview_file = fixtures_path / "5.1.2" / "dbcan_overview.tsv"
    input_gff = fixtures_path / "5.1.2" / "annotations.gff"

    substrates = load_substrates(hmm_file)
    genome_gff_lines = load_gff(input_gff)

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".gff") as outfile:
        output_path = outfile.name

    try:
        print_gff(overview_file, output_path, "5.1.2", substrates, genome_gff_lines)

        with open(output_path, "r") as f:
            content = f.read()
            lines = content.strip().split("\n")

        # Check header
        assert lines[0] == "##gff-version 3"

        # Find lines for specific genes and check substrate fields
        # ERZ23872251_126_1 should have substrates: cellulose,chitin,xylan
        cbm2_line = [line for line in lines if "ERZ23872251_126_1" in line and "CBM2_e104" in line]
        assert len(cbm2_line) == 1
        assert "substrate_dbcan-sub=cellulose,chitin,xylan" in cbm2_line[0]

        # ERZ23872251_132_1 should have substrates: beta-glucan,plant glycoside
        gh1_line = [line for line in lines if "ERZ23872251_132_1" in line and "GH1_e43" in line]
        assert len(gh1_line) == 1
        assert "substrate_dbcan-sub=beta-glucan,plant glycoside" in gh1_line[0]

        # ERZ23872251_29_1 should have substrate: lignin
        aa3_line = [line for line in lines if "ERZ23872251_29_1" in line and "AA3_4" in line]
        assert len(aa3_line) == 1
        assert "substrate_dbcan-sub=lignin" in aa3_line[0]

        # ERZ23872251_451_1 should have substrates: chitin,peptidoglycan (deduplicated from 2 hits)
        cbm50_line = [line for line in lines if "ERZ23872251_451_1" in line and "CBM50_e958" in line]
        assert len(cbm50_line) == 1
        assert "substrate_dbcan-sub=chitin,peptidoglycan" in cbm50_line[0]

        # ERZ23872251_422_1 and ERZ23872251_46_1 have no substrate (-)
        gt20_line = [line for line in lines if "ERZ23872251_422_1" in line and "GT20_e6" in line]
        assert len(gt20_line) == 1
        assert "substrate_dbcan-sub=N/A" in gt20_line[0]

        # ERZ23872251_495_1 should have substrate: lignin
        aa2_line = [line for line in lines if "ERZ23872251_495_1" in line and "AA2_e11" in line]
        assert len(aa2_line) == 1
        assert "substrate_dbcan-sub=lignin" in aa2_line[0]

    finally:
        Path(output_path).unlink(missing_ok=True)


def test_legacy_fixtures_checksum(fixtures_path: Path) -> None:
    """Test that legacy fixtures still produce expected output.

    This ensures backwards compatibility with older dbcan versions.

    :param fixtures_path: Path to the fixtures directory
    :type fixtures_path: Path
    """
    hmm_file = fixtures_path / "dbCANsub_hmm_results.tsv"
    overview_file = fixtures_path / "overview.tsv"
    input_gff = fixtures_path / "input_cgc.gff"

    substrates = load_substrates(hmm_file)
    genome_gff_lines = load_gff(input_gff)

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".gff") as outfile:
        output_path = outfile.name

    try:
        print_gff(overview_file, output_path, "4.6", substrates, genome_gff_lines)

        with open(output_path, "rb") as f:
            content = f.read()
            md5 = hashlib.md5(content).hexdigest()

        assert md5 == "00b2c1d3f1024163099086d70ff12d50"

    finally:
        Path(output_path).unlink(missing_ok=True)


@pytest.mark.parametrize(
    "raw_ec, expected_ec",
    [
        ("3.2.1.24:1|3.2.1.113:1|3.2.1.-:8", "3.2.1.24|3.2.1.113|3.2.1.-"),
        ("3.2.1.24:1;3.2.1.113:1;3.2.1.-:8", "3.2.1.24|3.2.1.113|3.2.1.-"),
        ("3.2.1.4:1;3.2.1.8:1|3.2.1.4:2;3.2.1.78:2", "3.2.1.4|3.2.1.8|3.2.1.4|3.2.1.78"),
        (" - | N/A ; ; 3.2.1.-:8 ; 3.2.1.24 ", "3.2.1.-|3.2.1.24"),
        ("-|N/A;; ", None),
    ],
)
def test_print_gff_ec_formats(tmp_path: Path, raw_ec: str, expected_ec: str | None) -> None:
    """Preserve complete EC lists from old, new and mixed dbCAN formats."""
    gene = "ERZ19210871_100019_1"
    overview = tmp_path / "overview.tsv"
    # Reorder columns and add the new Substrate column to exercise header-based reading.
    overview.write_text(f"Substrate\tRecommend Results\t#ofTools\tEC#\tGene ID\nhost glycan\tGH92_e0\t2\t{raw_ec}\t{gene}\n")
    output = tmp_path / "output.gff"
    print_gff(overview, str(output), "5.2.9", {}, {gene: [f"contig\t.\tCDS\t1\t501\t.\t+\t.\tID={gene}\n"]})
    content = output.read_text()
    if expected_ec is None:
        assert "eC_number=" not in content
    else:
        assert f"eC_number={expected_ec};" in content
    assert "protein_family=GH92_e0;" in content
    assert "num_tools=2;" in content
    # The overview's substrate column must not replace the HMM-derived substrate.
    assert "substrate_dbcan-sub=N/A;" in content


@pytest.mark.parametrize(
    "content, message",
    [
        ("Gene ID\tEC#\t#ofTools\n", "Missing required dbCAN overview columns: Recommend Results"),
        ("", "Missing required dbCAN overview columns:"),
        ("Gene ID\tEC#\t#ofTools\tRecommend Results\nERZ1\t-\t2\n", "Incomplete dbCAN overview row 2: Recommend Results"),
    ],
)
def test_print_gff_rejects_incomplete_overview(tmp_path: Path, content: str, message: str) -> None:
    """Reject unsupported schemas and truncated rows rather than misassigning columns."""
    overview = tmp_path / "overview.tsv"
    overview.write_text(content)
    with pytest.raises(ValueError, match=message):
        print_gff(overview, str(tmp_path / "output.gff"), "5.2.9", {}, {})


@pytest.mark.parametrize(
    "raw_substrate, expected",
    [
        ("alpha-glucan;alpha-glucan;-", {"alpha-glucan"}),
        (" host glycan ; - , N/A , chitin ", {"host glycan", "chitin"}),
        ("-;N/A, ; -", None),
        ("", None),
    ],
)
def test_load_substrates_filters_missing_tokens(tmp_path: Path, raw_substrate: str, expected: set[str] | None) -> None:
    """Filter missing values inside mixed lists and preserve the N/A fallback."""
    gene = "ERZ19210871_61249_1"
    hmm = tmp_path / "hmm.tsv"
    hmm.write_text(f"Substrate\tTarget Name\ti-Evalue\n{raw_substrate}\t{gene}\t1e-29\n")
    substrates = load_substrates(hmm)
    assert substrates.get(gene) == expected
    overview = tmp_path / "overview.tsv"
    overview.write_text(f"Gene ID\tEC#\t#ofTools\tRecommend Results\n{gene}\t-\t2\tCBM48_e3\n")
    output = tmp_path / "output.gff"
    print_gff(overview, str(output), "5.2.9", substrates, {gene: [f"contig\t.\tCDS\t1\t78\t.\t+\t.\tID={gene}\n"]})
    expected_text = ",".join(sorted(expected)) if expected else "N/A"
    assert f"substrate_dbcan-sub={expected_text};" in output.read_text()


@pytest.mark.parametrize(
    "run, dbcan_version, expected_ec",
    [
        ("6.0", "5.1.2", "3.2.1.-|3.2.1.113|3.2.1.114|3.2.1.24"),
        ("6.1", "5.2.9", "3.2.1.24|3.2.1.113|3.2.1.114|3.2.1.-"),
    ],
)
def test_pipeline_comparison_preserves_all_ecs(fixtures_path: Path, tmp_path: Path, run: str, dbcan_version: str, expected_ec: str) -> None:
    """Preserve every EC in real dbCAN records from the 6.0 and 6.1 runs."""
    example = fixtures_path / "pipeline_comparison" / run
    output = tmp_path / "cazys.gff"
    print_gff(
        example / "overview.tsv",
        str(output),
        dbcan_version,
        load_substrates(example / "sub_hmm.tsv"),
        load_gff(example / "annotation_summary.gff"),
    )
    assert output.read_text() == (
        "##gff-version 3\n"
        f"ERZ19210871_100019\tdbCAN:{dbcan_version}\tCAZyme\t1\t501\t.\t+\t.\t"
        "ID=ERZ19210871_100019_1;protein_family=GH92_e0;substrate_dbcan-sub=host glycan;"
        f"eC_number={expected_ec};num_tools=2;\n"
    )
