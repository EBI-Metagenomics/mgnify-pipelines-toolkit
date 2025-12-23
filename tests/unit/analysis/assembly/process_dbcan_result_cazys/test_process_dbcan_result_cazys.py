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
