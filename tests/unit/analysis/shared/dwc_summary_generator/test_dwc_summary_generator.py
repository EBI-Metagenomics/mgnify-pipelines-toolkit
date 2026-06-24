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

import pandas as pd
import shutil
from pathlib import Path
from mgnify_pipelines_toolkit.analysis.shared.dwc_summary_generator import get_asv_dict, get_closedref_dict
from click.testing import CliRunner
from mgnify_pipelines_toolkit.analysis.shared.dwc_summary_generator import cli


def test_get_asv_dict_empty_files(tmp_path: Path):
    """
    Test get_asv_dict to ensure it correctly skips runs with empty critical files.
    """
    # Setup
    run_acc = "ERR123456"
    runs_df = pd.DataFrame([{"run": run_acc, "status": "all_results"}])
    db = "DADA2-SILVA"
    otu_dir = tmp_path / "otu"
    otu_dir.mkdir()

    # Use the fixture OTU file instead of creating a dummy one
    fixture_otu = Path("tests/fixtures/refdb_otus/SILVA-SSU.otu")
    shutil.copy(fixture_otu, otu_dir / "SILVA-SSU.otu")

    root_path = tmp_path / "results"
    run_dir = root_path / run_acc
    tax_summary_dir = run_dir / "taxonomy-summary" / db
    asv_dir = run_dir / "asv"

    tax_summary_dir.mkdir(parents=True)
    asv_dir.mkdir(parents=True)

    # Files to be tested
    mapseq_file = tax_summary_dir / f"test_{db}.mseq"
    tax_file = asv_dir / f"test_{db}_asv_tax.tsv"
    asv_fasta_file = asv_dir / "test_asv_seqs.fasta"
    count_dir = asv_dir / "16S-V3-V4"
    count_dir.mkdir()
    count_file = count_dir / "test_16S-V3-V4_counts.tsv"

    # 1. Test mapseq_file empty
    mapseq_file.touch()
    tax_file.write_text("dummy")
    asv_fasta_file.write_text("dummy")
    count_file.write_text("dummy")

    res = get_asv_dict(runs_df, root_path, db, otu_dir)
    assert run_acc not in res

    # 2. Test tax_file empty
    mapseq_file.write_text(
        "asv\tdbhit\tdbhitEval\tdbhitIdentity\tdbhitQueryStart\tdbhitQueryEnd\tdbhitTargetStart\tdbhitTargetEnd\tdbhitScore\tdbhitStart\tdbhitEnd\nASV1\thit1\t0.0\t100.0\t1\t100\t1\t100\t100\t1\t100\n"
    )
    tax_file.write_text("")  # Empty
    res = get_asv_dict(runs_df, root_path, db, otu_dir)
    assert run_acc not in res

    # 3. Test asv_fasta_file empty
    tax_file.write_text(
        "ASV\tSuperkingdom\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\nASV1\td__Bacteria\tk__Kingdom\tp__Phylum\tc__Class\to__Order\tf__Family\tg__Genus\ts__Species\n"
    )
    asv_fasta_file.write_text("")  # Empty
    res = get_asv_dict(runs_df, root_path, db, otu_dir)
    assert run_acc not in res

    # 4. Test count_file empty
    asv_fasta_file.write_text(">ASV1\nATGC\n")
    count_file.write_text("")  # Empty
    res = get_asv_dict(runs_df, root_path, db, otu_dir)
    assert run_acc not in res

    # 5. Test with mixed count files (one empty, one not)
    count_file.write_text("asv\tcount\nASV1\t10\n")  # non-empty
    count_file_2 = count_dir / "test2_18S-V4_counts.tsv"
    count_file_2.touch()  # empty

    # It should NOT skip the run because count_file is non-empty
    res = get_asv_dict(runs_df, root_path, db, otu_dir)
    assert run_acc in res
    assert len(res[run_acc]) == 1
    assert res[run_acc].iloc[0]["asv"] == "ASV1"


def test_get_closedref_dict_empty_files(tmp_path: Path):
    """
    Test get_closedref_dict to ensure it correctly skips runs with empty krona files.
    """
    # Setup
    run_acc = "ERR123456"
    runs_df = pd.DataFrame([{"run": run_acc, "status": "all_results"}])
    db = "SILVA-SSU"
    otu_dir = tmp_path / "otu"
    otu_dir.mkdir()

    root_path = tmp_path / "results"
    run_dir = root_path / run_acc
    tax_summary_dir = run_dir / "taxonomy-summary" / db
    tax_summary_dir.mkdir(parents=True)

    kronatxt_file = tax_summary_dir / "test.krona.txt"
    kronatxt_file.touch()  # Empty

    res = get_closedref_dict(runs_df, root_path, db, otu_dir)
    assert run_acc not in res

    # Test non-empty krona file
    kronatxt_file.write_text("10\td__Bacteria\tp__Phylum\tc__Class\to__Order\tf__Family\tg__Genus\ts__Species\n")

    # Use the fixture OTU file
    fixture_otu = Path("tests/fixtures/refdb_otus/SILVA-SSU.otu")
    shutil.copy(fixture_otu, otu_dir / "SILVA-SSU.otu")

    res = get_closedref_dict(runs_df, root_path, db, otu_dir)
    assert run_acc in res
    assert res[run_acc].iloc[0]["count"] == 10


def test_generate_dwcready_summaries_single_run_empty_asv(tmp_path: Path):
    """
    Test generate_dwcready_summaries when there is a single run but ASV parsing
    returns an empty dict due to empty critical files.

    Mimics real pipeline structure with minimal required files only.
    """
    run_acc = "SRR3612607"
    # Runs table
    runs_df = pd.DataFrame([{"run": run_acc, "status": "all_results"}])
    runs_csv = tmp_path / "runs.csv"
    runs_df.to_csv(runs_csv, index=False, header=False)
    # Root structure
    root_path = tmp_path / "results"
    run_dir = root_path / run_acc
    db = "DADA2-SILVA"
    tax_summary_dir = run_dir / "taxonomy-summary" / db
    asv_dir = run_dir / "asv"
    count_dir = asv_dir / "16S-V4"
    tax_summary_dir.mkdir(parents=True)
    count_dir.mkdir(parents=True)
    # Empty mseq file → should cause get_asv_dict to skip run
    (tax_summary_dir / f"{run_acc}_{db}.mseq").touch()
    # Other required files exist but are valid
    (asv_dir / f"{run_acc}_{db}_asv_tax.tsv").write_text("ASV\tKingdom\nASV1\tBacteria\n")
    (asv_dir / f"{run_acc}_asv_seqs.fasta").write_text(">ASV1\nATGC\n")
    (count_dir / f"{run_acc}_16S-V4_asv_read_counts.tsv").write_text("asv\tcount\nASV1\t10\n")
    # OTU dir (required by function)
    otu_dir = tmp_path / "otu"
    otu_dir.mkdir()
    fixture_otu = Path("tests/fixtures/refdb_otus/SILVA-SSU.otu")
    shutil.copy(fixture_otu, otu_dir / "SILVA-SSU.otu")
    # Closed-reference files
    db_closed = "SILVA-SSU"
    # Closedref minimal structure (empty krona file → should be skipped)
    closedref_dir = root_path / run_acc / "taxonomy-summary" / db_closed
    closedref_dir.mkdir(parents=True, exist_ok=True)
    krona_file = closedref_dir / "empty.txt"
    krona_file.touch()  # empty → triggers skip in get_closedref_dict
    # Output dir
    output_dir = tmp_path / "output"
    # Run script
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "summarise",
            "-r",
            str(runs_csv),
            "-a",
            str(root_path),
            "-otu",
            str(otu_dir),
            "-p",
            str(output_dir),
        ],
    )
    # --- Assertions ---
    # Since get_asv_dict returns empty dict, no ASV-based outputs should be created
    assert result.exit_code == 0
    output_files = list(output_dir.glob("*"))
    assert output_files == [], "Expected no output summaries when ASV dict and closed-reference dict are empty"
