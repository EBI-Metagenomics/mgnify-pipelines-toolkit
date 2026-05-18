import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from mgnify_pipelines_toolkit.ena.webin_cli_handler import (
    ManifestValidationError,
    REPORT_FILE,
    check_manifest,
    check_report,
    check_submission_status_live,
    check_submission_status_test,
    main,
    parse_accession_from_report,
    parse_manifest,
)


class TestParseManifest:
    def test_single_fastq_stored_as_list(self, tmp_path):
        manifest = tmp_path / "run.manifest"
        manifest.write_text("NAME\tRUN1\nFASTQ\tfile1.fq.gz\n")
        result = parse_manifest(manifest)
        assert result["NAME"] == "RUN1"
        assert result["FASTQ"] == ["file1.fq.gz"]

    def test_two_fastq_lines_accumulate(self, tmp_path):
        manifest = tmp_path / "run.manifest"
        manifest.write_text("NAME\tRUN1\nFASTQ\tfile1.fq.gz\nFASTQ\tfile2.fq.gz\n")
        result = parse_manifest(manifest)
        assert result["FASTQ"] == ["file1.fq.gz", "file2.fq.gz"]

    def test_genome_fasta_is_string_not_list(self, tmp_path):
        manifest = tmp_path / "genome.manifest"
        manifest.write_text("ASSEMBLYNAME\tASM1\nFASTA\tasm.fasta.gz\n")
        result = parse_manifest(manifest)
        assert result["ASSEMBLYNAME"] == "ASM1"
        assert result["FASTA"] == "asm.fasta.gz"
        assert isinstance(result["FASTA"], str)

    def test_empty_lines_are_skipped(self, tmp_path):
        manifest = tmp_path / "blank.manifest"
        manifest.write_text("STUDY\tPRJ1\n\n\nSAMPLE\tSAMEA1\n")
        result = parse_manifest(manifest)
        assert result == {"STUDY": "PRJ1", "SAMPLE": "SAMEA1"}

    def test_malformed_line_without_whitespace_is_skipped(self, tmp_path):
        manifest = tmp_path / "bad.manifest"
        manifest.write_text("BADLINE\nSTUDY\tPRJ1\n")
        result = parse_manifest(manifest)
        assert result == {"STUDY": "PRJ1"}

    def test_empty_manifest_returns_empty_dict(self, tmp_path):
        manifest = tmp_path / "empty.manifest"
        manifest.write_text("")
        assert parse_manifest(manifest) == {}


class TestCheckManifest:
    # --- Genome path (regression) ---

    def test_genome_valid_absolute_fasta(self, tmp_path):
        fasta = tmp_path / "asm.fasta.gz"
        fasta.touch()
        manifest = tmp_path / "genome.manifest"
        manifest.write_text(f"ASSEMBLYNAME\tASM1\nFASTA\t{fasta}\n")
        alias, path = check_manifest(str(manifest), None, "genome")
        assert alias == "ASM1"
        assert path == str(manifest)

    def test_genome_missing_fasta_field_raises(self, tmp_path):
        manifest = tmp_path / "genome.manifest"
        manifest.write_text("ASSEMBLYNAME\tASM1\nSTUDY\tPRJ1\n")
        with pytest.raises(ManifestValidationError, match="FASTA field not found"):
            check_manifest(str(manifest), None, "genome")

    def test_genome_fasta_file_not_found_raises(self, tmp_path):
        manifest = tmp_path / "genome.manifest"
        manifest.write_text(f"ASSEMBLYNAME\tASM1\nFASTA\t{tmp_path / 'missing.fasta.gz'}\n")
        with pytest.raises(ManifestValidationError, match="Data file not found"):
            check_manifest(str(manifest), None, "genome")

    def test_genome_relative_fasta_without_fasta_root_skips_check(self, tmp_path):
        manifest = tmp_path / "genome.manifest"
        manifest.write_text("ASSEMBLYNAME\tASM1\nFASTA\trelative.fasta.gz\n")
        # Relative path + no fasta_root → file existence not checked
        alias, _ = check_manifest(str(manifest), None, "genome")
        assert alias == "ASM1"

    def test_genome_relative_fasta_with_fasta_root_validates_file(self, tmp_path):
        fasta_dir = tmp_path / "fastas"
        fasta_dir.mkdir()
        (fasta_dir / "asm.fasta.gz").touch()
        manifest = tmp_path / "genome.manifest"
        manifest.write_text("ASSEMBLYNAME\tASM1\nFASTA\tasm.fasta.gz\n")
        alias, _ = check_manifest(str(manifest), str(fasta_dir), "genome")
        assert alias == "ASM1"

    def test_genome_relative_fasta_with_fasta_root_missing_file_raises(self, tmp_path):
        fasta_dir = tmp_path / "fastas"
        fasta_dir.mkdir()
        manifest = tmp_path / "genome.manifest"
        manifest.write_text("ASSEMBLYNAME\tASM1\nFASTA\tmissing.fasta.gz\n")
        with pytest.raises(ManifestValidationError, match="Data file not found"):
            check_manifest(str(manifest), str(fasta_dir), "genome")

    # --- Reads path (new) ---

    def test_reads_valid_single_fastq(self, tmp_path):
        fq = tmp_path / "reads.fq.gz"
        fq.touch()
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\tRUN1\nFASTQ\t{fq}\n")
        alias, path = check_manifest(str(manifest), None, "reads")
        assert alias == "RUN1"
        assert path == str(manifest)

    def test_reads_valid_paired_fastq_both_files_exist(self, tmp_path):
        fq1 = tmp_path / "reads_1.fq.gz"
        fq2 = tmp_path / "reads_2.fq.gz"
        fq1.touch()
        fq2.touch()
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\tRUN1\nFASTQ\t{fq1}\nFASTQ\t{fq2}\n")
        alias, _ = check_manifest(str(manifest), None, "reads")
        assert alias == "RUN1"

    def test_reads_missing_fastq_field_raises(self, tmp_path):
        manifest = tmp_path / "run.manifest"
        manifest.write_text("NAME\tRUN1\nSTUDY\tPRJ1\n")
        with pytest.raises(ManifestValidationError, match="FASTQ field not found"):
            check_manifest(str(manifest), None, "reads")

    def test_reads_fastq_file_not_found_raises(self, tmp_path):
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\tRUN1\nFASTQ\t{tmp_path / 'missing.fq.gz'}\n")
        with pytest.raises(ManifestValidationError, match="Data file not found"):
            check_manifest(str(manifest), None, "reads")

    def test_reads_one_of_paired_fastq_missing_raises(self, tmp_path):
        fq1 = tmp_path / "reads_1.fq.gz"
        fq1.touch()
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\tRUN1\nFASTQ\t{fq1}\nFASTQ\t{tmp_path / 'missing_2.fq.gz'}\n")
        with pytest.raises(ManifestValidationError, match="Data file not found"):
            check_manifest(str(manifest), None, "reads")

    def test_reads_name_used_when_assemblyname_absent(self, tmp_path):
        fq = tmp_path / "reads.fq.gz"
        fq.touch()
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\tRUN1\nFASTQ\t{fq}\n")
        alias, _ = check_manifest(str(manifest), None, "reads")
        assert alias == "RUN1"

    def test_reads_assemblyname_takes_priority_over_name(self, tmp_path):
        fq = tmp_path / "reads.fq.gz"
        fq.touch()
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"ASSEMBLYNAME\tASM1\nNAME\tRUN1\nFASTQ\t{fq}\n")
        alias, _ = check_manifest(str(manifest), None, "reads")
        assert alias == "ASM1"

    # --- Shared ---

    def test_missing_name_and_assemblyname_raises(self, tmp_path):
        fq = tmp_path / "reads.fq.gz"
        fq.touch()
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"STUDY\tPRJ1\nFASTQ\t{fq}\n")
        with pytest.raises(ManifestValidationError, match="Neither ASSEMBLYNAME, nor NAME"):
            check_manifest(str(manifest), None, "reads")

    def test_nonexistent_manifest_raises(self, tmp_path):
        with pytest.raises(ManifestValidationError, match="does not exist"):
            check_manifest(str(tmp_path / "nonexistent.manifest"), None, "genome")


class TestParseAccessionFromReport:
    def _write_report(self, tmp_path, text):
        report = tmp_path / REPORT_FILE
        report.write_text(text)
        return report

    def test_assembly_accession(self, tmp_path):
        report = self._write_report(tmp_path, "Assembly accession: ERZ123456\nDone.")
        result = parse_accession_from_report(report)
        assert result == {"assembly": "ERZ123456", "experiment": None, "run": None}

    def test_run_accession(self, tmp_path):
        report = self._write_report(tmp_path, "Run accession: ERR654321\nDone.")
        result = parse_accession_from_report(report)
        assert result == {"assembly": None, "experiment": None, "run": "ERR654321"}

    def test_experiment_accession(self, tmp_path):
        report = self._write_report(tmp_path, "Experiment accession: ERX111111\nDone.")
        result = parse_accession_from_report(report)
        assert result == {"assembly": None, "experiment": "ERX111111", "run": None}

    def test_reads_submission_experiment_and_run_both_present(self, tmp_path):
        report = self._write_report(tmp_path, "ERX111111 ERR654321 submission done.")
        result = parse_accession_from_report(report)
        assert result["experiment"] == "ERX111111"
        assert result["run"] == "ERR654321"
        assert result["assembly"] is None

    def test_no_accessions_returns_all_none(self, tmp_path):
        report = self._write_report(tmp_path, "No accession assigned.")
        assert parse_accession_from_report(report) == {"assembly": None, "experiment": None, "run": None}

    def test_repeated_accession_is_not_ambiguous(self, tmp_path):
        report = self._write_report(tmp_path, "ERZ123456 ERZ123456 ERZ123456")
        result = parse_accession_from_report(report)
        assert result["assembly"] == "ERZ123456"

    def test_multiple_unique_assembly_accessions_yields_none(self, tmp_path):
        report = self._write_report(tmp_path, "ERZ111111 ERZ222222")
        result = parse_accession_from_report(report)
        assert result["assembly"] is None

    def test_ddbj_prefix_assembly(self, tmp_path):
        report = self._write_report(tmp_path, "DRZ123456")
        result = parse_accession_from_report(report)
        assert result["assembly"] == "DRZ123456"

    def test_ncbi_prefix_assembly(self, tmp_path):
        report = self._write_report(tmp_path, "SRZ123456")
        result = parse_accession_from_report(report)
        assert result["assembly"] == "SRZ123456"

    def test_short_accession_below_six_digits_not_matched(self, tmp_path):
        report = self._write_report(tmp_path, "ERZ12345")  # 5 digits — regex requires {6,}
        result = parse_accession_from_report(report)
        assert result == {"assembly": None, "experiment": None, "run": None}


class TestCheckReport:
    def _write_report(self, tmp_path, text):
        (tmp_path / REPORT_FILE).write_text(text)

    # --- Validate mode ---

    def test_validate_success(self, tmp_path):
        self._write_report(tmp_path, "Submission(s) validated successfully.")
        success, is_resub = check_report(str(tmp_path), "validate", False, "genome")
        assert success is True
        assert is_resub is False

    def test_validate_failure(self, tmp_path):
        self._write_report(tmp_path, "Error: validation failed.")
        success, is_resub = check_report(str(tmp_path), "validate", False, "genome")
        assert success is False
        assert is_resub is False

    # --- Submit + live + genome (regression) ---

    def test_submit_live_genome_first_time(self, tmp_path):
        self._write_report(tmp_path, "ERZ123456\nsubmission has been completed successfully")
        success, is_resub = check_report(str(tmp_path), "submit", False, "genome")
        assert success is True
        assert is_resub is False

    def test_submit_live_genome_resubmission(self, tmp_path):
        self._write_report(
            tmp_path,
            "object being added already exists in the submission account with accession ERZ123456",
        )
        success, is_resub = check_report(str(tmp_path), "submit", False, "genome")
        assert success is True
        assert is_resub is True

    def test_submit_live_genome_no_accession_fails(self, tmp_path):
        self._write_report(tmp_path, "submission has been completed successfully")
        success, _ = check_report(str(tmp_path), "submit", False, "genome")
        assert success is False

    # --- Submit + live + reads ---

    def test_submit_live_reads_first_time(self, tmp_path):
        self._write_report(
            tmp_path,
            "ERX111111 ERR654321\nsubmission has been completed successfully",
        )
        success, is_resub = check_report(str(tmp_path), "submit", False, "reads")
        assert success is True
        assert is_resub is False

    def test_submit_live_reads_missing_both_accessions(self, tmp_path):
        self._write_report(tmp_path, "submission has been completed successfully")
        success, _ = check_report(str(tmp_path), "submit", False, "reads")
        assert success is False

    def test_submit_live_reads_only_run_accession_fails(self, tmp_path):
        # Needs both experiment (ERX) and run (ERR)
        self._write_report(tmp_path, "ERR654321\nsubmission has been completed successfully")
        success, _ = check_report(str(tmp_path), "submit", False, "reads")
        assert success is False

    def test_submit_live_reads_only_experiment_accession_fails(self, tmp_path):
        self._write_report(tmp_path, "ERX111111\nsubmission has been completed successfully")
        success, _ = check_report(str(tmp_path), "submit", False, "reads")
        assert success is False

    # --- Submit + test + genome ---

    def test_submit_test_genome_first_time(self, tmp_path):
        self._write_report(
            tmp_path,
            "ERZ123456\nThis was a TEST submission(s).\nsubmission has been completed successfully",
        )
        success, is_resub = check_report(str(tmp_path), "submit", True, "genome")
        assert success is True
        assert is_resub is False

    # --- Submit + test + genome resubmission (regression: accession check must not block) ---

    def test_submit_test_genome_resubmission_minimal_report(self, tmp_path):
        """Test-server genome resubmission: minimal report with no accession must still succeed."""
        self._write_report(tmp_path, "This was a TEST submission(s).")
        success, is_resub = check_report(str(tmp_path), "submit", True, "genome")
        assert success is True
        assert is_resub is True

    def test_submit_test_reads_resubmission_minimal_report(self, tmp_path):
        """Test-server reads resubmission: minimal report with no ERX/ERR must still succeed."""
        self._write_report(tmp_path, "This was a TEST submission(s).")
        success, is_resub = check_report(str(tmp_path), "submit", True, "reads")
        assert success is True
        assert is_resub is True

    # --- Report file absent ---

    def test_report_file_missing_returns_false(self, tmp_path):
        success, is_resub = check_report(str(tmp_path), "submit", False, "genome")
        assert success is False
        assert is_resub is False


@pytest.mark.parametrize(
    "report_text, expected_success, expected_resub",
    [
        (
            "This was a TEST submission(s).\nsubmission has been completed successfully",
            True,
            False,
        ),
        ("This was a TEST submission(s).", True, True),
        ("This was a TEST submission(s).\n", True, True),
        (
            "This was a TEST submission(s).\nobject being added already exists in the submission account with accession ERZ123",
            True,
            True,
        ),
        ("Error: validation failed", False, False),
        ("", False, False),
    ],
    ids=[
        "first_time_success",
        "resubmission_exact_minimal",
        "resubmission_one_trailing_newline",
        "resubmission_explicit_exists_message",
        "failure_no_test_marker",
        "empty_report",
    ],
)
def test_check_submission_status_test(report_text, expected_success, expected_resub):
    success, is_resub = check_submission_status_test(report_text)
    assert success is expected_success
    assert is_resub is expected_resub


@pytest.mark.parametrize(
    "report_text, expected_success, expected_resub",
    [
        ("submission has been completed successfully", True, False),
        (
            "object being added already exists in the submission account with accession ERZ123456",
            True,
            True,
        ),
        ("Error: validation failed", False, False),
        ("", False, False),
    ],
    ids=[
        "first_time_success",
        "resubmission",
        "failure",
        "empty_report",
    ],
)
def test_check_submission_status_live(report_text, expected_success, expected_resub):
    success, is_resub = check_submission_status_live(report_text)
    assert success is expected_success
    assert is_resub is expected_resub


# ---------------------------------------------------------------------------
# Integration tests for main()
# ---------------------------------------------------------------------------

_FAKE_JAR = "/fake/webin-cli.jar"
_FAKE_WEBIN = "Webin-99999"
_FAKE_PASSWORD = "testpassword"


def _make_subprocess_side_effect(report_text, context=None, alias=None, mode=None):
    """Side-effect for subprocess.run that writes webin-cli artifacts to outputDir.

    For validate mode or resubmission reports only the report file is needed.
    For first-time submit, pass context/alias/mode to also create the result
    directory structure that check_result expects.
    """

    def side_effect(cmd, **_):
        output_dir = next(
            (Path(a.split("=", 1)[1]) for a in cmd if a.startswith("-outputDir=")),
            None,
        )
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / REPORT_FILE).write_text(report_text)
            if context and alias and mode:
                result_dir = output_dir / context / alias / mode
                result_dir.mkdir(parents=True, exist_ok=True)
                for fname in ["receipt.xml", "submission.xml", "webin-submission.xml"]:
                    (result_dir / fname).touch()
        result = MagicMock(spec=subprocess.CompletedProcess)
        result.returncode = 0
        result.stdout = ""
        result.stderr = ""
        return result

    return side_effect


class TestMain:
    def _argv(self, manifest_path, context, mode, outdir, extra=None):
        argv = [
            "webin_cli_handler",
            "-m",
            str(manifest_path),
            "-c",
            context,
            "--mode",
            mode,
            "--webin-cli-jar",
            _FAKE_JAR,
            "--outdir",
            str(outdir),
        ]
        if extra:
            argv.extend(extra)
        return argv

    def _genome_manifest(self, tmp_path, alias="ASM1"):
        fasta = tmp_path / "asm.fasta.gz"
        fasta.touch()
        manifest = tmp_path / "genome.manifest"
        manifest.write_text(f"ASSEMBLYNAME\t{alias}\nFASTA\t{fasta}\n")
        return manifest

    def _reads_manifest(self, tmp_path, alias="RUN1"):
        fq = tmp_path / "reads.fq.gz"
        fq.touch()
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\t{alias}\nFASTQ\t{fq}\n")
        return manifest

    # --- validate mode ---

    def test_validate_genome_success(self, tmp_path, monkeypatch):
        manifest = self._genome_manifest(tmp_path)
        monkeypatch.setenv("ENA_WEBIN", _FAKE_WEBIN)
        monkeypatch.setenv("ENA_WEBIN_PASSWORD", _FAKE_PASSWORD)
        argv = self._argv(manifest, "genome", "validate", tmp_path / "out")
        report = "Submission(s) validated successfully."
        with patch("sys.argv", argv), patch("subprocess.run", side_effect=_make_subprocess_side_effect(report)):
            assert main() == 0

    def test_validate_reads_success(self, tmp_path, monkeypatch):
        manifest = self._reads_manifest(tmp_path)
        monkeypatch.setenv("ENA_WEBIN", _FAKE_WEBIN)
        monkeypatch.setenv("ENA_WEBIN_PASSWORD", _FAKE_PASSWORD)
        argv = self._argv(manifest, "reads", "validate", tmp_path / "out")
        report = "Submission(s) validated successfully."
        with patch("sys.argv", argv), patch("subprocess.run", side_effect=_make_subprocess_side_effect(report)):
            assert main() == 0

    def test_validate_failure_returns_1(self, tmp_path, monkeypatch):
        manifest = self._genome_manifest(tmp_path)
        monkeypatch.setenv("ENA_WEBIN", _FAKE_WEBIN)
        monkeypatch.setenv("ENA_WEBIN_PASSWORD", _FAKE_PASSWORD)
        argv = self._argv(manifest, "genome", "validate", tmp_path / "out")
        report = "Error: validation failed."
        with patch("sys.argv", argv), patch("subprocess.run", side_effect=_make_subprocess_side_effect(report)):
            assert main() == 1

    # --- submit mode ---

    def test_submit_genome_resubmission(self, tmp_path, monkeypatch):
        manifest = self._genome_manifest(tmp_path)
        accessions_file = tmp_path / "accessions.tsv"
        monkeypatch.setenv("ENA_WEBIN", _FAKE_WEBIN)
        monkeypatch.setenv("ENA_WEBIN_PASSWORD", _FAKE_PASSWORD)
        argv = self._argv(
            manifest,
            "genome",
            "submit",
            tmp_path / "out",
            extra=["--output-accessions", str(accessions_file)],
        )
        # Resubmission report: accession in text + already-exists message
        report = "object being added already exists in the submission account with accession ERZ123456"
        with patch("sys.argv", argv), patch("subprocess.run", side_effect=_make_subprocess_side_effect(report)):
            assert main() == 0
        content = accessions_file.read_text()
        assert "ASM1" in content
        assert "ERZ123456" in content

    def test_submit_reads_first_time(self, tmp_path, monkeypatch):
        manifest = self._reads_manifest(tmp_path, alias="RUN1")
        accessions_file = tmp_path / "accessions.tsv"
        monkeypatch.setenv("ENA_WEBIN", _FAKE_WEBIN)
        monkeypatch.setenv("ENA_WEBIN_PASSWORD", _FAKE_PASSWORD)
        argv = self._argv(
            manifest,
            "reads",
            "submit",
            tmp_path / "out",
            extra=["--output-accessions", str(accessions_file)],
        )
        report = "ERX111111 ERR654321\nsubmission has been completed successfully"
        side_effect = _make_subprocess_side_effect(
            report,
            context="reads",
            alias="RUN1",
            mode="submit",
        )
        with patch("sys.argv", argv), patch("subprocess.run", side_effect=side_effect):
            assert main() == 0
        content = accessions_file.read_text()
        assert "RUN1" in content
        assert "ERX111111" in content
        assert "ERR654321" in content

    def test_submit_failure_returns_1(self, tmp_path, monkeypatch):
        manifest = self._genome_manifest(tmp_path)
        monkeypatch.setenv("ENA_WEBIN", _FAKE_WEBIN)
        monkeypatch.setenv("ENA_WEBIN_PASSWORD", _FAKE_PASSWORD)
        argv = self._argv(manifest, "genome", "submit", tmp_path / "out")
        report = "Error: submission failed."
        with patch("sys.argv", argv), patch("subprocess.run", side_effect=_make_subprocess_side_effect(report)):
            assert main() == 1

    # --- error handling ---

    def test_missing_credentials_returns_2(self, tmp_path, monkeypatch):
        manifest = self._genome_manifest(tmp_path)
        monkeypatch.delenv("ENA_WEBIN", raising=False)
        monkeypatch.delenv("ENA_WEBIN_PASSWORD", raising=False)
        argv = self._argv(manifest, "genome", "validate", tmp_path / "out")
        with patch("sys.argv", argv):
            assert main() == 2

    def test_invalid_manifest_returns_3(self, tmp_path, monkeypatch):
        manifest = tmp_path / "bad.manifest"
        manifest.write_text("STUDY\tPRJ1\n")  # no NAME or ASSEMBLYNAME
        monkeypatch.setenv("ENA_WEBIN", _FAKE_WEBIN)
        monkeypatch.setenv("ENA_WEBIN_PASSWORD", _FAKE_PASSWORD)
        argv = self._argv(manifest, "genome", "validate", tmp_path / "out")
        with patch("sys.argv", argv):
            assert main() == 3
