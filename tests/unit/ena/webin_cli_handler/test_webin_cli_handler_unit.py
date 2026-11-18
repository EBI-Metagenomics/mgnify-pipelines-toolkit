import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from mgnify_pipelines_toolkit.ena.webin_cli_handler import (
    ManifestValidationError,
    REPORT_FILE,
    check_manifest,
    WebinAbsConsentError,
    WebinCLIExecutionError,
    WebinCredentialsError,
    check_report,
    check_submission_status,
    main,
    parse_accession_from_report,
    parse_manifest,
    run_webin_cli,
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

    def test_reads_valid_bams(self, tmp_path):
        bam = tmp_path / "reads.bam"
        bam.touch()
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\tRUN1\nBAM\t{bam}\n")
        alias, _ = check_manifest(str(manifest), None, "reads")
        assert alias == "RUN1"

    def test_reads_missing_fastq_field_raises(self, tmp_path):
        manifest = tmp_path / "run.manifest"
        manifest.write_text("NAME\tRUN1\nSTUDY\tPRJ1\n")
        with pytest.raises(ManifestValidationError, match="Neither FASTQ nor BAM field found"):
            check_manifest(str(manifest), None, "reads")

    def test_reads_fastq_file_not_found_raises(self, tmp_path):
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\tRUN1\nFASTQ\t{tmp_path / 'missing.fq.gz'}\n")
        with pytest.raises(ManifestValidationError, match="Data file not found"):
            check_manifest(str(manifest), None, "reads")

    def test_reads_bam_file_not_found_raises(self, tmp_path):
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\tRUN1\nBAM\t{tmp_path / 'missing.bam'}\n")
        with pytest.raises(ManifestValidationError, match="Data file not found"):
            check_manifest(str(manifest), None, "reads")

    def test_reads_one_of_paired_fastq_missing_raises(self, tmp_path):
        fq1 = tmp_path / "reads_1.fq.gz"
        fq1.touch()
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\tRUN1\nFASTQ\t{fq1}\nFASTQ\t{tmp_path / 'missing_2.fq.gz'}\n")
        with pytest.raises(ManifestValidationError, match="Data file not found"):
            check_manifest(str(manifest), None, "reads")

    def test_bam_and_fastq_both_present_raises(self, tmp_path):
        fq = tmp_path / "reads.fq.gz"
        fq.touch()
        bam = tmp_path / "reads.bam"
        bam.touch()
        manifest = tmp_path / "run.manifest"
        manifest.write_text(f"NAME\tRUN1\nBAM\t{bam}\nFASTQ\t{fq}\n")
        with pytest.raises(ManifestValidationError, match="Both FASTQ and BAM fields found"):
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


# Reports produced by webin-cli 11.0.0 (paths shortened). Validate and genome reports were captured from
# ENA test server runs; reads reports are constructed in the same format.
VALIDATE_REPORT = "2026-11-18T11:36:33 INFO : Submission(s) validated successfully.\n"

VALIDATE_FAILED_REPORT = (
    "2026-11-18T11:35:49 ERROR: Submission validation failed because of a user error. "
    "Manifest name : ASM1. See reports for details : out/genome/ASM1/validate.\n"
)

GENOME_FIRST_SUBMIT_TEST_REPORT = """\
2026-11-18T11:36:43 INFO : Submission has not been validated previously.
2026-11-18T11:36:44 INFO : Submission(s) validated successfully.
2026-11-18T11:36:44 INFO : Connecting to FTP server : webin2.ebi.ac.uk
2026-11-18T11:36:46 INFO : Uploading file: /data/bin.fa.gz
2026-11-18T11:36:46 INFO : Files have been successfully uploaded.
2026-11-18T11:36:47 INFO : The submission has been completed successfully. The following analysis accession was assigned to the submission: ERZ30663524
2026-11-18T11:36:47 INFO : This submission is a TEST submission and will be discarded within 24 hours
"""

GENOME_RESUBMIT_TEST_REPORT = """\
2026-11-18T11:36:49 INFO : Submission has not been validated previously.
2026-11-18T11:36:50 INFO : Submission(s) validated successfully.
2026-11-18T11:36:50 INFO : Connecting to FTP server : webin2.ebi.ac.uk
2026-11-18T11:36:51 INFO : Uploading file: /data/bin.fa.gz
2026-11-18T11:36:51 INFO : Files have been successfully uploaded.
2026-11-18T11:36:51 INFO : This submission is a TEST submission and will be discarded within 24 hours
2026-11-18T11:36:51 ERROR: In analysis, alias: "webin-genome-ASM1". The object being added already exists in the submission account with accession: "ERZ30663524". The submission has failed because of a system error.
uk.ac.ebi.ena.webin.cli.WebinCliException: In analysis, alias: "webin-genome-ASM1". The object being added already exists in the submission account with accession: "ERZ30663524". The submission has failed because of a system error.
Caused by: uk.ac.ebi.ena.webin.cli.WebinCliException: In analysis, alias: "webin-genome-ASM1". The object being added already exists in the submission account with accession: "ERZ30663524".
\t... 11 common frames omitted
2026-11-18T11:36:51 ERROR: Some or all submissions failed. Please see application logs.
uk.ac.ebi.ena.webin.cli.WebinCliException: Some or all submissions failed. Please see application logs.
"""

# Live server reports have the same format, without the TEST submission line
GENOME_FIRST_SUBMIT_LIVE_REPORT = GENOME_FIRST_SUBMIT_TEST_REPORT.replace(
    "2026-11-18T11:36:47 INFO : This submission is a TEST submission and will be discarded within 24 hours\n", ""
)
GENOME_RESUBMIT_LIVE_REPORT = GENOME_RESUBMIT_TEST_REPORT.replace(
    "2026-11-18T11:36:51 INFO : This submission is a TEST submission and will be discarded within 24 hours\n", ""
)

READS_FIRST_SUBMIT_TEST_REPORT = """\
2026-11-18T11:36:44 INFO : Submission(s) validated successfully.
2026-11-18T11:36:46 INFO : Files have been successfully uploaded.
2026-11-18T11:36:47 INFO : The submission has been completed successfully. The following experiment accession was assigned to the submission: ERX111111
2026-11-18T11:36:47 INFO : The submission has been completed successfully. The following run accession was assigned to the submission: ERR654321
2026-11-18T11:36:47 INFO : This submission is a TEST submission and will be discarded within 24 hours
"""

READS_RESUBMIT_TEST_REPORT = """\
2026-11-18T11:36:50 INFO : Submission(s) validated successfully.
2026-11-18T11:36:51 INFO : This submission is a TEST submission and will be discarded within 24 hours
2026-11-18T11:36:51 ERROR: In experiment, alias: "webin-reads-RUN1". The object being added already exists in the submission account with accession: "ERX111111".In run, alias: "webin-reads-RUN1". The object being added already exists in the submission account with accession: "ERR654321". The submission has failed because of a system error.
2026-11-18T11:36:51 ERROR: Some or all submissions failed. Please see application logs.
"""

SUBMIT_NO_ACCESSION_REPORT = (
    "2026-11-18T11:36:47 INFO : The submission has been completed successfully. "
    "No accession was assigned to the analysis submission. Please contact the helpdesk.\n"
)


class TestCheckReport:
    def _write_report(self, tmp_path, text):
        (tmp_path / REPORT_FILE).write_text(text)

    # --- Validate mode ---

    def test_validate_success(self, tmp_path):
        self._write_report(tmp_path, VALIDATE_REPORT)
        success, is_resub = check_report(str(tmp_path), "validate", True, "genome")
        assert success is True
        assert is_resub is False

    def test_validate_failure(self, tmp_path):
        self._write_report(tmp_path, VALIDATE_FAILED_REPORT)
        success, is_resub = check_report(str(tmp_path), "validate", True, "genome")
        assert success is False
        assert is_resub is False

    # --- Submit + genome (test and live servers report the same way since webin-cli 10.0.0) ---

    @pytest.mark.parametrize(
        "report, test",
        [(GENOME_FIRST_SUBMIT_TEST_REPORT, True), (GENOME_FIRST_SUBMIT_LIVE_REPORT, False)],
        ids=["test_server", "live_server"],
    )
    def test_submit_genome_first_time(self, tmp_path, report, test):
        self._write_report(tmp_path, report)
        success, is_resub = check_report(str(tmp_path), "submit", test, "genome")
        assert success is True
        assert is_resub is False

    @pytest.mark.parametrize(
        "report, test",
        [(GENOME_RESUBMIT_TEST_REPORT, True), (GENOME_RESUBMIT_LIVE_REPORT, False)],
        ids=["test_server", "live_server"],
    )
    def test_submit_genome_resubmission(self, tmp_path, report, test):
        self._write_report(tmp_path, report)
        success, is_resub = check_report(str(tmp_path), "submit", test, "genome")
        assert success is True
        assert is_resub is True

    def test_submit_genome_no_accession_fails(self, tmp_path):
        self._write_report(tmp_path, SUBMIT_NO_ACCESSION_REPORT)
        success, _ = check_report(str(tmp_path), "submit", True, "genome")
        assert success is False

    # --- Submit + reads ---

    def test_submit_reads_first_time(self, tmp_path):
        self._write_report(tmp_path, READS_FIRST_SUBMIT_TEST_REPORT)
        success, is_resub = check_report(str(tmp_path), "submit", True, "reads")
        assert success is True
        assert is_resub is False

    def test_submit_reads_resubmission(self, tmp_path):
        self._write_report(tmp_path, READS_RESUBMIT_TEST_REPORT)
        success, is_resub = check_report(str(tmp_path), "submit", True, "reads")
        assert success is True
        assert is_resub is True

    def test_submit_reads_missing_both_accessions(self, tmp_path):
        self._write_report(tmp_path, SUBMIT_NO_ACCESSION_REPORT)
        success, _ = check_report(str(tmp_path), "submit", True, "reads")
        assert success is False

    def test_submit_reads_only_run_accession_fails(self, tmp_path):
        # Needs both experiment (ERX) and run (ERR)
        self._write_report(tmp_path, READS_FIRST_SUBMIT_TEST_REPORT.replace("ERX111111", "-"))
        success, _ = check_report(str(tmp_path), "submit", True, "reads")
        assert success is False

    def test_submit_reads_only_experiment_accession_fails(self, tmp_path):
        self._write_report(tmp_path, READS_FIRST_SUBMIT_TEST_REPORT.replace("ERR654321", "-"))
        success, _ = check_report(str(tmp_path), "submit", True, "reads")
        assert success is False

    # --- Report file absent ---

    def test_report_file_missing_returns_false(self, tmp_path):
        success, is_resub = check_report(str(tmp_path), "submit", False, "genome")
        assert success is False
        assert is_resub is False


@pytest.mark.parametrize("test", [True, False], ids=["test_server", "live_server"])
@pytest.mark.parametrize(
    "report_text, expected_success, expected_resub",
    [
        (GENOME_FIRST_SUBMIT_TEST_REPORT, True, False),
        (GENOME_RESUBMIT_TEST_REPORT, True, True),
        (READS_RESUBMIT_TEST_REPORT, True, True),
        (VALIDATE_FAILED_REPORT, False, False),
        ("", False, False),
    ],
    ids=[
        "first_time_success",
        "genome_resubmission",
        "reads_resubmission",
        "failure",
        "empty_report",
    ],
)
def test_check_submission_status(report_text, expected_success, expected_resub, test):
    success, is_resub = check_submission_status(report_text, test)
    assert success is expected_success
    assert is_resub is expected_resub


# ---------------------------------------------------------------------------
# run_webin_cli retry and error handling
# ---------------------------------------------------------------------------


def _completed(returncode, stdout=""):
    return subprocess.CompletedProcess(args=[], returncode=returncode, stdout=stdout, stderr="")


class TestRunWebinCli:
    def _run(self, tmp_path, results, retries=3):
        with patch("subprocess.run", side_effect=results) as mock_run, patch("time.sleep") as mock_sleep:
            outcome = run_webin_cli(
                manifest="x.manifest",
                context="genome",
                webin="Webin-99999",
                mode="submit",
                test=True,
                output_dir=str(tmp_path),
                jar="/fake/webin-cli.jar",
                retries=retries,
                retry_delay=1,
            )
        return outcome, mock_run, mock_sleep

    def test_success_first_attempt(self, tmp_path):
        result, mock_run, _ = self._run(tmp_path, [_completed(0)])
        assert result.returncode == 0
        assert mock_run.call_count == 1

    def test_system_error_is_retried_then_succeeds(self, tmp_path):
        result, mock_run, mock_sleep = self._run(tmp_path, [_completed(1, "ERROR: Connection reset"), _completed(0)])
        assert result.returncode == 0
        assert mock_run.call_count == 2
        mock_sleep.assert_called_once_with(1)

    def test_system_error_fails_after_all_retries(self, tmp_path):
        with pytest.raises(WebinCLIExecutionError, match="after 3 attempts"):
            self._run(tmp_path, [_completed(1, "ERROR: Connection reset")] * 3)

    @pytest.mark.parametrize(
        "returncode, output",
        [
            (2, "ERROR: Unknown sample ERS32871415 or the sample cannot be referenced by your submission account."),
            (3, VALIDATE_FAILED_REPORT),
        ],
        ids=["user_error", "validation_error"],
    )
    def test_user_and_validation_errors_are_not_retried(self, tmp_path, returncode, output):
        with patch("subprocess.run", return_value=_completed(returncode, output)) as mock_run, patch("time.sleep") as mock_sleep:
            with pytest.raises(WebinCLIExecutionError, match=f"exit code {returncode}"):
                run_webin_cli("x.manifest", "genome", "Webin-99999", "submit", True, str(tmp_path), jar="/fake/webin-cli.jar")
        assert mock_run.call_count == 1
        mock_sleep.assert_not_called()

    def test_invalid_credentials_not_retried(self, tmp_path):
        output = "ERROR: Invalid submission account user name or password. Please try enclosing your password in single quotes."
        with patch("subprocess.run", return_value=_completed(2, output)) as mock_run, patch("time.sleep"):
            with pytest.raises(WebinCredentialsError):
                run_webin_cli("x.manifest", "genome", "Webin-99999", "submit", True, str(tmp_path), jar="/fake/webin-cli.jar")
        assert mock_run.call_count == 1

    def test_abs_consent_not_retried(self, tmp_path):
        output = (
            "ERROR: Failed submission: Incomplete Access and Benefit Sharing (ABS) acknowledgement. "
            "Please login to Webin Submissions Portal to complete."
        )
        with patch("subprocess.run", return_value=_completed(1, output)) as mock_run, patch("time.sleep"):
            with pytest.raises(WebinAbsConsentError):
                run_webin_cli("x.manifest", "genome", "Webin-99999", "submit", True, str(tmp_path), jar="/fake/webin-cli.jar")
        assert mock_run.call_count == 1

    def test_already_exists_returns_result(self, tmp_path):
        result, mock_run, _ = self._run(tmp_path, [_completed(1, GENOME_RESUBMIT_TEST_REPORT)])
        assert result.returncode == 1
        assert mock_run.call_count == 1


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
        report = VALIDATE_REPORT
        with patch("sys.argv", argv), patch("subprocess.run", side_effect=_make_subprocess_side_effect(report)):
            assert main() == 0

    def test_validate_reads_success(self, tmp_path, monkeypatch):
        manifest = self._reads_manifest(tmp_path)
        monkeypatch.setenv("ENA_WEBIN", _FAKE_WEBIN)
        monkeypatch.setenv("ENA_WEBIN_PASSWORD", _FAKE_PASSWORD)
        argv = self._argv(manifest, "reads", "validate", tmp_path / "out")
        report = VALIDATE_REPORT
        with patch("sys.argv", argv), patch("subprocess.run", side_effect=_make_subprocess_side_effect(report)):
            assert main() == 0

    def test_validate_failure_returns_1(self, tmp_path, monkeypatch):
        manifest = self._genome_manifest(tmp_path)
        monkeypatch.setenv("ENA_WEBIN", _FAKE_WEBIN)
        monkeypatch.setenv("ENA_WEBIN_PASSWORD", _FAKE_PASSWORD)
        argv = self._argv(manifest, "genome", "validate", tmp_path / "out")
        report = VALIDATE_FAILED_REPORT
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
        report = GENOME_RESUBMIT_TEST_REPORT
        with patch("sys.argv", argv), patch("subprocess.run", side_effect=_make_subprocess_side_effect(report)):
            assert main() == 0
        content = accessions_file.read_text()
        assert "ASM1" in content
        assert "ERZ30663524" in content

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
        report = READS_FIRST_SUBMIT_TEST_REPORT
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
        report = VALIDATE_FAILED_REPORT
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
