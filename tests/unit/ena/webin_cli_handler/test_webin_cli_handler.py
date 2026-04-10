import os
import shutil
import subprocess
import time
from pathlib import Path

import pytest


timestamp = int(time.time())
timestamp_genomes = int(time.time())
webin_version = os.getenv("WEBIN_CLI_VERSION")
test_output_dir = Path("test_output")


@pytest.mark.webin_cli
class TestWebinCliHandler:
    def test_validate_assembly_upload(self, tmp_path):
        command = [
            "webin_cli_handler",
            "-c",
            "genome",
            "-m",
            "tests/fixtures/webin_cli_handler/assembly.manifest",
            "--mode",
            "validate",
            "--outdir",
            test_output_dir,
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Submission validation succeeded" in result.stderr
        assert "Submission/validation done for" in result.stderr

    def test_validate_genome_upload(self, tmp_path):
        command = [
            "webin_cli_handler",
            "-c",
            "genome",
            "-m",
            "tests/fixtures/webin_cli_handler/genome.manifest",
            "--mode",
            "validate",
            "--outdir",
            test_output_dir,
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Submission validation succeeded" in result.stderr
        assert "Submission/validation done for" in result.stderr

    def test_validate_genome_upload_relative_fasta_without_fasta_dir_fails(self, tmp_path):
        test_manifest = "tests/fixtures/webin_cli_handler/genome.manifest"
        fasta_dir = tmp_path / "fasta"
        fasta_dir.mkdir()
        shutil.copyfile("tests/fixtures/webin_cli_handler/test.fasta.gz", fasta_dir / "test.fasta.gz")

        manifest_path = tmp_path / "relative_no_fasta_dir.manifest"
        with open(manifest_path, "w") as file_out, open(test_manifest, "r") as file_in:
            for line in file_in:
                if "ASSEMBLYNAME" in line:
                    line = f"ASSEMBLYNAME\trelative_fail_{timestamp_genomes}\n"
                if "FASTA" in line:
                    line = "FASTA\ttest.fasta.gz\n"
                file_out.write(line)

        command = [
            "webin_cli_handler",
            "-c",
            "genome",
            "-m",
            str(manifest_path),
            "--mode",
            "validate",
            "--outdir",
            test_output_dir,
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode != 0, "Validation unexpectedly succeeded without --fasta-dir for a relative FASTA path"
        assert 'Could not read data file: "test.fasta.gz"' in result.stderr
        assert "test.fasta.gz" in result.stderr

    def test_validate_genome_upload_relative_fasta_with_fasta_dir_succeeds(self, tmp_path):
        test_manifest = "tests/fixtures/webin_cli_handler/genome.manifest"
        fasta_dir = tmp_path / "fasta"
        fasta_dir.mkdir()
        shutil.copyfile("tests/fixtures/webin_cli_handler/test.fasta.gz", fasta_dir / "test.fasta.gz")

        manifest_path = tmp_path / "relative_with_fasta_dir.manifest"
        with open(manifest_path, "w") as file_out, open(test_manifest, "r") as file_in:
            for line in file_in:
                if "ASSEMBLYNAME" in line:
                    line = f"ASSEMBLYNAME\trelative_success_{timestamp_genomes}\n"
                if "FASTA" in line:
                    line = "FASTA\ttest.fasta.gz\n"
                file_out.write(line)

        command = [
            "webin_cli_handler",
            "-c",
            "genome",
            "-m",
            str(manifest_path),
            "--mode",
            "validate",
            "--fasta-dir",
            str(fasta_dir),
            "--outdir",
            test_output_dir,
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Submission validation succeeded" in result.stderr
        assert "Submission/validation done for" in result.stderr

    def test_submit_assembly_upload_first_time_test_server(self, tmp_path):
        test_manifest = "tests/fixtures/webin_cli_handler/assembly.manifest"
        # create a new manifest with another alias for unique submission
        with open("new_assembly.manifest", "w") as file_out, open(test_manifest, "r") as file_in:
            for line in file_in:
                if "ASSEMBLYNAME" in line:
                    line = f"ASSEMBLYNAME\ttest_{timestamp}_a\n"
                file_out.write(line)
        command = [
            "webin_cli_handler",
            "-c",
            "genome",
            "-m",
            "new_assembly.manifest",
            "--mode",
            "submit",
            "--test",
            "--outdir",
            test_output_dir,
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Assigned accessions for" in result.stderr
        assert "Successfully submitted object for the first time on TEST server" in result.stderr
        assert "Submission/validation done" in result.stderr

    def test_submit_assembly_upload_second_time_test_server(self, tmp_path):
        test_manifest = "tests/fixtures/webin_cli_handler/assembly.manifest"
        # create a new manifest with same alias for second submission
        with open("repeat_assembly.manifest", "w") as file_out, open(test_manifest, "r") as file_in:
            for line in file_in:
                if "ASSEMBLYNAME" in line:
                    line = f"ASSEMBLYNAME\ttest_{timestamp}_a\n"
                file_out.write(line)
        command = [
            "webin_cli_handler",
            "-c",
            "genome",
            "-m",
            "repeat_assembly.manifest",
            "--mode",
            "submit",
            "--test",
            "--outdir",
            test_output_dir,
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
            "--debug",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Command completed successfully" in result.stderr
        assert "Submitted object already exists on TEST server" in result.stderr

    def test_submit_genome_upload_first_time_test_server(self, tmp_path):
        test_manifest = "tests/fixtures/webin_cli_handler/genome.manifest"
        # create a new manifest with timestamp alias for the first submission
        with open("new_genome.manifest", "w") as file_out, open(test_manifest, "r") as file_in:
            for line in file_in:
                if "ASSEMBLYNAME" in line:
                    line = f"ASSEMBLYNAME\ttest_{timestamp_genomes}_g\n"
                file_out.write(line)
        command = [
            "webin_cli_handler",
            "-c",
            "genome",
            "-m",
            "new_genome.manifest",
            "--mode",
            "submit",
            "--test",
            "--outdir",
            test_output_dir,
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Assigned accessions for" in result.stderr
        assert "Successfully submitted object for the first time on TEST server" in result.stderr
        assert "Submission/validation done for new_genome.manifest" in result.stderr

    def test_submit_genome_upload_second_time_test_server(self, tmp_path):
        test_manifest = "tests/fixtures/webin_cli_handler/genome.manifest"
        # create a second manifest with same timestamp alias as first submission
        with open("repeat_genome.manifest", "w") as file_out, open(test_manifest, "r") as file_in:
            for line in file_in:
                if "ASSEMBLYNAME" in line:
                    line = f"ASSEMBLYNAME\ttest_{timestamp_genomes}_g\n"
                file_out.write(line)
        command = [
            "webin_cli_handler",
            "-c",
            "genome",
            "-m",
            "repeat_genome.manifest",
            "--mode",
            "submit",
            "--test",
            "--outdir",
            test_output_dir,
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
            "--debug",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        print(result.stderr)
        assert "Submitted object already exists on TEST server" in result.stderr
        assert "Submission/validation done for repeat_genome.manifest" in result.stderr

    def test_submit_genome_upload_multiple_manifests_test_server(self, tmp_path):
        test_manifest = "tests/fixtures/webin_cli_handler/genome.manifest"
        manifest_dir = tmp_path / "manifests"
        manifest_dir.mkdir()
        # create 2 manifests in a directory for submission
        with open(manifest_dir / "first_genome.manifest", "w") as file_out, open(test_manifest, "r") as file_in:
            for line in file_in:
                if "ASSEMBLYNAME" in line:
                    line = f"ASSEMBLYNAME\ttest_{timestamp_genomes}_g1\n"
                file_out.write(line)
        with open(manifest_dir / "second_genome.manifest", "w") as file_out, open(test_manifest, "r") as file_in:
            for line in file_in:
                if "ASSEMBLYNAME" in line:
                    line = f"ASSEMBLYNAME\ttest_{timestamp_genomes}_g2\n"
                file_out.write(line)

        command = [
            "webin_cli_handler",
            "-c",
            "genome",
            "-m",
            str(manifest_dir),
            "--mode",
            "submit",
            "--test",
            "--outdir",
            test_output_dir,
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Assigned accessions for" in result.stderr
        assert "Successfully submitted object for the first time on TEST server" in result.stderr
        assert "first_genome.manifest" in result.stderr
        assert "second_genome.manifest" in result.stderr
        first_output_dirs = list(test_output_dir.glob(f"test_{timestamp_genomes}_g1_*"))
        second_output_dirs = list(test_output_dir.glob(f"test_{timestamp_genomes}_g2_*"))
        assert first_output_dirs, f"No timestamped output directory found for g1 in {test_output_dir}"
        assert second_output_dirs, f"No timestamped output directory found for g2 in {test_output_dir}"

        first_report = first_output_dirs[0] / "webin-cli.report"
        second_report = second_output_dirs[0] / "webin-cli.report"

        assert first_report.exists(), f"webin-cli.report for g1 not found after submission: {first_report}"
        assert second_report.exists(), f"webin-cli.report for g2 not found after submission: {second_report}"
        with open(first_report) as f:
            webin_report_lines = "\n".join(f.readlines())
        assert "ERZ" in webin_report_lines, f"No ERZ found in {first_report}, got {webin_report_lines}"
        assert os.path.exists("ena_accessions.tsv"), "ena_accessions.tsv not found after submission"
        with open("ena_accessions.tsv") as f:
            accession_lines = [line for line in f if line.strip()]
        assert len(accession_lines) == 3, (
            f"Expected 3 lines (header and 2 genomes) in ena_accessions.tsv, got {len(accession_lines)}: {accession_lines}"
        )

    def test_submit_genome_upload_resume_skips_existing_aliases(self, tmp_path):
        test_manifest = "tests/fixtures/webin_cli_handler/genome.manifest"
        manifest_dir = tmp_path / "resume_manifests"
        manifest_dir.mkdir()
        output_accessions = tmp_path / "resume_accessions.tsv"
        outdir = tmp_path / "resume_output"
        skipped_alias = f"test_{timestamp_genomes}_resume_existing"
        new_alias = f"test_{timestamp_genomes}_resume_new"

        with open(manifest_dir / "existing_genome.manifest", "w") as file_out, open(test_manifest, "r") as file_in:
            for line in file_in:
                if "ASSEMBLYNAME" in line:
                    line = f"ASSEMBLYNAME\t{skipped_alias}\n"
                file_out.write(line)

        with open(manifest_dir / "new_genome.manifest", "w") as file_out, open(test_manifest, "r") as file_in:
            for line in file_in:
                if "ASSEMBLYNAME" in line:
                    line = f"ASSEMBLYNAME\t{new_alias}\n"
                file_out.write(line)

        with open(output_accessions, "w") as file_out:
            file_out.write("alias\taccession\n")
            file_out.write(f"{skipped_alias}\tERZ123456\n")

        command = [
            "webin_cli_handler",
            "-c",
            "genome",
            "-m",
            str(manifest_dir),
            "--mode",
            "submit",
            "--test",
            "--resume",
            "--output-accessions",
            str(output_accessions),
            "--outdir",
            str(outdir),
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)

        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert f"Skipping {skipped_alias}: already present in accessions file with accession ERZ123456" in result.stderr
        assert "Submission/validation done for" in result.stderr

        skipped_output_dir = outdir / skipped_alias
        new_output_dirs = list(outdir.glob(f"{new_alias}_*"))
        skipped_output_dirs = list(outdir.glob(f"{skipped_alias}_*"))
        assert not skipped_output_dir.exists(), f"Skipped alias unexpectedly created output directory: {skipped_output_dir}"
        assert not skipped_output_dirs, f"Skipped alias unexpectedly created timestamped output directory: {skipped_output_dirs}"
        assert new_output_dirs, f"Expected timestamped output directory for resumed submission in: {outdir}"

        with open(output_accessions) as f:
            accession_lines = [line.strip() for line in f if line.strip()]

        assert len(accession_lines) == 3, f"Expected header and 2 genome lines in resume accessions file, got {accession_lines}"
        assert f"{skipped_alias}\tERZ123456" in accession_lines
        assert any(line.startswith(f"{new_alias}\tERZ") for line in accession_lines), accession_lines
