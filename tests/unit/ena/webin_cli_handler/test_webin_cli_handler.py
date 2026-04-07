import os
import subprocess
import time


timestamp = int(time.time())
timestamp_genomes = int(time.time())
webin_version = os.getenv("WEBIN_CLI_VERSION")


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
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Assigned accessions: ERZ" in result.stderr
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
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
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
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Assigned accessions: ERZ" in result.stderr
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
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
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
            "--webin-cli-jar",
            f"webin-cli-{webin_version}.jar",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Assigned accessions: ERZ" in result.stderr
        assert "Successfully submitted object for the first time on TEST server" in result.stderr
        assert "first_genome.manifest" in result.stderr
        assert "second_genome.manifest" in result.stderr
        assert (manifest_dir / f"webin-cli.test_{timestamp_genomes}_g1.report").exists(), "webin-cli.report for g1 not found after submission"
        assert (manifest_dir / f"webin-cli.test_{timestamp_genomes}_g2.report").exists(), "webin-cli.report for g2 not found after submission"
        with open(manifest_dir / f"webin-cli.test_{timestamp_genomes}_g1.report") as f:
            webin_report_lines = "\n".join(f.readlines())
        assert "ERZ" in webin_report_lines, f"No ERZ found in webin-cli.test_{timestamp_genomes}_g1.report, got {webin_report_lines}"
        assert os.path.exists("ena_accessions.tsv"), "ena_accessions.tsv not found after submission"
        with open("ena_accessions.tsv") as f:
            accession_lines = [line for line in f if line.strip()]
        assert len(accession_lines) == 3, (
            f"Expected 3 lines (header and 2 genomes) in ena_accessions.tsv, got {len(accession_lines)}: {accession_lines}"
        )
