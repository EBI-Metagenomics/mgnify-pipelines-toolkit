import os
import subprocess
import time

timestamp = int(time.time())
timestamp_genomes = int(time.time())
webin_version = os.getenv("WEBIN_CLI_VERSION")

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
