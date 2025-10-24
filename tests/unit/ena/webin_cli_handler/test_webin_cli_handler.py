from datetime import datetime as dt
import subprocess

timestamp = str(int(dt.timestamp(dt.now())))

class Tests:
    def test_validate_assembly_upload(tmp_path):
        command = [
            "python",
            "mgnify_pipelines_toolkit/ena/webin_cli_handler.py",
            "-c",
            "genome",
            "-m",
            "tests/fixtures/webin_cli_handler/assembly.manifest",
            "--mode",
            "validate",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        print(result.stdout)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Submission validation succeeded" in result.stdout
        assert "Submission/validation done for" in result.stdout

    def test_validate_genome_upload(tmp_path):
        command = [
            "python",
            "mgnify_pipelines_toolkit/ena/webin_cli_handler.py",
            "-c",
            "genome",
            "-m",
            "tests/fixtures/webin_cli_handler/genome.manifest",
            "--mode",
            "validate",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        print(result.stdout)
        assert result.returncode == 0, f"Run failed: {result.stderr}"
        assert "Submission validation succeeded" in result.stdout
        assert "Submission/validation done for" in result.stdout

    def test_submit_assembly_upload_first_time_test_server(tmp_path):
        print(tmp_path)
        test_manifest = "tests/fixtures/webin_cli_handler/assembly_test.manifest"
        # create a new manifest with another alias for unique submission
        with open("new_assembly.manifest", 'w') as file_out, open(test_manifest, 'r') as file_in:
            for line in file_in:
                if 'ASSEMBLYNAME' in line:
                    line = f"ASSEMBLYNAME\ttest_{timestamp}\n"
                file_out.write(line)
        command = [
            "python",
            "mgnify_pipelines_toolkit/ena/webin_cli_handler.py",
            "-c",
            "genome",
            "-m",
            "new_assembly.manifest",
            "--mode",
            "submit",
            "--test"
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"

    def test_submit_assembly_upload_second_time(tmp_path):
        command = [
            "python",
            "mgnify_pipelines_toolkit/ena/webin_cli_handler.py",
            "-c",
            "genome",
            "-m",
            "tests/fixtures/webin_cli_handler/assembly.manifest",
            "--mode",
            "submit",
            "--test"
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"

    def test_submit_genome_upload_first_time(tmp_path):
        command = [
            "python",
            "mgnify_pipelines_toolkit/ena/webin_cli_handler.py",
            "-c",
            "genome",
            "-m",
            "tests/fixtures/webin_cli_handler/genome.manifest",
            "--mode",
            "submit",
            "--test"
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"

    def test_submit_genome_upload_second_time(tmp_path):
        command = [
            "python",
            "mgnify_pipelines_toolkit/ena/webin_cli_handler.py",
            "-c",
            "genome",
            "-m",
            "tests/fixtures/webin_cli_handler/genome.manifest",
            "--mode",
            "submit",
            "--test"
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"