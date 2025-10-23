import subprocess


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

    def test_submit_assembly_upload_first_time(tmp_path):
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