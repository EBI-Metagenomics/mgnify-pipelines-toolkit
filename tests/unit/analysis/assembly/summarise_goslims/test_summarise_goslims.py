import gzip
from pathlib import Path

from click.testing import CliRunner
from mgnify_pipelines_toolkit.analysis.assembly.summarise_goslims import main


FIXTURE_DIR = Path("tests/fixtures/goslim")


def gzip_fixture(source: Path, destination: Path) -> None:
    with open(source, "rb") as source_fh, gzip.open(destination, "wb") as destination_fh:
        destination_fh.write(source_fh.read())


def run_summarise_goslims(runner: CliRunner, output: Path, **inputs: Path) -> None:
    result = runner.invoke(
        main,
        [
            "-go",
            str(inputs["go_obo"]),
            "-gb",
            str(inputs["go_banding"]),
            "-gaf",
            str(inputs["gaf_input"]),
            "-i",
            str(inputs["ips_input"]),
            "-o",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output


def test_summarise_goslims_accepts_gzipped_inputs(tmp_path):
    plain_inputs = {
        "go_obo": FIXTURE_DIR / "go-dummy.obo",
        "go_banding": FIXTURE_DIR / "goslim_banding_2024.txt",
        "gaf_input": FIXTURE_DIR / "ERRTESTING_ips_annotations.gaf",
        "ips_input": FIXTURE_DIR / "ips_out.tsv",
    }
    runner = CliRunner()
    plain_output = tmp_path / "plain_terms"
    run_summarise_goslims(runner, plain_output, **plain_inputs)

    gzipped_inputs = {}
    for name, source in plain_inputs.items():
        destination = tmp_path / f"{source.name}.gz"
        gzip_fixture(source, destination)
        gzipped_inputs[name] = destination

    gzipped_output = tmp_path / "gzipped_terms"
    run_summarise_goslims(runner, gzipped_output, **gzipped_inputs)

    assert gzipped_output.read_text() == plain_output.read_text()
    assert Path(f"{gzipped_output}_slim").read_text() == Path(f"{plain_output}_slim").read_text()
