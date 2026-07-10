import gzip
from pathlib import Path

from click.testing import CliRunner
from mgnify_pipelines_toolkit.analysis.assembly.generate_gaf import main


FIXTURE_DIR = Path("tests/fixtures/goslim")


def test_generate_gaf_accepts_gzipped_interproscan_input(tmp_path):
    source = FIXTURE_DIR / "ips_out.tsv"
    gzipped_input = tmp_path / "ips_out.tsv.gz"
    with open(source, "rb") as source_fh, gzip.open(gzipped_input, "wb") as destination_fh:
        destination_fh.write(source_fh.read())

    output_prefix = tmp_path / "ERRTESTING"
    result = CliRunner().invoke(main, ["-i", str(gzipped_input), "-o", str(output_prefix)])
    assert result.exit_code == 0, result.output

    output = Path(f"{output_prefix}_ips_annotations.gaf")
    assert output.exists()
    assert "EMG\tGO:0042773\tGO\t\tGO:0042773" in output.read_text()
