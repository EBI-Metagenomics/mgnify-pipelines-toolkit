import gzip
import sys
from pathlib import Path

from mgnify_pipelines_toolkit.analysis.assembly.generate_gaf import main


FIXTURE_DIR = Path("tests/fixtures/goslim")


def test_generate_gaf_accepts_gzipped_interproscan_input(tmp_path, monkeypatch):
    source = FIXTURE_DIR / "ips_out.tsv"
    gzipped_input = tmp_path / "ips_out.tsv.gz"
    with open(source, "rb") as source_fh, gzip.open(gzipped_input, "wb") as destination_fh:
        destination_fh.write(source_fh.read())

    output_prefix = tmp_path / "ERRTESTING"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "generate_gaf",
            "-i",
            str(gzipped_input),
            "-o",
            str(output_prefix),
        ],
    )

    main()

    output = Path(f"{output_prefix}_ips_annotations.gaf")
    assert output.exists()
    assert "EMG\tGO:0042773\tGO\t\tGO:0042773" in output.read_text()
