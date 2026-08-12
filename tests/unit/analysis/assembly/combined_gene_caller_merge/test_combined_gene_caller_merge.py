import gzip
from pathlib import Path

import pytest
from Bio import SeqIO

from mgnify_pipelines_toolkit.analysis.assembly.combined_gene_caller_merge import (
    output_fasta_files,
    parse_cmsearch_output,
    parse_fgsrs_output,
    parse_gff,
    parse_pyrodigal_output,
)


FIXTURES = Path("tests/fixtures/combined_gene_caller_merge")


@pytest.fixture(scope="session")
def compressed_fixtures(tmp_path_factory):
    compressed_dir = tmp_path_factory.mktemp("compressed-fixtures")
    names = [
        "input_prodigal.out",
        "input_fgs.out",
        "input_pyrodigal.gff",
        "cmsearch.all.tblout.deoverlapped",
        "input_prodigal.faa",
        "input_prodigal.ffn",
    ]
    compressed = {}
    for name in names:
        compressed_path = compressed_dir / f"{name}.gz"
        with gzip.open(compressed_path, "wt") as output:
            output.write((FIXTURES / name).read_text())
        compressed[name] = compressed_path
    return compressed


def test_annotation_inputs_can_be_gzipped(compressed_fixtures):
    assert parse_pyrodigal_output(compressed_fixtures["input_prodigal.out"])
    assert parse_fgsrs_output(compressed_fixtures["input_fgs.out"])
    assert parse_gff(compressed_fixtures["input_pyrodigal.gff"])
    assert parse_cmsearch_output(
        compressed_fixtures["cmsearch.all.tblout.deoverlapped"]
    )


def test_fasta_inputs_can_be_gzipped(tmp_path, compressed_fixtures):
    predictions = parse_pyrodigal_output(FIXTURES / "input_prodigal.out")
    output_faa = tmp_path / "output.faa"
    output_ffn = tmp_path / "output.ffn"

    output_fasta_files(
        {"Pyrodigal": predictions},
        {
            "Pyrodigal": {
                "proteins": compressed_fixtures["input_prodigal.faa"],
                "transcripts": compressed_fixtures["input_prodigal.ffn"],
            }
        },
        output_faa,
        output_ffn,
    )

    assert list(SeqIO.parse(output_faa, "fasta"))
    assert list(SeqIO.parse(output_ffn, "fasta"))
