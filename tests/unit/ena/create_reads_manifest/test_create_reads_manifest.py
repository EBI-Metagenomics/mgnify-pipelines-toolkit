import re
from pathlib import Path

from mgnify_pipelines_toolkit.ena.create_reads_manifest import (
    generate_reads_manifest,
    main,
)


def _read_fields(manifest_path: Path) -> dict[str, list[str]]:
    """Parse a manifest into {KEY: [value, ...]} (all values as lists for uniformity)."""
    result: dict[str, list[str]] = {}
    for line in manifest_path.read_text().splitlines():
        if "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        result.setdefault(key, []).append(value)
    return result


class TestGenerateReadsManifest:
    def test_required_fields_present(self, tmp_path):
        out = tmp_path / "run.manifest"
        generate_reads_manifest(
            output_path=out,
            study="PRJEB1",
            sample="SAMEA1",
            name="sample1",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2000",
            library_source="METAGENOMIC",
            library_selection="RANDOM",
            library_strategy="WGS",
            fastq_files=["reads_1.fq.gz"],
        )
        fields = _read_fields(out)
        assert fields["STUDY"] == ["PRJEB1"]
        assert fields["SAMPLE"] == ["SAMEA1"]
        assert fields["PLATFORM"] == ["ILLUMINA"]
        assert fields["INSTRUMENT"] == ["Illumina HiSeq 2000"]
        assert fields["LIBRARY_SOURCE"] == ["METAGENOMIC"]
        assert fields["LIBRARY_SELECTION"] == ["RANDOM"]
        assert fields["LIBRARY_STRATEGY"] == ["WGS"]
        assert fields["FASTQ"] == ["reads_1.fq.gz"]

    def test_name_has_timestamp_suffix(self, tmp_path):
        out = tmp_path / "run.manifest"
        generate_reads_manifest(
            output_path=out,
            study="PRJEB1",
            sample="SAMEA1",
            name="mysample",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2000",
            library_source="METAGENOMIC",
            library_selection="RANDOM",
            library_strategy="WGS",
            fastq_files=["r.fq.gz"],
        )
        fields = _read_fields(out)
        name_value = fields["NAME"][0]
        assert name_value.startswith("mysample_")
        assert re.match(r"mysample_\d{14}$", name_value), f"Unexpected NAME: {name_value}"

    def test_paired_end_two_fastq_lines(self, tmp_path):
        out = tmp_path / "run.manifest"
        generate_reads_manifest(
            output_path=out,
            study="PRJEB1",
            sample="SAMEA1",
            name="pe",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2000",
            library_source="METAGENOMIC",
            library_selection="RANDOM",
            library_strategy="WGS",
            fastq_files=["reads_1.fq.gz", "reads_2.fq.gz"],
        )
        fields = _read_fields(out)
        assert fields["FASTQ"] == ["reads_1.fq.gz", "reads_2.fq.gz"]

    def test_fastq_paths_written_as_is(self, tmp_path):
        """Paths must not be stripped to basename — caller controls what is written."""
        out = tmp_path / "run.manifest"
        full_path = "/data/project/reads/sample_R1.fq.gz"
        generate_reads_manifest(
            output_path=out,
            study="PRJEB1",
            sample="SAMEA1",
            name="s",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2000",
            library_source="METAGENOMIC",
            library_selection="RANDOM",
            library_strategy="WGS",
            fastq_files=[full_path],
        )
        fields = _read_fields(out)
        assert fields["FASTQ"] == [full_path]

    def test_optional_insert_size(self, tmp_path):
        out = tmp_path / "run.manifest"
        generate_reads_manifest(
            output_path=out,
            study="PRJEB1",
            sample="SAMEA1",
            name="s",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2000",
            library_source="METAGENOMIC",
            library_selection="RANDOM",
            library_strategy="WGS",
            fastq_files=["r.fq.gz"],
            insert_size=300,
        )
        fields = _read_fields(out)
        assert fields["INSERT_SIZE"] == ["300"]

    def test_optional_library_name(self, tmp_path):
        out = tmp_path / "run.manifest"
        generate_reads_manifest(
            output_path=out,
            study="PRJEB1",
            sample="SAMEA1",
            name="s",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2000",
            library_source="METAGENOMIC",
            library_selection="RANDOM",
            library_strategy="WGS",
            fastq_files=["r.fq.gz"],
            library_name="MyLib",
        )
        fields = _read_fields(out)
        assert fields["LIBRARY_NAME"] == ["MyLib"]

    def test_optional_description(self, tmp_path):
        out = tmp_path / "run.manifest"
        generate_reads_manifest(
            output_path=out,
            study="PRJEB1",
            sample="SAMEA1",
            name="s",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2000",
            library_source="METAGENOMIC",
            library_selection="RANDOM",
            library_strategy="WGS",
            fastq_files=["r.fq.gz"],
            description="test run",
        )
        fields = _read_fields(out)
        assert fields["DESCRIPTION"] == ["test run"]

    def test_optional_fields_absent_when_not_provided(self, tmp_path):
        out = tmp_path / "run.manifest"
        generate_reads_manifest(
            output_path=out,
            study="PRJEB1",
            sample="SAMEA1",
            name="s",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2000",
            library_source="METAGENOMIC",
            library_selection="RANDOM",
            library_strategy="WGS",
            fastq_files=["r.fq.gz"],
        )
        fields = _read_fields(out)
        assert "INSERT_SIZE" not in fields
        assert "LIBRARY_NAME" not in fields
        assert "DESCRIPTION" not in fields

    def test_lines_are_tab_delimited(self, tmp_path):
        out = tmp_path / "run.manifest"
        generate_reads_manifest(
            output_path=out,
            study="PRJEB1",
            sample="SAMEA1",
            name="s",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2000",
            library_source="METAGENOMIC",
            library_selection="RANDOM",
            library_strategy="WGS",
            fastq_files=["r.fq.gz"],
        )
        for line in out.read_text().splitlines():
            if line:
                assert "\t" in line, f"Line is not tab-delimited: {line!r}"
                assert line.count("\t") == 1, f"Line has multiple tabs: {line!r}"


class TestCLI:
    def _base_args(self, tmp_path) -> list[str]:
        return [
            "--study",
            "PRJEB1",
            "--sample",
            "SAMEA1",
            "--name",
            "cli_test",
            "--platform",
            "ILLUMINA",
            "--instrument",
            "Illumina HiSeq 2000",
            "--library-source",
            "METAGENOMIC",
            "--library-selection",
            "RANDOM",
            "--library-strategy",
            "WGS",
            "--fastq",
            "reads_1.fq.gz",
            "--output",
            str(tmp_path / "out.manifest"),
        ]

    def test_cli_creates_manifest(self, tmp_path):
        argv = self._base_args(tmp_path)
        assert main(argv) == 0
        assert (tmp_path / "out.manifest").exists()

    def test_cli_single_end(self, tmp_path):
        argv = self._base_args(tmp_path)
        assert main(argv) == 0
        fields = _read_fields(tmp_path / "out.manifest")
        assert len(fields["FASTQ"]) == 1

    def test_cli_paired_end(self, tmp_path):
        argv = self._base_args(tmp_path) + ["--fastq", "reads_2.fq.gz"]
        assert main(argv) == 0
        fields = _read_fields(tmp_path / "out.manifest")
        assert len(fields["FASTQ"]) == 2

    def test_cli_optional_insert_size(self, tmp_path):
        argv = self._base_args(tmp_path) + ["--insert-size", "400"]
        assert main(argv) == 0
        fields = _read_fields(tmp_path / "out.manifest")
        assert fields["INSERT_SIZE"] == ["400"]

    def test_cli_optional_library_name(self, tmp_path):
        argv = self._base_args(tmp_path) + ["--library-name", "LibA"]
        assert main(argv) == 0
        fields = _read_fields(tmp_path / "out.manifest")
        assert fields["LIBRARY_NAME"] == ["LibA"]

    def test_cli_optional_description(self, tmp_path):
        argv = self._base_args(tmp_path) + ["--description", "some description"]
        assert main(argv) == 0
        fields = _read_fields(tmp_path / "out.manifest")
        assert fields["DESCRIPTION"] == ["some description"]

    def test_cli_full_path_fastq_written_verbatim(self, tmp_path):
        """Absolute paths must not be stripped to basename."""
        argv = [
            "--study",
            "PRJEB1",
            "--sample",
            "SAMEA1",
            "--name",
            "s",
            "--platform",
            "ILLUMINA",
            "--instrument",
            "Illumina HiSeq 2000",
            "--library-source",
            "METAGENOMIC",
            "--library-selection",
            "RANDOM",
            "--library-strategy",
            "WGS",
            "--fastq",
            "/full/path/to/reads.fq.gz",
            "--output",
            str(tmp_path / "out.manifest"),
        ]
        assert main(argv) == 0
        fields = _read_fields(tmp_path / "out.manifest")
        assert fields["FASTQ"] == ["/full/path/to/reads.fq.gz"]
