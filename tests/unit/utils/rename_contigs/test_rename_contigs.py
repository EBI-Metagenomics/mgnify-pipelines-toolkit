#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest
import tempfile
import shutil
from pathlib import Path


from mgnify_pipelines_toolkit.utils.rename_contigs import (
    _parse_fasta_header,
    read_mapping_file,
    rename_fasta,
    restore_fasta,
    rename_gff,
    rename_genbank,
    write_mapping_file,
)


@pytest.fixture
def fixtures_dir():
    """Return the path to the fixtures directory."""
    return Path(__file__).parent.parent.parent.parent / "fixtures" / "rename_contigs"


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    temp_path = tempfile.mkdtemp()
    yield Path(temp_path)
    shutil.rmtree(temp_path)


class TestParseFastaHeader:
    """Test cases for _parse_fasta_header function."""

    def test_simple_header(self):
        """Test parsing a simple FASTA header."""
        seq_name, short_name, metadata, viral_id = _parse_fasta_header(">simple_seq")
        assert seq_name == "simple_seq"
        assert short_name == "simple_seq"
        assert metadata == []
        assert viral_id is None

    def test_header_with_description(self):
        """Test parsing a header with description."""
        seq_name, short_name, metadata, viral_id = _parse_fasta_header(
            ">NODE_1_length_1000 extra description"
        )
        assert seq_name == "NODE_1_length_1000 extra description"
        assert short_name == "NODE_1_length_1000"
        assert metadata == []
        assert viral_id is None

    def test_header_with_viral_id(self):
        """Test parsing a header with viral identifier."""
        seq_name, short_name, metadata, viral_id = _parse_fasta_header(
            ">viral_seq|viral_id_001"
        )
        assert seq_name == "viral_seq"
        assert short_name == "viral_seq"
        assert "viral_id_001" in metadata
        assert viral_id == "viral_id_001"

    def test_header_with_phage_circular(self):
        """Test parsing a header with phage-circular metadata."""
        seq_name, short_name, metadata, viral_id = _parse_fasta_header(
            ">viral_seq|viral_id_002|phage-circular"
        )
        assert seq_name == "viral_seq"
        assert "phage-circular" in metadata
        assert "viral_id_002" in metadata

    def test_header_with_prophage(self):
        """Test parsing a header with prophage metadata."""
        seq_name, short_name, metadata, viral_id = _parse_fasta_header(
            ">viral_seq|prophage-100:500"
        )
        assert seq_name == "viral_seq"
        assert "prophage-100:500" in metadata

    def test_header_without_prefix(self):
        """Test parsing a header without '>' prefix."""
        seq_name, short_name, metadata, viral_id = _parse_fasta_header("simple_seq")
        assert seq_name == "simple_seq"
        assert short_name == "simple_seq"


class TestMappingFile:
    """Test cases for mapping file operations."""

    def test_write_and_read_mapping(self, temp_dir):
        """Test writing and reading a mapping file."""
        mapping_file = temp_dir / "test_mapping.tsv"

        # Create test mapping
        test_mapping = {
            "contig1": ("NODE_1_length_1000", "NODE_1_length_1000"),
            "contig2": ("NODE_2_length_500 description", "NODE_2_length_500"),
        }

        # Write mapping
        write_mapping_file(test_mapping, str(mapping_file))

        # Read mapping back
        read_mapping = read_mapping_file(
            str(mapping_file), from_col="original", to_col="renamed"
        )

        assert "NODE_1_length_1000" in read_mapping
        assert read_mapping["NODE_1_length_1000"] == "contig1"
        assert read_mapping["NODE_2_length_500 description"] == "contig2"

    def test_read_mapping_with_custom_columns(self, temp_dir):
        """Test reading mapping file with custom columns."""
        mapping_file = temp_dir / "test_mapping.tsv"

        # Create test mapping
        test_mapping = {
            "contig1": ("NODE_1_length_1000", "NODE_1"),
            "contig2": ("NODE_2_length_500", "NODE_2"),
        }

        write_mapping_file(test_mapping, str(mapping_file))

        # Read with different columns
        read_mapping = read_mapping_file(
            str(mapping_file), from_col="renamed", to_col="short"
        )

        assert read_mapping["contig1"] == "NODE_1"
        assert read_mapping["contig2"] == "NODE_2"


class TestRenameFasta:
    """Test cases for FASTA renaming."""

    def test_rename_fasta_simple(self, fixtures_dir, temp_dir):
        """Test basic FASTA renaming."""
        input_file = fixtures_dir / "test_input.fasta"
        output_file = temp_dir / "output.fasta"

        mapping = rename_fasta(
            str(input_file),
            str(output_file),
            prefix="seq",
            use_pyfastx=False
        )

        # Check mapping was created
        assert len(mapping) == 3
        assert "seq1" in mapping
        assert "seq2" in mapping
        assert "seq3" in mapping

        # Check output file exists
        assert output_file.exists()

        # Check output content
        with open(output_file) as f:
            content = f.read()
            assert ">seq1" in content
            assert ">seq2" in content
            assert ">seq3" in content
            assert "NODE_1" not in content

    def test_rename_fasta_with_custom_prefix(self, fixtures_dir, temp_dir):
        """Test FASTA renaming with custom prefix."""
        input_file = fixtures_dir / "test_input.fasta"
        output_file = temp_dir / "output.fasta"

        mapping = rename_fasta(
            str(input_file),
            str(output_file),
            prefix="contig_",
            use_pyfastx=False
        )

        # Check mapping uses custom prefix
        assert "contig_1" in mapping or "contig1" in mapping

    def test_rename_fasta_preserve_metadata(self, fixtures_dir, temp_dir):
        """Test FASTA renaming with metadata preservation."""
        input_file = fixtures_dir / "test_viral.fasta"
        output_file = temp_dir / "output.fasta"

        mapping = rename_fasta(
            str(input_file),
            str(output_file),
            prefix="seq",
            preserve_metadata=True,
            use_pyfastx=False
        )

        # Check output preserves metadata
        with open(output_file) as f:
            content = f.read()
            assert "viral_id" in content
            assert "phage-circular" in content
            assert "prophage-100:500" in content

    def test_rename_fasta_with_existing_mapping(self, fixtures_dir, temp_dir):
        """Test FASTA renaming with existing mapping."""
        input_file = fixtures_dir / "test_input.fasta"
        output_file = temp_dir / "output.fasta"

        # Create custom mapping
        custom_mapping = {
            "NODE_1_length_1000_cov_10.5": "custom1",
            "NODE_2_length_500_cov_5.2": "custom2",
            "NODE_3_length_750_cov_8.1": "custom3",
        }

        mapping = rename_fasta(
            str(input_file),
            str(output_file),
            mapping=custom_mapping,
            use_pyfastx=False
        )

        # Check custom mapping was used
        with open(output_file) as f:
            content = f.read()
            assert ">custom1" in content
            assert ">custom2" in content
            assert ">custom3" in content


class TestRestoreFasta:
    """Test cases for FASTA restoration."""

    def test_restore_fasta(self, fixtures_dir, temp_dir):
        """Test restoring original FASTA names."""
        input_file = fixtures_dir / "test_input.fasta"
        renamed_file = temp_dir / "renamed.fasta"
        restored_file = temp_dir / "restored.fasta"

        # First rename
        mapping = rename_fasta(
            str(input_file),
            str(renamed_file),
            prefix="seq",
            use_pyfastx=False
        )

        # Create reverse mapping
        reverse_mapping = {new: old for new, (old, _) in mapping.items()}

        # Restore
        restore_fasta(str(renamed_file), str(restored_file), reverse_mapping)

        # Check restored content matches original
        with open(input_file) as f1, open(restored_file) as f2:
            original_headers = [line for line in f1 if line.startswith(">")]
            restored_headers = [line for line in f2 if line.startswith(">")]

            # Check we have the same number of sequences
            assert len(original_headers) == len(restored_headers)

    def test_restore_fasta_preserve_metadata(self, fixtures_dir, temp_dir):
        """Test restoring FASTA with metadata preservation."""
        input_file = fixtures_dir / "test_viral.fasta"
        renamed_file = temp_dir / "renamed.fasta"
        restored_file = temp_dir / "restored.fasta"

        # First rename with metadata
        mapping = rename_fasta(
            str(input_file),
            str(renamed_file),
            prefix="seq",
            preserve_metadata=True,
            use_pyfastx=False
        )

        # Create reverse mapping
        reverse_mapping = {new: old for new, (old, _) in mapping.items()}

        # Restore with metadata preservation
        restore_fasta(
            str(renamed_file),
            str(restored_file),
            reverse_mapping,
            preserve_metadata=True
        )

        # Check metadata is preserved
        with open(restored_file) as f:
            content = f.read()
            assert "viral_id" in content or "viral_seq" in content

    def test_restore_fasta_with_virsorter_prefix(self, fixtures_dir, temp_dir):
        """Test FASTA renaming with metadata preservation."""
        input_file = fixtures_dir / "virsorter_renamed.fasta"
        mapfile = fixtures_dir / "virsorter_map.tsv"
        output_file = temp_dir / "restored.fasta"
        mapping = read_mapping_file(mapfile, from_col='temporary', to_col='short')

        restore_fasta(
            str(input_file),
            str(output_file),
            mapping,
            preserve_metadata=True
        )
        expected = [
            "ERZ1.1",
            "ERZ1.2|phage-circular",
            "ERZ1.3|prophage-21696:135184",
            "ERZ1.3|prophage-100:500",
            "ERZ1.4|prophage-100:500|prophage-500:700",
            "ERZ1.5|phage-circular|prophage-100:500|prophage-500:700",
        ]

        obtained = []
        with open(output_file, "r") as restored_file:
            lines = 0
            records = 0
            for line in restored_file:
                if line.startswith(">"):
                    obtained = line.replace(">", "").strip()
                    records += 1
                assert expected[records-1] == obtained
                lines += 1
            assert lines == 24
            assert records == 6


class TestRenameGff:
    """Test cases for GFF renaming."""

    def test_rename_gff(self, fixtures_dir, temp_dir):
        """Test GFF renaming."""
        input_file = fixtures_dir / "test_input.gff"
        output_file = temp_dir / "output.gff"

        # Create mapping
        mapping = {
            "NODE_1_length_1000_cov_10.5": "contig1",
            "NODE_2_length_500_cov_5.2": "contig2",
            "NODE_3_length_750_cov_8.1": "contig3",
        }

        rename_gff(str(input_file), str(output_file), mapping)

        # Check output
        with open(output_file) as f:
            content = f.read()
            assert "contig1" in content
            assert "contig2" in content
            assert "contig3" in content
            assert "NODE_1" not in content

            # Check sequence-region directives were updated
            assert "##sequence-region contig1" in content

            # Check FASTA section was updated
            assert ">contig1" in content


class TestRenameGenbank:
    """Test cases for GenBank renaming."""

    def test_rename_genbank(self, fixtures_dir, temp_dir):
        """Test GenBank renaming."""
        input_file = fixtures_dir / "test_input.gbk"
        output_file = temp_dir / "output.gbk"

        # Create mapping
        mapping = {
            "NODE_1_length_1000_cov_10.5": "contig1",
            "NODE_2_length_500_cov_5.2": "contig2",
        }

        rename_genbank(str(input_file), str(output_file), mapping)

        # Check output
        with open(output_file) as f:
            content = f.read()
            assert "LOCUS       contig1" in content
            assert "LOCUS       contig2" in content
            assert "NODE_1" not in content


class TestIntegration:
    """Integration tests for the rename_contigs module."""

    def test_rename_multiple_file_types(self, fixtures_dir, temp_dir):
        """Test renaming FASTA, GFF together."""
        fasta_in = fixtures_dir / "test_input.fasta"
        gff_in = fixtures_dir / "test_input.gff"

        fasta_out = temp_dir / "test_input.fasta"
        gff_out = temp_dir / "test_input.gff"

        # Rename FASTA first
        mapping = rename_fasta(
            str(fasta_in),
            str(fasta_out),
            prefix="seq",
            use_pyfastx=False
        )

        # Convert to old -> new mapping
        old_to_new = {old: new for new, (old, _) in mapping.items()}

        # Rename GFF
        rename_gff(str(gff_in), str(gff_out), old_to_new)

        # Check both files exist and have consistent naming
        assert fasta_out.exists()
        assert gff_out.exists()

        with open(fasta_out) as f:
            fasta_content = f.read()

        with open(gff_out) as f:
            gff_content = f.read()

        # Check that sequence names match between files
        assert ">seq1" in fasta_content
        assert "seq1" in gff_content

    def test_full_rename_restore_cycle(self, fixtures_dir, temp_dir):
        """Test a full rename and restore cycle."""
        input_file = fixtures_dir / "test_input.fasta"
        renamed_file = temp_dir / "renamed.fasta"
        restored_file = temp_dir / "restored.fasta"
        mapping_file = temp_dir / "mapping.tsv"

        # Rename
        mapping = rename_fasta(
            str(input_file),
            str(renamed_file),
            prefix="seq",
            use_pyfastx=False
        )

        # Write mapping
        write_mapping_file(mapping, str(mapping_file))

        # Read mapping for restoration
        restore_mapping = read_mapping_file(
            str(mapping_file), from_col="renamed", to_col="original"
        )

        # Restore
        restore_fasta(str(renamed_file), str(restored_file), restore_mapping)

        # Compare original and restored (sequence content should match)
        with open(input_file) as f:
            original_lines = [line for line in f if not line.startswith(">")]

        with open(restored_file) as f:
            restored_lines = [line for line in f if not line.startswith(">")]

        # Sequences should be identical
        assert original_lines == restored_lines


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
