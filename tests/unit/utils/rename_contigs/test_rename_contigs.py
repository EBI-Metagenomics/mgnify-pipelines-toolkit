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


from mgnify_pipelines_toolkit.utils.rename_contig import (
    _parse_fasta_header,
    read_mapping_file,
    rename_fasta,
    rename_fasta_by_size,
    rename_genbank,
    rename_gff,
    restore_fasta,
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
            prefix="seq"
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
            prefix="contig_"
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
            formatter_mode="generic"
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

        rename_fasta(
            str(input_file),
            str(output_file),
            mapping=custom_mapping
        )

        # Check custom mapping was used
        with open(output_file) as f:
            content = f.read()
            assert ">custom1" in content
            assert ">custom2" in content
            assert ">custom3" in content

    def test_rename_fasta_separate_by_size(self, fixtures_dir, temp_dir):
        """Test FASTA to create 1k.fasta, 5k.fasta, etc (used in mobilome pipeline)"""
        input_file = fixtures_dir / "test_input.fasta"

        rename_fasta_by_size(
            str(input_file),
            output_prefix='assembly',
            prefix='contig_'
        )
        assert Path(temp_dir / "assembly_1kb_contigs.fasta").exists()


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
            prefix="seq"
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
            formatter_mode="generic"
        )

        # Create reverse mapping
        reverse_mapping = {new: old for new, (old, _) in mapping.items()}

        # Restore with metadata preservation
        restore_fasta(
            str(renamed_file),
            str(restored_file),
            reverse_mapping,
            formatter_mode="generic"
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
            formatter_mode="virify"
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
            prefix="seq"
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
            prefix="seq"
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


class TestMetaviraversePipeline:
    """Test cases specific to Metaviraverse pipeline requirements."""

    def test_simple_viral_identifier_preservation(self, temp_dir):
        """Test that viral identifiers are preserved in metaviraverse format."""
        # Create test input
        input_file = temp_dir / "metaviraverse_input.fasta"
        output_file = temp_dir / "metaviraverse_output.fasta"

        with open(input_file, "w") as f:
            f.write(">ERZ123|viral_id_001\n")
            f.write("ATCGATCGATCGATCG\n")
            f.write(">ERZ456|viral_id_002\n")
            f.write("GCTAGCTAGCTAGCTA\n")

        # Rename with virify mode (handles metaviraverse format)
        from mgnify_pipelines_toolkit.utils.rename_contig.parsers import parse_virify_header

        mapping = rename_fasta(
            str(input_file),
            str(output_file),
            prefix="contig_",
            parser_fn=parse_virify_header,
            formatter_mode="virify",
        )

        # Verify output
        with open(output_file) as f:
            content = f.read()
            assert ">contig_1|viral_id_001" in content
            assert ">contig_2|viral_id_002" in content
            assert "ERZ123" not in content
            assert "ERZ456" not in content

        # Verify mapping
        assert len(mapping) == 2

    def test_viral_id_with_phage_circular(self, temp_dir):
        """Test viral ID with phage-circular metadata."""
        input_file = temp_dir / "viral_circular.fasta"
        output_file = temp_dir / "viral_circular_out.fasta"

        with open(input_file, "w") as f:
            f.write(">seq1|viral_id_003|phage-circular\n")
            f.write("ATCGATCG\n")

        from mgnify_pipelines_toolkit.utils.rename_contig.parsers import parse_virify_header

        rename_fasta(
            str(input_file),
            str(output_file),
            prefix="contig_",
            parser_fn=parse_virify_header,
            formatter_mode="virify",
        )

        with open(output_file) as f:
            content = f.read()
            assert ">contig_1|viral_id_003|phage-circular" in content

    def test_viral_id_with_prophage(self, temp_dir):
        """Test viral ID with prophage coordinates."""
        input_file = temp_dir / "viral_prophage.fasta"
        output_file = temp_dir / "viral_prophage_out.fasta"

        with open(input_file, "w") as f:
            f.write(">seq1|viral_id_004|prophage-100:500\n")
            f.write("ATCGATCG\n")

        from mgnify_pipelines_toolkit.utils.rename_contig.parsers import parse_virify_header

        rename_fasta(
            str(input_file),
            str(output_file),
            prefix="contig_",
            parser_fn=parse_virify_header,
            formatter_mode="virify",
        )

        with open(output_file) as f:
            content = f.read()
            assert ">contig_1|viral_id_004|prophage-100:500" in content

    def test_metaviraverse_with_gff(self, temp_dir):
        """Test renaming FASTA and GFF together for metaviraverse."""
        fasta_in = temp_dir / "meta_viral.fasta"
        gff_in = temp_dir / "meta_viral.gff"
        fasta_out = temp_dir / "meta_viral_out.fasta"
        gff_out = temp_dir / "meta_viral_out.gff"

        # Create FASTA
        with open(fasta_in, "w") as f:
            f.write(">ERZ1.1|viral_001\n")
            f.write("ATCGATCG\n")

        # Create GFF
        with open(gff_in, "w") as f:
            f.write("##gff-version 3\n")
            f.write("##sequence-region ERZ1.1 1 8\n")
            f.write("ERZ1.1\t.\tgene\t1\t8\t.\t+\t.\tID=gene1\n")

        # Rename FASTA
        from mgnify_pipelines_toolkit.utils.rename_contig.parsers import parse_virify_header

        mapping = rename_fasta(
            str(fasta_in),
            str(fasta_out),
            prefix="contig_",
            parser_fn=parse_virify_header,
            formatter_mode="virify",
        )

        # Create mapping for GFF
        old_to_new = {old: new for new, (old, short) in mapping.items()}
        # Add short names
        for new, (old, short) in mapping.items():
            if short != old:
                old_to_new[short] = new

        # Rename GFF
        rename_gff(str(gff_in), str(gff_out), old_to_new)

        # Verify FASTA has viral ID
        with open(fasta_out) as f:
            content = f.read()
            assert ">contig_1|viral_001" in content

        # Verify GFF was renamed
        with open(gff_out) as f:
            content = f.read()
            assert "contig_1" in content
            assert "ERZ1.1" not in content


class TestVirifyPipeline:
    """Test cases specific to Virify pipeline requirements."""

    def test_virsorter_metadata_preservation(self, temp_dir):
        """Test VirSorter metadata is preserved correctly."""
        input_file = temp_dir / "virsorter.fasta"
        output_file = temp_dir / "virsorter_out.fasta"

        with open(input_file, "w") as f:
            f.write(">contig1|phage-circular\n")
            f.write("ATCGATCG\n")
            f.write(">contig2|prophage-100:500\n")
            f.write("GCTAGCTA\n")
            f.write(">contig3|prophage-100:500|prophage-600:900\n")
            f.write("TTAATTAA\n")

        from mgnify_pipelines_toolkit.utils.rename_contig.parsers import parse_virify_header

        rename_fasta(
            str(input_file),
            str(output_file),
            prefix="seq",
            parser_fn=parse_virify_header,
            formatter_mode="virify",
        )

        with open(output_file) as f:
            lines = [line.strip() for line in f if line.startswith(">")]
            assert ">seq1|phage-circular" in lines
            assert ">seq2|prophage-100:500" in lines
            assert ">seq3|prophage-100:500|prophage-600:900" in lines

    def test_virify_restore_with_metadata(self, temp_dir):
        """Test restore functionality preserves viral metadata."""
        input_file = temp_dir / "virify_in.fasta"
        renamed_file = temp_dir / "virify_renamed.fasta"
        restored_file = temp_dir / "virify_restored.fasta"

        # Create input with viral metadata
        with open(input_file, "w") as f:
            f.write(">original_seq|viral_id|phage-circular\n")
            f.write("ATCGATCG\n")

        from mgnify_pipelines_toolkit.utils.rename_contig.parsers import parse_virify_header

        # Rename
        mapping = rename_fasta(
            str(input_file),
            str(renamed_file),
            prefix="temp_",
            parser_fn=parse_virify_header,
            formatter_mode="virify",
        )

        # Restore
        reverse_mapping = {new: old for new, (old, _) in mapping.items()}
        restore_fasta(
            str(renamed_file),
            str(restored_file),
            reverse_mapping,
            parser_fn=parse_virify_header,
            formatter_mode="virify",
        )

        # Verify metadata preserved
        with open(restored_file) as f:
            content = f.read()
            assert "viral_id" in content
            assert "phage-circular" in content


class TestMETTPipeline:
    """Test cases specific to METT pipeline requirements."""

    def test_mett_basic_renaming(self, temp_dir):
        """Test basic METT renaming without metadata."""
        input_file = temp_dir / "mett_input.fasta"
        output_file = temp_dir / "mett_output.fasta"

        with open(input_file, "w") as f:
            f.write(">NODE_1_length_1000_cov_10.5\n")
            f.write("ATCGATCGATCGATCG\n")
            f.write(">NODE_2_length_500_cov_5.2\n")
            f.write("GCTAGCTAGCTAGCTA\n")

        # METT uses simple mode (no metadata preservation)
        mapping = rename_fasta(
            str(input_file),
            str(output_file),
            prefix="contig_",
            formatter_mode="simple",
        )

        with open(output_file) as f:
            lines = [line.strip() for line in f if line.startswith(">")]
            assert ">contig_1" in lines
            assert ">contig_2" in lines
            assert len(lines) == 2

    def test_mett_with_genbank(self, temp_dir):
        """Test METT GenBank file renaming."""
        fasta_in = temp_dir / "mett.fasta"
        gbk_in = temp_dir / "mett.gbk"
        fasta_out = temp_dir / "mett_out.fasta"
        gbk_out = temp_dir / "mett_out.gbk"

        # Create FASTA
        with open(fasta_in, "w") as f:
            f.write(">contig1\n")
            f.write("ATCGATCG\n")

        # Create GenBank
        with open(gbk_in, "w") as f:
            f.write("LOCUS       contig1                    8 bp    DNA     linear   01-JAN-1980\n")
            f.write("DEFINITION  Test sequence.\n")
            f.write("//\n")

        # Rename FASTA
        mapping = rename_fasta(
            str(fasta_in),
            str(fasta_out),
            prefix="renamed_",
            formatter_mode="simple",
        )

        # Rename GenBank
        old_to_new = {old: new for new, (old, _) in mapping.items()}
        rename_genbank(str(gbk_in), str(gbk_out), old_to_new)

        # Verify GenBank renamed
        with open(gbk_out) as f:
            content = f.read()
            assert "renamed_1" in content
            assert "contig1" not in content

    def test_mett_multi_format(self, temp_dir):
        """Test METT with FASTA, GFF, and GenBank together."""
        fasta_in = temp_dir / "mett_multi.fasta"
        gff_in = temp_dir / "mett_multi.gff"
        gbk_in = temp_dir / "mett_multi.gbk"

        fasta_out = temp_dir / "mett_multi_out.fasta"
        gff_out = temp_dir / "mett_multi_out.gff"
        gbk_out = temp_dir / "mett_multi_out.gbk"

        # Create test files
        with open(fasta_in, "w") as f:
            f.write(">seq1\n")
            f.write("ATCGATCG\n")

        with open(gff_in, "w") as f:
            f.write("##gff-version 3\n")
            f.write("seq1\t.\tgene\t1\t8\t.\t+\t.\tID=gene1\n")

        with open(gbk_in, "w") as f:
            f.write("LOCUS       seq1                       8 bp    DNA     linear   01-JAN-1980\n")
            f.write("//\n")

        # Rename all
        mapping = rename_fasta(str(fasta_in), str(fasta_out), prefix="c", formatter_mode="simple")
        old_to_new = {old: new for new, (old, _) in mapping.items()}

        rename_gff(str(gff_in), str(gff_out), old_to_new)
        rename_genbank(str(gbk_in), str(gbk_out), old_to_new)

        # Verify all renamed consistently
        with open(fasta_out) as f:
            assert ">c1" in f.read()
        with open(gff_out) as f:
            assert "c1" in f.read()
        with open(gbk_out) as f:
            assert "c1" in f.read()


class TestASAPipeline:
    """Test cases specific to ASA pipeline requirements."""

    def test_asa_simple_renaming(self, temp_dir):
        """Test simple ASA renaming."""
        input_file = temp_dir / "asa_input.fasta"
        output_file = temp_dir / "asa_output.fasta"

        with open(input_file, "w") as f:
            f.write(">sequence_1\n")
            f.write("ATCG\n")
            f.write(">sequence_2\n")
            f.write("GCTA\n")

        mapping = rename_fasta(
            str(input_file),
            str(output_file),
            prefix="seq",
            formatter_mode="simple",
        )

        with open(output_file) as f:
            content = f.read()
            assert ">seq1" in content
            assert ">seq2" in content

    def test_asa_mapping_file_format(self, temp_dir):
        """Test ASA mapping file has correct format."""
        input_file = temp_dir / "asa.fasta"
        output_file = temp_dir / "asa_out.fasta"
        mapping_file = temp_dir / "asa_mapping.tsv"

        with open(input_file, "w") as f:
            f.write(">original_name_1\n")
            f.write("ATCG\n")

        mapping = rename_fasta(str(input_file), str(output_file), prefix="renamed_")
        write_mapping_file(mapping, str(mapping_file))

        # Verify mapping file has correct columns
        with open(mapping_file) as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            assert header == ["original", "renamed", "short"]

            data = lines[1].strip().split("\t")
            assert data[0] == "original_name_1"  # original
            assert data[1] == "renamed_1"  # renamed
            assert data[2] == "original_name_1"  # short


class TestMAPPipeline:
    """Test cases specific to MAP (Mobilome) pipeline requirements."""

    def test_map_size_based_separation(self, temp_dir):
        """Test MAP pipeline size-based file separation."""
        input_file = temp_dir / "assembly.fasta"
        output_prefix = str(temp_dir / "assembly")

        # Create sequences of different sizes
        with open(input_file, "w") as f:
            # 500 bp - should not appear in any output
            f.write(">seq1\n")
            f.write("A" * 500 + "\n")

            # 1500 bp - should appear in 1kb only
            f.write(">seq2\n")
            f.write("C" * 1500 + "\n")

            # 6000 bp - should appear in 1kb and 5kb
            f.write(">seq3\n")
            f.write("G" * 6000 + "\n")

            # 150000 bp - should appear in all three
            f.write(">seq4\n")
            f.write("T" * 150000 + "\n")

        # Use rename_fasta_by_size function
        mapping = rename_fasta_by_size(
            str(input_file),
            output_prefix,
            prefix="contig_",
        )

        # Check all three output files exist
        file_1kb = temp_dir / "assembly_1kb_contigs.fasta"
        file_5kb = temp_dir / "assembly_5kb_contigs.fasta"
        file_100kb = temp_dir / "assembly_100kb_contigs.fasta"

        assert file_1kb.exists()
        assert file_5kb.exists()
        assert file_100kb.exists()

        # Verify 1kb file (sequences > 1000 bp)
        with open(file_1kb) as f:
            content_1kb = f.read()
            assert ">contig_1" not in content_1kb  # 500 bp
            assert ">contig_2" in content_1kb  # 1500 bp
            assert ">contig_3" in content_1kb  # 6000 bp
            assert ">contig_4" in content_1kb  # 150000 bp

        # Verify 5kb file (sequences > 5000 bp)
        with open(file_5kb) as f:
            content_5kb = f.read()
            assert ">contig_1" not in content_5kb  # 500 bp
            assert ">contig_2" not in content_5kb  # 1500 bp
            assert ">contig_3" in content_5kb  # 6000 bp
            assert ">contig_4" in content_5kb  # 150000 bp

        # Verify 100kb file (sequences >= 100000 bp)
        with open(file_100kb) as f:
            content_100kb = f.read()
            assert ">contig_1" not in content_100kb  # 500 bp
            assert ">contig_2" not in content_100kb  # 1500 bp
            assert ">contig_3" not in content_100kb  # 6000 bp
            assert ">contig_4" in content_100kb  # 150000 bp

    def test_map_uppercase_sequences(self, temp_dir):
        """Test MAP pipeline converts sequences to uppercase."""
        input_file = temp_dir / "mixed_case.fasta"
        output_prefix = str(temp_dir / "upper")

        with open(input_file, "w") as f:
            f.write(">seq1\n")
            f.write("atcgATCG" * 200 + "\n")  # 1600 bp

        rename_fasta_by_size(str(input_file), output_prefix, prefix="contig_")

        # Check 1kb file has uppercase
        file_1kb = temp_dir / "upper_1kb_contigs.fasta"
        with open(file_1kb) as f:
            content = f.read()
            assert "ATCGATCG" in content
            assert "atcg" not in content.split("\n")[1]  # Check sequence line, not header

    def test_map_mapping_file(self, temp_dir):
        """Test MAP pipeline creates correct mapping."""
        input_file = temp_dir / "map_test.fasta"
        output_prefix = str(temp_dir / "map_out")

        with open(input_file, "w") as f:
            f.write(">original_seq_1\n")
            f.write("A" * 2000 + "\n")

        mapping = rename_fasta_by_size(str(input_file), output_prefix, prefix="contig_")

        # Verify mapping structure
        assert "contig_1" in mapping
        original, short = mapping["contig_1"]
        assert original == "original_seq_1"
        assert short == "original_seq_1"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
