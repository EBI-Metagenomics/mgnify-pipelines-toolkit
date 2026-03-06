#!/usr/bin/env python3
# Copyright 2024-2026 EMBL - European Bioinformatics Institute
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

"""Tests for the sideload_json module."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pytest

from mgnify_pipelines_toolkit.analysis.shared.bgc.models import BGCRegion, MergedRegion
from mgnify_pipelines_toolkit.analysis.shared.bgc.sideload_json import (
    _choose_subregion_label,
    _derive_output_json_path,
    _sanitize_details_value,
    build_sideload_json_payload,
    validate_sideload_json,
    write_sideload_json,
)


class TestDeriveOutputJsonPath:
    """Tests for _derive_output_json_path() function."""

    def test_simple_gff_extension(self):
        """Test simple .gff extension replacement."""
        result = _derive_output_json_path(Path("output.gff"))
        assert result == Path("output.json")

    def test_gff3_extension(self):
        """Test .gff3 extension replacement."""
        result = _derive_output_json_path(Path("output.gff3"))
        assert result == Path("output.json")

    def test_gff_gz_extension(self):
        """Test .gff.gz extension (two suffixes)."""
        result = _derive_output_json_path(Path("output.gff.gz"))
        assert result == Path("output.json")

    def test_gff3_gz_extension(self):
        """Test .gff3.gz extension (two suffixes)."""
        result = _derive_output_json_path(Path("output.gff3.gz"))
        assert result == Path("output.json")

    def test_path_with_directory(self):
        """Test paths with directories."""
        result = _derive_output_json_path(Path("/tmp/data/output.gff"))
        assert result == Path("/tmp/data/output.json")

    def test_fallback_other_extension(self):
        """Test fallback for non-gff extensions."""
        result = _derive_output_json_path(Path("output.txt"))
        assert result == Path("output.json")


class TestSanitizeDetailsValue:
    """Tests for _sanitize_details_value() function."""

    def test_normal_string(self):
        """Test normal string passes through."""
        result = _sanitize_details_value("normal_value")
        assert result == "normal_value"

    def test_leading_whitespace(self):
        """Test leading whitespace is stripped."""
        result = _sanitize_details_value("  value")
        assert result == "value"

    def test_trailing_whitespace(self):
        """Test trailing whitespace is stripped."""
        result = _sanitize_details_value("value  ")
        assert result == "value"

    def test_leading_underscore(self):
        """Test leading underscore is removed."""
        result = _sanitize_details_value("_value")
        assert result == "value"

    def test_leading_equals(self):
        """Test leading equals is removed."""
        result = _sanitize_details_value("=value")
        assert result == "value"

    def test_leading_comma(self):
        """Test leading comma is removed."""
        result = _sanitize_details_value(",value")
        assert result == "value"

    def test_multiple_leading_special_chars(self):
        """Test multiple leading special characters are removed."""
        result = _sanitize_details_value("_=,value")
        assert result == "value"

    def test_empty_string_returns_na(self):
        """Test empty string returns 'NA'."""
        result = _sanitize_details_value("")
        assert result == "NA"

    def test_whitespace_only_returns_na(self):
        """Test whitespace-only string returns 'NA'."""
        result = _sanitize_details_value("   ")
        assert result == "NA"

    def test_special_chars_only_returns_na(self):
        """Test string with only special chars returns 'NA'."""
        result = _sanitize_details_value("_=,")
        assert result == "NA"

    def test_numeric_string(self):
        """Test numeric string is preserved."""
        result = _sanitize_details_value("12345")
        assert result == "12345"


class TestChooseSubregionLabel:
    """Tests for _choose_subregion_label() function."""

    def test_merged_region_multiple_members(self):
        """Test merged region with multiple members returns 'MergedBGC'."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[
                BGCRegion(contig="contig1", start=100, end=300, tool="gecco", source="GECCO", label="NRP"),
                BGCRegion(contig="contig1", start=250, end=500, tool="antismash", source="antiSMASH", label="PKS"),
            ],
        )
        result = _choose_subregion_label(mr)
        assert result == "MergedBGC"

    def test_single_member_with_label(self):
        """Test single member region uses its label."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=300,
            members=[BGCRegion(contig="contig1", start=100, end=300, tool="gecco", source="GECCO", label="NRP,Polyketide")],
        )
        result = _choose_subregion_label(mr)
        assert result == "NRP_Polyketide"

    def test_single_member_without_label_uses_tool(self):
        """Test single member without label uses tool name."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=300,
            members=[BGCRegion(contig="contig1", start=100, end=300, tool="gecco", source="GECCO", label="")],
        )
        result = _choose_subregion_label(mr)
        assert result == "gecco"

    def test_label_whitespace_replaced(self):
        """Test whitespace in label is replaced with underscores."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=300,
            members=[BGCRegion(contig="contig1", start=100, end=300, tool="gecco", source="GECCO", label="NRP Polyketide")],
        )
        result = _choose_subregion_label(mr)
        assert result == "NRP_Polyketide"

    def test_label_special_chars_replaced(self):
        """Test special characters in label are replaced with underscores."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=300,
            members=[BGCRegion(contig="contig1", start=100, end=300, tool="gecco", source="GECCO", label="NRP@Polyketide!")],
        )
        result = _choose_subregion_label(mr)
        assert result == "NRP_Polyketide_"

    def test_label_truncation(self):
        """Test label is truncated to 20 characters."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=300,
            members=[
                BGCRegion(
                    contig="contig1",
                    start=100,
                    end=300,
                    tool="gecco",
                    source="GECCO",
                    label="VeryLongLabelThatExceedsTwentyCharacters",
                )
            ],
        )
        result = _choose_subregion_label(mr)
        assert len(result) == 20
        assert result == "VeryLongLabelThatExc"

    def test_empty_label_falls_back_to_bgc(self):
        """Test empty label falls back to 'BGC'."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=300,
            members=[BGCRegion(contig="contig1", start=100, end=300, tool="", source="", label="")],
        )
        result = _choose_subregion_label(mr)
        assert result == "BGC"


class TestBuildSideloadJsonPayload:
    """Tests for build_sideload_json_payload() function."""

    def test_basic_payload_structure(self):
        """Test basic payload structure with one region."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[
                BGCRegion(
                    contig="contig1",
                    start=100,
                    end=500,
                    tool="gecco",
                    source="GECCO v0.9.8",
                    label="NRP",
                    attrs={"Type": "NRP"},
                )
            ],
        )

        payload = build_sideload_json_payload(
            [mr],
            tool_name="Test Tool",
            tool_version="1.0.0",
            tool_description="Test Description",
        )

        # Check top-level structure
        assert "tool" in payload
        assert "records" in payload
        assert "timestamp" in payload

        # Check tool metadata
        assert payload["tool"]["name"] == "Test Tool"
        assert payload["tool"]["version"] == "1.0.0"
        assert payload["tool"]["description"] == "Test Description"
        assert payload["tool"]["configuration"] == {}

        # Check records
        assert len(payload["records"]) == 1
        record = payload["records"][0]
        assert record["name"] == "contig1"
        assert "subregions" in record

        # Check subregions
        assert len(record["subregions"]) == 1
        subregion = record["subregions"][0]
        assert subregion["start"] == 99  # GFF 1-based to 0-based
        assert subregion["end"] == 500  # end-exclusive
        assert subregion["label"] == "NRP"
        assert "details" in subregion

    def test_coordinate_conversion(self):
        """Test GFF 1-based inclusive to JSON 0-based exclusive conversion."""
        mr = MergedRegion(
            contig="contig1",
            start=1,
            end=1000,
            members=[BGCRegion(contig="contig1", start=1, end=1000, tool="gecco", source="GECCO", label="NRP")],
        )

        payload = build_sideload_json_payload(
            [mr],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        subregion = payload["records"][0]["subregions"][0]
        assert subregion["start"] == 0  # 1 - 1 = 0
        assert subregion["end"] == 1000  # stays the same (end-exclusive)

    def test_multiple_contigs(self):
        """Test payload with multiple contigs."""
        regions = [
            MergedRegion(
                contig="contig1",
                start=100,
                end=500,
                members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
            ),
            MergedRegion(
                contig="contig2",
                start=200,
                end=600,
                members=[BGCRegion(contig="contig2", start=200, end=600, tool="antismash", source="antiSMASH", label="PKS")],
            ),
        ]

        payload = build_sideload_json_payload(
            regions,
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        assert len(payload["records"]) == 2
        contig_names = [r["name"] for r in payload["records"]]
        assert "contig1" in contig_names
        assert "contig2" in contig_names

    def test_multiple_subregions_same_contig(self):
        """Test multiple subregions on the same contig."""
        regions = [
            MergedRegion(
                contig="contig1",
                start=100,
                end=500,
                members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
            ),
            MergedRegion(
                contig="contig1",
                start=1000,
                end=1500,
                members=[BGCRegion(contig="contig1", start=1000, end=1500, tool="antismash", source="antiSMASH", label="PKS")],
            ),
        ]

        payload = build_sideload_json_payload(
            regions,
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        assert len(payload["records"]) == 1
        assert len(payload["records"][0]["subregions"]) == 2

    def test_subregions_sorted_by_start(self):
        """Test subregions are sorted by start coordinate."""
        regions = [
            MergedRegion(
                contig="contig1",
                start=1000,
                end=1500,
                members=[BGCRegion(contig="contig1", start=1000, end=1500, tool="gecco", source="GECCO", label="PKS")],
            ),
            MergedRegion(
                contig="contig1",
                start=100,
                end=500,
                members=[BGCRegion(contig="contig1", start=100, end=500, tool="antismash", source="antiSMASH", label="NRP")],
            ),
        ]

        payload = build_sideload_json_payload(
            regions,
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        subregions = payload["records"][0]["subregions"]
        assert subregions[0]["start"] == 99  # First region starts at 100 (1-based)
        assert subregions[1]["start"] == 999  # Second region starts at 1000 (1-based)

    def test_merged_region_details(self):
        """Test details for merged region with multiple members."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[
                BGCRegion(
                    contig="contig1",
                    start=100,
                    end=300,
                    tool="gecco",
                    source="GECCO v0.9.8",
                    label="NRP",
                    attrs={"Type": "NRP"},
                ),
                BGCRegion(
                    contig="contig1",
                    start=250,
                    end=500,
                    tool="antismash",
                    source="antiSMASH 8.0",
                    label="PKS",
                    attrs={"product": "T1PKS"},
                ),
            ],
        )

        payload = build_sideload_json_payload(
            [mr],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        details = payload["records"][0]["subregions"][0]["details"]

        # Check basic details
        assert details["ID"] == "contig1|bgc:100-500"
        assert "gecco" in details["bgc_tools"]
        assert "antismash" in details["bgc_tools"]
        assert details["member_bgcs"] == "2"

        # Check sources
        assert "sources" in details
        assert "GECCO v0.9.8" in details["sources"]
        assert "antiSMASH 8.0" in details["sources"]

        # Check attributes
        assert "Type" in details
        assert "NRP" in details["Type"]
        assert "product" in details
        assert "T1PKS" in details["product"]

    def test_attribute_union(self):
        """Test that attributes from multiple members are merged."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[
                BGCRegion(
                    contig="contig1",
                    start=100,
                    end=300,
                    tool="gecco",
                    source="GECCO",
                    label="NRP",
                    attrs={"Type": "NRP", "confidence": "high"},
                ),
                BGCRegion(
                    contig="contig1",
                    start=250,
                    end=500,
                    tool="antismash",
                    source="antiSMASH",
                    label="PKS",
                    attrs={"Type": "PKS", "product": "T1PKS"},
                ),
            ],
        )

        payload = build_sideload_json_payload(
            [mr],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        details = payload["records"][0]["subregions"][0]["details"]

        # Check that both Type values are present
        assert "Type" in details
        assert isinstance(details["Type"], list)
        assert "NRP" in details["Type"]
        assert "PKS" in details["Type"]

        # Check other attributes
        assert "confidence" in details
        assert "high" in details["confidence"]
        assert "product" in details
        assert "T1PKS" in details["product"]

    def test_empty_merged_regions_list(self):
        """Test empty merged regions list produces valid payload."""
        payload = build_sideload_json_payload(
            [],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        assert payload["tool"]["name"] == "Test"
        assert payload["records"] == []

    def test_merged_region_with_no_members(self):
        """Test merged region with empty members list is skipped."""
        mr = MergedRegion(contig="contig1", start=100, end=500, members=[])

        payload = build_sideload_json_payload(
            [mr],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        assert payload["records"] == []

    def test_timestamp_is_valid_iso_format(self):
        """Test that timestamp is in valid ISO format."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        payload = build_sideload_json_payload(
            [mr],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        # Should be able to parse the timestamp
        timestamp = payload["timestamp"]
        datetime.fromisoformat(timestamp.replace("Z", "+00:00"))

    def test_empty_tool_version_uses_unknown(self):
        """Test empty tool version defaults to 'unknown'."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        payload = build_sideload_json_payload(
            [mr],
            tool_name="Test",
            tool_version="",
            tool_description="Test",
        )

        assert payload["tool"]["version"] == "unknown"

    def test_whitespace_tool_version_uses_unknown(self):
        """Test whitespace-only tool version defaults to 'unknown'."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        payload = build_sideload_json_payload(
            [mr],
            tool_name="Test",
            tool_version="   ",
            tool_description="Test",
        )

        assert payload["tool"]["version"] == "unknown"


class TestValidateSideloadJson:
    """Tests for validate_sideload_json() function."""

    def test_valid_payload_passes(self):
        """Test that a valid payload passes validation."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[
                BGCRegion(
                    contig="contig1",
                    start=100,
                    end=500,
                    tool="gecco",
                    source="GECCO",
                    label="NRP",
                    attrs={"Type": "NRP"},
                )
            ],
        )

        payload = build_sideload_json_payload(
            [mr],
            tool_name="Test Tool",
            tool_version="1.0.0",
            tool_description="Test Description",
        )

        # Should not raise an exception
        validate_sideload_json(payload)

    def test_invalid_payload_raises_valueerror(self):
        """Test that an invalid payload raises ValueError."""
        invalid_payload = {
            "tool": {
                # Missing required 'name' field
                "version": "1.0.0",
                "description": "Test",
                "configuration": {},
            },
            "records": [],
            "timestamp": datetime.now().isoformat(),
        }

        with pytest.raises(ValueError, match="Schema validation failed"):
            validate_sideload_json(invalid_payload)

    def test_schema_dir_parameter_accepted_but_ignored(self):
        """Test that schema_dir parameter is accepted for backward compatibility but not used."""
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        payload = build_sideload_json_payload(
            [mr],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        # Should work with schema_dir=None (default)
        validate_sideload_json(payload, schema_dir=None)

        # Should also work with a Path value (ignored)
        validate_sideload_json(payload, schema_dir=Path("/fake/path"))

    def test_validation_with_complex_payload(self):
        """Test validation with a complex payload containing multiple regions."""
        regions = [
            MergedRegion(
                contig="contig1",
                start=100,
                end=500,
                members=[
                    BGCRegion(contig="contig1", start=100, end=300, tool="gecco", source="GECCO", label="NRP"),
                    BGCRegion(contig="contig1", start=250, end=500, tool="antismash", source="antiSMASH", label="PKS"),
                ],
            ),
            MergedRegion(
                contig="contig2",
                start=200,
                end=600,
                members=[
                    BGCRegion(
                        contig="contig2",
                        start=200,
                        end=600,
                        tool="sanntis",
                        source="SanntiS",
                        label="Terpene",
                        attrs={"nearest_MiBIG": "BGC0001234"},
                    )
                ],
            ),
        ]

        payload = build_sideload_json_payload(
            regions,
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        # Should not raise an exception
        validate_sideload_json(payload)


class TestWriteSideloadJson:
    """Tests for write_sideload_json() function."""

    def test_file_is_written(self, tmp_path):
        """Test that JSON file is written to disk."""
        out_json = tmp_path / "output.json"
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        write_sideload_json(
            out_json,
            [mr],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        assert out_json.exists()

    def test_file_content_is_valid_json(self, tmp_path):
        """Test that written file contains valid JSON."""
        out_json = tmp_path / "output.json"
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        write_sideload_json(
            out_json,
            [mr],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        # Should be able to parse the JSON
        content = json.loads(out_json.read_text())
        assert "tool" in content
        assert "records" in content
        assert "timestamp" in content

    def test_validation_enabled(self, tmp_path):
        """Test that validation runs when enabled."""
        out_json = tmp_path / "output.json"
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        # Should not raise an exception for valid data
        write_sideload_json(
            out_json,
            [mr],
            validate=True,
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        assert out_json.exists()

    def test_validation_disabled_by_default(self, tmp_path):
        """Test that validation is disabled by default."""
        out_json = tmp_path / "output.json"
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        # Should not run validation
        write_sideload_json(
            out_json,
            [mr],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        assert out_json.exists()

    def test_parent_directory_created(self, tmp_path):
        """Test that parent directories are created if they don't exist."""
        out_json = tmp_path / "subdir" / "output.json"

        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        write_sideload_json(
            out_json,
            [mr],
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        assert out_json.exists()
        assert out_json.parent.exists()

    def test_custom_tool_metadata(self, tmp_path):
        """Test custom tool metadata is written."""
        out_json = tmp_path / "output.json"
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        write_sideload_json(
            out_json,
            [mr],
            tool_name="Custom Tool",
            tool_version="2.5.1",
            tool_description="Custom Description",
        )

        content = json.loads(out_json.read_text())
        assert content["tool"]["name"] == "Custom Tool"
        assert content["tool"]["version"] == "2.5.1"
        assert content["tool"]["description"] == "Custom Description"

    def test_default_tool_metadata(self, tmp_path):
        """Test default tool metadata is used when not specified."""
        out_json = tmp_path / "output.json"
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        write_sideload_json(out_json, [mr])

        content = json.loads(out_json.read_text())
        assert content["tool"]["name"] == "MGnify BGC mapper"
        assert "mgnify_pipelines_toolkit" in content["tool"]["version"]
        assert "BGC subregions integrated" in content["tool"]["description"]

    def test_schema_dir_parameter_backward_compatibility(self, tmp_path):
        """Test that schema_dir parameter is accepted for backward compatibility."""
        out_json = tmp_path / "output.json"
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        # Should accept schema_dir parameter without error
        write_sideload_json(
            out_json,
            [mr],
            schema_dir=Path("/fake/path"),
            validate=False,
            tool_name="Test",
            tool_version="1.0",
            tool_description="Test",
        )

        assert out_json.exists()

    def test_file_written_even_if_validation_fails(self, tmp_path, monkeypatch):
        """Test that file is written even if validation fails."""
        out_json = tmp_path / "output.json"
        mr = MergedRegion(
            contig="contig1",
            start=100,
            end=500,
            members=[BGCRegion(contig="contig1", start=100, end=500, tool="gecco", source="GECCO", label="NRP")],
        )

        # Mock validate_sideload_json to always raise an error
        def mock_validate(*args, **kwargs):
            raise ValueError("Validation failed")

        import mgnify_pipelines_toolkit.analysis.shared.bgc.sideload_json as module

        monkeypatch.setattr(module, "validate_sideload_json", mock_validate)

        # Should raise error but file should still be written
        with pytest.raises(ValueError, match="Validation failed"):
            write_sideload_json(
                out_json,
                [mr],
                validate=True,
                tool_name="Test",
                tool_version="1.0",
                tool_description="Test",
            )

        # File should exist even though validation failed
        assert out_json.exists()
