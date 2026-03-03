#!/usr/bin/env python3
"""
test_io.py

Tests for mgnify_pipelines_toolkit.analysis.shared.gff.io – generic GFF3 primitives.
"""

import pytest

from mgnify_pipelines_toolkit.analysis.shared.gff.io import (
    attrs_to_str_with_id_first,
    extract_id_from_attr,
    iter_gff_rows,
    parse_attr_str,
    replace_or_append_attr,
)

# ──────────────────────────────────────────────────────────────────────────────
# iter_gff_rows
# ──────────────────────────────────────────────────────────────────────────────


def test_iter_gff_rows_basic(tmp_path):
    """iter_gff_rows yields data rows and skips comments/blank lines."""
    gff = tmp_path / "test.gff"
    gff.write_text("##gff-version 3\n\n# a comment\nctg1\tsrc\tCDS\t100\t200\t.\t+\t0\tID=gene1\n")
    rows = list(iter_gff_rows(gff))
    assert len(rows) == 1
    line, cols = rows[0]
    assert cols[0] == "ctg1"
    assert cols[2] == "CDS"
    assert cols[8] == "ID=gene1"


def test_iter_gff_rows_skips_short_rows(tmp_path):
    """iter_gff_rows skips rows with fewer than 9 tab-separated columns."""
    gff = tmp_path / "test.gff"
    gff.write_text("ctg1\tsrc\tCDS\t100\t200\t.\t+\n")  # only 7 columns
    rows = list(iter_gff_rows(gff))
    assert rows == []


def test_iter_gff_rows_preserves_raw_line(tmp_path):
    """iter_gff_rows returns the raw line without trailing newline."""
    raw = "ctg1\tsrc\tCDS\t100\t200\t.\t+\t0\tID=gene1"
    gff = tmp_path / "test.gff"
    gff.write_text(raw + "\n")
    rows = list(iter_gff_rows(gff))
    assert rows[0][0] == raw


# ──────────────────────────────────────────────────────────────────────────────
# parse_attr_str
# ──────────────────────────────────────────────────────────────────────────────


def test_parse_attr_str_basic():
    """parse_attr_str handles a normal key=value string."""
    result = parse_attr_str("ID=gene1;Name=foo;product=bar")
    assert result == {"ID": "gene1", "Name": "foo", "product": "bar"}


def test_parse_attr_str_dot():
    """parse_attr_str returns empty dict for '.' attribute column."""
    assert parse_attr_str(".") == {}


def test_parse_attr_str_empty():
    """parse_attr_str returns empty dict for blank string."""
    assert parse_attr_str("") == {}
    assert parse_attr_str("   ") == {}


def test_parse_attr_str_no_value():
    """parse_attr_str skips parts without '='."""
    result = parse_attr_str("ID=gene1;badpart;Name=foo")
    assert "ID" in result
    assert "Name" in result
    assert "badpart" not in result


def test_parse_attr_str_value_with_equals():
    """parse_attr_str splits only on the first '=' in each part."""
    result = parse_attr_str("ID=a=b")
    assert result["ID"] == "a=b"


# ──────────────────────────────────────────────────────────────────────────────
# extract_id_from_attr
# ──────────────────────────────────────────────────────────────────────────────


def test_extract_id_from_attr_present():
    """extract_id_from_attr returns the ID value."""
    assert extract_id_from_attr("ID=gene1;Name=foo") == "gene1"


def test_extract_id_from_attr_id_not_first():
    """extract_id_from_attr finds ID even when it is not the first attribute."""
    assert extract_id_from_attr("Name=foo;ID=gene1;product=bar") == "gene1"


def test_extract_id_from_attr_missing():
    """extract_id_from_attr returns None when ID attribute is absent."""
    assert extract_id_from_attr("Name=foo;product=bar") is None


def test_extract_id_from_attr_empty():
    """extract_id_from_attr returns None for empty string."""
    assert extract_id_from_attr("") is None


# ──────────────────────────────────────────────────────────────────────────────
# replace_or_append_attr
# ──────────────────────────────────────────────────────────────────────────────


def test_replace_or_append_attr_append():
    """replace_or_append_attr appends a new key=value pair."""
    result = replace_or_append_attr("ID=gene1;Name=foo", "bgc_tools", "gecco")
    assert result == "ID=gene1;Name=foo;bgc_tools=gecco"


def test_replace_or_append_attr_replace():
    """replace_or_append_attr replaces an existing key without moving it."""
    result = replace_or_append_attr("ID=gene1;bgc_tools=old;Name=foo", "bgc_tools", "gecco,sanntis")
    assert "bgc_tools=gecco,sanntis" in result
    assert "bgc_tools=old" not in result


def test_replace_or_append_attr_dot():
    """replace_or_append_attr treats '.' as empty and returns key=value."""
    assert replace_or_append_attr(".", "bgc_support", "1.00") == "bgc_support=1.00"


def test_replace_or_append_attr_empty_string():
    """replace_or_append_attr treats empty string as empty."""
    assert replace_or_append_attr("", "bgc_support", "0.50") == "bgc_support=0.50"


def test_replace_or_append_attr_sanitizes_semicolons():
    """replace_or_append_attr replaces semicolons in value with commas."""
    result = replace_or_append_attr("ID=gene1", "note", "a;b")
    assert "note=a,b" in result


# ──────────────────────────────────────────────────────────────────────────────
# attrs_to_str_with_id_first
# ──────────────────────────────────────────────────────────────────────────────


def test_attrs_to_str_with_id_first_id_first():
    """attrs_to_str_with_id_first places ID before other keys."""
    attrs = {"bgc_tools": "gecco", "ID": "ctg1|bgc:100-200"}
    result = attrs_to_str_with_id_first(attrs)
    assert result.startswith("ID=ctg1|bgc:100-200")
    assert "bgc_tools=gecco" in result


def test_attrs_to_str_with_id_first_no_id():
    """attrs_to_str_with_id_first works when no ID key is present."""
    attrs = {"bgc_tools": "gecco", "note": "test"}
    result = attrs_to_str_with_id_first(attrs)
    assert "bgc_tools=gecco" in result
    assert "note=test" in result
    assert not result.startswith("ID=")


def test_attrs_to_str_with_id_first_empty():
    """attrs_to_str_with_id_first returns '.' for empty dict."""
    assert attrs_to_str_with_id_first({}) == "."


def test_attrs_to_str_with_id_first_sanitizes_semicolons():
    """attrs_to_str_with_id_first replaces semicolons in values."""
    attrs = {"ID": "x;y"}
    result = attrs_to_str_with_id_first(attrs)
    assert "ID=x,y" in result


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
