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

from __future__ import annotations

import fileinput
import logging
import re
from collections.abc import Iterable
from pathlib import Path

logger = logging.getLogger(__name__)


def iter_gff_rows(path: Path) -> Iterable[tuple[str, list[str]]]:
    """Yield (raw_line_no_newline, cols) for non-comment GFF rows with 9 columns.

    :param path: Path to a GFF file (may be gzip-compressed).
    :type path: Path
    :returns: Iterable of (line, cols) tuples for valid 9-column rows.
    :rtype: Iterable[tuple[str, list[str]]]
    """
    with fileinput.input(files=[str(path)], openhook=fileinput.hook_compressed) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                logger.debug(f"Skipping non-9col line in {path}: {line}")
                continue
            yield line, cols


def parse_attr_str(attr: str) -> dict[str, str]:
    """Parse a GFF3 attribute string into a key-value dictionary.

    :param attr: GFF3 column 9 attribute string (semicolon-separated key=value pairs).
    :type attr: str
    :returns: Parsed attributes as a dictionary.
    :rtype: dict[str, str]
    """
    if attr == "." or attr.strip() == "":
        return {}
    out: dict[str, str] = {}
    for part in attr.split(";"):
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k] = v
    return out


def _sanitize_attr_value(v: str) -> str:
    """Sanitize an attribute value, replacing raw semicolons to avoid field breakage.

    :param v: Raw attribute value string.
    :type v: str
    :returns: Sanitized value with semicolons replaced by commas.
    :rtype: str
    """
    return str(v).replace(";", ",")


def join_unique(values: Iterable[str]) -> str:
    """Join unique, non-empty sanitized values into a comma-separated sorted string.

    :param values: Iterable of string values.
    :type values: Iterable[str]
    :returns: Sorted unique values joined by commas.
    :rtype: str
    """
    uniq = sorted({_sanitize_attr_value(v) for v in values if v is not None and str(v).strip() != ""})
    return ",".join(uniq)


def attrs_to_str_with_id_first(attrs: dict[str, str]) -> str:
    """Render a GFF3 attributes dict as a string with ID first.

    Used only for synthetic BGC region rows.

    :param attrs: Attributes dictionary; ID key will be placed first.
    :type attrs: dict[str, str]
    :returns: Semicolon-separated attribute string, or "." if empty.
    :rtype: str
    """
    if not attrs:
        return "."
    parts: list[str] = []
    if "ID" in attrs:
        parts.append(f"ID={_sanitize_attr_value(attrs['ID'])}")
    for k in sorted(attrs.keys()):
        if k == "ID":
            continue
        parts.append(f"{k}={_sanitize_attr_value(attrs[k])}")
    return ";".join(parts) if parts else "."


def extract_id_from_attr(attr: str) -> str | None:
    """Extract the ID value from a GFF3 attribute string.

    :param attr: GFF3 column 9 attribute string.
    :type attr: str
    :returns: The value of the ID attribute, or None if not present.
    :rtype: str | None
    """
    m = re.search(r"(?:(?<=;)|^)ID=([^;]+)", attr)
    return m.group(1) if m else None


def replace_or_append_attr(attr: str, key: str, value: str) -> str:
    """Replace key=... in a GFF3 attribute string if present, else append ;key=value.

    Does NOT reformat or reorder anything else in the attribute string.

    :param attr: Existing GFF3 attribute string.
    :type attr: str
    :param key: Attribute key to replace or append.
    :type key: str
    :param value: Attribute value to set.
    :type value: str
    :returns: Updated attribute string.
    :rtype: str
    """
    value = _sanitize_attr_value(value)
    if attr == "." or attr.strip() == "":
        return f"{key}={value}"
    if f"{key}=" in attr:
        pattern = rf"(?:(?<=;)|^){re.escape(key)}=[^;]*"
        return re.sub(pattern, f"{key}={value}", attr)
    return attr + f";{key}={value}"


def write_gff(output: Path, lines: list[tuple[str, int, str]]) -> None:
    """Write a GFF3 file from a flat list of (contig, start, line) tuples.

    Rows are sorted by (contig, start) before writing. The caller is
    responsible for building the unified list from all feature types.

    :param output: Output GFF3 file path.
    :type output: Path
    :param lines: (contig, start, line) tuples for all feature rows.
    :type lines: list[tuple[str, int, str]]
    """
    sorted_lines = sorted(lines, key=lambda x: (x[0], x[1]))
    with output.open("w") as out:
        out.write("##gff-version 3\n")
        for _, _, ln in sorted_lines:
            out.write(ln + "\n")
