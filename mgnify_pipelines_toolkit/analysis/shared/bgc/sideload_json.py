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

import json
import logging
import re
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path

from mgnify_pipelines_toolkit.analysis.shared.bgc.models import MergedRegion

logger = logging.getLogger(__name__)


def _derive_output_json_path(output_gff: Path) -> Path:
    """Derive a JSON output path with the same basename as output_gff.

    Examples:
      - out.gff    -> out.json
      - out.gff3   -> out.json
      - out.gff.gz -> out.json

    :param output_gff: Path to the output GFF file.
    :type output_gff: Path
    :returns: Derived JSON output path.
    :rtype: Path
    """
    suffixes = output_gff.suffixes

    # Handle .gff.gz, .gff3.gz cases (two suffixes)
    if len(suffixes) >= 2 and suffixes[-1] == ".gz" and suffixes[-2] in (".gff", ".gff3"):
        return output_gff.with_suffix("").with_suffix(".json")

    # Handle .gff, .gff3 cases (single suffix)
    if suffixes and suffixes[-1] in (".gff", ".gff3"):
        return output_gff.with_suffix(".json")

    # Fallback for other extensions
    return output_gff.with_suffix(".json")


def validate_sideload_json(payload: dict, schema_dir: Path | None = None) -> None:
    """Validate payload using Pydantic models instead of JSON schemas.

    :param payload: JSON payload to validate.
    :type payload: dict
    :param schema_dir: Kept for backward compatibility but not used (Pydantic doesn't need it).
    :type schema_dir: Path | None
    :raises ValueError: If the payload fails Pydantic validation.
    """
    print("Validating json format")

    from mgnify_pipelines_toolkit.analysis.shared.bgc.antismash_sideload_schemas.antismash_schema import (
        Model,
    )

    try:
        Model.model_validate(payload)
        logger.info("✓ Validation successful")
    except Exception as e:
        raise ValueError(f"Schema validation failed:\n{e}") from e


def _sanitize_details_value(value: str) -> str:
    """Make sure details values are schema-friendly (no leading whitespace or special chars).

    :param value: Raw value string.
    :type value: str
    :returns: Sanitized value, falling back to "NA" if empty after stripping.
    :rtype: str
    """
    v = str(value).strip()
    while v and v[0] in ("_", "=", ","):
        v = v[1:]
    return v if v else "NA"


def _choose_subregion_label(mr: MergedRegion) -> str:
    """Choose a conservative <=20 char label for a subregion.

    Uses the BGCRegion.label field set by the parser, falling back to the
    tool name if no label was provided.

    :param mr: MergedRegion to label.
    :type mr: MergedRegion
    :returns: Short alphanumeric label string (max 20 chars).
    :rtype: str
    """
    if len(mr.members) > 1:
        return "MergedBGC"
    m = mr.members[0]
    lbl = m.label or m.tool
    lbl = re.sub(r"\s+", "_", lbl)
    lbl = re.sub(r"[^A-Za-z0-9_\-\.]", "_", lbl)
    return (lbl[:20] if len(lbl) > 20 else lbl) or "BGC"


def build_sideload_json_payload(
    merged_regions: Sequence[MergedRegion],
    *,
    tool_name: str,
    tool_version: str,
    tool_description: str,
) -> dict:
    """Build an antiSMASH sideloader JSON payload using only records/subregions/details.

    :param merged_regions: Merged BGC regions to include in the payload.
    :type merged_regions: Sequence[MergedRegion]
    :param tool_name: Name of the tool producing the sideloader output.
    :type tool_name: str
    :param tool_version: Version string of the tool.
    :type tool_version: str
    :param tool_description: Description of the tool.
    :type tool_description: str
    :returns: antiSMASH sideloader JSON payload as a dictionary.
    :rtype: dict
    """
    records_by_name: dict[str, list[dict]] = {}

    for mr in merged_regions:
        if not mr.members:
            continue

        # GFF is 1-based inclusive; sideloader JSON is 0-based start, end-exclusive
        start0 = max(0, int(mr.start) - 1)
        end_excl = int(mr.end)

        tools = sorted({m.tool for m in mr.members})
        sources = sorted({m.source for m in mr.members if m.source})

        union: dict[str, set[str]] = {}
        for m in mr.members:
            for k, v in m.attrs.items():
                if v is None or str(v).strip() == "":
                    continue
                union.setdefault(k, set()).add(_sanitize_details_value(v))

        details: dict[str, object] = {
            "ID": _sanitize_details_value(f"{mr.contig}|bgc:{mr.start}-{mr.end}"),
            "bgc_tools": [_sanitize_details_value(t) for t in tools],
        }
        if len(mr.members) > 1:
            details["member_bgcs"] = _sanitize_details_value(str(len(mr.members)))
        if sources:
            details["sources"] = [_sanitize_details_value(s) for s in sources]

        for k, vals in union.items():
            if k in details:
                k = f"member_{k}"
            details[k] = sorted(vals)

        records_by_name.setdefault(mr.contig, []).append(
            {
                "start": start0,
                "end": end_excl,
                "label": _choose_subregion_label(mr),
                "details": details,
            }
        )

    records: list[dict] = []
    for name in sorted(records_by_name.keys()):
        subregions = sorted(records_by_name[name], key=lambda d: int(d["start"]))
        records.append({"name": name, "subregions": subregions})

    return {
        "tool": {
            "name": tool_name,
            "version": tool_version if str(tool_version).strip() else "unknown",
            "description": tool_description,
            "configuration": {},
        },
        "records": records,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


def write_sideload_json(
    out_json: Path,
    merged_regions: Sequence[MergedRegion],
    *,
    schema_dir: Path | None = None,
    validate: bool = False,
    tool_name: str = "MGnify BGC mapper",
    tool_version: str = "mgnify_pipelines_toolkit_v1.4.18",
    tool_description: str = "BGC subregions integrated from GECCO/antiSMASH/SanntiS into a base GFF.",
) -> None:
    """Build and write the antiSMASH sideloader JSON file. Optionally validate against the official antiSMASH schema files.

    :param out_json: Path where the JSON output should be written.
    :type out_json: Path
    :param merged_regions: Merged BGC regions to serialise.
    :type merged_regions: Sequence[MergedRegion]
    :param schema_dir: Optional schema directory override.
    :type schema_dir: Path | None
    :param validate: Whether to validate the JSON against the schema (default: False).
    :type validate: bool
    :param tool_name: Name of the tool for the JSON payload.
    :type tool_name: str
    :param tool_version: Version string for the JSON payload.
    :type tool_version: str
    :param tool_description: Description for the JSON payload.
    :type tool_description: str
    """
    print("writing json output")

    payload = build_sideload_json_payload(
        merged_regions,
        tool_name=tool_name,
        tool_version=tool_version,
        tool_description=tool_description,
    )

    # Write first (so you still get an artifact even if validation is slow/fails when enabled)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    logger.info("Wrote sideloader JSON: %s", out_json)

    # Validate only if requested
    if validate:
        validate_sideload_json(payload, schema_dir=schema_dir)
