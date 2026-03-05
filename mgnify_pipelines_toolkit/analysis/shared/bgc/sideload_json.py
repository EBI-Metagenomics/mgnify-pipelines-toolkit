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


def _default_antismash_schema_dir() -> Path:
    """Expected location for vendored antiSMASH sideloader schemas.

    :returns: Path to the antismash_sideload_schemas/general directory.
    :rtype: Path
    """
    # __file__ is analysis/shared/bgc/sideload_json.py → parent.parent = analysis/shared/
    return (Path(__file__).resolve().parent.parent / "antismash_sideload_schemas" / "general").resolve()


def _fix_refs(node):
    """Recursively strip 'file:' from $ref values in a JSON schema node.

    The vendored antiSMASH schemas use 'file:subschemas/foo.json' as $ref
    values. jsonschema's RefResolver looks up refs in the store by URI/path
    key; the store is keyed by plain paths (e.g. 'subschemas/foo.json').
    Stripping 'file:' makes the $ref values match the store keys.

    :param node: JSON-like node (dict, list, or scalar).
    :returns: Node with file: prefixes stripped from $ref values.
    """
    if isinstance(node, dict):
        out = {}
        for k, v in node.items():
            if k == "$ref" and isinstance(v, str):
                out[k] = v.replace("file:", "")
            else:
                out[k] = _fix_refs(v)
        return out
    if isinstance(node, list):
        return [_fix_refs(x) for x in node]
    return node


def _load_schema_files(schema_dir: Path) -> tuple[dict, dict[str, dict]]:
    """Load root schema.json and subschemas into a resolver store.

    jsonschema's RefResolver resolves $ref values by looking them up in a
    store dict. The antiSMASH schemas can appear under several different
    reference formats depending on how the resolver constructs the URI, so
    each document is registered under multiple keys to cover all cases.

    :param schema_dir: Directory containing schema.json and subschemas/.
    :type schema_dir: Path
    :returns: Tuple of (root_schema, store) where store maps URI/path keys to schema dicts.
    :rtype: tuple[dict, dict[str, dict]]
    :raises FileNotFoundError: If schema.json or subschemas/ directory is missing.
    """
    schema_path = schema_dir / "schema.json"
    subs_dir = schema_dir / "subschemas"
    if not schema_path.exists():
        raise FileNotFoundError(f"schema.json not found at: {schema_path}")
    if not subs_dir.exists():
        raise FileNotFoundError(f"subschemas/ dir not found at: {subs_dir}")

    def _read(p: Path) -> dict:
        return _fix_refs(json.loads(p.read_text(encoding="utf-8")))

    root = _read(schema_path)
    store: dict[str, dict] = {}

    def _register_keys_for_doc(p: Path, doc: dict) -> None:
        store[p.resolve().as_uri()] = doc  # absolute file:// URI
        store[p.name] = doc  # bare filename (e.g. "schema.json")
        store[str(p)] = doc  # absolute POSIX path string
        if p.parent.name == "subschemas":
            store[f"subschemas/{p.name}"] = doc  # relative path used in some $ref values
            # Also register under a URI as if the subschema were at the schema root,
            # because the resolver may construct a URI relative to schema.json.
            fake_root_path = (schema_dir / p.name).resolve()
            store[fake_root_path.as_uri()] = doc
            store[str(fake_root_path)] = doc

    _register_keys_for_doc(schema_path, root)
    for p in sorted(subs_dir.glob("*.json")):
        _register_keys_for_doc(p, _read(p))

    return root, store


def validate_sideload_json(payload: dict, schema_dir: Path | None = None) -> None:
    """Validate payload against the official antiSMASH sideloader schema (local copies).

    Uses jsonschema Draft7Validator with a referencing.Registry to resolve $ref
    from in-memory schemas.

    :param payload: JSON payload to validate.
    :type payload: dict
    :param schema_dir: Optional path to schema directory; defaults to vendored schemas.
    :type schema_dir: Path | None
    :raises RuntimeError: If jsonschema is not available.
    :raises jsonschema.ValidationError: If the payload fails schema validation.
    """
    print("Validating json format")

    try:
        from jsonschema import Draft7Validator  # type: ignore
        from referencing import Registry, Resource  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "Validation requires 'jsonschema' (and its 'referencing' dependency). "
            "Install it (e.g. pip install jsonschema) or disable validation."
        ) from e

    schema_dir = _default_antismash_schema_dir() if schema_dir is None else schema_dir
    root, store = _load_schema_files(schema_dir)

    # Build a registry that can resolve $ref values seen in antiSMASH schemas.
    # We register *every* key from store as a Resource.
    registry = Registry()
    for key, doc in store.items():
        registry = registry.with_resource(key, Resource.from_contents(doc))

    # Validate (raises jsonschema.ValidationError on failure)
    validator = Draft7Validator(schema=root, registry=registry)
    validator.validate(payload)


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
