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
# /// script
# dependencies = [
#   "requests",
#   "datamodel-code-generator",
# ]
# ///

"""Generate Pydantic v2 models from the antiSMASH sideloader JSON schemas.

The schemas use ``file:`` URI references and are split across a subschemas/
subdirectory, so this script downloads everything, flattens the layout, strips
the ``file:`` prefix from ``$ref`` values, and then runs datamodel-codegen.

Usage::

    uvx --with datamodel-code-generator build_models.py [OUTPUT]

``OUTPUT`` defaults to ``models.py``.
"""

import json
import re
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

BASE_URL = "https://raw.githubusercontent.com/antismash/antismash/master/antismash/detection/sideloader/schemas/general"
ANTISMASH_REPO = "https://api.github.com/repos/antismash/antismash"

# Maps local filename -> remote URL (flattened into a single directory).
SCHEMAS: dict[str, str] = {
    "schema.json": f"{BASE_URL}/schema.json",
    "tool.json": f"{BASE_URL}/subschemas/tool.json",
    "record.json": f"{BASE_URL}/subschemas/record.json",
    "details.json": f"{BASE_URL}/subschemas/details.json",
    "protocluster.json": f"{BASE_URL}/subschemas/protocluster.json",
    "subregion.json": f"{BASE_URL}/subschemas/subregion.json",
}


def fix_refs(node: Any) -> Any:
    """Recursively strip the ``file:`` prefix from all ``$ref`` values.

    datamodel-codegen resolves ``$ref`` as relative file paths, so
    ``file:tool.json`` must become ``tool.json``.

    :param node: A parsed JSON value (dict, list, or scalar).
    :return: The same structure with fixed ``$ref`` strings.
    """
    if isinstance(node, dict):
        return {k: (v.replace("file:", "") if k == "$ref" and isinstance(v, str) else fix_refs(v)) for k, v in node.items()}
    if isinstance(node, list):
        return [fix_refs(item) for item in node]
    return node


def get_commit_info() -> dict[str, str]:
    """Fetch the latest commit information from the antiSMASH repository.

    :return: Dictionary with 'sha', 'short_sha', 'date', and 'message' keys.
    :raises requests.RequestException: If the GitHub API request fails.
    """
    print("  Fetching commit information from antiSMASH repository ...")
    response = requests.get(
        f"{ANTISMASH_REPO}/commits?path=antismash/detection/sideloader/schemas/general&per_page=1",
        timeout=30,
    )
    response.raise_for_status()
    commits = response.json()
    if not commits:
        return {
            "sha": "unknown",
            "short_sha": "unknown",
            "date": "unknown",
            "message": "unknown",
        }

    commit = commits[0]
    sha = commit["sha"]
    short_sha = sha[:7]
    date = commit["commit"]["committer"]["date"]
    message = commit["commit"]["message"].split("\n")[0]  # First line only

    return {
        "sha": sha,
        "short_sha": short_sha,
        "date": date,
        "message": message,
    }


def download_schemas(dest: Path) -> None:
    """Download all schema files into *dest*, fixing ``$ref`` paths.

    :param dest: Directory to write the processed schema files into.
    """
    for filename, url in SCHEMAS.items():
        print(f"  Downloading {filename} ...")
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        schema = fix_refs(response.json())
        (dest / filename).write_text(json.dumps(schema, indent=2))


def generate_models(schema_path: Path, output: Path, commit_info: dict[str, str]) -> None:
    """Generate Pydantic v2 models from JSON schema using datamodel-code-generator as a module.

    :param schema_path: Path to the root ``schema.json`` file.
    :param output: Destination ``.py`` file for the generated models.
    :param commit_info: Dictionary containing commit metadata from GitHub.
    :raises Exception: If model generation fails (with full traceback).
    """
    from datamodel_code_generator import DataModelType, InputFileType, generate

    print(f"  Generating models from {schema_path} -> {output}")
    
    # Use datamodel-code-generator as a module instead of CLI
    # This provides better error handling and integrated tracebacks
    generate(
        input_=schema_path,
        input_file_type=InputFileType.JsonSchema,
        output=output,
        output_model_type=DataModelType.PydanticV2BaseModel,
        use_annotated=True,
        field_constraints=True,
        use_double_quotes=True,
    )

    # Add commit information to the generated file header
    add_commit_metadata(output, commit_info)


def add_commit_metadata(output: Path, commit_info: dict[str, str]) -> None:
    """Append commit metadata to the generated models file header.

    :param output: Path to the generated models.py file.
    :param commit_info: Dictionary containing commit metadata.
    """
    content = output.read_text()

    # Find the end of the auto-generated header (after timestamp line)
    # and insert our metadata before the first blank line
    lines = content.split("\n")
    insert_idx = 0

    for i, line in enumerate(lines):
        if line.startswith("#   timestamp:"):
            insert_idx = i + 1
            break

    if insert_idx > 0:
        metadata_lines = [
            f"#   source: antiSMASH repository (master branch)",
            f"#   commit: {commit_info['sha']}",
            f"#   commit_date: {commit_info['date']}",
            f"#   commit_msg: {commit_info['message']}",
        ]
        lines[insert_idx:insert_idx] = metadata_lines
        output.write_text("\n".join(lines))
        print(f"  Added commit metadata (commit {commit_info['short_sha']})")


def main() -> None:
    """Entry point: download schemas, fix refs, and generate models.

    :raises SystemExit: On download or code-generation failure.
    """
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("models.py")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        print("Downloading schemas ...")
        commit_info = get_commit_info()
        download_schemas(tmppath)

        print("Generating models ...")
        generate_models(tmppath / "schema.json", output.resolve(), commit_info)

    print(f"Done — models written to {output}")


if __name__ == "__main__":
    main()
