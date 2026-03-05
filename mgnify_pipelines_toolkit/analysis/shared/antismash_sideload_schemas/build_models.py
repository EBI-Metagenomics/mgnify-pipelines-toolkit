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
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import requests

BASE_URL = "https://raw.githubusercontent.com/antismash/antismash/master/antismash/detection/sideloader/schemas/general"

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


def generate_models(schema_path: Path, output: Path) -> None:
    """Run datamodel-codegen to produce Pydantic v2 models.

    :param schema_path: Path to the root ``schema.json`` file.
    :param output: Destination ``.py`` file for the generated models.
    :raises SystemExit: If datamodel-codegen returns a non-zero exit code.
    """
    cmd = [
        "datamodel-codegen",
        "--input",
        str(schema_path),
        "--input-file-type",
        "jsonschema",
        "--output",
        str(output),
        "--output-model-type",
        "pydantic_v2.BaseModel",
        "--use-annotated",
        "--field-constraints",
        "--use-double-quotes",
    ]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)


def main() -> None:
    """Entry point: download schemas, fix refs, and generate models.

    :raises SystemExit: On download or code-generation failure.
    """
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("models.py")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        print("Downloading schemas ...")
        download_schemas(tmppath)

        print("Generating models ...")
        generate_models(tmppath / "schema.json", output.resolve())

    print(f"Done — models written to {output}")


if __name__ == "__main__":
    main()
