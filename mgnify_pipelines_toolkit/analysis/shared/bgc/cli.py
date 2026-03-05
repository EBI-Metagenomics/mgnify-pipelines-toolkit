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

import logging
from pathlib import Path

import click

from mgnify_pipelines_toolkit.analysis.shared.bgc.gff_output import build_region_lines
from mgnify_pipelines_toolkit.analysis.shared.bgc.merge import (
    merge_overlaps,
    support_and_filter_cds,
)
from mgnify_pipelines_toolkit.analysis.shared.bgc.models import (
    BGCRegion,
    CDSRec,
    MergedRegion,
)
from mgnify_pipelines_toolkit.analysis.shared.bgc.sideload_json import (
    _derive_output_json_path,
    write_sideload_json,
)
from mgnify_pipelines_toolkit.analysis.shared.bgc.tool_parsers import (
    AntiSMASHParser,
    BGCToolParser,
    GECCOParser,
    SanntiSParser,
)
from mgnify_pipelines_toolkit.analysis.shared.gff.io import iter_gff_rows, write_gff

logger = logging.getLogger(__name__)

_PARSERS: dict[str, BGCToolParser] = {
    "gecco": GECCOParser(),
    "sanntis": SanntiSParser(),
    "antismash": AntiSMASHParser(),
}


def validate_inputs(
    base_gff: Path,
    gecco_gff: Path | None,
    antismash_gff: Path | None,
    sanntis_gff: Path | None,
) -> list[tuple[str, Path]]:
    """Validate input paths and return list of (tool, path) pairs.

    :param base_gff: Path to the mandatory base GFF file.
    :type base_gff: Path
    :param gecco_gff: Optional path to GECCO output GFF.
    :type gecco_gff: Path | None
    :param antismash_gff: Optional path to antiSMASH output GFF.
    :type antismash_gff: Path | None
    :param sanntis_gff: Optional path to SanntiS output GFF.
    :type sanntis_gff: Path | None
    :returns: List of (tool_name, path) tuples for provided optional inputs.
    :rtype: list[tuple[str, Path]]
    :raises FileNotFoundError: If base_gff or any optional GFF does not exist.
    :raises ValueError: If no optional predictor GFF is provided.
    """
    if not base_gff.exists():
        raise FileNotFoundError(f"Base GFF not found: {base_gff}")

    optional: list[tuple[str, Path]] = []
    if gecco_gff:
        optional.append(("gecco", gecco_gff))
    if antismash_gff:
        optional.append(("antismash", antismash_gff))
    if sanntis_gff:
        optional.append(("sanntis", sanntis_gff))

    if not optional:
        raise ValueError("At least one optional predictor GFF must be provided (gecco/antismash/sanntis).")

    for tool, path in optional:
        if not path.exists():
            raise FileNotFoundError(f"{tool} GFF not found: {path}")

    return optional


def load_base_cds(base_gff: Path) -> dict[str, list[CDSRec]]:
    """Load CDS records from the base GFF, preserving original lines.

    :param base_gff: Path to the base GFF file (may be compressed).
    :type base_gff: Path
    :returns: Dictionary mapping contig ID to sorted list of CDSRec objects.
    :rtype: dict[str, list[CDSRec]]
    """
    contig_to_cds: dict[str, list[CDSRec]] = {}
    for line, cols in iter_gff_rows(base_gff):
        if cols[2] != "CDS":
            continue
        contig = cols[0]
        contig_to_cds.setdefault(contig, []).append(CDSRec(contig=contig, start=int(cols[3]), end=int(cols[4]), line=line))

    for contig in contig_to_cds:
        contig_to_cds[contig].sort(key=lambda x: x.start)

    logger.info(
        "Loaded base CDS: contigs=%d cds=%d",
        len(contig_to_cds),
        sum(len(v) for v in contig_to_cds.values()),
    )
    return contig_to_cds


@click.command()
@click.option(
    "--base_gff",
    required=True,
    type=click.Path(path_type=Path),
    help="Mandatory base GFF containing the coordinates of the original CDS used for BGC prediction (may be compressed).",
)
@click.option("--gecco_gff", type=click.Path(path_type=Path), default=None, help="Optional GECCO output GFF (may be compressed).")
@click.option("--antismash_gff", type=click.Path(path_type=Path), default=None, help="Optional antiSMASH output GFF (may be compressed).")
@click.option("--sanntis_gff", type=click.Path(path_type=Path), default=None, help="Optional SanntiS output GFF (may be compressed).")
@click.option("--output_gff", required=True, type=click.Path(path_type=Path), help="Output integrated GFF3.")
@click.option(
    "--validate_json",
    is_flag=True,
    default=False,
    help="Validate the produced antiSMASH sideloader JSON against the local schema copies (dev/debug). Disabled by default.",
)
@click.option(
    "--verbose",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    help="Logging verbosity level.",
)
def main(
    base_gff: Path,
    gecco_gff: Path | None,
    antismash_gff: Path | None,
    sanntis_gff: Path | None,
    output_gff: Path,
    validate_json: bool,
    verbose: str,
) -> None:
    """Integrate optional GECCO/antiSMASH/SanntiS BGC calls into a base GFF3 file."""
    logging.basicConfig(
        level=getattr(logging, verbose),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    try:
        optional_inputs = validate_inputs(base_gff, gecco_gff, antismash_gff, sanntis_gff)
    except Exception as e:
        raise click.ClickException(str(e))

    n_tools = len(optional_inputs)
    contig_to_cds = load_base_cds(base_gff)

    regions: list[BGCRegion] = []
    antismash_gene_ann_by_id: dict[str, dict[str, str]] = {}

    for tool, path in optional_inputs:
        parser = _PARSERS.get(tool)
        if parser is None:
            logger.warning("Unknown tool %s (skipping)", tool)
            continue
        tool_regions, gene_ann = parser.parse_regions(path)
        regions.extend(tool_regions)
        antismash_gene_ann_by_id.update(gene_ann)

    if not regions:
        logger.warning("No BGC regions parsed from optional inputs; writing header only.")
        write_gff(output_gff, [])
        return

    merged_regions: list[MergedRegion] = merge_overlaps(regions)

    cds_lines_by_contig = support_and_filter_cds(
        contig_to_cds=contig_to_cds,
        merged_regions=merged_regions,
        antismash_gene_ann_by_id=antismash_gene_ann_by_id,
        n_tools=n_tools,
    )

    all_lines = list(build_region_lines(merged_regions))
    for contig, cds_lines in cds_lines_by_contig.items():
        for ln in cds_lines:
            all_lines.append((contig, int(ln.split("\t")[3]), ln))
    write_gff(output_gff, all_lines)
    logger.info("Wrote integrated GFF: %s", output_gff)

    # Write antiSMASH-compatible sideloader JSON (subregions only), same basename as output_gff
    out_json = _derive_output_json_path(output_gff)
    try:
        write_sideload_json(out_json, merged_regions, validate=validate_json)
    except Exception as e:
        raise click.ClickException(f"Failed to write/validate sideloader JSON: {e}")

if __name__ == "__main__":
    main()
