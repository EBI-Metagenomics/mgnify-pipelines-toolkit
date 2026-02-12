#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

"""
gbk_generator.py

Unified converter from:
  - contigs FASTA + CDS GFF3 + proteins FASTA -> GenBank (default mode)
  - contigs FASTA + Prodigal-style proteins FAA (headers carry coords) -> GenBank (Prodigal mode)

Key features:
  - Uses BCBio.GFF for GFF parsing (mandatory requirement in default mode).
  - Adds GenBank-conventional record metadata (molecule_type/topology).
  - Adds a contig-wide "source" feature.
  - Writes CDS features with qualifiers:
        /protein_id /locus_tag /product [/gene optional] [/translation optional]
  - Tags per-CDS provenance from GFF column 2: /note="gene_caller=<source>"
  - Can filter CDS by GFF column 2 (source) via --include-sources
  - Prodigal FAA header parsing with strict validation + reporting for malformed headers.
  - Translation fallback from contig nucleotides; drops terminal stop '*'.
"""

from __future__ import annotations

import argparse
import fileinput
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

PRODIGAL_SPLIT_RE = re.compile(r"\s+#\s+")
TRAILING_INDEX_RE = re.compile(r"^(?P<contig>.+?)_(?P<idx>\d+)$")


# ──────────────────────────────────────────────────────────────────────────────
# IO helpers
# ──────────────────────────────────────────────────────────────────────────────
def _file_is_readable(path: str | None) -> bool:
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


# ──────────────────────────────────────────────────────────────────────────────
# FASTA loaders
# ──────────────────────────────────────────────────────────────────────────────
def read_fasta_records_dict(fasta_path: str) -> dict[str, SeqRecord]:
    records: dict[str, SeqRecord] = {}
    with fileinput.hook_compressed(fasta_path, "rt") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            if rec.id in records:
                raise ValueError(f"Duplicate FASTA record id: {rec.id}")
            rec.seq = Seq(str(rec.seq).replace(" ", "").replace("\t", "").upper())
            records[rec.id] = rec
    return records


def read_protein_seqs_dict(proteins_path: str, drop_terminal_stop: bool) -> dict[str, str]:
    proteins: dict[str, str] = {}
    with fileinput.hook_compressed(proteins_path, "rt") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            prot = str(rec.seq).replace(" ", "").replace("\t", "").strip()
            if drop_terminal_stop and prot.endswith("*"):
                prot = prot[:-1]
            if not prot:
                continue
            proteins[rec.id] = prot
    return proteins


# ──────────────────────────────────────────────────────────────────────────────
# Qualifier helpers
# ──────────────────────────────────────────────────────────────────────────────
def _first_qualifier(feature: SeqFeature, keys: Iterable[str]) -> str | None:
    """
    Extract the first available qualifier value from a feature.
    
    Example qualifiers dict structure:
    {
        'ID': ['gene_001'],
        'locus_tag': ['SAMPLE_001'],
        'product': ['hypothetical protein'],
        'note': ['some annotation note']
    }
    """
    q = getattr(feature, "qualifiers", {}) or {}
    for k in keys:
        if k in q:
            v = q[k]
            if isinstance(v, list):
                return v[0] if v else None
            return v
    return None


def _safe_feature_id(feature: SeqFeature) -> str | None:
    # Prefer explicit GFF3 ID, then locus_tag, then Name.
    return _first_qualifier(feature, ["ID", "locus_tag", "Name", "protein_id"])


def _safe_feature_product(feature: SeqFeature) -> str | None:
    return _first_qualifier(feature, ["product", "product_name", "gene_product", "description"])


def _gff_source(feature: SeqFeature) -> str | None:
    # BCBio maps GFF column 2 into qualifier "source" in most cases.
    return _first_qualifier(feature, ["source"])


def _append_note(feature: SeqFeature, note: str) -> None:
    q = feature.qualifiers or {}
    existing = q.get("note", [])
    if isinstance(existing, str):
        existing = [existing]
    if not isinstance(existing, list):
        existing = [str(existing)]
    existing.append(note)
    q["note"] = existing
    feature.qualifiers = q


# ──────────────────────────────────────────────────────────────────────────────
# Prodigal FAA parsing + strict validation
# ──────────────────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class ProdigalCall:
    gene_id: str
    contig_id: str
    start_1: int  # 1-based inclusive
    end_1: int  # 1-based inclusive
    strand: int  # 1 or -1


def parse_prodigal_header(record_id: str, record_description: str) -> ProdigalCall:
    """
    Expect Prodigal style header with coords and strand:
        <gene_id> # <start> # <end> # <strand> ...

    strand must be 1 or -1. start/end must be >=1.
    If start > end, they are swapped.
    contig_id is inferred from gene_id by stripping a trailing _<idx> if present.
    """
    line = record_description or record_id
    parts = PRODIGAL_SPLIT_RE.split(line)
    if len(parts) < 4:
        raise ValueError(f"FAA header not parseable as Prodigal: {line!r}")

    gene_id_full = parts[0].strip().split()[0]
    start = int(parts[1].strip())
    end = int(parts[2].strip())
    strand = int(parts[3].strip())

    if strand not in (1, -1):
        raise ValueError(f"FAA header has invalid strand (expected 1 or -1): {line!r}")
    if start < 1 or end < 1:
        raise ValueError(f"FAA header has invalid coordinates (<1): {line!r}")
    if start > end:
        start, end = end, start

    m = TRAILING_INDEX_RE.match(gene_id_full)
    contig_id = m.group("contig") if m else gene_id_full

    return ProdigalCall(gene_id=gene_id_full, contig_id=contig_id, start_1=start, end_1=end, strand=strand)


def validate_prodigal_faa_headers(faa_path: str) -> tuple[int, int]:
    """
    Strict validator:
      - counts total records
      - counts records that do NOT follow Prodigal header convention

    Returns (total, bad).
    """
    total = 0
    bad = 0
    with fileinput.hook_compressed(faa_path, "rt") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            total += 1
            try:
                parse_prodigal_header(rec.id, rec.description or rec.id)
            except Exception:
                bad += 1
    return total, bad


def read_calls_from_prodigal_faa(faa_path: str) -> list[ProdigalCall]:
    calls: list[ProdigalCall] = []
    with fileinput.hook_compressed(faa_path, "rt") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            call = parse_prodigal_header(rec.id, rec.description or rec.id)
            calls.append(call)
    return calls


# ──────────────────────────────────────────────────────────────────────────────
# GenBank record construction utilities
# ──────────────────────────────────────────────────────────────────────────────
def make_source_feature(rec: SeqRecord, note: str) -> SeqFeature:
    return SeqFeature(
        location=FeatureLocation(0, len(rec.seq), strand=1),
        type="source",
        qualifiers={"note": [note]},
    )


def ensure_record_annotations(rec: SeqRecord, *, prefix: str) -> None:
    rec.id = rec.id
    rec.name = rec.id
    rec.description = f"{prefix} {rec.id}".strip() if prefix else rec.id
    rec.annotations.setdefault("molecule_type", "DNA")
    rec.annotations.setdefault("topology", "linear")


def drop_terminal_stop_if_requested(prot: str, drop_terminal_stop: bool) -> str:
    if drop_terminal_stop and prot.endswith("*"):
        return prot[:-1]
    return prot


def translate_feature_from_contig(
    rec: SeqRecord,
    feature: SeqFeature,
    *,
    transl_table: int,
    drop_terminal_stop: bool,
) -> str | None:
    # Extract nucleotide sequence of feature and translate.
    nuc = feature.extract(rec.seq)
    nuc = nuc[: (len(nuc) // 3) * 3]
    if len(nuc) == 0:
        return None
    prot = str(nuc.translate(table=transl_table))
    prot = drop_terminal_stop_if_requested(prot, drop_terminal_stop=drop_terminal_stop)
    return prot or None


def should_keep_feature_by_source(
    feature: SeqFeature,
    include_sources: set[str] | None,
) -> bool:
    """
    Determine if a feature should be kept based on its source.
    If include_sources is None, keep all features.
    If include_sources is provided, only keep features whose source is in the set.
    """
    if include_sources is None:
        return True
    
    src = _gff_source(feature)
    return src in include_sources


# ──────────────────────────────────────────────────────────────────────────────
# Default mode (GFF + FAA)
# ──────────────────────────────────────────────────────────────────────────────
def build_records_from_gff_and_faa(
    *,
    contigs_path: str,
    gff_path: str,
    faa_path: str,
    prefix: str,
    default_product: str,
    locus_tag_prefix: str,
    gene_from_locus_tag: bool,
    include_sources: set[str] | None,
    gene_caller_version: str,
    transl_table: int,
    drop_terminal_stop: bool,
) -> list[SeqRecord]:
    contigs = read_fasta_records_dict(contigs_path)

    proteins = read_protein_seqs_dict(faa_path, drop_terminal_stop=drop_terminal_stop)
    logger.info("Loaded proteins: %d", len(proteins))

    records: list[SeqRecord] = []
    with fileinput.hook_compressed(gff_path, "rt") as handle:
        for rec in GFF.parse(handle, base_dict=contigs):
            ensure_record_annotations(rec, prefix=prefix)

            # Add GenBank-conventional source feature first
            features_out: list[SeqFeature] = [make_source_feature(rec, note="generated by gbk_generator.py (GFF+FAA mode)")]

            for feat in getattr(rec, "features", []) or []:
                if feat.type != "CDS":
                    # Keep non-CDS features as-is (optional). If you prefer to drop,
                    # change this branch.
                    features_out.append(feat)
                    continue

                if not should_keep_feature_by_source(feat, include_sources):
                    continue

                src = _gff_source(feat)
                if src:
                    _append_note(feat, f"gene_caller={gene_caller_version}")

                cds_id = _safe_feature_id(feat)
                if not cds_id:
                    # If no usable ID, skip (can also choose to synthesize one)
                    continue

                product = _safe_feature_product(feat) or default_product
                protein_id_val = cds_id
                locus_tag_val = f"{locus_tag_prefix},{cds_id}" if locus_tag_prefix else cds_id

                # Build a new CDS feature with GenBank-conventional qualifiers
                qualifiers: dict[str, list[str]] = {
                    "protein_id": [protein_id_val],
                    "locus_tag": [locus_tag_val],
                    "product": [product],
                }
                if gene_from_locus_tag:
                    qualifiers["gene"] = [locus_tag_val]

                # Carry over provenance notes (and any existing note) to new feature
                note = feat.qualifiers.get("note")
                if note:
                    if isinstance(note, str):
                        qualifiers["note"] = [note]
                    else:
                        qualifiers["note"] = [str(x) for x in note]

                # Prefer protein FASTA translation, else translate from contig
                translation: str | None = None
                if cds_id in proteins:
                    translation = proteins[cds_id]
                else:
                    translation = translate_feature_from_contig(
                        rec,
                        feat,
                        transl_table=transl_table,
                        drop_terminal_stop=drop_terminal_stop,
                    )

                if translation:
                    qualifiers["translation"] = [translation]

                # Use the original location from BCBio feature
                new_feat = SeqFeature(
                    location=feat.location,
                    type="CDS",
                    qualifiers=qualifiers,
                )
                features_out.append(new_feat)

            rec.features = features_out
            records.append(rec)

    if not records:
        raise ValueError("No records parsed from GFF/contigs (check contig IDs and inputs).")

    return records


# ──────────────────────────────────────────────────────────────────────────────
# Prodigal mode (Prodigal FAA only)
# ──────────────────────────────────────────────────────────────────────────────
def build_records_from_prodigal_faa(
    *,
    contigs_path: str,
    prodigal_headers_faa_path: str,
    prefix: str,
    default_product: str,
    locus_tag_prefix: str,
    gene_from_locus_tag: bool,
    skip_missing_contigs: bool,
    require_translation: bool,
    gene_caller_version: str,
    transl_table: int,
    drop_terminal_stop: bool,
    strict_headers: bool,
) -> list[SeqRecord]:
    total, bad = validate_prodigal_faa_headers(prodigal_headers_faa_path)
    logger.info("Prodigal header validation: total=%d bad=%d", total, bad)
    if strict_headers and bad > 0:
        raise ValueError(f"{bad} of {total} FAA records do not follow Prodigal header convention")

    contigs = read_fasta_records_dict(contigs_path)

    # Use prodigal FAA sequences as proteins
    proteins_by_id = read_protein_seqs_dict(prodigal_headers_faa_path, drop_terminal_stop=drop_terminal_stop)
    logger.info("Loaded proteins (from Prodigal FAA): %d", len(proteins_by_id))

    calls = []
    with fileinput.hook_compressed(prodigal_headers_faa_path, "rt") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            try:
                call = parse_prodigal_header(rec.id, rec.description or rec.id)
            except Exception:
                continue
            calls.append(call)

    by_contig: dict[str, list[ProdigalCall]] = {}
    for c in calls:
        if c.contig_id not in contigs:
            if skip_missing_contigs:
                continue
            raise ValueError(f"Prodigal FAA calls reference contig not found in contigs FASTA: {c.contig_id}")
        by_contig.setdefault(c.contig_id, []).append(c)

    out_records: list[SeqRecord] = []
    emitted_any_cds = False

    for contig_id, base_rec in contigs.items():
        rec = base_rec[:]
        ensure_record_annotations(rec, prefix=prefix)

        features_out: list[SeqFeature] = [make_source_feature(rec, note="generated by gbk_generator.py (Prodigal mode)")]

        feats = sorted(by_contig.get(contig_id, []), key=lambda x: (x.start_1, x.end_1, x.gene_id))
        if feats:
            emitted_any_cds = True

        for cds in feats:
            # Convert 1-based inclusive to 0-based half-open
            start0 = cds.start_1 - 1
            end0 = cds.end_1
            strand = 1 if cds.strand == 1 else -1

            protein_id_val = cds.gene_id
            locus_tag_val = f"{locus_tag_prefix},{protein_id_val}" if locus_tag_prefix else protein_id_val
            product_val = default_product

            qualifiers: dict[str, list[str]] = {
                "protein_id": [protein_id_val],
                "locus_tag": [locus_tag_val],
                "product": [product_val],
                "note": [f"gene_caller={gene_caller_version}"],  # provenance in Prodigal mode
            }
            if gene_from_locus_tag:
                qualifiers["gene"] = [locus_tag_val]

            translation: str | None = proteins_by_id.get(cds.gene_id)
            if not translation:
                # Fallback translation from contig
                tmp_feat = SeqFeature(
                    location=FeatureLocation(start0, end0, strand=strand),
                    type="CDS",
                    qualifiers={},
                )
                translation = translate_feature_from_contig(
                    rec,
                    tmp_feat,
                    transl_table=transl_table,
                    drop_terminal_stop=drop_terminal_stop,
                )

            if require_translation and not translation:
                raise ValueError(f"Missing translation for {cds.gene_id} on {cds.contig_id}")

            if translation:
                qualifiers["translation"] = [translation]

            features_out.append(
                SeqFeature(
                    location=FeatureLocation(start0, end0, strand=strand),
                    type="CDS",
                    qualifiers=qualifiers,
                )
            )

        rec.features = features_out
        out_records.append(rec)

    if not emitted_any_cds:
        raise ValueError("No CDS features were emitted (check inputs/headers/contig IDs).")

    return out_records


# ──────────────────────────────────────────────────────────────────────────────
# Writer
# ──────────────────────────────────────────────────────────────────────────────
def write_genbank(records: list[SeqRecord], out_path: str) -> None:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    mode = "wt"
    with fileinput.hook_compressed(out_path, mode) as out_handle:
        SeqIO.write(records, out_handle, "genbank")


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate GenBank from contigs + (GFF+FAA or Prodigal FAA)",
    )

    p.add_argument("--contigs", required=True, help="Contigs FASTA (.fa/.fna/.fasta, optionally .gz)")
    
    # Default mode arguments
    p.add_argument("--gff", required=False, help="CDS GFF3 (optionally .gz). Uses BCBio.GFF. Required for default mode.")
    p.add_argument("--faa", required=False, help="Proteins FASTA for translations (optionally .gz). Required for default mode.")
    
    # Prodigal mode arguments
    p.add_argument("--prodigal-headers-faa", required=False, help="Prodigal proteins FAA with coords in headers (optionally .gz). Enables Prodigal mode.")
    
    p.add_argument("--output_gbk", required=True, help="Output GenBank path (.gb/.gbk, optionally .gz)")

    p.add_argument(
        "--include-sources",
        default=None,
        help="Comma-separated list of allowed GFF column-2 sources (e.g. Pyrodigal,FragGeneScanRS). Applies to CDS only in default mode.",
    )

    p.add_argument("--prefix", default="", help="Prefix used in record description (GenBank DEFINITION-like text)")
    p.add_argument("--default-product", default="hypothetical protein", help="Default /product when absent")
    p.add_argument("--locus-tag-prefix", default="", help="Prefix prepended to locus_tag values")
    p.add_argument("--gene-from-locus-tag", action="store_true", help="Also emit /gene=<locus_tag>")
    p.add_argument("--gene-caller-version", default="unknown", help="Gene caller tool and version for provenance notes")

    # Prodigal mode specific options
    p.add_argument("--skip-missing-contigs", action="store_true", help="In Prodigal mode, drop CDS calls whose contig is missing")
    p.add_argument("--require-translation", action="store_true", help="In Prodigal mode, fail if translation cannot be produced")
    p.add_argument(
        "--strict-prodigal-headers",
        action="store_true",
        help="In Prodigal mode, fail if any FAA headers do not match Prodigal convention. Always logs counts.",
    )

    p.add_argument("--transl-table", type=int, default=11, help="NCBI translation table for nucleotide translation fallback")
    p.add_argument("--drop-terminal-stop", action="store_true", default=True, help="Drop trailing '*' in translations")
    p.add_argument(
        "--no-drop-terminal-stop",
        action="store_false",
        dest="drop_terminal_stop",
        help="Keep trailing '*' in translations",
    )

    p.add_argument("--log-level", default="INFO", help="Logging level (INFO, DEBUG, etc.)")
    return p.parse_args()


def _parse_csv_set(v: str | None) -> set[str] | None:
    if not v:
        return None
    items = [x.strip() for x in v.split(",") if x.strip()]
    return set(items) if items else None


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, str(args.log_level).upper(), logging.INFO))

    include_sources = _parse_csv_set(args.include_sources)

    # Determine processing mode
    prodigal_mode = _file_is_readable(args.prodigal_headers_faa)
    default_mode = _file_is_readable(args.gff) and _file_is_readable(args.faa)

    # Validate input combinations
    if prodigal_mode and default_mode:
        raise ValueError("Cannot use both default mode (--gff + --faa) and Prodigal mode (--prodigal-headers-faa). Choose one mode.")
    
    if prodigal_mode:
        # Prodigal mode: only allow prodigal-headers-faa
        if _file_is_readable(args.gff) or _file_is_readable(args.faa):
            raise ValueError("Prodigal mode: cannot provide --gff or --faa when using --prodigal-headers-faa")
        
        logger.info("Running in Prodigal mode")
        records = build_records_from_prodigal_faa(
            contigs_path=args.contigs,
            prodigal_headers_faa_path=args.prodigal_headers_faa,
            prefix=args.prefix,
            default_product=args.default_product,
            locus_tag_prefix=args.locus_tag_prefix,
            gene_from_locus_tag=args.gene_from_locus_tag,
            skip_missing_contigs=args.skip_missing_contigs,
            require_translation=args.require_translation,
            gene_caller_version=args.gene_caller_version,
            transl_table=args.transl_table,
            drop_terminal_stop=args.drop_terminal_stop,
            strict_headers=args.strict_prodigal_headers,
        )
    
    elif default_mode:
        # Default mode: require both GFF and FAA
        if args.prodigal_headers_faa:
            raise ValueError("Default mode: cannot provide --prodigal-headers-faa when using --gff + --faa")
        
        logger.info("Running in default mode (GFF + FAA)")
        records = build_records_from_gff_and_faa(
            contigs_path=args.contigs,
            gff_path=args.gff,
            faa_path=args.faa,
            prefix=args.prefix,
            default_product=args.default_product,
            locus_tag_prefix=args.locus_tag_prefix,
            gene_from_locus_tag=args.gene_from_locus_tag,
            include_sources=include_sources,
            gene_caller_version=args.gene_caller_version,
            transl_table=args.transl_table,
            drop_terminal_stop=args.drop_terminal_stop,
        )
    
    else:
        raise ValueError("Must provide either:\n"
                        "  Default mode: --gff + --faa\n"
                        "  Prodigal mode: --prodigal-headers-faa")

    write_genbank(records, args.output_gbk)
    logger.info("Wrote GenBank: %s (records=%d)", args.output_gbk, len(records))


if __name__ == "__main__":
    main()
