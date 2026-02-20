# tests/unit/analysis/shared/test_bgc_mapper.py
#
# End-to-end + targeted unit checks for bgc_mapper.py
#
# Key things we validate:
# - Non-overlapping BGC regions keep EXACT original predictor "source" (column 2), including versions.
# - Overlapping regions are merged and get source="bgc_merged".
# - bgc_support uses N = number of tools provided as input (NOT number of merged predictions).
# - bgc_tools on CDS are DISTINCT tools covering the CDS, sorted, comma-separated.
# - antiSMASH gene-level metadata is applied ONLY when CDS ID matches antiSMASH gene ID.
# - Region-level metadata union is propagated to merged region lines (e.g. antismash_product, gecco_bgc_type, nearest_MiBIG*).
#
# Import module under test:
from __future__ import annotations

from pathlib import Path

from mgnify_pipelines_toolkit.analysis.shared.bgc_mapper import (
    main,
)


def _write_gff(path: Path, rows: list[dict[str, str]]) -> None:
    """
    Write a small GFF3 file from a list of dict rows.
    Each dict should contain:
      seqid, source, type, start, end, score, strand, phase, attributes
    """
    with path.open("w") as fh:
        fh.write("##gff-version 3\n")
        for r in rows:
            fh.write(
                "\t".join(
                    [
                        r["seqid"],
                        r["source"],
                        r["type"],
                        str(r["start"]),
                        str(r["end"]),
                        r.get("score", "."),
                        r.get("strand", "."),
                        r.get("phase", "."),
                        r.get("attributes", "."),
                    ]
                )
                + "\n"
            )


def _read_non_header_lines(path: Path) -> list[str]:
    return [ln.rstrip("\n") for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]


def _find_lines_containing(lines: list[str], needle: str) -> list[str]:
    return [ln for ln in lines if needle in ln]


def test_bgc_mapper_end_to_end_sources_support_and_metadata(tmp_path: Path):
    """
    Build minimal-but-representative inputs inspired by the real files, then run main()
    and assert critical properties of the output.
    """
    # ── Files
    base_gff = tmp_path / "base.gff"
    antismash_gff = tmp_path / "antismash.gff"
    gecco_gff = tmp_path / "gecco.gff"
    sanntis_gff = tmp_path / "sanntis.gff"
    out_gff = tmp_path / "result.gff"

    contig = "b1_04_2.circ"

    # ── Minimal base CDS set: include a few CDS that cover all logic branches
    # CDS IDs chosen from your real result.gff excerpt.
    base_rows = [
        # SanntiS-only region (non-overlapping with others in this minimal set)
        {
            "seqid": contig,
            "source": "Prodigal:002006",
            "type": "CDS",
            "start": 142110,
            "end": 142544,
            "score": ".",
            "strand": "+",
            "phase": "0",
            "attributes": "ID=b1_04_05961;inference=ab initio prediction:Prodigal:002006;locus_tag=b1_04_05961;product=hypothetical protein",
        },
        # antiSMASH-only CDS within region1 but outside SanntiS overlap segment
        {
            "seqid": contig,
            "source": "Prodigal:002006",
            "type": "CDS",
            "start": 211322,
            "end": 211720,
            "score": ".",
            "strand": "-",
            "phase": "0",
            "attributes": "ID=b1_04_06024;inference=ab initio prediction:Prodigal:002006;locus_tag=b1_04_06024;product=hypothetical protein",
        },
        # CDS in overlap between antiSMASH region1 and SanntiS cluster_3 => support 2/3
        {
            "seqid": contig,
            "source": "Prodigal:002006",
            "type": "CDS",
            "start": 217314,
            "end": 218099,
            "score": ".",
            "strand": "+",
            "phase": "0",
            "attributes": "ID=b1_04_06033;inference=ab initio prediction:Prodigal:002006;locus_tag=b1_04_06033;product=Glucose 1-dehydrogenase",
        },
        # CDS in overlap between antiSMASH region2 + SanntiS cluster_4 + GECCO cluster_1 => support 3/3
        {
            "seqid": contig,
            "source": "Prodigal:002006",
            "type": "CDS",
            "start": 241259,
            "end": 242098,
            "score": ".",
            "strand": "-",
            "phase": "0",
            "attributes": "ID=b1_04_06058;inference=ab initio prediction:Prodigal:002006;locus_tag=b1_04_06058;product=IS3 family transposase",
        },
        # A CDS outside any region should be filtered out (not present in output)
        {
            "seqid": contig,
            "source": "Prodigal:002006",
            "type": "CDS",
            "start": 1000,
            "end": 1200,
            "score": ".",
            "strand": "+",
            "phase": "0",
            "attributes": "ID=cds_outside;product=outside",
        },
    ]
    _write_gff(base_gff, base_rows)

    # ── antiSMASH: include region1 + gene b1_04_06033 with gene_functions/as_type/as_gene_clusters
    antismash_rows = [
        {
            "seqid": contig,
            "source": "antiSMASH:8.0.1",
            "type": "region",
            "start": 210661,
            "end": 232040,
            "score": ".",
            "strand": ".",
            "phase": ".",
            "attributes": "ID=b1_04_2.circ_region1;product=sactipeptide",
        },
        {
            "seqid": contig,
            "source": "antiSMASH:8.0.1",
            "type": "gene",
            "start": 211322,
            "end": 211720,
            "score": ".",
            "strand": "-",
            "phase": ".",
            "attributes": "ID=b1_04_06024;as_type=other;Parent=b1_04_2.circ_region1",
        },
        {
            "seqid": contig,
            "source": "antiSMASH:8.0.1",
            "type": "gene",
            "start": 217314,
            "end": 218099,
            "score": ".",
            "strand": "+",
            "phase": ".",
            "attributes": (
                "ID=b1_04_06033;as_type=biosynthetic-additional;"
                "gene_functions=biosynthetic-additional_(rule-based-clusters)_adh_short;"
                "as_gene_clusters=sactipeptide;Parent=b1_04_2.circ_region1"
            ),
        },
        # region2 (for the big merged region with GECCO+SanntiS)
        {
            "seqid": contig,
            "source": "antiSMASH:8.0.1",
            "type": "region",
            "start": 237180,
            "end": 342929,
            "score": ".",
            "strand": ".",
            "phase": ".",
            "attributes": "ID=b1_04_2.circ_region2;product=NRPS,lanthipeptide-class-ii,T1PKS",
        },
    ]
    _write_gff(antismash_gff, antismash_rows)

    # ── GECCO: in your mapper you parse any 9-col row with Type as a "region"
    gecco_rows = [
        {
            "seqid": contig,
            "source": "GECCO v0.9.8",
            "type": "BGC",
            "start": 241259,
            "end": 308608,
            "score": "0.99",
            "strand": ".",
            "phase": ".",
            "attributes": "ID=b1_04_2.circ_cluster_1;Type=NRP,Polyketide",
        }
    ]
    _write_gff(gecco_gff, gecco_rows)

    # ── SanntiS: in your mapper you parse any 9-col row with nearest_MiBIG* attrs as a "region"
    # Include:
    #  - a SanntiS-only region (142110-148437) to check non-overlapping source preservation
    #  - a region overlapping antiSMASH region1 (217314-226605) to form a merged region
    #  - a region overlapping antiSMASH region2 and GECCO (234258-330050) to form a big merged region
    sanntis_rows = [
        {
            "seqid": contig,
            "source": "SanntiSv0.9.3.3",
            "type": "CLUSTER",
            "start": 142110,
            "end": 148437,
            "score": ".",
            "strand": ".",
            "phase": ".",
            "attributes": "ID=b1_04_2.circ_sanntis_1;nearest_MiBIG=BGC0000849;nearest_MiBIG_class=Other",
        },
        {
            "seqid": contig,
            "source": "SanntiSv0.9.3.3",
            "type": "CLUSTER",
            "start": 217314,
            "end": 226605,
            "score": ".",
            "strand": ".",
            "phase": ".",
            "attributes": "ID=b1_04_2.circ_sanntis_3;nearest_MiBIG=BGC0000600;nearest_MiBIG_class=RiPP",
        },
        {
            "seqid": contig,
            "source": "SanntiSv0.9.3.3",
            "type": "CLUSTER",
            "start": 234258,
            "end": 330050,
            "score": ".",
            "strand": ".",
            "phase": ".",
            "attributes": "ID=b1_04_2.circ_sanntis_4;nearest_MiBIG=BGC0001059;nearest_MiBIG_class=NRP Polyketide",
        },
    ]
    _write_gff(sanntis_gff, sanntis_rows)

    # ── Run main with ALL THREE tools present => N=3 for bgc_support
    argv = [
        "--base_gff",
        str(base_gff),
        "--antismash_gff",
        str(antismash_gff),
        "--gecco_gff",
        str(gecco_gff),
        "--sanntis_gff",
        str(sanntis_gff),
        "--output_gff",
        str(out_gff),
        "--log_level",
        "ERROR",
    ]
    rc = main(argv)
    assert rc == 0
    assert out_gff.exists()

    lines = _read_non_header_lines(out_gff)

    # ──────────────────────────────────────────────────────────────────────────
    # 1) Non-overlapping region keeps EXACT source string
    # ──────────────────────────────────────────────────────────────────────────
    # Expect the SanntiS-only region line to have column 2 == "SanntiSv0.9.3.3" (not "sanntis")
    sanntis_only_region = [ln for ln in lines if ln.startswith(f"{contig}\tSanntiSv0.9.3.3\tbgc_region\t142110\t148437\t")]
    assert len(sanntis_only_region) == 1
    assert "bgc_tools=sanntis" in sanntis_only_region[0]
    assert "nearest_MiBIG=BGC0000849" in sanntis_only_region[0]

    # ──────────────────────────────────────────────────────────────────────────
    # 2) Overlapping regions merge => source "bgc_merged"
    # ──────────────────────────────────────────────────────────────────────────
    merged_region1 = [ln for ln in lines if ln.startswith(f"{contig}\tbgc_merged\tbgc_region\t210661\t232040\t")]
    assert len(merged_region1) == 1
    # should mention both tools in region-level tools set
    assert "bgc_tools=antismash,sanntis" in merged_region1[0]
    assert "member_bgcs=2" in merged_region1[0]
    assert "antismash_product=sactipeptide" in merged_region1[0]
    assert "nearest_MiBIG=BGC0000600" in merged_region1[0]

    merged_region2 = [ln for ln in lines if ln.startswith(f"{contig}\tbgc_merged\tbgc_region\t234258\t342929\t")]
    assert len(merged_region2) == 1
    assert "bgc_tools=antismash,gecco,sanntis" in merged_region2[0]
    # union of metadata keys
    assert "antismash_product=NRPS,lanthipeptide-class-ii,T1PKS" in merged_region2[0]
    assert "gecco_bgc_type=NRP,Polyketide" in merged_region2[0]
    assert "nearest_MiBIG=BGC0001059" in merged_region2[0]

    # ──────────────────────────────────────────────────────────────────────────
    # 3) bgc_support values use N=3 (tools provided), and tool coverage is distinct tools
    # ──────────────────────────────────────────────────────────────────────────
    # SanntiS-only CDS => 1/3 => 0.33
    cds_05961 = _find_lines_containing(lines, "ID=b1_04_05961;")
    assert len(cds_05961) == 1
    assert "bgc_tools=sanntis" in cds_05961[0]
    assert "bgc_support=0.33" in cds_05961[0]

    # antiSMASH-only CDS (even though region is merged with SanntiS, this CDS is outside the SanntiS interval)
    # => tools covering CDS == {antismash} => 1/3 => 0.33
    cds_06024 = _find_lines_containing(lines, "ID=b1_04_06024;")
    assert len(cds_06024) == 1
    assert "bgc_tools=antismash" in cds_06024[0]
    assert "bgc_support=0.33" in cds_06024[0]
    # and should inherit region product via antiSMASH gene->parent region mapping
    assert "antismash_product=sactipeptide" in cds_06024[0]
    assert "antismash_as_type=other" in cds_06024[0]

    # antiSMASH + SanntiS covered CDS => 2/3 => 0.67 and gene-level metadata present
    cds_06033 = _find_lines_containing(lines, "ID=b1_04_06033;")
    assert len(cds_06033) == 1
    assert "bgc_tools=antismash,sanntis" in cds_06033[0]
    assert "bgc_support=0.67" in cds_06033[0]
    assert "antismash_gene_function=biosynthetic-additional_(rule-based-clusters)_adh_short" in cds_06033[0]
    assert "antismash_as_type=biosynthetic-additional" in cds_06033[0]
    assert "antismash_as_gene_clusters=sactipeptide" in cds_06033[0]
    assert "antismash_product=sactipeptide" in cds_06033[0]
    # SanntiS meta should also be present
    assert "nearest_MiBIG=BGC0000600" in cds_06033[0]
    assert "nearest_MiBIG_class=RiPP" in cds_06033[0]

    # antiSMASH + GECCO + SanntiS covered CDS => 3/3 => 1.00 and GECCO meta present
    cds_06058 = _find_lines_containing(lines, "ID=b1_04_06058;")
    assert len(cds_06058) == 1
    assert "bgc_tools=antismash,gecco,sanntis" in cds_06058[0]
    assert "bgc_support=1.00" in cds_06058[0]
    assert "gecco_bgc_type=NRP,Polyketide" in cds_06058[0]
    assert "antismash_product=NRPS,lanthipeptide-class-ii,T1PKS" in cds_06058[0]
    assert "nearest_MiBIG=BGC0001059" in cds_06058[0]

    # ──────────────────────────────────────────────────────────────────────────
    # 4) Filter rule: CDS outside any BGC region must not appear
    # ──────────────────────────────────────────────────────────────────────────
    assert not _find_lines_containing(lines, "ID=cds_outside;")


def test_bgc_support_denominator_is_number_of_tools_provided(tmp_path: Path):
    """
    Explicitly confirm N changes when fewer tool inputs are provided.

    With only SanntiS provided (N=1), a SanntiS-covered CDS must get bgc_support=1.00.
    """
    base_gff = tmp_path / "base.gff"
    sanntis_gff = tmp_path / "sanntis.gff"
    out_gff = tmp_path / "result.gff"

    contig = "b1_04_2.circ"

    _write_gff(
        base_gff,
        [
            {
                "seqid": contig,
                "source": "Prodigal:002006",
                "type": "CDS",
                "start": 142110,
                "end": 142544,
                "score": ".",
                "strand": "+",
                "phase": "0",
                "attributes": "ID=b1_04_05961;product=hypothetical protein",
            }
        ],
    )

    _write_gff(
        sanntis_gff,
        [
            {
                "seqid": contig,
                "source": "SanntiSv0.9.3.3",
                "type": "CLUSTER",
                "start": 142110,
                "end": 148437,
                "score": ".",
                "strand": ".",
                "phase": ".",
                "attributes": "ID=b1_04_2.circ_sanntis_1;nearest_MiBIG=BGC0000849;nearest_MiBIG_class=Other",
            }
        ],
    )

    argv = [
        "--base_gff",
        str(base_gff),
        "--sanntis_gff",
        str(sanntis_gff),
        "--output_gff",
        str(out_gff),
        "--log_level",
        "ERROR",
    ]
    rc = main(argv)
    assert rc == 0

    lines = _read_non_header_lines(out_gff)
    cds_05961 = _find_lines_containing(lines, "ID=b1_04_05961;")
    assert len(cds_05961) == 1
    # N=1 => 1/1
    assert "bgc_support=1.00" in cds_05961[0]
    assert "bgc_tools=sanntis" in cds_05961[0]


def test_fails_if_no_optional_predictor_is_provided(tmp_path: Path):
    """
    validate_inputs() is called by main; if no predictor is provided, main should return 2.
    """
    base_gff = tmp_path / "base.gff"
    out_gff = tmp_path / "result.gff"

    _write_gff(
        base_gff,
        [
            {
                "seqid": "b1_04_2.circ",
                "source": "Prodigal:002006",
                "type": "CDS",
                "start": 1,
                "end": 100,
                "score": ".",
                "strand": "+",
                "phase": "0",
                "attributes": "ID=cds1",
            }
        ],
    )

    argv = [
        "--base_gff",
        str(base_gff),
        "--output_gff",
        str(out_gff),
        "--log_level",
        "ERROR",
    ]
    rc = main(argv)
    assert rc == 2
    assert not out_gff.exists()
