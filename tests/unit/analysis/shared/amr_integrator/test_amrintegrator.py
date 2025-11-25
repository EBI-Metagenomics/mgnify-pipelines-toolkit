#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024-2025 EMBL - European Bioinformatics Institute
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


import argparse
import pytest
import pandas as pd
import gzip
from pathlib import Path
from unittest.mock import patch
from typing import Dict, List
from collections import defaultdict

from mgnify_pipelines_toolkit.analysis.shared.amrintegrator import (
    validate_inputs,
    normalize_drug_class,
    parse_hamronized,
    parse_amrfinderplus,
    parse_amr_dict,
    parse_gff,
    main,
)


# Test data structures for DeepARG hamronized output
DEEPARG_TEST_RECORDS: List[Dict[str, str]] = [
    {
        "input_file_name": "b20_05",
        "gene_symbol": "ROSA",
        "gene_name": "YP_001401993|FEATURES|rosA|fosmidomycin|rosA",
        "reference_database_name": "deeparg_db",
        "reference_database_version": "2",
        "reference_accession": "YP_001401993",
        "analysis_software_name": "deeparg",
        "analysis_software_version": "1.0.4",
        "genetic_variation_type": "gene_presence_detected",
        "antimicrobial_agent": "",
        "coverage_percentage": "",
        "coverage_depth": "",
        "coverage_ratio": "",
        "drug_class": "fosmidomycin",
        "input_gene_length": "",
        "input_gene_start": "10",
        "input_gene_stop": "401",
        "input_protein_length": "",
        "input_protein_start": "",
        "input_protein_stop": "",
        "input_sequence_id": "b20_05_1.circ_3209",
        "nucleotide_mutation": "",
        "nucleotide_mutation_interpretation": "",
        "predicted_phenotype": "",
        "predicted_phenotype_confidence_level": "",
        "amino_acid_mutation": "",
        "amino_acid_mutation_interpretation": "",
        "reference_gene_length": "",
        "reference_gene_start": "",
        "reference_gene_stop": "",
        "reference_protein_length": "",
        "reference_protein_start": "",
        "reference_protein_stop": "",
        "resistance_mechanism": "",
        "strand_orientation": "",
        "sequence_identity": "54.3",
    },
    {
        "input_file_name": "b20_05",
        "gene_symbol": "ERMF",
        "gene_name": "AAA27431|FEATURES|ermF|MLS|ermF",
        "reference_database_name": "deeparg_db",
        "reference_database_version": "2",
        "reference_accession": "AAA27431",
        "analysis_software_name": "deeparg",
        "analysis_software_version": "1.0.4",
        "genetic_variation_type": "gene_presence_detected",
        "antimicrobial_agent": "",
        "coverage_percentage": "",
        "coverage_depth": "",
        "coverage_ratio": "",
        "drug_class": "MLS",
        "input_gene_length": "",
        "input_gene_start": "1",
        "input_gene_stop": "266",
        "input_protein_length": "",
        "input_protein_start": "",
        "input_protein_stop": "",
        "input_sequence_id": "b20_05_1.circ_2868",
        "nucleotide_mutation": "",
        "nucleotide_mutation_interpretation": "",
        "predicted_phenotype": "",
        "predicted_phenotype_confidence_level": "",
        "amino_acid_mutation": "",
        "amino_acid_mutation_interpretation": "",
        "reference_gene_length": "",
        "reference_gene_start": "",
        "reference_gene_stop": "",
        "reference_protein_length": "",
        "reference_protein_start": "",
        "reference_protein_stop": "",
        "resistance_mechanism": "",
        "strand_orientation": "",
        "sequence_identity": "100.0",
    },
    {
        "input_file_name": "b20_05",
        "gene_symbol": "TETX",
        "gene_name": "gi:1004705391:gb:AMP52595.1:|FEATURES|tetX|tetracycline|tetX",
        "reference_database_name": "deeparg_db",
        "reference_database_version": "2",
        "reference_accession": "gi:1004705391:gb:AMP52595.1:",
        "analysis_software_name": "deeparg",
        "analysis_software_version": "1.0.4",
        "genetic_variation_type": "gene_presence_detected",
        "antimicrobial_agent": "",
        "coverage_percentage": "",
        "coverage_depth": "",
        "coverage_ratio": "",
        "drug_class": "tetracycline",
        "input_gene_length": "",
        "input_gene_start": "1",
        "input_gene_stop": "388",
        "input_protein_length": "",
        "input_protein_start": "",
        "input_protein_stop": "",
        "input_sequence_id": "b20_05_1.circ_2148",
        "nucleotide_mutation": "",
        "nucleotide_mutation_interpretation": "",
        "predicted_phenotype": "",
        "predicted_phenotype_confidence_level": "",
        "amino_acid_mutation": "",
        "amino_acid_mutation_interpretation": "",
        "reference_gene_length": "",
        "reference_gene_start": "",
        "reference_gene_stop": "",
        "reference_protein_length": "",
        "reference_protein_start": "",
        "reference_protein_stop": "",
        "resistance_mechanism": "",
        "strand_orientation": "",
        "sequence_identity": "99.7",
    },
]

# Test data structures for RGI hamronized output
RGI_TEST_RECORDS: List[Dict[str, str]] = [
    {
        "input_file_name": "b20_05",
        "gene_symbol": "tet(Q)",
        "gene_name": "tetracycline-resistant ribosomal protection protein",
        "reference_database_name": "CARD",
        "reference_database_version": "4.0.1",
        "reference_accession": "3000191",
        "analysis_software_name": "rgi",
        "analysis_software_version": "6.0.5",
        "genetic_variation_type": "gene_presence_detected",
        "antimicrobial_agent": "tetracycline; doxycycline; minocycline; chlortetracycline; demeclocycline; oxytetracycline",
        "coverage_percentage": "97.56",
        "coverage_depth": "",
        "coverage_ratio": "",
        "drug_class": "tetracycline antibiotic",
        "input_gene_length": "",
        "input_gene_start": "",
        "input_gene_stop": "",
        "input_protein_length": "",
        "input_protein_start": "",
        "input_protein_stop": "",
        "input_sequence_id": "b20_05_1.circ_2009 # 2925164 # 2927089 # -1 # ID=1_2009;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=5bp;gc_cont=0.398",
        "nucleotide_mutation": "",
        "nucleotide_mutation_interpretation": "",
        "predicted_phenotype": "",
        "predicted_phenotype_confidence_level": "",
        "amino_acid_mutation": "",
        "amino_acid_mutation_interpretation": "",
        "reference_gene_length": "",
        "reference_gene_start": "",
        "reference_gene_stop": "",
        "reference_protein_length": "",
        "reference_protein_start": "",
        "reference_protein_stop": "",
        "resistance_mechanism": "antibiotic target protection",
        "strand_orientation": "",
        "sequence_identity": "96.41",
    },
    {
        "input_file_name": "b20_05",
        "gene_symbol": "ErmF",
        "gene_name": "Erm 23S ribosomal RNA methyltransferase",
        "reference_database_name": "CARD",
        "reference_database_version": "4.0.1",
        "reference_accession": "3000498",
        "analysis_software_name": "rgi",
        "analysis_software_version": "6.0.5",
        "genetic_variation_type": "gene_presence_detected",
        "antimicrobial_agent": "erythromycin; roxithromycin; lincomycin; telithromycin; clarithromycin; clindamycin; tylosin; spiramycin; azithromycin; dirithromycin; pristinamycin IA; quinupristin; virginiamycin M1; madumycin II; griseoviridin; dalfopristin; pristinamycin IB; virginiamycin S2; pristinamycin IC; vernamycin C; patricin A; patricin B; ostreogrycin B3; oleandomycin",
        "coverage_percentage": "100.0",
        "coverage_depth": "",
        "coverage_ratio": "",
        "drug_class": "macrolide antibiotic; lincosamide antibiotic; streptogramin antibiotic; streptogramin A antibiotic; streptogramin B antibiotic",
        "input_gene_length": "",
        "input_gene_start": "",
        "input_gene_stop": "",
        "input_protein_length": "",
        "input_protein_start": "",
        "input_protein_stop": "",
        "input_sequence_id": "b20_05_1.circ_2149 # 3109927 # 3110727 # -1 # ID=1_2149;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=10bp;gc_cont=0.341",
        "nucleotide_mutation": "",
        "nucleotide_mutation_interpretation": "",
        "predicted_phenotype": "",
        "predicted_phenotype_confidence_level": "",
        "amino_acid_mutation": "",
        "amino_acid_mutation_interpretation": "",
        "reference_gene_length": "",
        "reference_gene_start": "",
        "reference_gene_stop": "",
        "reference_protein_length": "",
        "reference_protein_start": "",
        "reference_protein_stop": "",
        "resistance_mechanism": "antibiotic target alteration",
        "strand_orientation": "",
        "sequence_identity": "99.62",
    },
    {
        "input_file_name": "b20_05",
        "gene_symbol": "tet(X)",
        "gene_name": "tetracycline inactivation enzyme",
        "reference_database_name": "CARD",
        "reference_database_version": "4.0.1",
        "reference_accession": "3000205",
        "analysis_software_name": "rgi",
        "analysis_software_version": "6.0.5",
        "genetic_variation_type": "gene_presence_detected",
        "antimicrobial_agent": "tigecycline; tetracycline; doxycycline; minocycline; chlortetracycline; demeclocycline; oxytetracycline",
        "coverage_percentage": "100.0",
        "coverage_depth": "",
        "coverage_ratio": "",
        "drug_class": "glycylcycline; tetracycline antibiotic",
        "input_gene_length": "",
        "input_gene_start": "",
        "input_gene_stop": "",
        "input_protein_length": "",
        "input_protein_start": "",
        "input_protein_stop": "",
        "input_sequence_id": "b20_05_1.circ_2148 # 3108548 # 3109714 # 1 # ID=1_2148;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.374",
        "nucleotide_mutation": "",
        "nucleotide_mutation_interpretation": "",
        "predicted_phenotype": "",
        "predicted_phenotype_confidence_level": "",
        "amino_acid_mutation": "",
        "amino_acid_mutation_interpretation": "",
        "reference_gene_length": "",
        "reference_gene_start": "",
        "reference_gene_stop": "",
        "reference_protein_length": "",
        "reference_protein_start": "",
        "reference_protein_stop": "",
        "resistance_mechanism": "antibiotic inactivation",
        "strand_orientation": "",
        "sequence_identity": "99.74",
    },
]

# Test data structures for AMRFinderPlus output
AMRFP_TEST_RECORDS: List[Dict[str, str]] = [
    {
        "Protein id": "b20_05_1.circ_2009",
        "Element symbol": "tet(Q)",
        "Element name": "tetracycline resistance ribosomal protection protein Tet(Q)",
        "Scope": "core",
        "Type": "AMR",
        "Subtype": "AMR",
        "Class": "TETRACYCLINE",
        "Subclass": "TETRACYCLINE",
        "Method": "BLASTP",
        "Target length": "641",
        "Reference sequence length": "641",
        "% Coverage of reference": "100.00",
        "% Identity to reference": "99.53",
        "Alignment length": "641",
        "Closest reference accession": "WP_063856407.1",
        "Closest reference name": "tetracycline resistance ribosomal protection protein Tet(Q)",
        "HMM accession": "NF012154.0",
        "HMM description": "tetracycline resistance ribosomal protection protein Tet(Q)",
    },
    {
        "Protein id": "b20_05_1.circ_2149",
        "Element symbol": "erm(F)",
        "Element name": "23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)",
        "Scope": "core",
        "Type": "AMR",
        "Subtype": "AMR",
        "Class": "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN",
        "Subclass": "CLINDAMYCIN/MACROLIDE/STREPTOGRAMIN",
        "Method": "EXACTP",
        "Target length": "266",
        "Reference sequence length": "266",
        "% Coverage of reference": "100.00",
        "% Identity to reference": "100.00",
        "Alignment length": "266",
        "Closest reference accession": "WP_002682030.1",
        "Closest reference name": "23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)",
        "HMM accession": "NF012223.0",
        "HMM description": "23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)",
    },
    {
        "Protein id": "b20_05_1.circ_3091",
        "Element symbol": "bla",
        "Element name": "class A beta-lactamase, subclass A2",
        "Scope": "core",
        "Type": "AMR",
        "Subtype": "AMR",
        "Class": "BETA-LACTAM",
        "Subclass": "BETA-LACTAM",
        "Method": "HMM",
        "Target length": "301",
        "Reference sequence length": "296",
        "% Coverage of reference": "99.66",
        "% Identity to reference": "59.53",
        "Alignment length": "299",
        "Closest reference accession": "WP_005827792.1",
        "Closest reference name": "class A beta-lactamase, subclass A2",
        "HMM accession": "NF012099.1",
        "HMM description": "class A beta-lactamase, subclass A2",
    },
    {
        "Protein id": "b20_05_1.circ_2148",
        "Element symbol": "tet(X2)",
        "Element name": "tetracycline-inactivating monooxygenase Tet(X2)",
        "Scope": "core",
        "Type": "AMR",
        "Subtype": "AMR",
        "Class": "TETRACYCLINE",
        "Subclass": "TIGECYCLINE",
        "Method": "BLASTP",
        "Target length": "388",
        "Reference sequence length": "388",
        "% Coverage of reference": "100.00",
        "% Identity to reference": "99.74",
        "Alignment length": "388",
        "Closest reference accession": "WP_008651082.1",
        "Closest reference name": "tetracycline-inactivating monooxygenase Tet(X2)",
        "HMM accession": "NF033111.2",
        "HMM description": "tetracycline-inactivating monooxygenase Tet(X)",
    },
]

# Test data structures for GFF files
GFF_TEST_RECORDS: List[Dict[str, str]] = [
    {
        "seqname": "b20_05_1.circ",
        "source": "Prodigal_v2.6.3",
        "feature": "CDS",
        "start": "2925164",
        "end": "2927089",
        "score": "152.2",
        "strand": "-",
        "frame": "0",
        "attribute": "ID=b20_05_1.circ_2009;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=5bp;gc_cont=0.398;conf=100.00;score=152.15;cscore=148.05;sscore=4.10;rscore=-0.03;uscore=0.76;tscore=4.02",
    },
    {
        "seqname": "b20_05_1.circ",
        "source": "Prodigal_v2.6.3",
        "feature": "CDS",
        "start": "3108548",
        "end": "3109714",
        "score": "78.3",
        "strand": "+",
        "frame": "0",
        "attribute": "ID=b20_05_1.circ_2148;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.374;conf=100.00;score=78.33;cscore=89.46;sscore=-11.13;rscore=-9.55;uscore=-5.60;tscore=4.02",
    },
    {
        "seqname": "b20_05_1.circ",
        "source": "Prodigal_v2.6.3",
        "feature": "CDS",
        "start": "3109927",
        "end": "3110727",
        "score": "38.1",
        "strand": "-",
        "frame": "0",
        "attribute": "ID=b20_05_1.circ_2149;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=10bp;gc_cont=0.341;conf=99.98;score=38.14;cscore=31.14;sscore=7.00;rscore=-0.03;uscore=3.01;tscore=4.02",
    },
    {
        "seqname": "b20_05_1.circ",
        "source": "Prodigal_v2.6.3",
        "feature": "CDS",
        "start": "3798391",
        "end": "3799731",
        "score": "191.3",
        "strand": "+",
        "frame": "0",
        "attribute": "ID=b20_05_1.circ_2656;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=11bp;gc_cont=0.398;conf=99.99;score=191.32;cscore=179.05;sscore=12.27;rscore=8.19;uscore=-0.65;tscore=4.02",
    },
]


def create_hamronized_tsv(output_file: Path, records: List[Dict[str, str]]) -> None:
    """Convert hamronized test records to TSV format."""
    df = pd.DataFrame(records)
    df.to_csv(output_file, sep="\t", index=False)


def create_amrfp_tsv(output_file: Path, records: List[Dict[str, str]]) -> None:
    """Convert AMRFinderPlus test records to TSV format."""
    df = pd.DataFrame(records)
    df.to_csv(output_file, sep="\t", index=False)


def create_gff_file(output_file: Path, records: List[Dict[str, str]]) -> None:
    """Convert GFF test records to GFF format."""
    with open(output_file, "w") as f:
        f.write("##gff-version 3\n")
        for record in records:
            line = "\t".join(
                [
                    record["seqname"],
                    record["source"],
                    record["feature"],
                    record["start"],
                    record["end"],
                    record["score"],
                    record["strand"],
                    record["frame"],
                    record["attribute"],
                ]
            )
            f.write(line + "\n")


class TestValidateInputs:
    """Test the validate_inputs function."""

    def test_validate_inputs_all_exist(self, tmp_path):
        """Test validation when all files exist."""
        # Create temporary files
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file1.write_text("header\ndata")
        file2.write_text("header\ndata")

        optional_inputs = {"deeparg": str(file1), "rgi": str(file2)}

        result = validate_inputs(optional_inputs)
        assert result == ["deeparg", "rgi"]

    def test_validate_inputs_some_exist(self, tmp_path):
        """Test validation when some files exist."""
        file1 = tmp_path / "file1.txt"
        file1.write_text("header\ndata")

        optional_inputs = {"deeparg": str(file1), "rgi": "/nonexistent/path.txt", "amrfinderplus": None}

        result = validate_inputs(optional_inputs)
        assert result == ["deeparg"]

    def test_validate_inputs_none_exist(self):
        """Test validation when no files exist."""
        optional_inputs = {"deeparg": "/nonexistent1.txt", "rgi": "/nonexistent2.txt", "amrfinderplus": None}

        result = validate_inputs(optional_inputs)
        assert result == []


class TestNormalizeDrugClass:
    """Test the normalize_drug_class function."""

    def test_mls_replacement(self):
        """Test MLS drug class replacement."""
        result = normalize_drug_class("macrolide antibiotic; lincosamide antibiotic; streptogramin antibiotic")
        assert result == "MLS"

    def test_antibiotic_suffix_removal(self):
        """Test antibiotic suffix removal."""
        result = normalize_drug_class("tetracycline antibiotic")
        assert result == "tetracycline"

    def test_no_change_needed(self):
        """Test drug class that doesn't need normalization."""
        result = normalize_drug_class("tetracycline")
        assert result == "tetracycline"


class TestParseHamronized:
    """Test the parse_hamronized function for parsing DeepARG and RGI outputs."""

    def test_parse_deeparg_hamronized(self, tmp_path):
        """Test parsing DeepARG hamronized output."""
        # Create test file using structured data
        deeparg_file = tmp_path / "deeparg.tsv"
        create_hamronized_tsv(deeparg_file, DEEPARG_TEST_RECORDS)

        amr_annotation = defaultdict(dict)
        result = parse_hamronized(amr_annotation, str(deeparg_file))

        # Verify results
        assert "b20_05_1.circ_3209" in result
        assert "b20_05_1.circ_2868" in result

        assert result["b20_05_1.circ_3209"]["deeparg"]["drug_class"] == ["fosmidomycin"]
        assert result["b20_05_1.circ_3209"]["deeparg"]["seq_identity"] == "54.3"

        assert result["b20_05_1.circ_2868"]["deeparg"]["drug_class"] == ["MLS"]
        assert result["b20_05_1.circ_2868"]["deeparg"]["seq_identity"] == "100.0"

    def test_parse_rgi_hamronized(self, tmp_path):
        """Test parsing RGI hamronized output."""
        # Create test file using structured data
        rgi_file = tmp_path / "rgi.tsv"
        create_hamronized_tsv(rgi_file, RGI_TEST_RECORDS)

        amr_annotation = defaultdict(dict)
        result = parse_hamronized(amr_annotation, str(rgi_file))

        # Verify results
        assert "b20_05_1.circ_2009" in result
        assert "b20_05_1.circ_2149" in result

        assert result["b20_05_1.circ_2009"]["rgi"]["drug_class"] == ["tetracycline"]
        assert result["b20_05_1.circ_2009"]["rgi"]["seq_identity"] == "96.41"

        # Test MLS replacement for complex drug classes
        assert result["b20_05_1.circ_2149"]["rgi"]["drug_class"] == ["MLS", "streptogramin A", "streptogramin B"]
        assert result["b20_05_1.circ_2149"]["rgi"]["seq_identity"] == "99.62"


class TestParseAmrfinderplus:
    """Test the parse_amrfinderplus function."""

    def test_parse_amrfinderplus(self, tmp_path):
        """Test parsing AMRFinderPlus output."""
        # Create test file using structured data
        amrfp_file = tmp_path / "amrfp.tsv"
        create_amrfp_tsv(amrfp_file, AMRFP_TEST_RECORDS)

        amr_annotation = defaultdict(dict)
        result = parse_amrfinderplus(amr_annotation, str(amrfp_file))

        # Verify results
        assert "b20_05_1.circ_2009" in result
        assert "b20_05_1.circ_2149" in result
        assert "b20_05_1.circ_3091" in result

        assert result["b20_05_1.circ_2009"]["amrfinderplus"]["drug_class"] == ["tetracycline"]
        assert result["b20_05_1.circ_2009"]["amrfinderplus"]["seq_identity"] == "99.53"

        # Test MLS replacement
        assert result["b20_05_1.circ_2149"]["amrfinderplus"]["drug_class"] == ["MLS"]
        assert result["b20_05_1.circ_2149"]["amrfinderplus"]["seq_identity"] == "100.00"

        assert result["b20_05_1.circ_3091"]["amrfinderplus"]["drug_class"] == ["beta-lactam"]
        assert result["b20_05_1.circ_3091"]["amrfinderplus"]["seq_identity"] == "59.53"


class TestParseAmrDict:
    """Test the parse_amr_dict function for processing annotations."""

    def test_parse_amr_dict_single_tool(self):
        """Test parsing annotations from a single tool."""
        amr_annotation = {"protein1": {"rgi": {"drug_class": ["tetracycline", "glycylcycline"], "seq_identity": "99.74"}}}

        result = parse_amr_dict(amr_annotation)

        expected = {"protein1": ["drug_class=tetracycline,glycylcycline", "amr_tool=rgi", "amr_tool_ident=99.74"]}

        assert result == expected

    def test_parse_amr_dict_multiple_tools(self):
        """Test parsing annotations from multiple tools."""
        amr_annotation = {
            "protein1": {
                "deeparg": {"drug_class": ["tetracycline"], "seq_identity": "99.7"},
                "rgi": {"drug_class": ["tetracycline", "glycylcycline"], "seq_identity": "99.74"},
                "amrfinderplus": {"drug_class": ["tetracycline"], "seq_identity": "99.74"},
            }
        }

        result = parse_amr_dict(amr_annotation)

        # Check that all unique drug classes are included
        drug_class_str = result["protein1"][0]
        assert "tetracycline" in drug_class_str
        assert "glycylcycline" in drug_class_str

        # Check tool names
        tools_str = result["protein1"][1]
        assert "deeparg" in tools_str
        assert "rgi" in tools_str
        assert "amrfinderplus" in tools_str

        # Check identities
        idents_str = result["protein1"][2]
        assert "99.7" in idents_str
        assert "99.74" in idents_str

    def test_parse_amr_dict_drug_class_cleanup(self):
        """Test drug class cleanup (spaces to underscores)."""
        amr_annotation = {"protein1": {"rgi": {"drug_class": ["streptogramin A", "beta lactam"], "seq_identity": "95.0"}}}

        result = parse_amr_dict(amr_annotation)

        assert "streptogramin_A" in result["protein1"][0]
        assert "beta_lactam" in result["protein1"][0]


class TestParseGff:
    """Test the parse_gff function for integrating annotations into GFF."""

    def test_parse_gff_integration(self, tmp_path):
        """Test integration of AMR annotations into GFF file."""
        # Create test GFF file using structured data
        gff_file = tmp_path / "input.gff"
        create_gff_file(gff_file, GFF_TEST_RECORDS)

        # Test protein attributes
        protein_attributes = {
            "b20_05_1.circ_2009": ["drug_class=tetracycline", "amr_tool=rgi,amrfinderplus", "amr_tool_ident=96.41,99.53"],
            "b20_05_1.circ_2148": ["drug_class=tetracycline,glycylcycline", "amr_tool=deeparg,rgi,amrfinderplus", "amr_tool_ident=99.7,99.74,99.74"],
        }

        output_file = tmp_path / "output.gff"

        # Run the function
        parse_gff(str(gff_file), str(output_file), protein_attributes)

        # Verify output
        output_content = output_file.read_text()
        lines = output_content.strip().split("\n")

        # Check header
        assert lines[0] == "##gff-version 3"

        # Check that AMR annotations are added
        assert "drug_class=tetracycline" in output_content
        assert "amr_tool=rgi,amrfinderplus" in output_content
        assert "amr_tool_ident=96.41,99.53" in output_content

        # Check that CDS without AMR annotations are not included
        assert "b20_05_1.circ_2656" not in output_content

    def test_parse_gff_with_gzipped_input(self, tmp_path):
        """Test parsing gzipped GFF files."""
        gff_content = """##gff-version 3
b20_05_1.circ	Prodigal_v2.6.3	CDS	2925164	2927089	152.2	-	0	ID=b20_05_1.circ_2009;partial=00;start_type=ATG"""

        gff_file = tmp_path / "input.gff.gz"
        with gzip.open(str(gff_file), "wt") as f:
            f.write(gff_content)

        protein_attributes = {"b20_05_1.circ_2009": ["drug_class=tetracycline", "amr_tool=rgi", "amr_tool_ident=96.41"]}

        output_file = tmp_path / "output.gff"

        # Should work with gzipped input
        parse_gff(str(gff_file), str(output_file), protein_attributes)

        output_content = output_file.read_text()
        assert "drug_class=tetracycline" in output_content


class TestMainFunction:
    """Test the main function and command-line interface."""

    @patch("mgnify_pipelines_toolkit.analysis.shared.amrintegrator.validate_inputs")
    @patch("mgnify_pipelines_toolkit.analysis.shared.amrintegrator.parse_hamronized")
    @patch("mgnify_pipelines_toolkit.analysis.shared.amrintegrator.parse_amrfinderplus")
    @patch("mgnify_pipelines_toolkit.analysis.shared.amrintegrator.parse_amr_dict")
    @patch("mgnify_pipelines_toolkit.analysis.shared.amrintegrator.parse_gff")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main_with_valid_inputs(
        self, mock_args, mock_parse_gff, mock_parse_amr_dict, mock_parse_amrfinderplus, mock_parse_hamronized, mock_validate
    ):
        """Test main function with valid inputs."""
        # Mock arguments
        mock_args.return_value = argparse.Namespace(
            deeparg_hamr="deeparg.tsv", rgi_hamr="rgi.tsv", amrfp_out="amrfp.tsv", cds_gff="input.gff", output="output.gff", verbose=False
        )

        # Mock validate_inputs to return all inputs as valid
        mock_validate.return_value = ["deeparg", "rgi", "amrfinderplus"]

        # Mock other functions
        mock_parse_hamronized.return_value = {}
        mock_parse_amrfinderplus.return_value = {}
        mock_parse_amr_dict.return_value = {}

        # Run main
        main()

        # Verify functions were called
        assert mock_validate.called
        assert mock_parse_hamronized.call_count == 2  # Called for deeparg and rgi
        assert mock_parse_amrfinderplus.called
        assert mock_parse_amr_dict.called
        assert mock_parse_gff.called

    @patch("mgnify_pipelines_toolkit.analysis.shared.amrintegrator.validate_inputs")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main_no_valid_inputs(self, mock_args, mock_validate, tmp_path):
        """Test main function when no valid inputs are provided."""
        output_file = tmp_path / "output.gff"

        # Mock arguments
        mock_args.return_value = argparse.Namespace(
            deeparg_hamr=None, rgi_hamr=None, amrfp_out=None, cds_gff="input.gff", output=str(output_file), verbose=False
        )

        # Mock validate_inputs to return no valid inputs
        mock_validate.return_value = []

        # Run main
        main()


class TestIntegrationTests:
    """Integration tests using real data similar to your examples."""

    def test_full_integration_workflow(self, tmp_path):
        """Test the complete workflow with realistic data."""
        # Create test input files using structured data
        deeparg_file = tmp_path / "deeparg.tsv"
        create_hamronized_tsv(deeparg_file, DEEPARG_TEST_RECORDS[-2:])  # Last 2 records

        rgi_file = tmp_path / "rgi.tsv"
        create_hamronized_tsv(rgi_file, RGI_TEST_RECORDS[-2:])  # Last 2 records

        amrfp_file = tmp_path / "amrfp.tsv"
        create_amrfp_tsv(amrfp_file, AMRFP_TEST_RECORDS[-2:])  # Last 2 records

        gff_file = tmp_path / "input.gff"
        create_gff_file(gff_file, GFF_TEST_RECORDS[1:3])  # Records for 2148 and 2149

        output_file = tmp_path / "output.gff"

        # Test the full workflow
        optional_inputs = {"deeparg": str(deeparg_file), "rgi": str(rgi_file), "amrfinderplus": str(amrfp_file)}

        valid_inputs = validate_inputs(optional_inputs)
        assert set(valid_inputs) == {"deeparg", "rgi", "amrfinderplus"}

        # Parse all input files
        amr_annotation = defaultdict(dict)
        for name in valid_inputs:
            if name == "amrfinderplus":
                amr_annotation = parse_amrfinderplus(amr_annotation, optional_inputs[name])
            else:
                amr_annotation = parse_hamronized(amr_annotation, optional_inputs[name])

        # Process annotations
        protein_attributes = parse_amr_dict(amr_annotation)

        # Verify we got the expected proteins
        expected_proteins = {"b20_05_1.circ_2868", "b20_05_1.circ_2148", "b20_05_1.circ_2149", "b20_05_1.circ_3091"}
        assert set(protein_attributes.keys()) == expected_proteins

        # Check b20_05_1.circ_2148 (should have all three tools)
        attrs_2148 = protein_attributes["b20_05_1.circ_2148"]
        assert "tetracycline" in attrs_2148[0]
        assert "glycylcycline" in attrs_2148[0]
        assert "deeparg,rgi,amrfinderplus" in attrs_2148[1]
        assert "99.7,99.74,99.74" in attrs_2148[2]

        # Check b20_05_1.circ_2149 (should only have RGI)
        attrs_2149 = protein_attributes["b20_05_1.circ_2149"]
        assert "MLS" in attrs_2149[0]
        assert "streptogramin_A" in attrs_2149[0] or "streptogramin_B" in attrs_2149[0]
        assert attrs_2149[1] == "amr_tool=rgi"  # Only RGI, not all three!
        assert "99.62" in attrs_2149[2]

        # Check b20_05_1.circ_2868 (should only have deeparg)
        attrs_2868 = protein_attributes["b20_05_1.circ_2868"]
        assert "MLS" in attrs_2868[0]
        assert attrs_2868[1] == "amr_tool=deeparg"
        assert "100.0" in attrs_2868[2]

        # Check b20_05_1.circ_3091 (should only have amrfinderplus)
        attrs_3091 = protein_attributes["b20_05_1.circ_3091"]
        assert "beta-lactam" in attrs_3091[0]
        assert attrs_3091[1] == "amr_tool=amrfinderplus"
        assert "59.53" in attrs_3091[2]

        # Generate final GFF
        parse_gff(str(gff_file), str(output_file), protein_attributes)

        # Verify output
        output_content = output_file.read_text()
        assert "##gff-version 3" in output_content
        assert "drug_class=tetracycline,glycylcycline" in output_content
        assert "drug_class=MLS" in output_content
        assert "amr_tool=deeparg,rgi,amrfinderplus" in output_content


if __name__ == "__main__":
    pytest.main([__file__])
