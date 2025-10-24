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
import gzip
import pytest
from unittest.mock import patch


from mgnify_pipelines_toolkit.analysis.shared.amrintegrator import (
    validate_inputs,
    open_file,
    parse_hamronized,
    parse_amrfinderplus,
    parse_amr_dict,
    parse_gff,
    main,
)


class TestValidateInputs:
    """Test the validate_inputs function."""

    def test_validate_inputs_all_exist(self, tmp_path):
        """Test validation when all files exist."""
        # Create temporary files
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file1.write_text("test")
        file2.write_text("test")

        optional_inputs = {"deeparg": str(file1), "rgi": str(file2)}

        result = validate_inputs(optional_inputs)
        assert result == ["deeparg", "rgi"]

    def test_validate_inputs_some_exist(self, tmp_path):
        """Test validation when some files exist."""
        file1 = tmp_path / "file1.txt"
        file1.write_text("test")

        optional_inputs = {"deeparg": str(file1), "rgi": "/nonexistent/path.txt", "amrfinderplus": None}

        result = validate_inputs(optional_inputs)
        assert result == ["deeparg"]

    def test_validate_inputs_none_exist(self):
        """Test validation when no files exist."""
        optional_inputs = {"deeparg": "/nonexistent1.txt", "rgi": "/nonexistent2.txt", "amrfinderplus": None}

        result = validate_inputs(optional_inputs)
        assert result == []


class TestOpenFile:
    """Test the open_file function for handling compressed and uncompressed files."""

    def test_open_regular_file(self, tmp_path):
        """Test opening regular text file."""
        test_file = tmp_path / "test.txt"
        test_content = "line1\nline2\n"
        test_file.write_text(test_content)

        with open_file(str(test_file)) as f:
            content = f.read()

        assert content == test_content

    def test_open_gzipped_file(self, tmp_path):
        """Test opening gzipped file."""
        test_file = tmp_path / "test.txt.gz"
        test_content = "line1\nline2\n"

        with gzip.open(str(test_file), "wt") as f:
            f.write(test_content)

        with open_file(str(test_file)) as f:
            content = f.read()

        assert content == test_content


class TestParseHamronized:
    """Test the parse_hamronized function for parsing DeepARG and RGI outputs."""

    def test_parse_deeparg_hamronized(self, tmp_path):
        """Test parsing DeepARG hamronized output."""
        # Create test data based on your example
        deeparg_content = """input_file_name	gene_symbol	gene_name	reference_database_name	reference_database_version	reference_accession	analysis_software_name	analysis_software_version	genetic_variation_type	antimicrobial_agent	coverage_percentage	coverage_depth	coverage_ratio	drug_class	input_gene_length	input_gene_start	input_gene_stop	input_protein_length	input_protein_start	input_protein_stop	input_sequence_id	nucleotide_mutation	nucleotide_mutation_interpretation	predicted_phenotype	predicted_phenotype_confidence_level	amino_acid_mutation	amino_acid_mutation_interpretation	reference_gene_length	reference_gene_start	reference_gene_stop	reference_protein_length	reference_protein_start	reference_protein_stop	resistance_mechanism	strand_orientation	sequence_identity
b20_05	ROSA	YP_001401993|FEATURES|rosA|fosmidomycin|rosA	deeparg_db	2	YP_001401993	deeparg	1.0.4	gene_presence_detected					fosmidomycin		10	401				b20_05_1.circ_3209															54.3
b20_05	ERMF	AAA27431|FEATURES|ermF|MLS|ermF	deeparg_db	2	AAA27431	deeparg	1.0.4	gene_presence_detected					MLS		1	266				b20_05_1.circ_2868															100.0
"""

        test_file = tmp_path / "deeparg.tsv"
        test_file.write_text(deeparg_content)

        amr_annotation = {}
        result = parse_hamronized(amr_annotation, str(test_file))

        # Verify results
        assert "b20_05_1.circ_3209" in result
        assert "b20_05_1.circ_2868" in result

        assert result["b20_05_1.circ_3209"]["deeparg"]["drug_class"] == ["fosmidomycin"]
        assert result["b20_05_1.circ_3209"]["deeparg"]["seq_identity"] == "54.3"

        assert result["b20_05_1.circ_2868"]["deeparg"]["drug_class"] == ["MLS"]
        assert result["b20_05_1.circ_2868"]["deeparg"]["seq_identity"] == "100.0"

    def test_parse_rgi_hamronized(self, tmp_path):
        """Test parsing RGI hamronized output."""
        rgi_content = """input_file_name	gene_symbol	gene_name	reference_database_name	reference_database_version	reference_accession	analysis_software_name	analysis_software_version	genetic_variation_type	antimicrobial_agent	coverage_percentage	coverage_depth	coverage_ratio	drug_class	input_gene_length	input_gene_start	input_gene_stop	input_protein_length	input_protein_start	input_protein_stop	input_sequence_id	nucleotide_mutation	nucleotide_mutation_interpretation	predicted_phenotype	predicted_phenotype_confidence_level	amino_acid_mutation	amino_acid_mutation_interpretation	reference_gene_length	reference_gene_start	reference_gene_stop	reference_protein_length	reference_protein_start	reference_protein_stop	resistance_mechanism	strand_orientation	sequence_identity
b20_05	tet(Q)	tetracycline-resistant ribosomal protection protein	CARD	4.0.1	3000191	rgi	6.0.5	gene_presence_detected	tetracycline; doxycycline; minocycline; chlortetracycline; demeclocycline; oxytetracycline	97.56			tetracycline antibiotic							b20_05_1.circ_2009 # 2925164 # 2927089 # -1 # ID=1_2009;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=5bp;gc_cont=0.398										antibiotic target protection					96.41
b20_05	ErmF	Erm 23S ribosomal RNA methyltransferase	CARD	4.0.1	3000498	rgi	6.0.5	gene_presence_detected	erythromycin; roxithromycin; lincomycin; telithromycin; clarithromycin; clindamycin; tylosin; spiramycin; azithromycin; dirithromycin; pristinamycin IA; quinupristin; virginiamycin M1; madumycin II; griseoviridin; dalfopristin; pristinamycin IB; virginiamycin S2; pristinamycin IC; vernamycin C; patricin A; patricin B; ostreogrycin B3; oleandomycin	100.0			macrolide antibiotic; lincosamide antibiotic; streptogramin antibiotic; streptogramin A antibiotic; streptogramin B antibiotic							b20_05_1.circ_2149 # 3109927 # 3110727 # -1 # ID=1_2149;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=10bp;gc_cont=0.341									antibiotic target alteration						99.62
"""

        test_file = tmp_path / "rgi.tsv"
        test_file.write_text(rgi_content)

        amr_annotation = {}
        result = parse_hamronized(amr_annotation, str(test_file))

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
        amrfp_content = """Protein id	Element symbol	Element name	Scope	Type	Subtype	Class	Subclass	Method	Target length	Reference sequence length	% Coverage of reference	% Identity to reference	Alignment length	Closest reference accession	Closest reference name	HMM accession	HMM description
b20_05_1.circ_2009	tet(Q)	tetracycline resistance ribosomal protection protein Tet(Q)	core	AMR	AMR	TETRACYCLINE	TETRACYCLINE	BLASTP	641	641	100.00	99.53	641	WP_063856407.1	tetracycline resistance ribosomal protection protein Tet(Q)	NF012154.0	tetracycline resistance ribosomal protection protein Tet(Q)
b20_05_1.circ_2149	erm(F)	23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)	core	AMR	AMR	LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN	CLINDAMYCIN/MACROLIDE/STREPTOGRAMIN	EXACTP	266	266	100.00	100.00	266	WP_002682030.1	23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)	NF012223.0	23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)
b20_05_1.circ_3091	bla	class A beta-lactamase, subclass A2	core	AMR	AMR	BETA-LACTAM	BETA-LACTAM	HMM	301	296	99.66	59.53	299	WP_005827792.1	class A beta-lactamase, subclass A2	NF012099.1	class A beta-lactamase, subclass A2"""

        test_file = tmp_path / "amrfp.tsv"
        test_file.write_text(amrfp_content)

        amr_annotation = {}
        result = parse_amrfinderplus(amr_annotation, str(test_file))

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
        # Create test GFF content
        gff_content = """##gff-version 3
b20_05_1.circ	Prodigal_v2.6.3	CDS	2925164	2927089	152.2	-	0	ID=b20_05_1.circ_2009;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=5bp;gc_cont=0.398;conf=100.00;score=152.15;cscore=148.05;sscore=4.10;rscore=-0.03;uscore=0.76;tscore=4.02
b20_05_1.circ	Prodigal_v2.6.3	CDS	3108548	3109714	78.3	+	0	ID=b20_05_1.circ_2148;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.374;conf=100.00;score=78.33;cscore=89.46;sscore=-11.13;rscore=-9.55;uscore=-5.60;tscore=4.02
b20_05_1.circ	Prodigal_v2.6.3	CDS	3798391	3799731	191.3	+	0	ID=b20_05_1.circ_2656;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=11bp;gc_cont=0.398;conf=99.99;score=191.32;cscore=179.05;sscore=12.27;rscore=8.19;uscore=-0.65;tscore=4.02
"""

        gff_file = tmp_path / "input.gff"
        gff_file.write_text(gff_content)

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
            deeparg_hamr="deeparg.tsv", rgi_hamr="rgi.tsv", amrfp_out="amrfp.tsv", cds_gff="input.gff", output="output.gff"
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
        mock_args.return_value = argparse.Namespace(deeparg_hamr=None, rgi_hamr=None, amrfp_out=None, cds_gff="input.gff", output=str(output_file))

        # Mock validate_inputs to return no valid inputs
        mock_validate.return_value = []

        # Run main
        main()

        # Check that empty output file is created
        assert output_file.exists()
        content = output_file.read_text()
        assert content == "##gff-version 3\n"


class TestIntegrationTests:
    """Integration tests using real data similar to your examples."""

    def test_full_integration_workflow(self, tmp_path):
        """Test the complete workflow with realistic data."""
        # Create test input files
        deeparg_file = tmp_path / "deeparg.tsv"
        deeparg_content = """input_file_name	gene_symbol	gene_name	reference_database_name	reference_database_version	reference_accession	analysis_software_name	analysis_software_version	genetic_variation_type	antimicrobial_agent	coverage_percentage	coverage_depth	coverage_ratio	drug_class	input_gene_length	input_gene_start	input_gene_stop	input_protein_length	input_protein_start	input_protein_stop	input_sequence_id	nucleotide_mutation	nucleotide_mutation_interpretation	predicted_phenotype	predicted_phenotype_confidence_level	amino_acid_mutation	amino_acid_mutation_interpretation	reference_gene_length	reference_gene_start	reference_gene_stop	reference_protein_length	reference_protein_start	reference_protein_stop	resistance_mechanism	strand_orientation	sequence_identity
b20_05	ERMF	AAA27431|FEATURES|ermF|MLS|ermF	deeparg_db	2	AAA27431	deeparg	1.0.4	gene_presence_detected					MLS		1	266				b20_05_1.circ_2149															100.0
b20_05	TETX	gi:1004705391:gb:AMP52595.1:|FEATURES|tetX|tetracycline|tetX	deeparg_db	2	gi:1004705391:gb:AMP52595.1:	deeparg	1.0.4	gene_presence_detected					tetracycline		1	388				b20_05_1.circ_2148															99.7
"""
        deeparg_file.write_text(deeparg_content)

        rgi_file = tmp_path / "rgi.tsv"
        rgi_content = """input_file_name	gene_symbol	gene_name	reference_database_name	reference_database_version	reference_accession	analysis_software_name	analysis_software_version	genetic_variation_type	antimicrobial_agent	coverage_percentage	coverage_depth	coverage_ratio	drug_class	input_gene_length	input_gene_start	input_gene_stop	input_protein_length	input_protein_start	input_protein_stop	input_sequence_id	nucleotide_mutation	nucleotide_mutation_interpretation	predicted_phenotype	predicted_phenotype_confidence_level	amino_acid_mutation	amino_acid_mutation_interpretation	reference_gene_length	reference_gene_start	reference_gene_stop	reference_protein_length	reference_protein_start	reference_protein_stop	resistance_mechanism	strand_orientation	sequence_identity
b20_05	tet(X)	tetracycline inactivation enzyme	CARD	4.0.1	3000205	rgi	6.0.5	gene_presence_detected	tigecycline; tetracycline; doxycycline; minocycline; chlortetracycline; demeclocycline; oxytetracycline	100.0			glycylcycline; tetracycline antibiotic							b20_05_1.circ_2148 # 3108548 # 3109714 # 1 # ID=1_2148;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.374									antibiotic inactivation						99.74
b20_05	ErmF	Erm 23S ribosomal RNA methyltransferase	CARD	4.0.1	3000498	rgi	6.0.5	gene_presence_detected	erythromycin; roxithromycin; lincomycin; telithromycin; clarithromycin; clindamycin; tylosin; spiramycin; azithromycin; dirithromycin; pristinamycin IA; quinupristin; virginiamycin M1; madumycin II; griseoviridin; dalfopristin; pristinamycin IB; virginiamycin S2; pristinamycin IC; vernamycin C; patricin A; patricin B; ostreogrycin B3; oleandomycin	100.0			macrolide antibiotic; lincosamide antibiotic; streptogramin antibiotic; streptogramin A antibiotic; streptogramin B antibiotic							b20_05_1.circ_2149 # 3109927 # 3110727 # -1 # ID=1_2149;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=10bp;gc_cont=0.341									antibiotic target alteration						99.62
"""
        rgi_file.write_text(rgi_content)

        amrfp_file = tmp_path / "amrfp.tsv"
        amrfp_content = """Protein id	Element symbol	Element name	Scope	Type	Subtype	Class	Subclass	Method	Target length	Reference sequence length	% Coverage of reference	% Identity to reference	Alignment length	Closest reference accession	Closest reference name	HMM accession	HMM description
b20_05_1.circ_2148	tet(X2)	tetracycline-inactivating monooxygenase Tet(X2)	core	AMR	AMR	TETRACYCLINE	TIGECYCLINE	BLASTP	388	388	100.00	99.74	388	WP_008651082.1	tetracycline-inactivating monooxygenase Tet(X2)	NF033111.2	tetracycline-inactivating monooxygenase Tet(X)
b20_05_1.circ_2149	erm(F)	23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)	core	AMR	AMR	LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN	CLINDAMYCIN/MACROLIDE/STREPTOGRAMIN	EXACTP	266	266	100.00	100.00	266	WP_002682030.1	23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)	NF012223.0	23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)
"""
        amrfp_file.write_text(amrfp_content)

        gff_file = tmp_path / "input.gff"
        gff_content = """##gff-version 3
b20_05_1.circ	Prodigal_v2.6.3	CDS	3108548	3109714	78.3	+	0	ID=b20_05_1.circ_2148;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.374;conf=100.00;score=78.33;cscore=89.46;sscore=-11.13;rscore=-9.55;uscore=-5.60;tscore=4.02
b20_05_1.circ	Prodigal_v2.6.3	CDS	3109927	3110727	38.1	-	0	ID=b20_05_1.circ_2149;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=10bp;gc_cont=0.341;conf=99.98;score=38.14;cscore=31.14;sscore=7.00;rscore=-0.03;uscore=3.01;tscore=4.02
"""
        gff_file.write_text(gff_content)

        output_file = tmp_path / "output.gff"

        # Test the full workflow
        optional_inputs = {"deeparg": str(deeparg_file), "rgi": str(rgi_file), "amrfinderplus": str(amrfp_file)}

        valid_inputs = validate_inputs(optional_inputs)
        assert set(valid_inputs) == {"deeparg", "rgi", "amrfinderplus"}

        # Parse all input files
        amr_annotation = {}
        for name in valid_inputs:
            if name == "amrfinderplus":
                amr_annotation = parse_amrfinderplus(amr_annotation, optional_inputs[name])
            else:
                amr_annotation = parse_hamronized(amr_annotation, optional_inputs[name])

        # Process annotations
        protein_attributes = parse_amr_dict(amr_annotation)

        # Verify processed annotations
        assert "b20_05_1.circ_2148" in protein_attributes
        assert "b20_05_1.circ_2149" in protein_attributes

        # Check b20_05_1.circ_2148 (should have all three tools)
        attrs_2148 = protein_attributes["b20_05_1.circ_2148"]
        assert "tetracycline" in attrs_2148[0]
        assert "glycylcycline" in attrs_2148[0]
        assert "deeparg,rgi,amrfinderplus" in attrs_2148[1]
        assert "99.7,99.74,99.74" in attrs_2148[2]

        # Check b20_05_1.circ_2149 (should have all three tools)
        attrs_2149 = protein_attributes["b20_05_1.circ_2149"]
        assert "MLS" in attrs_2149[0]
        assert "streptogramin_A" in attrs_2149[0] or "streptogramin_B" in attrs_2149[0]
        assert "deeparg,rgi,amrfinderplus" in attrs_2149[1]
        assert "100.0,99.62,100.00" in attrs_2149[2]

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
