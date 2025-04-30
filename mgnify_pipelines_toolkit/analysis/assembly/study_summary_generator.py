#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
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

import click
from functools import reduce
import glob
import logging
from pathlib import Path
from typing import Literal

import pandas as pd

from mgnify_pipelines_toolkit.schemas.schemas import (
    CompletedAnalysisSchema,
    TaxonSchema,
    GOSummarySchema,
    InterProSummarySchema,
    KOSummarySchema,
    SanntisSummarySchema,
    AntismashSummarySchema,
    PFAMSummarySchema,
    KEGGModulesSummarySchema,
    GOStudySummarySchema,
    InterProStudySummarySchema,
    TaxonomyStudySummarySchema,
    KOStudySummarySchema,
    SanntisStudySummarySchema,
    AntismashStudySummarySchema,
    PFAMStudySummarySchema,
    KEGGModulesStudySummarySchema,
    validate_dataframe,
)

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)

# Keys are the original column names in the input files,
# values are the standardised column names used in the generated study summary files
# Note: "Count" or "count" column should be excluded
GO_COLUMN_NAMES = {
    "go": "GO",
    "term": "description",
    "category": "category",
}

INTERPRO_COLUMN_NAMES = {
    "interpro_accession": "IPR",
    "description": "description",
}

SANNTIS_COLUMN_NAMES = {
    "nearest_mibig": "nearest_mibig",
    "nearest_mibig_class": "nearest_mibig_class",
    "description": "description",
}

ANTISMASH_COLUMN_NAMES = {
    "label": "label",
    "description": "description",
}

KEGG_COLUMN_NAMES = {
    "ko": "KO",
    "description": "description",
}

PFAM_COLUMN_NAMES = {
    "pfam": "PFAM",
    "description": "description",
}

KEGG_MODULES_COLUMN_NAMES = {
    "module_accession": "module_accession",
    "pathway_name": "pathway_name",
    "pathway_class": "pathway_class",
}

# this mapping allows using 'for' cycle later to process all summary types in one way
SUMMARY_TYPES_MAP = {
    "go": {
        "folder": "functional-annotation/go",
        "column_names": GO_COLUMN_NAMES,
        "schema": GOSummarySchema,
        "study_schema": GOStudySummarySchema,
    },
    "goslim": {
        "folder": "functional-annotation/go",
        "column_names": GO_COLUMN_NAMES,
        "schema": GOSummarySchema,
        "study_schema": GOStudySummarySchema,
    },
    "interpro": {
        "folder": "functional-annotation/interpro",
        "column_names": INTERPRO_COLUMN_NAMES,
        "schema": InterProSummarySchema,
        "study_schema": InterProStudySummarySchema,
    },
    "ko": {
        "folder": "functional-annotation/kegg",
        "column_names": KEGG_COLUMN_NAMES,
        "schema": KOSummarySchema,
        "study_schema": KOStudySummarySchema,
    },
    "sanntis": {
        "folder": "pathways-and-systems/sanntis",
        "column_names": SANNTIS_COLUMN_NAMES,
        "schema": SanntisSummarySchema,
        "study_schema": SanntisStudySummarySchema,
    },
    "antismash": {
        "folder": "pathways-and-systems/antismash",
        "column_names": ANTISMASH_COLUMN_NAMES,
        "schema": AntismashSummarySchema,
        "study_schema": AntismashStudySummarySchema,
    },
    "pfam": {
        "folder": "functional-annotation/pfam",
        "column_names": PFAM_COLUMN_NAMES,
        "schema": PFAMSummarySchema,
        "study_schema": PFAMStudySummarySchema,
    },
    "kegg_modules": {
        "folder": "pathways-and-systems/kegg-modules",
        "column_names": KEGG_MODULES_COLUMN_NAMES,
        "schema": KEGGModulesSummarySchema,
        "study_schema": KEGGModulesStudySummarySchema,
    },
}

# The taxonomy file is a tab-separated file without any header
# containing of following columns:
TAXONOMY_COLUMN_NAMES = [
    "Count",
    "Superkingdom",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
]

OUTPUT_SUFFIX = "summary.tsv"


@click.group()
def cli():
    pass


def check_files_exist(file_list: list[Path]) -> None:
    """
    Check that all files in the given list exist on disk.

    :param file_list: List of file paths to check.
    :raises FileNotFoundError: If any file does not exist.
    """
    missing_files = [str(path) for path in file_list if not path.is_file()]
    if missing_files:
        raise FileNotFoundError(
            f"The following required files are missing: {', '.join(missing_files)}"
        )


def generate_taxonomy_summary(
    file_dict: dict[str, Path],
    output_file_name: str,
    outdir: Path = None,
) -> None:
    """
    Generate a combined study-level taxonomic classification summary from multiple input
    assembly-level summary files.

    :param file_dict: Dictionary mapping assembly accession to its taxonomy file.
    :param output_file_name: Output path for the output summary file.
    :param outdir: Optional output directory for the results.

    Example of the taxonomy file:
    23651	sk__Bacteria
    4985	sk__Archaea	k__Thermoproteati	p__Nitrososphaerota
    882	sk__Archaea	k__Nanobdellati	p__	c__	o__	f__	g__	s__Candidatus Pacearchaeota archaeon
    """
    check_files_exist(list(file_dict.values()))

    tax_dfs = []
    for assembly_acc, path in file_dict.items():
        df = pd.read_csv(path, sep="\t", names=TAXONOMY_COLUMN_NAMES).fillna("")

        # Note: schema validation will fail if the taxonomy file is empty
        df = validate_dataframe(df, TaxonSchema, str(path))

        # Combine all taxonomic ranks in the classification into a single string
        df["full_taxon"] = (
            df[TAXONOMY_COLUMN_NAMES[1:]].agg(";".join, axis=1).str.strip(";")
        )

        # Create a new DataFrame with taxonomy as index and count as the only column
        result = df[["Count", "full_taxon"]].set_index("full_taxon")
        result.columns = [assembly_acc]
        tax_dfs.append(result)

    summary_df = pd.concat(tax_dfs, axis=1)
    summary_df = summary_df.fillna(0).astype(int).sort_index()

    outfile = output_file_name
    if outdir:
        outfile = outdir / output_file_name

    summary_df.to_csv(outfile, sep="\t", index_label="taxonomy")


def generate_functional_summary(
    file_dict: dict[str, Path],
    column_names: dict[str, str],
    output_prefix: str,
    label: Literal[
        "go", "goslim", "interpro", "ko", "sanntis", "antismash", "pfam", "kegg_modules"
    ],
    outdir: Path = None,
) -> None:
    """
    Generate a combined study-level functional annotation summary from multiple input
    assembly-level summary files.

    :param file_dict: Dictionary mapping assembly accession to its summary file path.
    :param column_names: Dictionary mapping original column names to standard column names.
    :param output_prefix: Prefix for the output summary file.
    :param label: Label for the functional annotation type
    (expected one of ["go", "goslim", "interpro", "ko", "sanntis", "antismash", "pfam", "kegg_modules"]).
    :param outdir: Optional output directory for the results.

    In the input files, column orders may vary, but the following columns are expected:
    GO summary input file:
    go	term	category	count
    GO:0016020	membrane	cellular_component	30626
    GO:0005524	ATP binding	molecular_function	30524

    InterPro summary input file:
    interpro_accession	description	count
    IPR036291	NAD(P)-binding domain superfamily	16503
    IPR019734	Tetratricopeptide repeat	14694

    KEGG summary input file:
    ko      description	count
    K01552  energy-coupling factor transport system ATP-binding protein [EC:7.-.-.-]	562
    K18889  ATP-binding cassette, subfamily B, multidrug efflux pump	537
    K15497  molybdate/tungstate transport system ATP-binding protein [EC:7.3.2.5 7.3.2.6]	517

    Sanntis summary input file:
    nearest_mibig	nearest_mibig_class	description	count
    BGC0000787	Saccharide	Carbohydrate-based natural products (e.g., aminoglycoside antibiotics)	1
    BGC0000248	Polyketide	Built from iterative condensation of acetate units derived from acetyl-CoA	3
    BGC0001327	NRP Polyketide	Nonribosomal Peptide Polyketide	2

    Antismash summary input file:
    label	description	count
    terpene	Terpene	16
    betalactone	Beta-lactone containing protease inhibitor	8
    T1PKS	Type I PKS (Polyketide synthase)	3

    PFAM summary input file:
    pfam	description	count
    PF00265	Thymidine kinase	457
    PF01852	START domain	368
    PF13756	Stimulus-sensing domain	397

    KEGG modules summary input file:
    module_accession	completeness	pathway_name	pathway_class	matching_ko	missing_ko
    M00986	100.0	Sulfur reduction, sulfur => sulfide	Pathway modules; Energy metabolism; Sulfur metabolism	K18367
    M00163	83.33	Photosystem I	Pathway modules; Energy metabolism; Photosynthesis	K02689,K02690,K02691,K02692,K02694	K02693
    M00615	50.0	Nitrate assimilation	Signature modules; Module set; Metabolic capacity	K02575	M00531
    """
    check_files_exist(list(file_dict.values()))

    output_file_name = f"{output_prefix}_{label}_{OUTPUT_SUFFIX}"

    original_col_names = list(column_names.keys())
    renamed_col_names = list(column_names.values())
    value_col_name = "completeness" if label == "kegg_modules" else "count"

    dfs = []
    for assembly_acc, filepath in file_dict.items():
        try:
            df = pd.read_csv(filepath, sep="\t")
        except pd.errors.EmptyDataError:
            logging.warning(f"File {filepath.resolve()} is empty. Skipping.")
            continue

        schema = SUMMARY_TYPES_MAP[label]["schema"]
        df = validate_dataframe(df, schema, str(filepath))

        # Extract only relevant columns
        df = df[original_col_names + [value_col_name]].copy()

        # Rename columns: metadata columns are renamed according to column_names dict, "count"/"completeness" -> assembly acc
        df.rename(columns={**column_names, value_col_name: assembly_acc}, inplace=True)
        dfs.append(df)

    if not dfs:
        logging.warning(
            f"No valid files with functional annotation summary were found. Skipping creation of {output_file_name}."
        )
        return

    # Merge all dataframes on the renamed metadata columns
    merged_df = reduce(
        lambda left, right: pd.merge(left, right, on=renamed_col_names, how="outer"),
        dfs,
    )

    # Fill missing values appropriately, convert completeness percentages to float, counts to integers
    value_columns = [col for col in merged_df.columns if col not in renamed_col_names]
    fill_value = 0.0 if label == "kegg_modules" else 0
    dtype = float if label == "kegg_modules" else int
    merged_df[value_columns] = merged_df[value_columns].fillna(fill_value).astype(dtype)

    # Reorder columns: merge keys first, then sorted assembly accessions
    merged_df = merged_df[renamed_col_names + sorted(value_columns)]

    outfile = output_file_name
    if outdir:
        outfile = outdir / output_file_name

    merged_df.to_csv(outfile, sep="\t", index=False)


@cli.command(
    "summarise",
    options_metavar="-a <assemblies> -s <study_dir> -p <output_prefix>",
    short_help="Generate study-level summaries for assembly analysis results.",
)
@click.option(
    "-a",
    "--assemblies",
    required=True,
    help="CSV file containing successful analyses generated by the pipeline",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-s",
    "--study_dir",
    required=True,
    help="Input directory to where all the individual analyses subdirectories for summarising",
    type=click.Path(exists=True, path_type=Path, file_okay=False),
)
@click.option(
    "-p",
    "--output_prefix",
    required=True,
    help="Prefix for generated summary files",
    type=str,
)
@click.option(
    "-o",
    "--outdir",
    required=False,
    help="Directory for the output files, by default it will use the current working directory.",
    type=click.Path(exists=True, path_type=Path, file_okay=False),
)
def summarise_analyses(
    assemblies: Path, study_dir: Path, output_prefix: str, outdir: Path
) -> None:
    """
    Generate study-level summaries for successfully proccessed assemblies.

    :param assemblies: Path to a file listing completed assembly accessions and their status.
    :param study_dir: Path to the directory containing analysis results for each assembly.
    :param output_prefix: Prefix for the generated summary files.
    """
    logging.info(f"Reading assembly list from {assemblies.resolve()}")
    assemblies_df = pd.read_csv(assemblies, names=["assembly", "status"])
    CompletedAnalysisSchema(assemblies_df)
    assembly_list = assemblies_df["assembly"].tolist()
    logging.info("Assembly list was read successfully.")

    def get_file_paths(subdir: str, filename_template: str) -> dict[str, Path]:
        """
        Construct file paths for each assembly given a subdirectory and filename template.
        Template must contain {acc} as a placeholder.
        """
        return {
            acc: study_dir / acc / subdir / filename_template.format(acc=acc)
            for acc in assembly_list
        }

    logging.info("Start processing of assembly-level summaries.")

    logging.info(
        "Generating taxonomy summary from assembly-level summaries <accession>.krona.txt"
    )
    generate_taxonomy_summary(
        get_file_paths("taxonomy", "{acc}.krona.txt.gz"),
        f"{output_prefix}_taxonomy_{OUTPUT_SUFFIX}",
        outdir=outdir,
    )

    for summary_type, config in SUMMARY_TYPES_MAP.items():
        logging.info(
            f"Generating study-level {summary_type.capitalize()} summary from file <accession>_{summary_type}_summary.tsv.gz"
        )
        generate_functional_summary(
            get_file_paths(config["folder"], f"{{acc}}_{summary_type}_summary.tsv.gz"),
            config["column_names"],
            output_prefix,
            summary_type,
            outdir=outdir,
        )
    logging.info("Assembly-level summaries were generated successfully.")
    logging.info("Done.")


@cli.command(
    "merge",
    options_metavar="-a <study_dir> -p <output_prefix>",
    short_help="Merge multiple study-level summaries of assembly analysis.",
)
@click.option(
    "-s",
    "--study_dir",
    required=True,
    help="Input directory to where all the individual analyses subdirectories for merging",
    type=click.Path(exists=True, file_okay=False),
)
@click.option(
    "-p",
    "--output_prefix",
    required=True,
    help="Prefix for generated merged summary files",
    type=str,
)
def merge_summaries(study_dir: str, output_prefix: str) -> None:
    """
    Merge multiple study-level summary files into combined summary files.

    :param study_dir: Path to the directory containing study-level summary files.
    :param output_prefix: Prefix for the output merged summary files.
    """

    def get_file_paths(summary_type: str) -> list[str]:
        return glob.glob(f"{study_dir}/*_{summary_type}_{OUTPUT_SUFFIX}")

    logging.info("Generating combined assembly-level summaries")
    logging.info("Parsing summary files for taxonomic classification")
    merge_taxonomy_summaries(
        get_file_paths("taxonomy"), f"{output_prefix}_taxonomy_{OUTPUT_SUFFIX}"
    )

    for summary_type, config in SUMMARY_TYPES_MAP.items():
        logging.info(f"Parsing summary files for {summary_type.capitalize()}.")
        column_names = config["column_names"]
        merge_functional_summaries(
            get_file_paths(summary_type),
            list(column_names.values()),
            output_prefix,
            summary_type,
        )
    logging.info("Merged assembly-level summaries were generated successfully.")
    logging.info("Done.")


def merge_taxonomy_summaries(summary_files: list[str], output_file_name: str) -> None:
    """
    Merge multiple taxonomy study-level summary files into a single study-level summary.

    :param summary_files: List of paths to taxonomy summary files, each containing
                        taxonomic classifications and counts for an individual analysis.
    :param output_file_name: Output path for the merged taxonomy summary.

    Example of input taxonomy summary file:
    taxonomy	ERZ1049444	ERZ1049446
    sk__Eukaryota;k__Metazoa;p__Chordata	2	10
    sk__Eukaryota;k__Metazoa;p__Chordata;c__Mammalia;o__Primates	118	94
    """
    if not summary_files:
        raise FileNotFoundError(
            "The required taxonomic classification summary files are missing. Exiting."
        )

    summary_dfs = []
    for file in summary_files:
        df = pd.read_csv(file, sep="\t", index_col=0)
        df = validate_dataframe(df, TaxonomyStudySummarySchema, file)
        summary_dfs.append(df)
    merged_df = pd.concat(summary_dfs, axis=1)
    merged_df = merged_df.fillna(0).astype(int)

    # Reorder columns: taxonomy first, then sorted assembly accessions
    merged_df = merged_df[sorted(merged_df.columns)]
    merged_df = merged_df.sort_index()

    merged_df.to_csv(
        output_file_name,
        sep="\t",
        index_label="taxonomy",
    )


def merge_functional_summaries(
    summary_files: list[str],
    merge_keys: list[str],
    output_prefix: str,
    label: Literal[
        "go", "goslim", "interpro", "ko", "sanntis", "antismash", "pfam", "kegg_modules"
    ],
) -> None:
    """
    Merge multiple functional study-level summary files into a single study-level summary.

    :param summary_files: List of paths to functional summary files, each containing
                        annotation terms and counts for an individual analysis.
    :param merge_keys: List of column names to merge on (e.g. term ID, description).
    :param output_prefix: Prefix for the generated output file.
    :param label: Label describing the functional annotation type
    (expected one of ["go", "goslim", "interpro", "ko", "sanntis", "antismash", "pfam", "kegg_modules"]).

    In the input files, column orders may vary, but the following columns are expected:
    GO summary input:
    GO	description	category	ERZ1049444	ERZ1049446
    GO:0016020	membrane	cellular_component	30626	673
    GO:0005524	ATP binding	molecular_function	30524	2873

    Example of InterPro summary input:
    IPR	description	ERZ1049444	ERZ1049446
    IPR036291	NAD(P)-binding domain superfamily	16503	13450
    IPR019734	Tetratricopeptide repeat	14694	11021

    KEGG summary input:
    GO	description	category	ERZ1049440	ERZ1049443
    GO:0003677	DNA binding	molecular_function	6125	16417
    GO:0055085	transmembrane transport	biological_process	144	13926

    Sanntis summary input:
    nearest_mibig	nearest_mibig_class	description	ERZ1049440	ERZ1049443
    BGC0001356	RiPP	Ribosomally synthesised and Post-translationally modified Peptide	230	185
    BGC0001432	NRP Polyketide	Nonribosomal Peptide Polyketide	0	8

    Antismash summary input:
    label	description	ERZ1049440	ERZ1049443
    NRPS	Non-ribosomal peptide synthetase	368	0
    arylpolyene	Aryl polyene	149	447

    PFAM summary input:
    PFAM	description	ERZ1049440	ERZ1049443
    PF24718	HTH-like domain	468	1
    PF06039	Malate:quinone oxidoreductase (Mqo)	490	21

    KEGG modules summary input:
    module_accession	pathway_name	pathway_class	ERZ1049440	ERZ1049443
    M00109	C21-Steroid hormone biosynthesis, progesterone => cortisol/cortisone	Pathway modules; Lipid metabolism; Sterol biosynthesis	38.9	0.0
    M00153	Cytochrome bd ubiquinol oxidase	Pathway modules; Energy metabolism; ATP synthesis	44.7	84.4
    """
    output_file_name = f"{output_prefix}_{label}_{OUTPUT_SUFFIX}"

    if not summary_files:
        logging.warning(
            f"Skipping creation of {output_file_name} because no summaries were found for this type of functional annotation."
        )
        return

    validation_schema = SUMMARY_TYPES_MAP[label]["study_schema"]

    dfs = []
    for filepath in summary_files:
        df = pd.read_csv(filepath, sep="\t")
        df = validate_dataframe(df, validation_schema, filepath)
        dfs.append(df)

    if len(dfs) == 1:
        merged_df = dfs[0]
    else:
        merged_df = reduce(
            lambda left, right: pd.merge(left, right, on=merge_keys, how="outer"), dfs
        )

    # Identify non-key columns (i.e. counts)
    value_columns = [col for col in merged_df.columns if col not in merge_keys]

    # Fill NaNs and set dtype accordingly
    fill_value = 0.0 if label == "kegg_modules" else 0
    dtype = float if label == "kegg_modules" else int
    merged_df[value_columns] = merged_df[value_columns].fillna(fill_value).astype(dtype)

    # Reorder columns
    merged_df = merged_df[merge_keys + sorted(value_columns)]

    merged_df.to_csv(output_file_name, sep="\t", index=False)


if __name__ == "__main__":
    cli()
