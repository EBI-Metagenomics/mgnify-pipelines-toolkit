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
import glob
import logging
from pathlib import Path

import pandas as pd

from mgnify_pipelines_toolkit.schemas.schemas import (
    CompletedAnalysisSchema,
    TaxonSchema,
    GOSummarySchema,
    InterProSummarySchema,
    validate_dataframe,
)

# TODO add logging to the script
logging.basicConfig(level=logging.DEBUG)


GO_COLUMN_NAMES = {
    "go": "GO",
    "term": "description",
    "category": "category",
}

INTERPRO_COLUMN_NAMES = {
    "interpro_accession": "IPR",
    "description": "description",
}


@click.group()
def cli():
    pass


def check_files_exist(file_list: list[Path]) -> None:
    missing_files = [str(path) for path in file_list if not path.exists()]
    if missing_files:
        raise FileNotFoundError(
            f"The following required files are missing: {', '.join(missing_files)}"
        )


def generate_taxonomy_summary(
    file_dict: dict[str, Path], output_file_name: str
) -> None:
    """
    Generates a combined taxonomy summary from multiple taxonomy summary_files.

    :param file_dict: Dictionary mapping assembly accession to its taxonomy file.
    :param output_prefix: Prefix for the output summary file.

    Example of the taxonomy file:
    23651	sk__Bacteria
    4985	sk__Archaea	k__Thermoproteati	p__Nitrososphaerota
    882	sk__Archaea	k__Nanobdellati	p__	c__	o__	f__	g__	s__Candidatus Pacearchaeota archaeon
    """
    check_files_exist(list(file_dict.values()))

    tax_columns = [
        "Count",
        "superkingdom",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    tax_dfs = []
    for assembly_acc, path in file_dict.items():
        df = pd.read_csv(path, sep="\t", names=tax_columns).fillna("")

        # TODO Now if the dataframe is empty, SchemaErrors will be raised, but message is not clear
        # Maybe we should fix it?
        validate_dataframe(df, TaxonSchema, str(path))

        # Combine all taxonomic levels into a single string per row
        df["full_taxon"] = df[tax_columns[1:]].agg(";".join, axis=1).str.strip(";")

        # Create a DataFrame with taxonomy as index and count as the only column
        result = df[["Count", "full_taxon"]].set_index("full_taxon")
        result.columns = [assembly_acc]
        tax_dfs.append(result)

    summary_df = pd.concat(tax_dfs, axis=1).fillna(0).astype(int).sort_index()

    summary_df.to_csv(output_file_name, sep="\t", index_label="taxonomy")


def generate_functional_summary(
    file_dict: dict[str, Path],
    column_names: dict[str, str],
    output_prefix: str,
    label: str,
) -> None:
    """
    Merge multiple summary files into a single summary table.

    Parameters:
        file_dict (dict): Dictionary with assembly IDs as keys and file paths as values.
        output_file (str): Path to write the merged summary file.

    Output:
        A TSV file.

    Example of the Interpro summary file:
    count	interpro_accession	description
    16503	IPR036291	NAD(P)-binding domain superfamily
    14694	IPR019734	Tetratricopeptide repeat

    Example of the GO summary file:
    go	term	category	count
    GO:0016020	membrane	cellular_component	30626
    GO:0005524	ATP binding	molecular_function	30524
    """
    check_files_exist(list(file_dict.values()))

    output_file_name = f"{output_prefix}_{label}_summary.tsv"

    merged_df = None
    for assembly_acc, filepath in file_dict.items():
        try:
            df = pd.read_csv(filepath, sep="\t")
        except pd.errors.EmptyDataError:
            logging.warning(f"File {filepath} is empty. Skipping.")
            continue

        if label in ["go", "goslim"]:
            validate_dataframe(df, GOSummarySchema, str(filepath))
        elif label == "interpro":
            validate_dataframe(df, InterProSummarySchema, str(filepath))

        df = df.rename(columns={**column_names, "count": assembly_acc})

        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(
                merged_df, df, on=list(column_names.values()), how="outer"
            )

    if merged_df is None:
        logging.warning(
            f"No valid files with functional annotation summary were found. Skipping creation of {output_file_name}."
        )
        return

    # Fill NaNs with 0 and make sure count columns are integers
    count_columns = [
        col for col in merged_df.columns if col not in column_names.values()
    ]
    merged_df[count_columns] = merged_df[count_columns].fillna(0).astype(int)

    cols = list(column_names.values()) + [
        col for col in merged_df.columns if col not in column_names.values()
    ]
    merged_df = merged_df[cols]

    merged_df.to_csv(output_file_name, sep="\t", index=False)


@cli.command(
    "summarise",
    options_metavar="-r <assemblies> -a <study_dir> -p <output_prefix>",
    short_help="Generate study-level analysis summaries.",
)
@click.option(
    "-r",
    "--assemblies",
    required=True,
    help="CSV file containing successful analyses generated by the pipeline",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-a",
    "--study_dir",
    required=True,
    help="Input directory to where all the individual analyses subdirectories for summarising",
    type=click.Path(exists=True, path_type=Path, file_okay=False),
)
@click.option(
    "-p",
    "--output_prefix",
    required=True,
    help="Prefix to summary summary_files",
    type=str,
)
def summarise_analyses(assemblies: Path, study_dir: Path, output_prefix: str) -> None:
    """Generate study-level summary summary_files for successfuly proccessed assemblies.

    :param assemblies: Path to a file listing successful assembly accessions and their status.
    :param study_dir: Path to the directory containing analysis results for each assembly.
    :param output_prefix: Prefix for the generated summary summary_files.
    """
    # TODO: this file with successful jobs is not yet published by pipeline
    assemblies_df = pd.read_csv(assemblies, names=["assembly", "status"])
    CompletedAnalysisSchema(assemblies_df)
    assembly_list = assemblies_df["assembly"].tolist()

    def get_file_paths(subdir: str, filename_template: str) -> dict[str, Path]:
        return {
            acc: study_dir / acc / subdir / filename_template.format(acc=acc)
            for acc in assembly_list
        }

    # TODO handle both gz and non-gz summary_files
    generate_taxonomy_summary(
        get_file_paths("taxonomy", "{acc}.krona.txt"),
        f"{output_prefix}_taxonomy_summary.tsv",
    )

    generate_functional_summary(
        get_file_paths("functional-annotation/interpro", "{acc}_interpro_summary.tsv"),
        INTERPRO_COLUMN_NAMES,
        output_prefix,
        "interpro",
    )

    generate_functional_summary(
        get_file_paths("functional-annotation/go", "{acc}_go_summary.tsv"),
        GO_COLUMN_NAMES,
        output_prefix,
        "go",
    )

    generate_functional_summary(
        get_file_paths("functional-annotation/go", "{acc}_goslim_summary.tsv"),
        GO_COLUMN_NAMES,
        output_prefix,
        "goslim",
    )


@cli.command(
    "merge",
    options_metavar="-a <study_dir> -p <output_prefix>",
    short_help="Merge multiple study-level analysis summaries.",
)
@click.option(
    "-a",
    "--study_dir",
    required=True,
    help="Input directory to where all the individual analyses subdirectories for merging",
    type=click.Path(exists=True, file_okay=False),
)
@click.option(
    "-p",
    "--output_prefix",
    required=True,
    help="Prefix to merged summary summary_files",
    type=str,
)
def merge_summaries(study_dir: str, output_prefix: str) -> None:
    """Function that will take a file path containing study-level
    summaries that should be merged together.
    \f

    :param study_dir: The filepath to the directory containing all of the analyses.
    :type study_dir: str
    :param output_prefix: Prefix to be added to the generated summary file.
    :type output_prefix: str
    """

    def get_file_paths(summary_type: str) -> list[str]:
        return glob.glob(f"{study_dir}/*_{summary_type}_summary.tsv")

    merge_taxonomy_summaries(
        get_file_paths("taxonomy"), f"{output_prefix}_taxonomy_summary.tsv"
    )

    merge_functional_summaries(
        get_file_paths("interpro"),
        list(INTERPRO_COLUMN_NAMES.values()),
        f"{output_prefix}_interpro_summary.tsv",
    )

    merge_functional_summaries(
        get_file_paths("go"),
        list(GO_COLUMN_NAMES.values()),
        f"{output_prefix}_go_summary.tsv",
    )

    merge_functional_summaries(
        get_file_paths("goslim"),
        list(GO_COLUMN_NAMES.values()),
        f"{output_prefix}_goslim_summary.tsv",
    )


def merge_taxonomy_summaries(summary_files: list[str], output_file_name: str) -> None:
    """
    TODO: Add docstring

    Input summary file example:
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
        # TODO: Add validation for the dataframe
        # and check if dataframe is empty, raise error in this case
        summary_dfs.append(df)
    merged_df = pd.concat(summary_dfs, axis=1).fillna(0).astype(int)
    merged_df = merged_df.reindex(sorted(merged_df.columns), axis=1)
    merged_df.to_csv(
        output_file_name,
        sep="\t",
        index_label="taxonomy",
    )


def merge_functional_summaries(
    summary_files: list[str], merge_keys: list[str], output_file_name: str
) -> None:
    """
    TODO: Add docstring

    Input summary file example:
    GO	description	category	ERZ1049444	ERZ1049446
    GO:0016020	membrane	cellular_component	30626	673
    GO:0005524	ATP binding	molecular_function	30524	2873
    """
    if not summary_files:
        logging.warning(
            f"Skipping creation of {output_file_name} because no summaries were found for this type of functional annotation."
        )
        return

    # TODO validate the summary_files, can't be empty

    merged_df = pd.read_csv(summary_files[0], sep="\t")
    if len(summary_files) > 1:
        for path in summary_files[1:]:
            df = pd.read_csv(path, sep="\t")
            merged_df = pd.merge(merged_df, df, on=merge_keys, how="outer")

        # Fill NaNs with 0 and make sure count columns are integers
        count_columns = [col for col in merged_df.columns if col not in merge_keys]

        # Reorder columns: merge_keys first, then sorted count columns
        sorted_columns = merge_keys + sorted(count_columns)
        merged_df = merged_df[sorted_columns]
        merged_df[count_columns] = merged_df[count_columns].fillna(0).astype(int)

    merged_df.to_csv(output_file_name, sep="\t", index=False)


if __name__ == "__main__":
    cli()
