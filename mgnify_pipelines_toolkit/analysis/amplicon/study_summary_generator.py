#!/usr/bin/env python
# -*- coding: utf-8 -*-
import shutil
from shutil import SameFileError

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

import click
from collections import defaultdict
import glob
import logging
from pathlib import Path
from typing import Union, List

import pandas as pd

from mgnify_pipelines_toolkit.constants.db_labels import TAXDB_LABELS, ASV_TAXDB_LABELS
from mgnify_pipelines_toolkit.constants.tax_ranks import (
    _SILVA_TAX_RANKS,
    _PR2_TAX_RANKS,
)
from mgnify_pipelines_toolkit.schemas.schemas import (
    AmpliconPassedRunsSchema,
    AmpliconNonINSDCPassedRunsSchema,
    TaxonSchema,
    PR2TaxonSchema,
    validate_dataframe,
)

logging.basicConfig(level=logging.DEBUG)


@click.group()
def cli():
    pass


def get_tax_file(
    run_acc: str, analyses_dir: Path, db_label: str
) -> Union[Path, List[Path]]:
    """Takes path information for a particular analysis and db_label combo, and returns any existing files.

    :param run_acc: Run accession for the tax file that should be retrieved.
    :type run_acc: str
    :param analyses_dir: The path to the directory containing all of the analyses,
            including the tax file corresponding to :param:`run_acc`.
    :type analyses_dir: Path
    :param db_label: One of the database labels that results might exist for,
            values of which come from the imported constants ``TAXDB_LABELS`` and ``ASV_TAXDB_LABELS``.
    :type db_label: str
    :return: Either a :class:`Path` object if :param:`db_label` comes from ``TAXDB_LABELS``,
            or a list of :class:`Path` objects if from ``ASV_TAXDB_LABELS``.
    :rtype: Union[Path, List[Path]]
    """

    tax_file = None

    db_path = Path(f"{analyses_dir}/{run_acc}/taxonomy-summary/{db_label}")

    if not db_path.exists():
        logging.debug(
            f"DB {db_path} doesn't exist for {run_acc}. Skipping"
        )  # or error?
        return

    if db_label in TAXDB_LABELS:
        tax_file = Path(
            f"{analyses_dir}/{run_acc}/taxonomy-summary/{db_label}/{run_acc}_{db_label}.txt"
        )
        if not tax_file.exists():
            logging.error(
                f"DB path exists but file doesn't - exiting. Path: {tax_file}"
            )
            exit(1)

        file_size = tax_file.stat().st_size
        if (
            file_size == 0
        ):  # Pipeline can generate files that are empty for ITS DBs (UNITE and ITSoneDB),
            # so need to skip those. Should probably fix that at some point
            logging.debug(
                f"File {tax_file} exists but is empty, so will be skipping it."
            )
            tax_file = None
    elif db_label in ASV_TAXDB_LABELS:
        # ASV tax files could have up to two files, one for each amplified region (maximum two from the pipeline).
        # So will need to handle this differently to closed-reference files
        asv_tax_files = glob.glob(
            f"{analyses_dir}/{run_acc}/taxonomy-summary/{db_label}/*.txt"
        )
        asv_tax_files = [
            Path(file) for file in asv_tax_files if "concat" not in file
        ]  # Have to filter out concatenated file if it exists

        tax_file = asv_tax_files

    return tax_file


def parse_one_tax_file(
    run_acc: str, tax_file: Path, long_tax_ranks: list
) -> pd.DataFrame:
    """Parses a taxonomy file, and returns it as a pandas DataFrame object.

    :param run_acc: Run accession of the taxonomy file that will be parsed.
    :type run_acc: str
    :param tax_file: Taxonomy file that will be parsed.
    :type tax_file: Path
    :param long_tax_ranks: Either the imported list _SILVA_TAX_RANKS or _PR2_TAX_RANKS
            to validate the taxonomic ranks of the file.
    :type tax_file: list
    :return: The parsed :param:`tax_file` as a :class:`pd.DataFrame` object
    :rtype: pd.DataFrame
    """

    res_df = pd.read_csv(tax_file, sep="\t", names=["Count"] + long_tax_ranks)
    res_df = res_df.fillna("")

    # Two different schemas used for validation depending on the database
    # because PR2 schema has different taxonomic ranks than the standard
    if len(long_tax_ranks) == 8:
        validate_dataframe(res_df, TaxonSchema, str(tax_file))
    elif len(long_tax_ranks) == 9:
        validate_dataframe(res_df, PR2TaxonSchema, str(tax_file))

    res_df["full_taxon"] = res_df.iloc[:, 1:].apply(
        lambda x: ";".join(x).strip(";"), axis=1
    )
    final_df = res_df.iloc[:, [0, -1]]
    final_df = final_df.set_index("full_taxon")
    final_df.columns = [run_acc]

    return final_df


def generate_db_summary(
    db_label: str, tax_dfs: defaultdict[Path], output_prefix: str
) -> None:
    """Takes paired run accessions taxonomy dataframes in the form of a dictionary,
    and respective db_label, joins them together, and generates a study-wide summary
    in the form of a .tsv file.

    :param db_label: One of the database labels that results might exist for,
            values of which come from the imported constants ``TAXDB_LABELS`` and ``ASV_TAXDB_LABELS``.
    :param tax_dfs: Dictionary where the key is a run accession,
        and values are either one parsed taxonomy dataframe if the :param:db_label comes from ``TAXDB_LABELS``,
        or a list of at least 1 and at most 2 dataframes if it comes from ``ASV_TAXDB_LABELS``.
        These dataframes are parsed by :func:`parse_one_tax_file`
    :type tax_dfs: defaultdict[Path]
    :param output_prefix: Prefix to be added to the generated summary file.
    :type output_prefix: str
    """

    if db_label in TAXDB_LABELS:
        df_list = []

        if "PR2" in db_label:
            long_tax_ranks = _PR2_TAX_RANKS
        else:
            long_tax_ranks = _SILVA_TAX_RANKS

        for run_acc, tax_df in tax_dfs.items():
            res_df = parse_one_tax_file(run_acc, tax_df, long_tax_ranks)
            df_list.append(res_df)

        res_df = pd.concat(df_list, axis=1).fillna(0)
        res_df = res_df.sort_index()
        res_df = res_df.astype(int)

        res_df.to_csv(
            f"{output_prefix}_{db_label}_study_summary.tsv",
            sep="\t",
            index_label="taxonomy",
        )

    elif db_label in ASV_TAXDB_LABELS:

        if "PR2" in db_label:
            long_tax_ranks = _PR2_TAX_RANKS
        else:
            long_tax_ranks = _SILVA_TAX_RANKS

        amp_region_dict = defaultdict(list)

        for (
            run_acc,
            tax_df_asv_lst,
        ) in (
            tax_dfs.items()
        ):  # each `tax_file` will be a list containing at most two files (one for each amp_region)
            for tax_df in tax_df_asv_lst:
                amp_region = str(tax_df).split("_")[
                    -5
                ]  # there are a lot of underscores in these names... but it is consistent
                # e.g. ERR4334351_16S-V3-V4_DADA2-SILVA_asv_krona_counts.txt
                amp_region_df = parse_one_tax_file(run_acc, tax_df, long_tax_ranks)
                amp_region_dict[amp_region].append(amp_region_df)

        for amp_region, amp_region_dfs in amp_region_dict.items():
            if amp_region_dfs:
                amp_res_df = amp_region_dfs[0]
                for amp_df in amp_region_dfs[1:]:
                    amp_res_df = amp_res_df.join(amp_df, how="outer")
                amp_res_df = amp_res_df.fillna(0)
                amp_res_df = amp_res_df.astype(int)

                amp_res_df.to_csv(
                    f"{output_prefix}_{db_label}_{amp_region}_asv_study_summary.tsv",
                    sep="\t",
                    index_label="taxonomy",
                )


def organise_study_summaries(all_study_summaries: List[str]) -> defaultdict[List]:
    """Matches different summary files of the same database label and analysis
    type (and amplified region for ASVs) into a dictionary to help merge
    the correct summaries.

    :param all_study_summaries: List of file paths to different summary files
    :type all_study_summaries: List[str]
    :return: Organised dictionary where each summary is paired to a specific
        database label key to be merged together.
    :rtype: defaultdict[List]
    """
    summaries_dict = defaultdict(list)

    for summary in all_study_summaries:
        summary_path = Path(summary)
        summary_filename = summary_path.stem

        temp_lst = summary_filename.split("_")
        if "asv_study_summary" in summary_filename:
            summary_db_label = "_".join(
                temp_lst[1:3]
            )  # For ASVs we need to include the amp_region in the label
        else:
            summary_db_label = temp_lst[
                1
            ]  # For closed reference, just the db_label is needed

        summaries_dict[summary_db_label].append(summary_path)

    return summaries_dict


@cli.command(
    "summarise",
    options_metavar="-r <runs> -a <analyses_dir> -p <output_prefix>",
    short_help="Generate study-level summaries of amplicon analysis results.",
)
@click.option(
    "-r",
    "--runs",
    required=True,
    help="CSV file containing successful analyses generated by the pipeline",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-a",
    "--analyses_dir",
    required=True,
    help="Input directory to where all the individual analyses subdirectories for summarising",
    type=click.Path(exists=True, path_type=Path, file_okay=False),
)
@click.option(
    "-p", "--output_prefix", required=True, help="Prefix to summary files", type=str
)
@click.option(
    "--non_insdc",
    default=False,
    is_flag=True,
    help="If run accessions aren't INSDC-formatted",
)
def summarise_analyses(
    runs: Path, analyses_dir: Path, output_prefix: str, non_insdc: bool
) -> None:
    """Function that will take a file of pipeline-successful run accessions
    that should be used for the generation of the relevant db-specific
    study-level summary files. For ASV results, these will also be on a
    per-amplified-region basis.
    \f

    :param runs: Path to a qc_passed_runs file from the pipeline execution.
        Contains the accessions of runs that should therefore be included in the generated
        summaries.
    :type runs: Path
    :param analyses_dir: The path to the directory containing all of the analyses.
    :type analyses_dir: Path
    :param output_prefix: Prefix to be added to the generated summary file.
    :type output_prefix: str
    """
    runs_df = pd.read_csv(runs, names=["run", "status"])

    if not non_insdc:
        AmpliconPassedRunsSchema(
            runs_df
        )  # Run validation on the successful_runs .csv file
    else:
        AmpliconNonINSDCPassedRunsSchema(runs_df)

    all_db_labels = TAXDB_LABELS + ASV_TAXDB_LABELS
    for db_label in all_db_labels:

        tax_files = defaultdict(Path)
        for i in range(0, len(runs_df)):
            run_acc = runs_df.loc[i, "run"]
            tax_file = get_tax_file(run_acc, analyses_dir, db_label)

            if tax_file:
                tax_files[run_acc] = tax_file

        if tax_files:
            generate_db_summary(db_label, tax_files, output_prefix)


@cli.command(
    "merge",
    options_metavar="-a <analyses_dir> -p <output_prefix>",
    short_help="Merge multiple study-level summaries of amplicon analysis.",
)
@click.option(
    "-a",
    "--analyses_dir",
    required=True,
    help="Input directory to where all the individual analyses subdirectories for merging",
    type=click.Path(exists=True, file_okay=False),
)
@click.option(
    "-p",
    "--output_prefix",
    required=True,
    help="Prefix to merged summary files",
    type=str,
)
def merge_summaries(analyses_dir: str, output_prefix: str) -> None:
    """Function that will take a file path containing study-level
    summaries that should be merged together on a per-db-per-amplified-region
    basis.
    \f

    :param analyses_dir: The filepath to the directory containing all of the analyses.
    :type analyses_dir: str
    :param output_prefix: Prefix to be added to the generated summary file.
    :type output_prefix: str
    """

    all_study_summaries = glob.glob(f"{analyses_dir}/*_study_summary.tsv")

    summaries_dict = organise_study_summaries(all_study_summaries)

    for db_label, summaries in summaries_dict.items():
        merged_summary_name = f"{output_prefix}_{db_label}_study_summary.tsv"
        if len(summaries) > 1:
            res_df = pd.read_csv(summaries[0], sep="\t", index_col=0)
            for summary in summaries[1:]:
                curr_df = pd.read_csv(summary, sep="\t", index_col=0)
                res_df = res_df.join(curr_df, how="outer")
                res_df = res_df.fillna(0)
                res_df = res_df.astype(int)

            res_df = res_df.reindex(sorted(res_df.columns), axis=1)
            res_df.to_csv(
                merged_summary_name,
                sep="\t",
                index_label="taxonomy",
            )
        elif len(summaries) == 1:
            logging.info(
                f"Only one summary ({summaries[0]}) so will use that as {merged_summary_name}"
            )
            try:
                shutil.copyfile(summaries[0], merged_summary_name)
            except SameFileError:
                pass


if __name__ == "__main__":
    cli()
