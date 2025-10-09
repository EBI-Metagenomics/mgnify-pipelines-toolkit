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

import logging
from typing import Type

import pandas as pd
import pandera as pa
from pandera.typing import Series
from pandera.typing.common import DataFrameBase

from mgnify_pipelines_toolkit.constants.tax_ranks import (
    SHORT_MOTUS_TAX_RANKS,
    SHORT_PR2_TAX_RANKS,
    SHORT_TAX_RANKS,
)


class CoerceBaseDataFrameSchema(pa.DataFrameModel):
    """Base schema for all dataframe models.

    Provides common configuration for automatic type coercion.
    """

    class Config:
        """Pandera configuration.

        coerce: Automatically convert column dtypes to match schema
        """

        coerce = True


# This is the schema for the whole DF
class AmpliconPassedRunsSchema(CoerceBaseDataFrameSchema):
    """Class modelling a Pandera dataframe schema for amplicon passed runs.
    Validates the generated dataframe when read by pandas.read_csv.
    """

    run: Series[str] = pa.Field(str_matches=r"(E|D|S)RR[0-9]{6,}", unique=True)
    status: Series[str] = pa.Field(isin=["all_results", "no_asvs", "dada2_stats_fail"])


class CompletedAnalysisSchema(CoerceBaseDataFrameSchema):
    """Class modelling a Pandera dataframe schema for completed assemblies.
    Validates the generated dataframe when read by pandas.read_csv.
    """

    assembly: Series[str] = pa.Field(str_matches=r"ERZ\d{6,}", unique=True)
    status: Series[str] = pa.Field(isin=["success"])


class BaseSummarySchema(CoerceBaseDataFrameSchema):
    """Base schema for summary files.

    All summary schemas inherit from this base and use coerce=True by default.
    """


class InterProSummarySchema(BaseSummarySchema):
    """Schema for InterPro summary file validation."""

    count: Series[int] = pa.Field(ge=0)
    interpro_accession: Series[str] = pa.Field(str_matches=r"IPR\d{6}", unique=True)
    description: Series[str]


class GOSummarySchema(BaseSummarySchema):
    """Schema for GO or GOslim summary file validation."""

    go: Series[str] = pa.Field(str_matches=r"GO:\d{7}", unique=True)
    term: Series[str]
    category: Series[str]
    count: Series[int] = pa.Field(ge=0)


class SanntisSummarySchema(BaseSummarySchema):
    """Schema for Sanntis summary file validation."""

    nearest_mibig: Series[str] = pa.Field(str_matches=r"BGC\d{7}", unique=True)
    nearest_mibig_class: Series[str]
    description: Series[str]
    count: Series[int] = pa.Field(ge=0)


class AntismashSummarySchema(BaseSummarySchema):
    """Schema for Antismash summary file validation."""

    label: Series[str] = pa.Field(unique=True)
    description: Series[str]
    count: Series[int] = pa.Field(ge=0)


class KOSummarySchema(BaseSummarySchema):
    """Schema for KEGG Orthology summary file validation."""

    ko: Series[str] = pa.Field(str_matches=r"K\d{5,}", unique=True)
    description: Series[str]
    count: Series[int] = pa.Field(ge=0)


class PFAMSummarySchema(BaseSummarySchema):
    """Schema for PFAM summary file validation."""

    pfam: Series[str] = pa.Field(str_matches=r"PF\d{5}", unique=True)
    description: Series[str]
    count: Series[int] = pa.Field(ge=0)


class KEGGModulesSummarySchema(BaseSummarySchema):
    """Schema for KEGG Modules summary file validation."""

    module_accession: Series[str] = pa.Field(str_matches=r"M\d{5}", unique=True)
    completeness: Series[float] = pa.Field(ge=0)
    pathway_name: Series[str]
    pathway_class: Series[str]


class GOStudySummarySchema(BaseSummarySchema):
    GO: Series[str] = pa.Field(str_matches=r"^GO:\d{7}$", unique=True)
    description: Series[str]
    category: Series[str]


class InterProStudySummarySchema(BaseSummarySchema):
    IPR: Series[str] = pa.Field(str_matches=r"^IPR\d{6}$", unique=True)
    description: Series[str]


class AntismashStudySummarySchema(BaseSummarySchema):
    label: Series[str] = pa.Field(unique=True)


class SanntisStudySummarySchema(BaseSummarySchema):
    # TODO: limit mibig to the avaiable mibig categories
    nearest_mibig: Series[str] = pa.Field(unique=True)


class KOStudySummarySchema(BaseSummarySchema):
    KO: Series[str] = pa.Field(unique=True)


class PFAMStudySummarySchema(BaseSummarySchema):
    PFAM: Series[str] = pa.Field(unique=True)


class KEGGModulesStudySummarySchema(BaseSummarySchema):
    module_accession: Series[str] = pa.Field(unique=True)


class TaxonomyStudySummarySchema(BaseSummarySchema):
    pass


class AmpliconNonINSDCPassedRunsSchema(CoerceBaseDataFrameSchema):
    """Class modelling the same dataframe schema as the preceding one, except with no INSDC validation."""

    run: Series[str]
    status: Series[str] = pa.Field(isin=["all_results", "no_asvs"])


# This is the schema for the whole DF
class TaxonSchema(CoerceBaseDataFrameSchema):
    """Class modelling a Pandera dataframe schema for taxonomy records.
    Validates the generated dataframe when read by pandas.read_csv.
    """

    Superkingdom: Series[str] = pa.Field(nullable=True)
    Kingdom: Series[str] = pa.Field(nullable=True)
    Phylum: Series[str] = pa.Field(nullable=True)
    Class: Series[str] = pa.Field(nullable=True)
    Order: Series[str] = pa.Field(nullable=True)
    Family: Series[str] = pa.Field(nullable=True)
    Genus: Series[str] = pa.Field(nullable=True)
    Species: Series[str] = pa.Field(nullable=True)
    Count: Series[int]

    @pa.check(r"Superkingdom|Kingdom|Phylum|Class|Order|Family|Genus|Species", regex=True)
    def validate_tax_rank_format(self, series: Series[str]) -> Series[bool]:
        """Validate that taxonomy rank values follow the format: ${rank}__${taxon}
        or are 'Unclassified' or empty/null.

        :param series: Column series to validate
        :return: Boolean series indicating valid rows
        """
        valid_ranks = ["sk", "k", "p", "c", "o", "f", "g", "s"]

        def check_format(val):
            if pd.isna(val) or val == "" or val.capitalize() == "Unclassified":
                return True
            if "__" not in val:
                return False
            rank = val.split("__")[0]
            return rank in valid_ranks or rank == ""

        return series.apply(check_format)


class PR2TaxonSchema(CoerceBaseDataFrameSchema):
    """Class modelling a Pandera dataframe schema for PR2 taxonomy records."""

    Domain: Series[str] = pa.Field(nullable=True)
    Supergroup: Series[str] = pa.Field(nullable=True)
    Division: Series[str] = pa.Field(nullable=True)
    Subdivision: Series[str] = pa.Field(nullable=True)
    Class: Series[str] = pa.Field(nullable=True)
    Order: Series[str] = pa.Field(nullable=True)
    Family: Series[str] = pa.Field(nullable=True)
    Genus: Series[str] = pa.Field(nullable=True)
    Species: Series[str] = pa.Field(nullable=True)
    Count: Series[int]

    @pa.check(r"Domain|Supergroup|Division|Subdivision|Class|Order|Family|Genus|Species", regex=True)
    def validate_pr2_tax_rank_format(self, series: Series[str]) -> Series[bool]:
        """Validate that PR2 taxonomy rank values follow the format: ${rank}__${taxon}
        or are 'Unclassified' or empty/null.

        :param series: Column series to validate
        :return: Boolean series indicating valid rows
        """
        valid_ranks = SHORT_TAX_RANKS + SHORT_PR2_TAX_RANKS

        def check_format(val):
            if pd.isna(val) or val == "" or val.capitalize() == "Unclassified":
                return True
            if "__" not in val:
                return False
            rank = val.split("__")[0]
            return rank in valid_ranks or rank == ""

        return series.apply(check_format)


# This is the schema for the whole DF
class RawReadsPassedRunsSchema(CoerceBaseDataFrameSchema):
    """Class modelling a Pandera dataframe schema for raw reads passed runs.
    Validates the generated dataframe when read by pandas.read_csv.
    """

    run: Series[str] = pa.Field(str_matches=r"(E|D|S)RR[0-9]{6,}", unique=True)
    status: Series[str] = pa.Field(isin=["all_results", "no_reads", "all_empty_results", "some_empty_results"])


class RawReadsNonINSDCPassedRunsSchema(CoerceBaseDataFrameSchema):
    """Class modelling the same dataframe schema as the preceding one, except with no INSDC validation."""

    run: Series[str]
    status: Series[str] = pa.Field(isin=["all_results", "no_reads", "all_empty_results", "some_empty_results"])


class MotusTaxonSchema(CoerceBaseDataFrameSchema):
    """Class for modelling a single Taxonomic Rank in mOTUs output.
    Essentially is just a special string with validation of the structure:
    `${rank}__${taxon}`
    Where `${rank}` is one of the allowed short ranks defined by the imported
    `SHORT_MOTUS_TAX_RANKS` variables.
    And `${taxon}` is the actual taxon for that rank (this isn't validated).
    It will also validate if the whole string is the permitted "unassigned" or "unclassified".
    """

    Kingdom: Series[str] = pa.Field(nullable=True)
    Phylum: Series[str] = pa.Field(nullable=True)
    Class: Series[str] = pa.Field(nullable=True)
    Order: Series[str] = pa.Field(nullable=True)
    Family: Series[str] = pa.Field(nullable=True)
    Genus: Series[str] = pa.Field(nullable=True)
    Species: Series[str] = pa.Field(nullable=True)
    Count: Series[int]

    @pa.check(r"Kingdom|Phylum|Class|Order|Family|Genus|Species", regex=True)
    def validate_motus_tax_rank_format(self, series: Series[str]) -> Series[bool]:
        """Validate that mOTUs taxonomy rank values follow the format: ${rank}__${taxon}
        or are 'Unclassified', 'Unassigned', or empty/null.

        :param series: Column series to validate
        :return: Boolean series indicating valid rows
        """
        valid_ranks = SHORT_MOTUS_TAX_RANKS

        def check_format(val):
            if pd.isna(val) or val == "":
                return True
            if val.capitalize() in {"Unclassified", "Unassigned"}:
                return True
            if "__" not in val:
                return False
            rank = val.split("__")[0]
            return rank in valid_ranks or rank == ""

        return series.apply(check_format)


class FunctionProfileSchema(CoerceBaseDataFrameSchema):
    """Class modelling a Pandera dataframe schema for functional profile data.
    This is what actually validates the generated dataframe when read by pandas.read_csv.
    """

    read_count: Series[int]
    coverage_depth: Series[float]
    coverage_breadth: Series[float]


def validate_dataframe(df: pd.DataFrame, schema: Type[pa.DataFrameModel], df_metadata: str) -> DataFrameBase:
    """
    Validate a pandas dataframe using a pandera schema.
    df_metadata will be shown in logs on failure: example, the TSV filename from which the df was read.
    """
    try:
        dfs = schema.validate(df, lazy=True)
    except pa.errors.SchemaError as e:
        logging.error(f"{schema.__name__} validation failure for {df_metadata}")
        raise e
    return dfs
