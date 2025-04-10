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
import re

from enum import Enum
from typing import ClassVar, Optional, Type

import pandas as pd
import pandera as pa
from pandera.typing import Series
from pandera.typing.common import DataFrameBase

from pydantic import (
    Field,
    BaseModel,
    field_validator,
    RootModel,
)
from pandera.engines.pandas_engine import PydanticModel

from mgnify_pipelines_toolkit.constants.tax_ranks import (
    SHORT_TAX_RANKS,
    SHORT_PR2_TAX_RANKS,
)


class INSDCRunAccession(RootModel):
    """Class for modelling for INSDC-specific run accessions.
    Essentially is just a special string with regex-based validation of the accession.
    """

    # RootModel example:
    # https://stackoverflow.com/questions/78393675/how-to-make-a-custom-type-inheriting-from-uuid-work-as-a-pydantic-model

    root: str = Field(
        unique=True,
        description="The run needs to be a valid ENA accession",
        examples=["ERR123456", "DRR789012", "SRR345678"],
    )

    @field_validator("root", mode="after")
    @classmethod
    def run_validity_check(cls, run: str) -> bool:
        """Checks that the run string matches the regex code of an INSDC run accession.
        Throws a `ValueError` exception if not, which is what Pydantic prefers for validation errors.
        """

        run_accession_regex = "(E|D|S)RR[0-9]{6,}"
        regex_res = re.match(run_accession_regex, run)

        if regex_res is None:
            raise ValueError(
                f"Accession `{run}` does not fit INSDC format [ERR*,SRR*,DRR*]."
            )

        return run


class AmpliconResultTypes(str, Enum):
    """Class that models the two allowed statuses for successful amplicon analysis runs.
    Pydantic validates Enums very simply without needing to declare a new function.
    """

    all_results = "all_results"
    no_asvs = "no_asvs"


class AmpliconPassedRunsRecord(BaseModel):
    """Class defining a Pydantic model for a single "row" of an amplicon passed runs file.
    Uses the previous two classes.
    """

    run: INSDCRunAccession
    status: AmpliconResultTypes


class AmpliconNonINSDCSPassedRunsRecord(BaseModel):
    """Class modeling a very similar model as the preceding one, but with no INSDC-validation.
    This is achieved by replacing the type of the runs with just a simple string so no validation
    happens.
    """

    run: str
    status: AmpliconResultTypes


# This is the schema for the whole DF
class AmpliconPassedRunsSchema(pa.DataFrameModel):
    """Class modelling a Pandera dataframe schema that uses the AmpliconPassedRunsRecord class as dtype.
    This is what actually validates the generated dataframe when read by pandas.read_csv.
    """

    class Config:
        """Config with dataframe-level data type."""

        dtype = PydanticModel(AmpliconPassedRunsRecord)
        coerce = True


class CompletedAnalysisRecord(BaseModel):
    """Class defining a Pydantic model for a single "row" of an successfully analysed assemblies file.
    Uses the previous two classes.
    """

    assembly: str = Field(..., description="Assembly accession", examples=["ERZ789012"])
    status: str = Field(..., description="")

    @field_validator("assembly")
    @classmethod
    def validate_format(cls, accession: str) -> str:
        """
        Checks that the analysis string matches the regex code of an INSDC analysis accession.
        Throws a `ValueError` exception if not, which is what Pydantic prefers for validation errors.
        """
        if not re.match(r"ERZ[0-9]{6,}", accession):
            raise ValueError("Invalid accession")
        return accession

    @field_validator("status")
    @classmethod
    def validate_status(cls, value: str) -> str:
        allowed = {"success"}
        if value not in allowed:
            raise ValueError(
                f"Invalid assembly job status: {value}. Allowed values: {', '.join(allowed)}"
            )
        return value


class CompletedAnalysisSchema(pa.DataFrameModel):
    """Class modelling a Pandera dataframe schema that uses the AnalysisPassedAssembliesRecord class as dtype.
    This is what actually validates the generated dataframe when read by pandas.read_csv.
    """

    assembly: Series[str]

    @pa.check("assembly")
    def accessions_unique(cls, series: Series[str]) -> Series[bool]:
        return ~series.duplicated()

    class Config:
        """Config with dataframe-level data type."""

        dtype = PydanticModel(CompletedAnalysisRecord)
        coerce = True


class InterProSummaryRecord(BaseModel):
    """Model of a row in the InterPro summary file."""

    count: int = Field(
        ..., ge=0, description="Number of hits for the InterPro accession"
    )
    interpro_accession: str = Field(
        ..., description="InterPro accession ID", examples=["IPR123456"]
    )
    description: str = Field(..., description="Description of the InterPro domain")

    @field_validator("interpro_accession")
    @classmethod
    def check_interpro_format(cls, id) -> str:
        # TODO: check if this is the correct regex for InterPro accessions
        if not re.fullmatch(r"IPR\d{6}", id):
            raise ValueError(f"Invalid InterPro accession: {id}")
        return id


class GOSummaryRecord(BaseModel):
    """Model of a row in the GO summary file."""

    go: str = Field(..., description="GO term identifier", examples=["GO:1234567"])
    term: str = Field(..., description="GO term name")
    category: str = Field(
        ...,
        description="GO category",
        examples=["biological_process", "molecular_function", "cellular_component"],
    )
    count: int = Field(..., ge=0, description="Number of times the GO term is observed")

    @field_validator("go")
    @classmethod
    def check_goterm_format(cls, id) -> str:
        # TODO: check if this is the correct regex
        if not re.fullmatch(r"GO:\d{7}", id):
            raise ValueError(f"Invalid GO term format: {id}")
        return id


class InterProSummarySchema(pa.DataFrameModel):
    """Schema for InterPro summary file validation."""

    interpro_accession: Series[str]

    @pa.check("interpro_accession")
    def interpro_ids_unique(cls, series: Series[str]) -> Series[bool]:
        return ~series.duplicated()

    class Config:
        dtype = PydanticModel(InterProSummaryRecord)
        coerce = True


class GOSummarySchema(pa.DataFrameModel):
    """Schema for GO or GOslim summary file validation."""

    go: Series[str]

    @pa.check("go")
    def go_ids_unique(cls, series: Series[str]) -> Series[bool]:
        return ~series.duplicated()

    class Config:
        dtype = PydanticModel(GOSummaryRecord)
        coerce = True


class GOStudySummarySchema(pa.DataFrameModel):
    """Schema for GO or GOslim study summary file validation."""

    GO: Series[str] = pa.Field(str_matches=r"^GO:\d{7}$")
    description: Series[str]
    category: Series[str]

    @pa.check(regex=r"^ERZ\d+")
    def count_columns_are_non_negative(cls, s: Series[int]) -> Series[bool]:
        return s >= 0

    class Config:
        strict = False  # allow extra ERZ* columns not declared above


class InterProStudySummarySchema(pa.DataFrameModel):
    """Schema for InterPro study summary file validation."""

    IPR: Series[str] = pa.Field(str_matches=r"^IPR\d{6}$")
    description: Series[str]

    @pa.check(regex=r"^ERZ\d+")
    def count_columns_are_non_negative(cls, s: Series[int]) -> Series[bool]:
        return s >= 0

    class Config:
        strict = False  # allow extra ERZ* columns not declared above


class TaxonomyStudySummarySchema(pa.DataFrameModel):
    """Schema for validation of study summary file with taxonomic classification."""

    @pa.check(regex=r"^ERZ\d+")
    def count_columns_are_non_negative(cls, s: Series[int]) -> Series[bool]:
        return s >= 0

    class Config:
        strict = False  # allow extra ERZ* columns not declared above


class AmpliconNonINSDCPassedRunsSchema(pa.DataFrameModel):
    """Class modelling the same dataframe schema as the preceding one, except with no INSDC validation.
    Uses the AmpliconNonINSDCSPassedRunsRecord as a dtype to achieve this.
    """

    class Config:
        """Config with dataframe-level data type."""

        dtype = PydanticModel(AmpliconNonINSDCSPassedRunsRecord)
        coerce = True


class TaxRank(RootModel):
    """Class for modelling a single Taxonomic Rank.
    Essentially is just a special string with validation of the structure:
    `${rank}__${taxon}`
    Where `${rank}` is one of the allowed short ranks defined by the imported
    `SHORT_TAX_RANKS` and `SHORT_PR2_TAX_RANKS` variables.
    And `${taxon}` is the actual taxon for that rank (this isn't validated).
    It will also validate if the whole string is the permitted "Unclassified".
    """

    valid_tax_ranks: ClassVar = SHORT_TAX_RANKS + SHORT_PR2_TAX_RANKS

    root: str = Field(
        unique=True,
        description="A single taxon in a taxonomy record",
        examples=["sk__Bacteria", "p__Bacillota", "g__Tundrisphaera"],
    )

    @field_validator("root", mode="after")
    @classmethod
    def rank_structure_validity_check(cls, taxrank: str) -> bool:
        taxrank_list = taxrank.split("__")
        rank = taxrank_list[0]
        if rank != "" and rank != "Unclassified" and rank not in cls.valid_tax_ranks:
            raise ValueError(f"Invalid taxonomy rank {rank}.")

        return taxrank


# TODO: see if we can simplify the declaration of two Taxon classes by using one of these solutions
# None of the solutions have a model-only way of doing it, but worth considering maybe
# https://stackoverflow.com/questions/76537360/initialize-one-of-two-pydantic-models-depending-on-an-init-parameter


class Taxon(BaseModel):
    """Class for modelling an entire Taxon or taxonomic assignment.
    All of the ranks are optional, to model for the taxon being "Unclassified".
    """

    Superkingdom: Optional[TaxRank] = None
    Kingdom: Optional[TaxRank] = None
    Phylum: Optional[TaxRank] = None
    Class: Optional[TaxRank] = None
    Order: Optional[TaxRank] = None
    Family: Optional[TaxRank] = None
    Genus: Optional[TaxRank] = None
    Species: Optional[TaxRank] = None


class PR2Taxon(Taxon):
    """Class for modelling the same thing as the preceding class, but for PR2 ranks."""

    Domain: Optional[TaxRank] = None
    Supergroup: Optional[TaxRank] = None
    Division: Optional[TaxRank] = None
    Subdivision: Optional[TaxRank] = None


class TaxonRecord(Taxon):
    """Class for modelling a single taxon record in a taxonomy file.
    It inherits the Taxon class, and simply adds a Count field, modelling the read counts
    for that particular Taxon record.
    """

    Count: int


class PR2TaxonRecord(PR2Taxon):
    """Class for modelling the same thing as the preceding class, but for PR2 ranks."""

    Count: int


# This is the schema for the whole DF
class TaxonSchema(pa.DataFrameModel):
    """Class modelling a Pandera dataframe schema that uses the TaxonRecord class as dtype.
    This is what actually validates the generated dataframe when read by pandas.read_csv.
    """

    class Config:
        """Config with dataframe-level data type."""

        dtype = PydanticModel(TaxonRecord)
        coerce = True


class PR2TaxonSchema(pa.DataFrameModel):
    """Class modelling the same dataframe schema as the preceding one, except for the PR2 taxonomy.
    Uses the PR2TaxonSchema as a dtype to achieve this.
    """

    class Config:
        """Config with dataframe-level data type."""

        dtype = PydanticModel(PR2TaxonRecord)
        coerce = True


def validate_dataframe(
    df: pd.DataFrame, schema: Type[pa.DataFrameModel], df_metadata: str
) -> DataFrameBase:
    """
    Validate a pandas dataframe using a pandera schema.
    df_metadata will be shown in logs on failure: example, the TSV filename from which the df was read.
    """
    try:
        dfs = schema.validate(df, lazy=True)
    except pa.errors.SchemaErrors as e:
        logging.error(f"{schema.__name__} validation failure for {df_metadata}")
        raise e
    return dfs
