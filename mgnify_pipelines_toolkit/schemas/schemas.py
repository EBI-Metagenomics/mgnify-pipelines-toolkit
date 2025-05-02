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
from typing import ClassVar, Optional, Type, Literal

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
    """Class defining a Pydantic model for a single "row" of an successfully analysed assemblies file."""

    assembly: str = Field(
        ...,
        description="Assembly accession",
        examples=["ERZ789012"],
        pattern=r"ERZ\d{6,}",
    )
    status: Literal["success"] = Field(
        ...,
        description="Pipeline output for whether this assembly's analysis succeeded or not",
    )


class CompletedAnalysisSchema(pa.DataFrameModel):
    """Class modelling a Pandera dataframe schema that uses the CompletedAnalysisSchema class as dtype.
    This is what actually validates the generated dataframe when read by pandas.read_csv.
    """

    assembly: Series[str]

    @pa.check("assembly")
    def accessions_unique(self, series: Series[str]) -> Series[bool]:
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
        ...,
        description="InterPro accession ID",
        examples=["IPR123456"],
        pattern=r"IPR\d{6}",
    )
    description: str = Field(..., description="Description of the InterPro domain")


class GOSummaryRecord(BaseModel):
    """Model of a row in the GO summary file."""

    go: str = Field(
        ...,
        description="GO term identifier",
        examples=["GO:1234567"],
        pattern=r"GO:\d{7}",
    )
    term: str = Field(..., description="GO term name")
    category: str = Field(
        ...,
        description="GO category",
        examples=["biological_process", "molecular_function", "cellular_component"],
    )
    count: int = Field(..., ge=0, description="Number of times the GO term is observed")


class BaseSummarySchema(pa.DataFrameModel):
    """Base schema for summary files."""

    @staticmethod
    def is_unique(series: Series[str]) -> Series[bool]:
        return ~series.duplicated()


class InterProSummarySchema(BaseSummarySchema):
    """Schema for InterPro summary file validation."""

    interpro_accession: Series[str]

    @pa.check("interpro_accession")
    def interpro_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)

    class Config:
        dtype = PydanticModel(InterProSummaryRecord)
        coerce = True


class GOSummarySchema(BaseSummarySchema):
    """Schema for GO or GOslim summary file validation."""

    go: Series[str]

    @pa.check("go")
    def go_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)

    class Config:
        dtype = PydanticModel(GOSummaryRecord)
        coerce = True


class SanntisSummaryRecord(BaseModel):
    """Model of a row in the Sanntis assembly-level summary file."""

    nearest_mibig: str = Field(
        ...,
        description="The accession ID of the closest matching biosynthetic gene cluster (BGC) in the MIBiG database",
        examples=["BGC0000073"],
        pattern=r"BGC\d{7}",
    )
    nearest_mibig_class: str = Field(
        ...,
        description="The biosynthetic class of the nearest MIBiG BGC",
        examples=["Polyketide"],
    )
    description: str = Field(
        ...,
        description="A brief summary of the biosynthetic process or type of metabolite associated with the nearest MIBiG cluster",
    )

    count: int = Field(
        ..., ge=0, description="Number of times the MIBiG entry is observed"
    )


class AntismashSummaryRecord(BaseModel):
    """Model of a row in the Antismash summary file."""

    label: str = Field(
        ...,
        description="Biosynthetic class or label assigned by Antismash based on sequence similarity to known biosynthetic gene clusters.",
        examples=["RiPP-like", "T1PKS", "terpene"],
    )
    description: str = Field(
        ...,
        description="Brief explanation of the biosynthetic class, often indicating compound type or functional characteristics.",
        examples=["Type I PKS (Polyketide synthase)", "Redox-cofactors such as PQQ"],
    )
    count: int = Field(
        ...,
        ge=0,
        description="Number of BGCs (biosynthetic gene clusters) in the dataset assigned to this label.",
    )


class KOSummaryRecord(BaseModel):
    """Model of a row in the KEGG summary file."""

    ko: str = Field(
        ...,
        description="KEGG Orthology (KO) identifier representing a functional gene or pathway component.",
        examples=["K07547", "K04874", "K19946"],
        pattern=r"K\d{5,}",
    )
    description: str = Field(
        ...,
        description="Name or function of the KO, sometimes including EC numbers and protein families.",
        examples=["optineurin", "MFS transporter, POT/PTR family"],
    )
    count: int = Field(
        ...,
        ge=0,
        description="Number of times this KO identifier is observed in the dataset.",
    )


class PFAMSummaryRecord(BaseModel):
    """Model of a row in the PFAM summary file."""

    pfam: str = Field(
        ...,
        description="PFAM accession identifier representing a protein domain or family.",
        examples=["PF00265", "PF01956", "PF00673"],
        pattern=r"PF\d{5}",
    )
    description: str = Field(
        ...,
        description="Description of the protein domain or family associated with the PFAM ID.",
        examples=["Thymidine kinase", "Integral membrane protein EMC3/TMCO1-like"],
    )
    count: int = Field(
        ...,
        ge=0,
        description="Number of times the PFAM domain is observed in the dataset.",
    )


class KEGGModulesSummaryRecord(BaseModel):
    """Model of a row in the KEGG Modules summary file."""

    module_accession: str = Field(
        ...,
        description="KEGG Module identifier representing a specific metabolic pathway or module.",
        examples=["M00123", "M00234"],
        pattern=r"M\d{5}",
    )
    completeness: float = Field(
        ...,
        ge=0,
        description="Completeness score of the KEGG Module, indicating the extent to which the module is present in the metagenome.",
    )
    pathway_name: str = Field(
        ...,
        description="Name of the metabolic pathway associated with the KEGG Module.",
        examples=["Sulfur reduction, sulfur => sulfide"],
    )
    pathway_class: str = Field(
        ...,
        description="Biosynthetic class or category associated with the KEGG Module, semi colon separated.",
        examples=["Pathway modules; Energy metabolism; Photosynthesis"],
    )


class SanntisSummarySchema(BaseSummarySchema):
    nearest_mibig: Series[str]

    @pa.check("nearest_mibig")
    def mibig_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)

    class Config:
        dtype = PydanticModel(SanntisSummaryRecord)
        coerce = True


class AntismashSummarySchema(BaseSummarySchema):
    label: Series[str]

    @pa.check("label")
    def class_names_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)

    class Config:
        dtype = PydanticModel(AntismashSummaryRecord)
        coerce = True


class KOSummarySchema(BaseSummarySchema):
    ko: Series[str]

    @pa.check("ko")
    def ko_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)

    class Config:
        dtype = PydanticModel(KOSummaryRecord)
        coerce = True


class PFAMSummarySchema(BaseSummarySchema):
    pfam: Series[str]

    @pa.check("pfam")
    def pfam_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)

    class Config:
        dtype = PydanticModel(PFAMSummaryRecord)
        coerce = True


class KEGGModulesSummarySchema(BaseSummarySchema):
    module_accession: Series[str]

    @pa.check("module_accession")
    def module_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)

    class Config:
        dtype = PydanticModel(KEGGModulesSummaryRecord)
        coerce = True


class BaseStudySummarySchema(BaseSummarySchema):
    """Base schema for study summary files with ERZ* columns and count checks."""

    @pa.check(regex=r"^ERZ\d+")
    def count_columns_are_non_negative(self, s: Series[int]) -> Series[bool]:
        return s >= 0

    class Config:
        strict = False  # allow extra ERZ* columns not declared above


class GOStudySummarySchema(BaseStudySummarySchema):
    GO: Series[str] = pa.Field(str_matches=r"^GO:\d{7}$")
    description: Series[str]
    category: Series[str]

    @pa.check("GO")
    def go_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)


class InterProStudySummarySchema(BaseStudySummarySchema):
    IPR: Series[str] = pa.Field(str_matches=r"^IPR\d{6}$")
    description: Series[str]

    @pa.check("IPR")
    def interpro_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)


class AntismashStudySummarySchema(BaseStudySummarySchema):
    label: Series[str]

    @pa.check("label")
    def class_names_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)


class SanntisStudySummarySchema(BaseStudySummarySchema):
    nearest_mibig: Series[str]

    @pa.check("nearest_mibig")
    def mibig_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)


class KOStudySummarySchema(BaseStudySummarySchema):
    KO: Series[str]

    @pa.check("KO")
    def ko_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)


class PFAMStudySummarySchema(BaseStudySummarySchema):
    PFAM: Series[str]

    @pa.check("PFAM")
    def pfam_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)


class KEGGModulesStudySummarySchema(BaseStudySummarySchema):
    module_accession: Series[str]

    @pa.check("module_accession")
    def module_ids_unique(self, series: Series[str]) -> Series[bool]:
        return self.is_unique(series)


class TaxonomyStudySummarySchema(BaseStudySummarySchema):
    pass


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
        if (
            rank != ""
            and rank.capitalize() != "Unclassified"
            and rank not in cls.valid_tax_ranks
        ):
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
