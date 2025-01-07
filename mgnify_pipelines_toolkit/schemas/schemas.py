#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024 EMBL - European Bioinformatics Institute
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

import re

from enum import Enum
import pandas as pd
import pandera as pa
import pandera.extensions as extensions

from pydantic import (
    Field,
    BaseModel,
    field_validator,
    RootModel,
)
from pandera.engines.pandas_engine import PydanticModel


@extensions.register_check_method(statistics=["short_tax_ranks"])
def is_valid_tax_hierarchy(pandas_obj, *, short_tax_ranks):

    bool_list = []
    short_tax_ranks.append(
        "Unclassified"
    )  # This is the only non-hierarchical value we can still accept

    for taxa in pandas_obj:
        taxa_lst = [rank.split("__")[0] for rank in taxa.split(";")]
        if len(set.intersection(set(taxa_lst), set(short_tax_ranks))) == len(taxa_lst):
            bool_list.append(True)
        else:
            bool_list.append(False)

    return pd.Series(bool_list)


def generate_dynamic_tax_df_schema(
    run_acc: str, short_tax_ranks: list
) -> pa.DataFrameSchema:

    tax_schema = pa.DataFrameSchema(
        {run_acc: pa.Column(int, checks=pa.Check.ge(0))},
        index=pa.Index(
            str,
            unique=True,
            checks=[pa.Check.is_valid_tax_hierarchy(short_tax_ranks=short_tax_ranks)],
        ),
        strict=True,
        coerce=True,
    )

    return tax_schema


class ASVResultTypes(str, Enum):
    all_results = "all_results"
    no_asvs = "no_asvs"


class INSDCRunAccession(RootModel):
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
        run_accession_regex = "(E|D|S)RR[0-9]{6,}"
        regex_res = re.match(run_accession_regex, run)

        if regex_res is None:
            raise ValueError(
                f"Accession `{run}` does not fit INSDC format [ERR*,SRR*,DRR*]."
            )

        return run


# This is one row
class AmpliconPassedRunsRecord(BaseModel):
    run: INSDCRunAccession
    status: ASVResultTypes


class AmpliconNonINSDCSPassedRunsRecord(BaseModel):
    run: str
    status: ASVResultTypes


# This is the schema for the whole DF
class AmpliconPassedRunsSchema(pa.DataFrameModel):
    """Pandera schema using the pydantic model."""

    class Config:
        """Config with dataframe-level data type."""

        dtype = PydanticModel(AmpliconPassedRunsRecord)
        coerce = True


class AmpliconNonINSDCPassedRunsSchema(pa.DataFrameModel):
    """Pandera schema using the pydantic model."""

    class Config:
        """Config with dataframe-level data type."""

        dtype = PydanticModel(AmpliconNonINSDCSPassedRunsRecord)
        coerce = True
