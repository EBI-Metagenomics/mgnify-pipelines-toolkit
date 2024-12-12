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

import pandas as pd
import pandera as pa
from pandera.typing import DataFrame, Series
import pandera.extensions as extensions

import pydantic


class SuccessfulRunsSchema(pa.DataFrameModel):
    run: Series[str] = pa.Field(unique=True)
    status: Series[str]

    @pa.check(
        "run",
        name="run_validity_check",
        raise_warning=True,
        error="One or more run accessions do not fit the INSDC format [ERR*,SRR*,DRR*]. This is only a warning, not an error.",
    )
    def run_validity_check(self, run: Series[str]) -> Series[bool]:
        # This will only produce a WARNING, not an ERROR. This is to allow flexibility of running this on non-ENA/INSDC data
        run_accession_regex = "(E|D|S)RR[0-9]{6,}"
        return run.str.contains(run_accession_regex)

    @pa.check(
        "status",
        name="status_vality_check",
        error='The status column can only have values ["all_results", "no_asvs"].',
    )
    def status_vality_check(self, status: Series[str]) -> Series[bool]:
        possible_statuses = ["all_results", "no_asvs"]
        return status.isin(possible_statuses)


class PydanticModel(pydantic.BaseModel):
    df: DataFrame[pa.DataFrameModel]


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
