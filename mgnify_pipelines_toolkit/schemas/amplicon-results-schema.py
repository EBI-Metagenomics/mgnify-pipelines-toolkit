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


from pydantic import BaseModel, field_validator, model_validator
from pathlib import Path
from base_schemas import *


class QCFolderModel(DirectoryModel):
    path: RequiredDirectoryModel
    run_id: str

    @model_validator(mode="after")
    @classmethod
    def validate_qc_folder(cls, values):
        """Validates the QC folder structure and returns a dictionary of validated files."""
        required_filesuffixes = [
            "_seqfu.tsv"
        ]
        optional_filesuffixes = [
            ".merged.fastq.gz",
            ".fastp.fastq.gz",
            ".fastp.json",
            "_suffix_header_err.json",
            "_multiqc_report.html",
        ]

        DirectoryModel.validate_folder(
            run_id=values.run_id,
            path=values.path,
            required_filesuffixes=required_filesuffixes,
            optional_filesuffixes=optional_filesuffixes
        )


qc_folder = QCFolderModel(
    path=RequiredDirectoryModel(path=Path("/Users/kates/Desktop/EBI/MGnify/mgnify-pipelines-toolkit/tests/fixtures/study_summary_inputs/amplicon/ERR4334351/qc")),
    run_id="ERR4334351")
#validated_files = qc_folder.validate_qc_folder()



