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


from pydantic import BaseModel, validator
from pathlib import Path
import re

class FileModel(BaseModel):
    path: Path

    @validator('path')
    def file_must_exist(cls, v):
        if not v.is_file():
            raise ValueError(f'File does not exist: {v}')
        return v

class DirectoryModel(BaseModel):
    path: Path
    files: list[FileModel] = []
    subdirectories: list['DirectoryModel'] = []

    @validator('path')
    def directory_must_exist(cls, v):
        if not v.is_dir():
            raise ValueError(f'Directory does not exist: {v}')
        return v

    @validator('files', each_item=True)
    def files_must_exist(cls, v):
        if not v.path.is_file():
            raise ValueError(f'File does not exist: {v.path}')
        return v

    @validator('subdirectories', each_item=True)
    def subdirectories_must_exist(cls, v):
        if not v.path.is_dir():
            raise ValueError(f'Subdirectory does not exist: {v.path}')
        return v

class QCFolderModel(BaseModel):
    path: Path
    run_id: str

    @validator('path')
    def directory_must_exist(cls, v):
        if not v.is_dir():
            raise ValueError(f"Directory does not exist: {v}")
        return v

    def validate_file(self, filename_pattern: str, required: bool = False) -> FileModel:
        """Validates if a file exists in the folder following the given pattern."""
        expected_file = self.path / filename_pattern.format(run_id=self.run_id)
        return FileModel(path=expected_file, required=required)

    def validate_qc_folder(self):
        """Validates the QC folder structure and returns a dictionary of validated files."""
        required_files = {f"{self.run_id}_seqfu.tsv"}
        optional_files = {
            f"{self.run_id}.merged.fastq.gz",
            f"{self.run_id}.fastp.fastq.gz",
            f"{self.run_id}.fastp.json",
            f"{self.run_id}_suffix_header_err.json",
            f"{self.run_id}_multiqc_report.html",
        }

        for filename in self.path.iterdir():
            if filename in required_files:
                self.validate_file(filename, required=True)
            elif filename in optional_files:
                self.validate_file(filename)
            else:
                raise ValueError(f"Unexpected file {filename}")


qc_folder = QCFolderModel(path=Path("qc"), run_id="ERR")
validated_files = qc_folder.validate_qc_folder()


