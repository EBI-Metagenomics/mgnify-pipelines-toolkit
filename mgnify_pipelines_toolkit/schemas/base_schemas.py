#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Union, Optional

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


class FileModel(BaseModel):
    path: Path


class AccessionedFileModel(FileModel):
    @classmethod
    def validate_accession_file(cls, path, run_id) -> bool:
        """
        Does file start with run_id ?
        """
        return path.name.startswith(run_id)


class SuffixedFileModel(FileModel):
    @classmethod
    def validate_suffix_file(cls, path, pattern) -> bool:
        """
        Does filename end with suffix ?
        """
        return path.name.endswith(pattern)


class RequiredFileModel(FileModel):
    @field_validator("path", mode="after")
    @classmethod
    def file_must_exist(cls, filename) -> Path:
        if not filename.is_file():
            raise ValueError(f'File does not exist: {filename}')
        return filename


class OptionalFileModel(FileModel):

    @classmethod
    def file_can_exist(cls, filename) -> bool:
        return filename.is_file()


class ResultFile(AccessionedFileModel, SuffixedFileModel):
    pass


# ----------- files in directory -----------


class DirectoryModel(BaseModel):
    path: Path

    @classmethod
    def check_reqired(cls, path, suffixes, run_id):
        """
        suffixes = [ Single[], OneOf[], ...]
        """
        summ = 0
        required_files = []
        for filename in path.path.iterdir():
            for suffix in suffixes:
                res = [ ResultFile.validate_suffix_file(path=filename, pattern=suff) and ResultFile.validate_accession_file(path=filename, run_id=run_id) for suff in suffix ]
                summ += sum(res)
                if sum(res) == 1:
                    required_files.append(filename)
        if summ != len(suffixes):
            raise ValueError("Error in required step")
        return len(suffixes), required_files

    @classmethod
    def check_optional(cls, path, suffixes, run_id):
        """
        suffixes = [ Single[], OneOf[], Single[], ...]
        """
        summ = 0
        optional_files = []
        for filename in path.path.iterdir():
            for suffix in suffixes:
                res = [
                    ResultFile.validate_suffix_file(path=filename, pattern=suff) and ResultFile.validate_accession_file(
                        path=filename, run_id=run_id) for suff in suffix]
                summ += sum(res)
                if sum(res) == 1:
                    optional_files.append(filename)

        return summ, optional_files

    @classmethod
    def validate_folder(cls, path, required_suffixes, optional_suffixes, run_id):
        """
        AMPLICON_QC_REQ = [Single("*seq_fu.tsv"])
        AMPLICON_QC_OPT = [[Single("json"), Single("html"), OneOf(["merged", "fastp"])], Single("err.json")]
        """
        result_required, required_files = cls.check_reqired(path=path, suffixes=required_suffixes, run_id=run_id)
        result_optional, optional_files = cls.check_optional(path=path, suffixes=optional_suffixes, run_id=run_id)
        if result_required + result_optional < len(list(path.path.iterdir())):
            print([ parsed_file for parsed_file in path.path.iterdir() if parsed_file not in required_files+optional_files] )
            raise ValueError("Unexpected file presented")
        return result_optional

"""
        required_files, optional_files = [], []
        for filename in path.path.iterdir():

            for pattern in required_filesuffixes:
                if RequiredDirectoryModel.validate_file(
                        run_id=run_id,
                        filename=RequiredFileModel(path=filename).path.name,
                        filename_pattern=pattern
                ):
                    required_files.append(filename)

            if filename not in required_files:
                for pattern in optional_filesuffixes:
                    if OptionalFileModel.file_can_exist(filename):
                        if DirectoryModel.validate_file(
                                run_id=run_id,
                                filename=FileModel(path=filename).path.name,
                                filename_pattern=pattern):
                            optional_files.append(filename)
        other_files = set(path.path.iterdir()).difference(set(optional_files + required_files))

        for filename in other_files:
            if filename not in required_files and filename not in optional_files:
                raise ValueError(f"Unexpected file {filename}")
"""

class RequiredDirectoryModel(DirectoryModel):
    path: Path

    @field_validator("path")
    def directory_must_exist(cls, directory_name):
        if not directory_name.is_dir():
            raise ValueError(f'Directory does not exist: {directory_name}')
        return directory_name
