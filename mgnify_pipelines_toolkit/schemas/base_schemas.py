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


class FileModel(BaseModel):
    path: Path
    # TODO: check path is not a directory


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


class DirectoryModel(BaseModel):
    path: Path

    @classmethod
    def validate_file(cls, run_id, filename, filename_pattern: str) -> bool:
        pattern = f"{run_id}{filename_pattern}"
        return filename == pattern

    @classmethod
    def validate_folder(cls, run_id, path, required_filesuffixes, optional_filesuffixes):
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


class RequiredDirectoryModel(DirectoryModel):
    path: Path

    @field_validator("path")
    def directory_must_exist(cls, directory_name):
        if not directory_name.is_dir():
            raise ValueError(f'Directory does not exist: {directory_name}')
        return directory_name
