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

from pathlib import Path
import pytest

from mgnify_pipelines_toolkit.analysis.amplicon.amplicon_utils import get_read_count

# This is not very nice, we will improve later
fixtures_directory = Path(__file__).resolve().parent.parent.parent.parent / "fixtures"


def test_get_read_count_from_fastq():
    read_count = get_read_count(f"{fixtures_directory}/sequences/ERR4674038_1.fastq.gz")

    assert read_count == 10


def test_get_read_count_from_fasta():
    read_count = get_read_count(f"{fixtures_directory}/sequences/rpxx.fasta", "fasta")

    assert read_count == 2

    read_count = get_read_count(
        f"{fixtures_directory}/sequences/rpxx.fasta.gz", "fasta"
    )

    assert read_count == 2


def test_get_read_count_should_break():
    with pytest.raises(RuntimeError):
        get_read_count(f"{fixtures_directory}/sequences/rpxx.fasta", "fastq")


def test_get_read_count_should_break_if_no_records(tmp_path):
    empty_file = tmp_path / "empty.fasta"
    empty_file.touch()  # This creates an empty file

    with pytest.raises(RuntimeError, match="is not plain or gzip compressed fasta"):
        get_read_count(str(empty_file), "fasta")
