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

from unittest.mock import patch

import pytest

from mgnify_pipelines_toolkit.analysis.amplicon.amplicon_utils import get_read_count


@pytest.mark.parametrize(
    "read_count_output, expected_result",
    [
        ("123\n", 123),
        ("5435", 5435),
        ("1000000000000", 1000000000000),
        ("999999999999999999999999999999", 999999999999999999999999999999),
        ("\n999999999999999999999999999999", 999999999999999999999999999999),
        ("\n999999999999999999999999999999\n", 999999999999999999999999999999),
        ("     999999999999999999999999999999      ", 999999999999999999999999999999),
    ],
)
@patch("mgnify_pipelines_toolkit.analysis.amplicon.amplicon_utils.subprocess.Popen")
def test_get_read_count(mock_popen, read_count_output, expected_result):
    mock_process = mock_popen.return_value
    mock_process.communicate.return_value = (read_count_output, expected_result)

    print(expected_result)

    read_count = get_read_count("/path/to/read_file")

    assert read_count == expected_result


@pytest.mark.parametrize(
    "read_count_output",
    [
        "9999.9999",
        "abc\n",
        "",
        None,
        "8765\t888",
        "888 888",
    ],
)
@patch("mgnify_pipelines_toolkit.analysis.amplicon.amplicon_utils.subprocess.Popen")
def test_get_read_count_error(mock_popen, read_count_output):
    mock_process = mock_popen.return_value
    mock_process.communicate.return_value = (read_count_output, None)
    # Simulating non-digit output
    with pytest.raises(SystemExit) as exc_info:
        get_read_count("/path/to/read_file")

    assert exc_info.type == SystemExit
    assert exc_info.value.code == 1
