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

import gzip


def open_file(file_path: str, mode: str = "r"):
    r"""
    Open a file, handling both compressed and uncompressed formats.

    Automatically detects gzip format by file extension (.gz) and opens
    accordingly. For text mode on gzip files, automatically adds 't' flag.

    :param file_path: Path to the file
    :type file_path: str
    :param mode: File opening mode (default: "r")
    :type mode: str
    :returns: File handle (gzip.GzipFile for .gz files, standard file handle otherwise)
    :rtype: file object

    Example usage::

        with open_file("data.fasta.gz", "r") as f:
            # Handle as text file, gzip transparent
            for line in f:
                print(line.strip())

        with open_file("data.fasta", "r") as f:
            # Standard file handling
            for line in f:
                print(line.strip())
    """
    if file_path.endswith(".gz"):
        # Ensure text mode is explicit for gzip
        if "b" not in mode and "t" not in mode:
            mode = mode + "t"
        return gzip.open(file_path, mode)
    return open(file_path, mode)
