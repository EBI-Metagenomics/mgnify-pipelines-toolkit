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

import argparse
import sys
import re
import csv
import hashlib
import gzip
from Bio import SeqIO
from mgnify_pipelines_toolkit.constants.regex_fasta_header import _FORMAT_REGEX_MAP

def is_gzipped(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def guess_header_format(header):
    matches = [(format, re.search(regex, header)) for format, regex in _FORMAT_REGEX_MAP.items()]
    guesses = [(format, match.groups()) for format, match in matches if match is not None]

    if not guesses:
        raise ValueError("Header format could not be determined")

    guessed_format, _ = max(guesses, key=lambda g: g[1])

    return guessed_format

def md5_hash(s):
    md5 = hashlib.md5()
    md5.update(s.encode('utf-8'))

    return md5.hexdigest()

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Path to (gzipped) Fasta file")
    parser.add_argument("-o", "--output", type=str, help="Output path")
    parser.add_argument("-f", "--format", type=str, choices=["auto", "uniprotkb", "rpxx"], default="auto", help="Format of the input Fasta header")
    parser.add_argument("-d", "--delimiter", type=str, default="\t", help="Output column delimiter")
    parser.add_argument("--with-hash", action="store_true", help="Add a MD5 hash of the sequence to the output")
    parser.add_argument("--no-header", action="store_true", help="Do not add header to output file")

    args = parser.parse_args()
    
    _PATH = args.input
    _OUTPUT = args.output
    _FORMAT = args.format
    _DELIMITER = args.delimiter
    _WITH_HASH = args.with_hash
    _NO_HEADER = args.no_header

    return _PATH, _OUTPUT, _FORMAT, _DELIMITER, _WITH_HASH, _NO_HEADER


def main():
    _PATH, _OUTPUT, _FORMAT, _DELIMITER, _WITH_HASH, _NO_HEADER = parse_args()

    if is_gzipped(_PATH):
        input_fh = gzip.open(_PATH, mode="rt")
    else:
        input_fh = open(_PATH, mode="rt")

    if _OUTPUT is None:
        output_fh = sys.stdout
    else:
        output_fh = open(_PATH, mode="w", newline="")

    if _FORMAT != "auto":
        header_regex = re.compile(_FORMAT_REGEX_MAP[_FORMAT])
    else:
        header, _ = next(SeqIO.FastaIO.SimpleFastaParser(input_fh))
        format = guess_header_format(header)
        header_regex = re.compile(_FORMAT_REGEX_MAP[format])
        input_fh.seek(0)

    fieldnames = list(header_regex.groupindex.keys())

    if _WITH_HASH:
        fieldnames += ["sequence_hash"]

    fieldnames += ["sequence"]

    csv_writer = csv.DictWriter(output_fh, fieldnames=fieldnames, delimiter=_DELIMITER, extrasaction="ignore")

    if not _NO_HEADER:
        csv_writer.writeheader()

    for header, sequence in SeqIO.FastaIO.SimpleFastaParser(input_fh):
        header_match = header_regex.match(header)

        row = {"sequence": sequence}

        if header_match:
            row.update(header_match.groupdict())

        if _WITH_HASH:
            row["sequence_hash"] = md5_hash(sequence)

        csv_writer.writerow(row)

    input_fh.close()
    output_fh.close()


if __name__ == "__main__":
    main()
