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


import re
from enum import Enum
from typing import List, Optional


class FormatterMode(Enum):
    """
    Enum for FASTA header formatter modes.

    :ivar SIMPLE: Only sequence name (no metadata)
    :type SIMPLE: str
    :ivar GENERIC: All metadata joined by pipes
    :type GENERIC: str
    :ivar VIRIFY: Virify-specific metadata formatting
    :type VIRIFY: str
    """

    SIMPLE = "simple"
    GENERIC = "generic"
    VIRIFY = "virify"


class ParsedHeader:
    """
    Represents parsed FASTA header information.

    :ivar seq_name: Full sequence name
    :type seq_name: str
    :ivar short_name: Short name (first token before space)
    :type short_name: str
    :ivar metadata: Additional metadata extracted from header
    :type metadata: list of str
    :ivar viral_id: Viral identifier if present
    :type viral_id: str or None
    """

    def __init__(
        self,
        seq_name: str,
        short_name: str,
        metadata: Optional[List[str]] = None,
        viral_id: Optional[str] = None,
    ):
        """
        Initialize parsed header.

        :param seq_name: Full sequence name
        :type seq_name: str
        :param short_name: Short name (first token before space)
        :type short_name: str
        :param metadata: Additional metadata (default: empty list)
        :type metadata: list of str or None
        :param viral_id: Viral identifier if present
        :type viral_id: str or None
        """
        self.seq_name = seq_name
        self.short_name = short_name
        self.metadata = metadata or []
        self.viral_id = viral_id


def _extract_short_name(seq_name: str) -> str:
    """
    Extract short name (first token before space).

    :param seq_name: Full sequence name
    :type seq_name: str
    :returns: Short name (first whitespace-delimited token)
    :rtype: str
    """
    return seq_name.split()[0]


def _clean_header(header: str) -> str:
    """
    Clean header line (remove ">", newlines, whitespace).

    :param header: Raw header line
    :type header: str
    :returns: Cleaned header string
    :rtype: str
    """
    return header.replace(">", "").replace("\n", "").strip()


def parse_header(header: str) -> ParsedHeader:
    """
    Parse FASTA header to extract sequence and short names (base parser).

    :param header: Raw FASTA header
    :type header: str
    :returns: ParsedHeader with parsed components
    :rtype: ParsedHeader
    """
    clean = _clean_header(header)
    seq_name = clean.split("|")[0]
    short_name = _extract_short_name(seq_name)

    return ParsedHeader(
        seq_name=seq_name,
        short_name=short_name,
        metadata=[],
        viral_id=None,
    )


def parse_virify_header(header: str) -> ParsedHeader:
    """
    Parse FASTA header extracting Virify-specific metadata.

    Handles:
    - Viral identifiers with "|" separator
    - VirSorter metadata (phage-circular, prophage-<start>:<end>)

    :param header: Raw FASTA header
    :type header: str
    :returns: ParsedHeader with viral metadata components
    :rtype: ParsedHeader
    """
    clean = _clean_header(header)
    parts = clean.split("|")
    seq_name = parts[0]
    short_name = _extract_short_name(seq_name)

    metadata = []
    viral_identifier = None

    # Check for phage-circular in header
    if "phage-circular" in clean:
        metadata.append("phage-circular")

    # Check for viral identifier (second part if not prophage metadata)
    if len(parts) > 1:
        if not parts[1].startswith("prophage-") and parts[1] != "phage-circular":
            viral_identifier = parts[1]
            if viral_identifier not in metadata:
                metadata.append(viral_identifier)

    # Parse prophage metadata patterns
    for part in parts:
        match = re.search(r"prophage-\d+:\d+", part)
        if match and match[0] not in metadata:
            metadata.append(match[0])

    return ParsedHeader(
        seq_name=seq_name,
        short_name=short_name,
        metadata=metadata,
        viral_id=viral_identifier,
    )
