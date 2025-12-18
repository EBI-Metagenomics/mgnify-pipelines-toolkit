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


from .parsers import FormatterMode, ParsedHeader


def format_header(new_name: str, parsed: ParsedHeader, mode: str = "simple") -> str:
    """
    Format a FASTA header with new name and optional metadata.

    :param new_name: New sequence name
    :type new_name: str
    :param parsed: Parsed header information
    :type parsed: ParsedHeader
    :param mode: Formatting mode ('simple', 'generic', or 'virify')
    :type mode: str
    :returns: Formatted FASTA header line (without newline)
    :rtype: str
    """
    # Convert string mode to enum if needed
    if isinstance(mode, str):
        mode = FormatterMode(mode)

    if mode == FormatterMode.SIMPLE or not parsed.metadata:
        return f">{new_name}"

    if mode == FormatterMode.GENERIC:
        return f">{new_name}|{'|'.join(parsed.metadata)}"

    if mode == FormatterMode.VIRIFY:
        # Start with viral ID if present
        if parsed.viral_id:
            header_parts = [parsed.viral_id]
            # Add other metadata (excluding viral_id)
            other_metadata = [m for m in parsed.metadata if m != parsed.viral_id]
            if other_metadata:
                header_parts.extend(other_metadata)
            return f">{new_name}|{'|'.join(header_parts)}"
        else:
            # No viral ID, just join all metadata
            return f">{new_name}|{'|'.join(parsed.metadata)}"

    # Default to simple if unknown mode
    return f">{new_name}"
