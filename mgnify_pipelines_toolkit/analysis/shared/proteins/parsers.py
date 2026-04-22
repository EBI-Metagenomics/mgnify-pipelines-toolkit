#!/usr/bin/env python3

# Copyright 2026Give examples  EMBL - European Bioinformatics Institute
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
from dataclasses import dataclass

# Prodigal: protein id format:
#   <contig_id>_<protein_index>
# Example:
#   ERZ23299386_269_1
PRODIGAL_PROTEIN_ID_REGEX = re.compile(r"^(?P<contig_id>.+)_(?P<idx>\d+)$")

# FragGeneScan protein id format:
#   <contig_id>_<start>_<end>_<strand>
# Example:
#   ERZ23299386_9851_98_1372_+
FGS_PROTEIN_ID_REGEX = re.compile(r"^(?P<gene_id>(?P<contig_id>.+)_(?P<start>\d+)_(?P<end>\d+)_(?P<strand>[+-]))$")

# Prodigal full header line split pattern:
#   <gene_id> # <start> # <end> # <strand> ...
PRODIGAL_SPLIT_RE = re.compile(r"\s+#\s+")


@dataclass(frozen=True)
class ProdigalCall:
    gene_id: str
    contig_id: str
    start_1: int  # 1-based inclusive
    end_1: int  # 1-based inclusive
    strand: int  # 1 or -1


@dataclass(frozen=True)
class FragGeneScanCall:
    gene_id: str
    contig_id: str
    start_1: int  # 1-based inclusive
    end_1: int  # 1-based inclusive
    strand: str  # + or -


def parse_prodigal_header(record_id: str, record_description: str) -> ProdigalCall:
    """Parse Prodigal FASTA header into coordinates and strand.

    Expect Prodigal style header:
        <gene_id> # <start> # <end> # <strand> ...
    """
    line = record_description or record_id
    parts = PRODIGAL_SPLIT_RE.split(line)
    if len(parts) < 4:
        raise ValueError(f"FAA header not parseable as Prodigal: {line!r}")

    gene_id_full = parts[0].strip().split()[0]
    start = int(parts[1].strip())
    end = int(parts[2].strip())
    strand = int(parts[3].strip())

    if strand not in (1, -1):
        raise ValueError(f"FAA header has invalid strand (expected 1 or -1): {line!r}")
    if start < 1 or end < 1:
        raise ValueError(f"FAA header has invalid coordinates (<1): {line!r}")
    if start > end:
        start, end = end, start

    m = PRODIGAL_PROTEIN_ID_REGEX.match(gene_id_full)
    contig_id = m.group("contig_id") if m else gene_id_full

    return ProdigalCall(
        gene_id=gene_id_full,
        contig_id=contig_id,
        start_1=start,
        end_1=end,
        strand=strand,
    )


def parse_fraggenescan_header(record_id: str) -> FragGeneScanCall:
    """Parse FragGeneScan-style protein ID into coordinates and strand.

    Expected token format:
        <contig_id>_<start>_<end>_<strand>

    Example:
        ERZ23299386_9851_98_1372_+
    """
    token = record_id.strip().split()[0]
    m = FGS_PROTEIN_ID_REGEX.match(token)
    if not m:
        raise ValueError(f"FAA header not parseable as FragGeneScan: {record_id!r}")

    start = int(m.group("start"))
    end = int(m.group("end"))
    strand = m.group("strand")

    if start < 1 or end < 1:
        raise ValueError(f"FAA header has invalid coordinates (<1): {record_id!r}")
    if start > end:
        start, end = end, start

    return FragGeneScanCall(
        gene_id=m.group("gene_id"),
        contig_id=m.group("contig_id"),
        start_1=start,
        end_1=end,
        strand=strand,
    )
