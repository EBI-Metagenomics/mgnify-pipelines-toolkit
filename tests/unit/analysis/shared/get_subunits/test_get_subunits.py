#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2024-2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from argparse import Namespace
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mgnify_pipelines_toolkit.analysis.shared.get_subunits import main


@pytest.mark.parametrize("separate_models", [False, True])
@pytest.mark.parametrize("sequence_length", [40, 12000])
@patch("argparse.ArgumentParser.parse_args")
def test_shared_output_preserves_records(
    mock_args: MagicMock,
    tmp_path: Path,
    separate_models: bool,
    sequence_length: int) -> None:
    """Preserve interleaved model records across output buffer boundaries.

    :param mock_args: Mock for parsed command arguments.
    :param tmp_path: Temporary input and output directory.
    :param separate_models: Whether to write separate ribosomal model files.
    :param sequence_length: Length of the first sequence.
    """
    models = ["Hatchet", "Cobalamin", "Hatchet", "Bacteria_large_SRP"]
    records = [
        SeqRecord(Seq("ACGT" * (sequence_length // 4 + index)), id=f"ERZ123_{index}-{model}/1-100", description="")
        for index, model in enumerate(models)
    ]
    source = tmp_path / "coords.fasta"

    SeqIO.write(records, source, "fasta")

    mock_args.return_value = Namespace(
        input=str(source),
        prefix=None,
        name="sample",
        separate_subunits_by_models=separate_models
    )

    output_dir = tmp_path / "sequence-categorisation"

    with patch("mgnify_pipelines_toolkit.analysis.shared.get_subunits.DIRECTORY_SEQ_CAT", str(output_dir)):
        main()

    output = output_dir / "sample_other_ncRNA.fasta"
    assert output.read_text() == source.read_text()
