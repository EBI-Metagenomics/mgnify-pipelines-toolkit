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

DIRECTORY_SEQ_CAT = "sequence-categorisation"

SSU = "SSU_rRNA"
LSU = "LSU_rRNA"
Seq5S = "mtPerm-5S"
Seq5_8S = "5_8S_rRNA"

SSU_rRNA_archaea = "SSU_rRNA_archaea"
SSU_rRNA_bacteria = "SSU_rRNA_bacteria"
SSU_rRNA_eukarya = "SSU_rRNA_eukarya"
SSU_rRNA_microsporidia = "SSU_rRNA_microsporidia"

LSU_rRNA_archaea = "LSU_rRNA_archaea"
LSU_rRNA_bacteria = "LSU_rRNA_bacteria"
LSU_rRNA_eukarya = "LSU_rRNA_eukarya"

NON_CODING_RNA = [
    SSU_rRNA_archaea,
    SSU_rRNA_bacteria,
    SSU_rRNA_eukarya,
    SSU_rRNA_microsporidia,
    LSU_rRNA_archaea,
    LSU_rRNA_bacteria,
    LSU_rRNA_eukarya,
    Seq5S,
    Seq5_8S,
]

SSU_MODELS = [
    SSU_rRNA_archaea,
    SSU_rRNA_bacteria,
    SSU_rRNA_eukarya,
    SSU_rRNA_microsporidia,
]

LSU_MODELS = [LSU_rRNA_archaea, LSU_rRNA_bacteria, LSU_rRNA_eukarya]

RFAM_MODELS = {
    SSU_rRNA_archaea: "RF01959",
    SSU_rRNA_bacteria: "RF00177",
    SSU_rRNA_eukarya: "RF01960",
    SSU_rRNA_microsporidia: "RF02542",
    LSU_rRNA_archaea: "RF02540",
    LSU_rRNA_bacteria: "RF02541",
    LSU_rRNA_eukarya: "RF02543",
}

TRNA = [
    "Ala",
    "Gly",
    "Pro",
    "Thr",
    "Val",
    "Ser",
    "Arg",
    "Leu",
    "Phe",
    "Asn",
    "Lys",
    "Asp",
    "Glu",
    "His",
    "Gln",
    "Ile",
    "Tyr",
    "Cys",
    "Trp",
]
