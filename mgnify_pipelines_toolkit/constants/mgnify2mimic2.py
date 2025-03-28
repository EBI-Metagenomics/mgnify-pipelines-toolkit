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


PFAM_VERSIONS = [
    "32.0",
    "33.0",
    "33.1",
    "34.0",
    "35.0",
    "36.0",
    "37.0",
    "37.1",
    "37.2",
]

BIOMES_DICT = {
    "chicken-gut": "v1.0.1",
    "cow-rumen": "v1.0.1",
    "honeybee-gut": "v1.0.1",
    "human-gut": "v2.0.2",
    "human-oral": "v1.0.1",
    "human-vaginal": "v1.0",
    "marine": "v2.0",
    "mouse-gut": "v1.0",
    "non-model-fish-gut": "v2.0",
    "pig-gut": "v1.0",
    "sheep-rumen": "v1.0",
    "zebrafish-fecal": "v1.0",
}


MAX_RETRIES = 3
RETRY_DELAY = 5

BASE_STUDY_API = "https://www.ebi.ac.uk/metagenomics/api/v1/studies/{mgys_id}/analyses"
BASE_MGYA_URL = "https://www.ebi.ac.uk/metagenomics/api/v1/analyses/{mgya_id}"
BASE_MGYA_DOWNLOADS_URL = (
    "https://www.ebi.ac.uk/metagenomics/api/v1/analyses/{mgya_id}/downloads"
)
