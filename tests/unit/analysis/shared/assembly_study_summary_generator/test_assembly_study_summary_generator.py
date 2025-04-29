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

import pytest
import gzip
import hashlib
import os
from typing import Dict, List, Any

from click.testing import CliRunner

from mgnify_pipelines_toolkit.analysis.assembly.study_summary_generator import (
    summarise_analyses,
)


def tsv_string(rows: List[List[Any]]) -> str:
    """
    Generate a TSV string from the provided rows.

    :param rows: A list of lists containing the data for the TSV.
    :returns: A tab-separated values (TSV) string.
    :rtype: str
    """
    output = ""
    for row in rows:
        output += "\t".join(str(column) for column in row) + "\n"
    return output


@pytest.fixture
def go_summary_tsv_rows_per_accession() -> Dict:
    """
    Fixture that provides a TSV string representation of GO (Gene Ontology) summary data.

    The TSV string includes the following columns:

    - **go**: The Gene Ontology identifier (e.g., GO:0003824).
    - **term**: The name of the GO term (e.g., "catalytic activity").
    - **category**: The GO category (e.g., "molecular_function", "biological_process").
    - **count**: The number of occurrences of the GO term.

    The data is organized as a dictionary where keys are accession IDs (e.g., ERZ1049440)
    and values are lists of lists. The first list in each value is the header,
    and subsequent lists are rows of data.

    :returns: A dictionary mapping accession IDs to lists of GO summary data rows.
    :rtype: dict
    """
    rows_per_accession = {}

    header = ["go", "term", "category", "count"]

    rows_per_accession["ERZ1049440"] = [
        header,
        ["GO:0003824", "catalytic activity", "molecular_function", 2145],
        ["GO:0003677", "DNA binding", "molecular_function", 6125],
        ["GO:0055085", "transmembrane transport", "biological_process", 144],
        ["GO:0016491", "oxidoreductase activity", "molecular_function", 1513],
    ]
    rows_per_accession["ERZ1049443"] = [
        header,
        ["GO:0003824", "catalytic activity", "molecular_function", 18626],
        ["GO:0003677", "DNA binding", "molecular_function", 16417],
        ["GO:0055085", "transmembrane transport", "biological_process", 13926],
        ["GO:0005515", "protein binding", "molecular_function", 11917],
        ["GO:0016491", "oxidoreductase activity", "molecular_function", 16064],
    ]

    rows_per_accession["ERZ1049444"] = [
        header,
        ["GO:0003824", "catalytic activity", "molecular_function", 18626],
        ["GO:0003677", "DNA binding", "molecular_function", 16417],
        ["GO:0016491", "oxidoreductase activity", "molecular_function", 16064],
        ["GO:0055085", "transmembrane transport", "biological_process", 13926],
        ["GO:0005515", "protein binding", "molecular_function", 11917],
    ]

    rows_per_accession["ERZ1049445"] = [
        header,
        ["GO:0003824", "catalytic activity", "molecular_function", 2145],
        ["GO:0003677", "DNA binding", "molecular_function", 6125],
        ["GO:0016491", "oxidoreductase activity", "molecular_function", 1513],
        ["GO:0055085", "transmembrane transport", "biological_process", 144],
        ["GO:0005515", "protein binding", "molecular_function", 0],
    ]
    return rows_per_accession


@pytest.fixture
def interpro_summary_tsv_rows_per_accession() -> Dict:
    """
    Fixture that provides a TSV string representation of InterPro summary data.

    The TSV string includes the following columns:

    - **count**: The number of occurrences of the InterPro domain.
    - **interpro_accession**: The InterPro domain identifier (e.g., IPR027417).
    - **description**: A description of the InterPro domain.

    The data is organized as a dictionary where keys are accession IDs (e.g., ERZ1049440)
    and values are lists of lists. The first list in each value is the header,
    and subsequent lists are rows of data.

    :returns: A dictionary mapping accession IDs to lists of InterPro summary data rows.
    :rtype: dict
    """
    rows_per_accession = {}

    header = ["count", "interpro_accession", "description"]

    rows_per_accession["ERZ1049440"] = [
        header,
        [3241, "IPR027417", "P-loop containing nucleoside triphosphate hydrolase"],
        [3245, "IPR002347", "Short-chain dehydrogenase/reductase SDR"],
        [12565, "IPR036291", "NAD(P)-binding domain superfamily"],
        [13412, "IPR036188", "FAD/NAD(P)-binding domain superfamily"],
        [
            5362,
            "IPR029063",
            "S-adenosyl-L-methionine-dependent methyltransferase superfamily",
        ],
        [536, "IPR014729", "Rossmann-like alpha/beta/alpha sandwich fold"],
    ]
    rows_per_accession["ERZ1049443"] = [
        header,
        [60350, "IPR027417", "P-loop containing nucleoside triphosphate hydrolase"],
        [23116, "IPR002347", "Short-chain dehydrogenase/reductase SDR"],
        [16503, "IPR036291", "NAD(P)-binding domain superfamily"],
        [14221, "IPR036188", "FAD/NAD(P)-binding domain superfamily"],
        [
            13933,
            "IPR029063",
            "S-adenosyl-L-methionine-dependent methyltransferase superfamily",
        ],
        [6414, "IPR014729", "Rossmann-like alpha/beta/alpha sandwich fold"],
    ]

    rows_per_accession["ERZ1049444"] = [
        header,
        [6532, "IPR027417", "P-loop containing nucleoside triphosphate hydrolase"],
        [1356, "IPR002347", "Short-chain dehydrogenase/reductase SDR"],
        [9253, "IPR036291", "NAD(P)-binding domain superfamily"],
        [131, "IPR019734", "Tetratricopeptide repeat"],
        [1354, "IPR036188", "FAD/NAD(P)-binding domain superfamily"],
        [
            36331,
            "IPR029063",
            "S-adenosyl-L-methionine-dependent methyltransferase superfamily",
        ],
    ]

    rows_per_accession["ERZ1049445"] = [
        header,
        [3241, "IPR027417", "P-loop containing nucleoside triphosphate hydrolase"],
        [3245, "IPR002347", "Short-chain dehydrogenase/reductase SDR"],
        [12565, "IPR036291", "NAD(P)-binding domain superfamily"],
        [1543, "IPR019734", "Tetratricopeptide repeat"],
        [13412, "IPR036188", "FAD/NAD(P)-binding domain superfamily"],
        [
            5362,
            "IPR029063",
            "S-adenosyl-L-methionine-dependent methyltransferase superfamily",
        ],
    ]
    return rows_per_accession


@pytest.fixture
def ko_summary_summary_tsv_rows_per_accession() -> Dict:
    """
    Fixture that provides a TSV string representation of KO (KEGG Orthology) summary data.

    The TSV string includes the following columns:

    - **count**: The number of occurrences of the KO term.
    - **ko**: The KEGG Orthology identifier (e.g., K03807).
    - **description**: A description of the KO term.

    The data is organized as a dictionary where keys are accession IDs (e.g., ERZ1049440)
    and values are lists of lists. The first list in each value is the header,
    and subsequent lists are rows of data.

    :returns: A dictionary mapping accession IDs to lists of KO summary data rows.
    :rtype: dict
    """
    rows_per_accession = {}
    header = ["count", "ko", "description"]
    rows_per_accession["ERZ1049440"] = [
        header,
        [441, "K03807", "AmpE protein"],
        [399, "K10093", "galectin-9"],
        [480, "K09539", "DnaJ homolog subfamily C member 19"],
        [219, "K14642", "golgi apyrase [EC:3.6.1.5]"],
        [390, "K24212", "complement C1q and tumor necrosis factor-related protein 5"],
        [
            219,
            "K13020",
            "UDP-N-acetyl-2-amino-2-deoxyglucuronate dehydrogenase [EC:1.1.1.335]",
        ],
        [392, "K10893", "fanconi anemia group F protein"],
        [222, "K24258", "UDP-N-acetyl-L-fucosamine synthase [EC:5.1.3.28]"],
    ]
    rows_per_accession["ERZ1049443"] = [
        header,
        [456, "K10093", "galectin-9"],
        [215, "K09539", "DnaJ homolog subfamily C member 19"],
        [45, "K14642", "golgi apyrase [EC:3.6.1.5]"],
        [278, "K24212", "complement C1q and tumor necrosis factor-related protein 5"],
        [
            434,
            "K13020",
            "UDP-N-acetyl-2-amino-2-deoxyglucuronate dehydrogenase [EC:1.1.1.335]",
        ],
        [425, "K10893", "fanconi anemia group F protein"],
        [74, "K24258", "UDP-N-acetyl-L-fucosamine synthase [EC:5.1.3.28]"],
        [274, "K23160", "GlcNAc transferase"],
    ]
    rows_per_accession["ERZ1049444"] = [
        header,
        [409, "K03807", "AmpE protein"],
        [129, "K10093", "galectin-9"],
        [381, "K09539", "DnaJ homolog subfamily C member 19"],
        [189, "K24212", "complement C1q and tumor necrosis factor-related protein 5"],
        [
            323,
            "K13020",
            "UDP-N-acetyl-2-amino-2-deoxyglucuronate dehydrogenase [EC:1.1.1.335]",
        ],
        [247, "K10893", "fanconi anemia group F protein"],
        [438, "K10641", "E3 ubiquitin-protein ligase LRSAM1 [EC:2.3.2.27]"],
        [334, "K23160", "GlcNAc transferase"],
    ]
    rows_per_accession["ERZ1049445"] = [
        header,
        [90, "K10093", "galectin-9"],
        [117, "K14642", "golgi apyrase [EC:3.6.1.5]"],
        [423, "K24212", "complement C1q and tumor necrosis factor-related protein 5"],
        [
            42,
            "K13020",
            "UDP-N-acetyl-2-amino-2-deoxyglucuronate dehydrogenase [EC:1.1.1.335]",
        ],
        [423, "K10893", "fanconi anemia group F protein"],
        [335, "K10641", "E3 ubiquitin-protein ligase LRSAM1 [EC:2.3.2.27]"],
        [144, "K24258", "UDP-N-acetyl-L-fucosamine synthase [EC:5.1.3.28]"],
        [363, "K23160", "GlcNAc transferase"],
    ]
    return rows_per_accession


@pytest.fixture
def pfam_summary_summary_tsv_rows_per_accession() -> Dict:
    """
    Fixture that provides a TSV string representation of Pfam summary data.

    The TSV string includes the following columns:

    - **count**: The number of occurrences of the Pfam domain.
    - **pfam**: The Pfam domain identifier (e.g., PF02265).
    - **description**: A description of the Pfam domain.

    The data is organized as a dictionary where keys are accession IDs (e.g., ERZ1049440)
    and values are lists of lists. The first list in each value is the header,
    and subsequent lists are rows of data.

    :returns: A dictionary mapping accession IDs to lists of Pfam summary data rows.
    :rtype: dict
    """
    rows_per_accession = {}
    header = ["count", "pfam", "description"]
    rows_per_accession["ERZ1049440"] = [
        header,
        [230, "PF02265", "S1/P1 Nuclease"],
        [45, "PF02113", "D-Ala-D-Ala carboxypeptidase 3 (S13) family"],
        [2, "PF13558", "SbcC/RAD50-like, Walker B motif"],
        [34, "PF11655", "Protein of unknown function (DUF2589)"],
        [490, "PF06039", "Malate:quinone oxidoreductase (Mqo)"],
        [200, "PF07987", "Domain of unkown function (DUF1775)"],
        [108, "PF10726", "Protein of function (DUF2518)"],
        [468, "PF24718", "HTH-like domain"],
    ]
    rows_per_accession["ERZ1049443"] = [
        header,
        [213, "PF23864", "Domain of unknown function (DUF7222)"],
        [383, "PF02113", "D-Ala-D-Ala carboxypeptidase 3 (S13) family"],
        [149, "PF14269", "Arylsulfotransferase (ASST)"],
        [472, "PF11655", "Protein of unknown function (DUF2589)"],
        [21, "PF06039", "Malate:quinone oxidoreductase (Mqo)"],
        [225, "PF07987", "Domain of unkown function (DUF1775)"],
        [278, "PF10726", "Protein of function (DUF2518)"],
        [1, "PF24718", "HTH-like domain"],
    ]
    rows_per_accession["ERZ1049444"] = [
        header,
        [200, "PF23864", "Domain of unknown function (DUF7222)"],
        [291, "PF02113", "D-Ala-D-Ala carboxypeptidase 3 (S13) family"],
        [211, "PF14269", "Arylsulfotransferase (ASST)"],
        [48, "PF13558", "SbcC/RAD50-like, Walker B motif"],
        [498, "PF11655", "Protein of unknown function (DUF2589)"],
        [297, "PF06039", "Malate:quinone oxidoreductase (Mqo)"],
        [380, "PF07987", "Domain of unkown function (DUF1775)"],
        [429, "PF24718", "HTH-like domain"],
    ]
    rows_per_accession["ERZ1049445"] = [
        header,
        [280, "PF23864", "Domain of unknown function (DUF7222)"],
        [198, "PF02265", "S1/P1 Nuclease"],
        [265, "PF02113", "D-Ala-D-Ala carboxypeptidase 3 (S13) family"],
        [179, "PF14269", "Arylsulfotransferase (ASST)"],
        [428, "PF13558", "SbcC/RAD50-like, Walker B motif"],
        [252, "PF07987", "Domain of unkown function (DUF1775)"],
        [130, "PF10726", "Protein of function (DUF2518)"],
        [132, "PF24718", "HTH-like domain"],
    ]
    return rows_per_accession


@pytest.fixture
def antismash_summary_tsv_rows_per_accession() -> Dict:
    """
    Fixture that provides a TSV string representation of antiSMASH summary data.

    The TSV string includes the following columns:

    - **count**: The number of occurrences of the biosynthetic gene cluster (BGC).
    - **label**: The label or type of the BGC (e.g., "NRPS-like", "terpene").
    - **description**: A description of the BGC.

    The data is organized as a dictionary where keys are accession IDs (e.g., ERZ1049440)
    and values are lists of lists. The first list in each value is the header,
    and subsequent lists are rows of data.

    :returns: A dictionary mapping accession IDs to lists of antiSMASH summary data rows.
    :rtype: dict
    """
    rows_per_accession = {}
    header = ["count", "label", "description"]
    rows_per_accession["ERZ1049440"] = [
        header,
        [
            357,
            "NRPS-like,T1PKS",
            "NRPS-like fragment,Type I PKS (Polyketide synthase)",
        ],
        [114, "arylpolyene,resorcinol", "Aryl polyene,Resorcinol"],
        [367, "terpene", "Terpene"],
        [149, "arylpolyene", "Aryl polyene"],
        [185, "NRPS-like", "NRPS-like fragment"],
        [171, "betalactone", "Beta-lactone containing protease inhibitor"],
        [368, "NRPS", "Non-ribosomal peptide synthetase"],
        [
            341,
            "redox-cofactor",
            "Redox-cofactors such as PQQ (NC_021985:1458906-1494876)",
        ],
    ]
    rows_per_accession["ERZ1049443"] = [
        header,
        [251, "resorcinol", "Resorcinol"],
        [
            168,
            "NRPS-like,T1PKS",
            "NRPS-like fragment,Type I PKS (Polyketide synthase)",
        ],
        [389, "arylpolyene,resorcinol", "Aryl polyene,Resorcinol"],
        [447, "arylpolyene", "Aryl polyene"],
        [255, "NRPS-like", "NRPS-like fragment"],
        [324, "betalactone", "Beta-lactone containing protease inhibitor"],
        [452, "T1PKS", "Type I PKS (Polyketide synthase)"],
        [
            325,
            "redox-cofactor",
            "Redox-cofactors such as PQQ (NC_021985:1458906-1494876)",
        ],
    ]
    rows_per_accession["ERZ1049444"] = [
        header,
        [35, "resorcinol", "Resorcinol"],
        [
            456,
            "NRPS-like,T1PKS",
            "NRPS-like fragment,Type I PKS (Polyketide synthase)",
        ],
        [306, "arylpolyene,resorcinol", "Aryl polyene,Resorcinol"],
        [417, "terpene", "Terpene"],
        [19, "arylpolyene", "Aryl polyene"],
        [19, "betalactone", "Beta-lactone containing protease inhibitor"],
        [444, "NRPS", "Non-ribosomal peptide synthetase"],
        [
            297,
            "redox-cofactor",
            "Redox-cofactors such as PQQ (NC_021985:1458906-1494876)",
        ],
    ]
    rows_per_accession["ERZ1049445"] = [
        header,
        [160, "resorcinol", "Resorcinol"],
        [383, "arylpolyene,resorcinol", "Aryl polyene,Resorcinol"],
        [40, "terpene", "Terpene"],
        [92, "arylpolyene", "Aryl polyene"],
        [295, "NRPS-like", "NRPS-like fragment"],
        [388, "T1PKS", "Type I PKS (Polyketide synthase)"],
        [298, "NRPS", "Non-ribosomal peptide synthetase"],
        [
            247,
            "redox-cofactor",
            "Redox-cofactors such as PQQ (NC_021985:1458906-1494876)",
        ],
    ]
    return rows_per_accession


@pytest.fixture
def kegg_modules_summary_tsv_rows_per_accession() -> Dict:
    """
    Fixture that provides a TSV string representation of pathway summary data.

    The TSV string includes the following columns:

    - **completeness**: The completeness percentage of the pathway.
    - **module_accession**: The unique identifier for the pathway module.
    - **pathway_name**: The name of the pathway.
    - **pathway_class**: The classification of the pathway.
    - **matching_ko**: The KEGG Orthology identifiers that match the pathway.
    - **missing_ko**: The KEGG Orthology identifiers that are missing from the pathway (if any).

    The data includes several pathway entries with their respective details.

    :returns: A tab-separated values (TSV) string containing the pathway summary data.
    :rtype: str
    """
    rows_per_accession = {}
    header = [
        "completeness",
        "module_accession",
        "pathway_name",
        "pathway_class",
        "matching_ko",
        "missing_ko",
    ]

    rows_per_accession["ERZ1049440"] = [
        header,
        [
            52.1,
            "M00539",
            "Cumate degradation, p-cumate => 2-oxopent-4-enoate + 2-methylpropanoate",
            "Pathway modules; Xenobiotics biodegradation; Aromatics degradation",
            "K10619,K10620,K10621,K10622,K10623,K16303,K16304,K18227",
        ],
        [
            51.5,
            "M00851",
            "Carbapenem resistance",
            "Signature modules; Gene set; Drug resistance",
            "K18768",
        ],
        [
            36.5,
            "M00370",
            "Glucosinolate biosynthesis, tryptophan => glucobrassicin",
            "Pathway modules; Biosynthesis of other secondary metabolites; Biosynthesis of phytochemical compounds",
            "K11812,K11818,K11819,K11820,K11821",
        ],
        [
            37.7,
            "M00137",
            "Flavanone biosynthesis, phenylalanine => naringenin",
            "Pathway modules; Biosynthesis of other secondary metabolites; Biosynthesis of phytochemical compounds",
            "K00487,K00660,K01904,K10775 K01859",
        ],
        [
            98.0,
            "M00881",
            "Lipoic acid biosynthesis, plants and bacteria, octanoyl-ACP => dihydrolipoyl-E2/H",
            "Pathway modules; Metabolism of cofactors and vitamins; Cofactor and vitamin metabolism",
            "K03644,K03801",
        ],
        [
            39.2,
            "M00057",
            "Glycosaminoglycan biosynthesis, linkage tetrasaccharide",
            "Pathway modules; Glycan metabolism; Glycosaminoglycan metabolism",
            "K00733,K00734,K10158 K00771",
        ],
        [
            82.1,
            "M00543",
            "Biphenyl degradation, biphenyl => 2-oxopent-4-enoate + benzoate",
            "Pathway modules; Xenobiotics biodegradation; Aromatics degradation",
            "K00462,K08689,K08690,K10222,K15750,K18087,K18088",
        ],
        [
            39.5,
            "M00045",
            "Histidine degradation, histidine => N-formiminoglutamate => glutamate",
            "Pathway modules; Amino acid metabolism; Histidine metabolism",
            "K01468,K01479,K01712,K01745",
        ],
    ]
    rows_per_accession["ERZ1049443"] = [
        header,
        [
            58.3,
            "M00539",
            "Cumate degradation, p-cumate => 2-oxopent-4-enoate + 2-methylpropanoate",
            "Pathway modules; Xenobiotics biodegradation; Aromatics degradation",
            "K10619,K10620,K10621,K10622,K10623,K16303,K16304,K18227",
        ],
        [
            10.5,
            "M00851",
            "Carbapenem resistance",
            "Signature modules; Gene set; Drug resistance",
            "K18768",
        ],
        [
            86.8,
            "M00137",
            "Flavanone biosynthesis, phenylalanine => naringenin",
            "Pathway modules; Biosynthesis of other secondary metabolites; Biosynthesis of phytochemical compounds",
            "K00487,K00660,K01904,K10775 K01859",
        ],
        [
            81.5,
            "M00881",
            "Lipoic acid biosynthesis, plants and bacteria, octanoyl-ACP => dihydrolipoyl-E2/H",
            "Pathway modules; Metabolism of cofactors and vitamins; Cofactor and vitamin metabolism",
            "K03644,K03801",
        ],
        [
            8.6,
            "M00057",
            "Glycosaminoglycan biosynthesis, linkage tetrasaccharide",
            "Pathway modules; Glycan metabolism; Glycosaminoglycan metabolism",
            "K00733,K00734,K10158 K00771",
        ],
        [
            81.7,
            "M00974",
            "Betaine metabolism, animals, betaine => glycine",
            "Pathway modules; Amino acid metabolism; Serine and threonine metabolism",
            "K00306,K00315,K00544",
        ],
        [
            35.6,
            "M00045",
            "Histidine degradation, histidine => N-formiminoglutamate => glutamate",
            "Pathway modules; Amino acid metabolism; Histidine metabolism",
            "K01468,K01479,K01712,K01745",
        ],
        [
            13.7,
            "M00897",
            "Thiamine biosynthesis, plants, AIR (+ NAD+) => TMP/thiamine/TPP",
            "Pathway modules; Metabolism of cofactors and vitamins; Cofactor and vitamin metabolism",
            "K00949,K03146,K03147,K14153,K22911",
        ],
    ]
    rows_per_accession["ERZ1049444"] = [
        header,
        [
            98.7,
            "M00539",
            "Cumate degradation, p-cumate => 2-oxopent-4-enoate + 2-methylpropanoate",
            "Pathway modules; Xenobiotics biodegradation; Aromatics degradation",
            "K10619,K10620,K10621,K10622,K10623,K16303,K16304,K18227",
        ],
        [
            24.9,
            "M00851",
            "Carbapenem resistance",
            "Signature modules; Gene set; Drug resistance",
            "K18768",
        ],
        [
            8.2,
            "M00137",
            "Flavanone biosynthesis, phenylalanine => naringenin",
            "Pathway modules; Biosynthesis of other secondary metabolites; Biosynthesis of phytochemical compounds",
            "K00487,K00660,K01904,K10775 K01859",
        ],
        [
            86.1,
            "M00881",
            "Lipoic acid biosynthesis, plants and bacteria, octanoyl-ACP => dihydrolipoyl-E2/H",
            "Pathway modules; Metabolism of cofactors and vitamins; Cofactor and vitamin metabolism",
            "K03644,K03801",
        ],
        [
            11.3,
            "M00057",
            "Glycosaminoglycan biosynthesis, linkage tetrasaccharide",
            "Pathway modules; Glycan metabolism; Glycosaminoglycan metabolism",
            "K00733,K00734,K10158 K00771",
        ],
        [
            69.2,
            "M00974",
            "Betaine metabolism, animals, betaine => glycine",
            "Pathway modules; Amino acid metabolism; Serine and threonine metabolism",
            "K00306,K00315,K00544",
        ],
        [
            95.2,
            "M00543",
            "Biphenyl degradation, biphenyl => 2-oxopent-4-enoate + benzoate",
            "Pathway modules; Xenobiotics biodegradation; Aromatics degradation",
            "K00462,K08689,K08690,K10222,K15750,K18087,K18088",
        ],
        [
            63.8,
            "M00897",
            "Thiamine biosynthesis, plants, AIR (+ NAD+) => TMP/thiamine/TPP",
            "Pathway modules; Metabolism of cofactors and vitamins; Cofactor and vitamin metabolism",
            "K00949,K03146,K03147,K14153,K22911",
        ],
    ]
    rows_per_accession["ERZ1049445"] = [
        header,
        [
            29.7,
            "M00851",
            "Carbapenem resistance",
            "Signature modules; Gene set; Drug resistance",
            "K18768",
        ],
        [
            92.0,
            "M00370",
            "Glucosinolate biosynthesis, tryptophan => glucobrassicin",
            "Pathway modules; Biosynthesis of other secondary metabolites; Biosynthesis of phytochemical compounds",
            "K11812,K11818,K11819,K11820,K11821",
        ],
        [
            50.6,
            "M00137",
            "Flavanone biosynthesis, phenylalanine => naringenin",
            "Pathway modules; Biosynthesis of other secondary metabolites; Biosynthesis of phytochemical compounds",
            "K00487,K00660,K01904,K10775 K01859",
        ],
        [
            25.5,
            "M00057",
            "Glycosaminoglycan biosynthesis, linkage tetrasaccharide",
            "Pathway modules; Glycan metabolism; Glycosaminoglycan metabolism",
            "K00733,K00734,K10158 K00771",
        ],
        [
            55.1,
            "M00974",
            "Betaine metabolism, animals, betaine => glycine",
            "Pathway modules; Amino acid metabolism; Serine and threonine metabolism",
            "K00306,K00315,K00544",
        ],
        [
            15.4,
            "M00543",
            "Biphenyl degradation, biphenyl => 2-oxopent-4-enoate + benzoate",
            "Pathway modules; Xenobiotics biodegradation; Aromatics degradation",
            "K00462,K08689,K08690,K10222,K15750,K18087,K18088",
        ],
        [
            52.6,
            "M00045",
            "Histidine degradation, histidine => N-formiminoglutamate => glutamate",
            "Pathway modules; Amino acid metabolism; Histidine metabolism",
            "K01468,K01479,K01712,K01745",
        ],
        [
            22.9,
            "M00897",
            "Thiamine biosynthesis, plants, AIR (+ NAD+) => TMP/thiamine/TPP",
            "Pathway modules; Metabolism of cofactors and vitamins; Cofactor and vitamin metabolism",
            "K00949,K03146,K03147,K14153,K22911",
        ],
    ]
    return rows_per_accession


@pytest.fixture
def sanntis_modules_summary_tsv_rows_per_accession() -> Dict:
    """
    Fixture that provides a TSV string representation of SanntiS summary data.

    The TSV string includes the following columns:

    - **count**: The number of occurrences.
    - **nearest_mibig**: The nearest MIBiG identifier.
    - **nearest_mibig_class**: The class of the nearest MIBiG entry.
    - **description**: A brief description of the BGC.

    The data includes several BGC entries with their respective details.

    :returns: A tab-separated values (TSV) string containing the BGC summary data.
    :rtype: str
    """
    rows_per_accession = {}
    header = ["count", "nearest_mibig", "nearest_mibig_class", "description"]
    rows_per_accession["ERZ1049440"] = [
        header,
        [
            328,
            "BGC0001651",
            "Other",
            "Catch-all class for clusters encoding metabolites outside main classes (cyclitols, indolocarbazoles, and phosphonates).",
        ],
        [
            131,
            "BGC0000866",
            "Other",
            "Catch-all class for clusters encoding metabolites outside main classes (cyclitols, indolocarbazoles, and phosphonates).",
        ],
        [
            255,
            "BGC0000790",
            "Saccharide",
            "Carbohydrate-based natural products (e.g., aminoglycoside antibiotics)",
        ],
        [
            357,
            "BGC0000648",
            "Terpene",
            "Composed of isoprene (C5) units derived from isopentenyl pyrophosphate",
        ],
        [
            306,
            "BGC0000073",
            "Polyketide",
            "Built from iterative condensation of acetate units derived from acetyl-CoA",
        ],
        [
            15,
            "BGC0000248",
            "Polyketide",
            "Built from iterative condensation of acetate units derived from acetyl-CoA",
        ],
        [
            230,
            "BGC0001356",
            "RiPP",
            "Ribosomally synthesised and Post-translationally modified Peptide",
        ],
        [
            419,
            "BGC0001191",
            "Other",
            "Catch-all class for clusters encoding metabolites outside main classes (cyclitols, indolocarbazoles, and phosphonates).",
        ],
    ]
    rows_per_accession["ERZ1049443"] = [
        header,
        [
            341,
            "BGC0001651",
            "Other",
            "Catch-all class for clusters encoding metabolites outside main classes (cyclitols, indolocarbazoles, and phosphonates).",
        ],
        [
            79,
            "BGC0000866",
            "Other",
            "Catch-all class for clusters encoding metabolites outside main classes (cyclitols, indolocarbazoles, and phosphonates).",
        ],
        [
            10,
            "BGC0000790",
            "Saccharide",
            "Carbohydrate-based natural products (e.g., aminoglycoside antibiotics)",
        ],
        [
            373,
            "BGC0000648",
            "Terpene",
            "Composed of isoprene (C5) units derived from isopentenyl pyrophosphate",
        ],
        [
            425,
            "BGC0000073",
            "Polyketide",
            "Built from iterative condensation of acetate units derived from acetyl-CoA",
        ],
        [
            63,
            "BGC0000838",
            "Polyketide",
            "Built from iterative condensation of acetate units derived from acetyl-CoA",
        ],
        [
            314,
            "BGC0000248",
            "Polyketide",
            "Built from iterative condensation of acetate units derived from acetyl-CoA",
        ],
        [
            185,
            "BGC0001356",
            "RiPP",
            "Ribosomally synthesised and Post-translationally modified Peptide",
        ],
    ]
    rows_per_accession["ERZ1049444"] = [
        header,
        [
            233,
            "BGC0001651",
            "Other",
            "Catch-all class for clusters encoding metabolites outside main classes (cyclitols, indolocarbazoles, and phosphonates).",
        ],
        [
            427,
            "BGC0000866",
            "Other",
            "Catch-all class for clusters encoding metabolites outside main classes (cyclitols, indolocarbazoles, and phosphonates).",
        ],
        [
            179,
            "BGC0000790",
            "Saccharide",
            "Carbohydrate-based natural products (e.g., aminoglycoside antibiotics)",
        ],
        [
            436,
            "BGC0000073",
            "Polyketide",
            "Built from iterative condensation of acetate units derived from acetyl-CoA",
        ],
        [
            190,
            "BGC0000248",
            "Polyketide",
            "Built from iterative condensation of acetate units derived from acetyl-CoA",
        ],
        [
            182,
            "BGC0001356",
            "RiPP",
            "Ribosomally synthesised and Post-translationally modified Peptide",
        ],
        [
            70,
            "BGC0001191",
            "Other",
            "Catch-all class for clusters encoding metabolites outside main classes (cyclitols, indolocarbazoles, and phosphonates).",
        ],
        [6, "BGC0001432", "NRP Polyketide", "Nonribosomal Peptide Polyketide"],
    ]
    rows_per_accession["ERZ1049445"] = [
        header,
        [
            178,
            "BGC0001651",
            "Other",
            "Catch-all class for clusters encoding metabolites outside main classes (cyclitols, indolocarbazoles, and phosphonates).",
        ],
        [
            19,
            "BGC0000790",
            "Saccharide",
            "Carbohydrate-based natural products (e.g., aminoglycoside antibiotics)",
        ],
        [
            432,
            "BGC0000648",
            "Terpene",
            "Composed of isoprene (C5) units derived from isopentenyl pyrophosphate",
        ],
        [
            405,
            "BGC0000838",
            "Polyketide",
            "Built from iterative condensation of acetate units derived from acetyl-CoA",
        ],
        [
            146,
            "BGC0000248",
            "Polyketide",
            "Built from iterative condensation of acetate units derived from acetyl-CoA",
        ],
        [
            410,
            "BGC0001356",
            "RiPP",
            "Ribosomally synthesised and Post-translationally modified Peptide",
        ],
        [
            471,
            "BGC0001191",
            "Other",
            "Catch-all class for clusters encoding metabolites outside main classes (cyclitols, indolocarbazoles, and phosphonates).",
        ],
        [287, "BGC0001432", "NRP Polyketide", "Nonribosomal Peptide Polyketide"],
    ]
    return rows_per_accession


@pytest.fixture
def taxonomy_tsv_rows_per_accession() -> Dict:
    """
    Fixture that provides a TSV string representation of taxonomy summary data.

    The TSV string includes the following columns:
    - Identifier: A unique identifier for each entry.
    - Superkingdom: The highest taxonomic rank (e.g., Archaea or Bacteria).
    - Kingdom: The next taxonomic rank (e.g., Thermoproteati, Bacillati).
    - Phylum: The phylum classification (e.g., Nitrososphaerota, Chloroflexota).
    - Class: The class classification (e.g., Nitrososphaeria, Chloroflexia).
    - Order: The order classification (e.g., Nitrosopumilales, Dehalococcoidia).
    - Family: The family classification (e.g., Nitrosopumilaceae, Nitrosotaleaceae).
    - Genus: The genus classification (e.g., Nitrosopumilus, Nitrosotalea).
    - Species: The species classification (e.g., Candidatus Nitrosopumilus koreensis, Nitrosotalea devaniterrae).

    The data includes several entries with their respective taxonomic details.

    :returns: A tab-separated values (TSV) string containing the taxonomy summary data.
    :rtype: str
    """
    rows_per_accession = {}

    rows_per_accession["ERZ1049440"] = [
        [
            53,
            "sk__Archaea",
            "k__Thermoproteati",
            "p__Nitrososphaerota",
            "c__Nitrososphaeria",
            "o__Nitrosopumilales",
            "f__Nitrosopumilaceae",
            "g__Nitrosopumilus",
            "s__Candidatus Nitrosopumilus koreensis",
        ],
        [
            786,
            "sk__Archaea",
            "k__Thermoproteati",
            "p__Nitrososphaerota",
            "c__Nitrososphaeria",
            "o__Nitrosotaleales",
            "f__Nitrosotaleaceae",
            "g__Nitrosotalea",
            "s__Nitrosotalea devaniterrae",
        ],
        [29, "sk__Bacteria"],
        [275, "sk__Bacteria", "k__Bacillati"],
        [
            62,
            "sk__Bacteria",
            "k__Bacillati",
            "p__Chloroflexota",
            "c__Chloroflexia",
            "o__",
            "f__",
            "g__",
            "s__Chloroflexia bacterium",
        ],
        [7, "sk__Bacteria", "k__Bacillati", "p__Chloroflexota", "c__Dehalococcoidia"],
        [
            7293,
            "sk__Bacteria",
            "k__Bacillati",
            "p__Actinomycetota",
            "c__Acidimicrobiia",
            "o__Acidimicrobiales",
        ],
        [742, "sk__Bacteria", "k__Bacillati", "p__Actinomycetota", "c__Actinomycetes"],
        [31, "sk__Bacteria", "k__Pseudomonadati", "p__Bacteroidota", "c__Bacteroidia"],
        [13, "unclassified"],
    ]
    rows_per_accession["ERZ1049443"] = [
        [
            7,
            "sk__Archaea",
            "k__Thermoproteati",
            "p__Nitrososphaerota",
            "c__Nitrososphaeria",
            "o__Nitrosopumilales",
            "f__Nitrosopumilaceae",
            "g__Nitrosopumilus",
            "s__Candidatus Nitrosopumilus koreensis",
        ],
        [
            3,
            "sk__Archaea",
            "k__Thermoproteati",
            "p__Nitrososphaerota",
            "c__Nitrososphaeria",
            "o__Nitrosotaleales",
            "f__Nitrosotaleaceae",
            "g__Nitrosotalea",
            "s__Nitrosotalea devaniterrae",
        ],
        [98, "sk__Bacteria"],
        [2, "sk__Bacteria", "k__Bacillati"],
        [
            3,
            "sk__Bacteria",
            "k__Bacillati",
            "p__Chloroflexota",
            "c__Chloroflexia",
            "o__",
            "f__",
            "g__",
            "s__Chloroflexia bacterium",
        ],
        [2, "sk__Bacteria", "k__Bacillati", "p__Chloroflexota", "c__Dehalococcoidia"],
        [
            8,
            "sk__Bacteria",
            "k__Bacillati",
            "p__Actinomycetota",
            "c__Acidimicrobiia",
            "o__Acidimicrobiales",
        ],
        [1, "sk__Bacteria", "k__Bacillati", "p__Actinomycetota", "c__Actinomycetes"],
        [1, "sk__Bacteria", "k__Pseudomonadati", "p__Bacteroidota", "c__Bacteroidia"],
        [46, "unclassified"],
    ]
    rows_per_accession["ERZ1049444"] = [
        [
            3,
            "sk__Archaea",
            "k__Thermoproteati",
            "p__Nitrososphaerota",
            "c__Nitrososphaeria",
            "o__Nitrosopumilales",
            "f__Nitrosopumilaceae",
            "g__Nitrosopumilus",
            "s__Candidatus Nitrosopumilus koreensis",
        ],
        [
            1,
            "sk__Archaea",
            "k__Thermoproteati",
            "p__Nitrososphaerota",
            "c__Nitrososphaeria",
            "o__Nitrosopumilales",
            "f__Nitrosopumilaceae",
            "g__Nitrosopumilus",
            "s__Nitrosopumilus adriaticus",
        ],
        [
            3,
            "sk__Archaea",
            "k__Thermoproteati",
            "p__Nitrososphaerota",
            "c__Nitrososphaeria",
            "o__Nitrosotaleales",
            "f__Nitrosotaleaceae",
            "g__Nitrosotalea",
            "s__Nitrosotalea devaniterrae",
        ],
        [23651, "sk__Bacteria"],
        [2, "sk__Bacteria", "k__Bacillati"],
        [
            8,
            "sk__Bacteria",
            "k__Bacillati",
            "p__Actinomycetota",
            "c__Acidimicrobiia",
            "o__Acidimicrobiales",
        ],
        [1, "sk__Bacteria", "k__Pseudomonadati", "p__Bacteroidota", "c__Bacteroidia"],
        [6577, "unclassified"],
    ]
    rows_per_accession["ERZ1049445"] = [
        [
            76,
            "sk__Archaea",
            "k__Thermoproteati",
            "p__Nitrososphaerota",
            "c__Nitrososphaeria",
            "o__Nitrosopumilales",
            "f__Nitrosopumilaceae",
            "g__Nitrosopumilus",
            "s__Candidatus Nitrosopumilus koreensis",
        ],
        [
            575,
            "sk__Archaea",
            "k__Thermoproteati",
            "p__Nitrososphaerota",
            "c__Nitrososphaeria",
            "o__Nitrosopumilales",
            "f__Nitrosopumilaceae",
            "g__Nitrosopumilus",
            "s__Nitrosopumilus adriaticus",
        ],
        [
            786,
            "sk__Archaea",
            "k__Thermoproteati",
            "p__Nitrososphaerota",
            "c__Nitrososphaeria",
            "o__Nitrosotaleales",
            "f__Nitrosotaleaceae",
            "g__Nitrosotalea",
            "s__Nitrosotalea devaniterrae",
        ],
        [9654, "sk__Bacteria"],
        [275, "sk__Bacteria", "k__Bacillati"],
        [
            8259,
            "sk__Bacteria",
            "k__Bacillati",
            "p__Actinomycetota",
            "c__Acidimicrobiia",
            "o__Acidimicrobiales",
        ],
        [742, "sk__Bacteria", "k__Bacillati", "p__Actinomycetota", "c__Actinomycetes"],
        [31, "sk__Bacteria", "k__Pseudomonadati", "p__Bacteroidota", "c__Bacteroidia"],
        [1521, "unclassified"],
    ]
    return rows_per_accession


def create_results_dirs(fs: Any, accession: str) -> None:
    """
    Create the necessary directory structure for the results.

    :param fs: The fake filesystem instance.
    :param accession: The accession ID for which to create directories.
    """
    fs.create_dir(f"/{accession}/functional-annotation/go")
    fs.create_dir(f"/{accession}/functional-annotation/interpro")
    fs.create_dir(f"/{accession}/functional-annotation/kegg")
    fs.create_dir(f"/{accession}/functional-annotation/pfam")
    fs.create_dir(f"/{accession}/pathways-and-systems/antismash")
    fs.create_dir(f"/{accession}/pathways-and-systems/kegg")
    fs.create_dir(f"/{accession}/pathways-and-systems/sanntis")
    fs.create_dir(f"/{accession}/taxonomy")


def write_gzipped_tsv(file_path: str, rows: List[List[Any]]) -> None:
    """
    Write the provided TSV rows to a gzipped file.

    :param file_path: The path to the file to write.
    :param rows: The TSV rows to write.
    """
    with gzip.open(file_path, "wt") as fh:
        fh.write(tsv_string(rows))


def calculate_md5(file_path: str) -> str:
    """
    Calculate the MD5 checksum of a file.

    :param file_path: The path to the file.
    :return: The MD5 checksum as a hexadecimal string.
    """
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)
    return md5_hash.hexdigest()


def test_assembly_study_summary_generator_correct_summarise(
    fs: Any,
    go_summary_tsv_rows_per_accession: Dict,
    interpro_summary_tsv_rows_per_accession: Dict,
    ko_summary_summary_tsv_rows_per_accession: Dict,
    pfam_summary_summary_tsv_rows_per_accession: Dict,
    antismash_summary_tsv_rows_per_accession: Dict,
    kegg_modules_summary_tsv_rows_per_accession: Dict,
    sanntis_modules_summary_tsv_rows_per_accession: Dict,
    taxonomy_tsv_rows_per_accession: Dict,
) -> None:
    """
    Summarize a correct set of results :).
    """
    accessions = [
        "ERZ1049440",
        "ERZ1049443",
        "ERZ1049444",
        "ERZ1049445",
    ]

    for accession in accessions:
        create_results_dirs(fs, accession)
        # functional-annotation #
        write_gzipped_tsv(
            f"/{accession}/functional-annotation/go/{accession}_go_summary.tsv.gz",
            go_summary_tsv_rows_per_accession[accession],
        )
        # We are reusing the GO results for GO slim, the structure of the results is the same
        write_gzipped_tsv(
            f"/{accession}/functional-annotation/go/{accession}_goslim_summary.tsv.gz",
            go_summary_tsv_rows_per_accession[accession],
        )
        write_gzipped_tsv(
            f"/{accession}/functional-annotation/interpro/{accession}_interpro_summary.tsv.gz",
            interpro_summary_tsv_rows_per_accession[accession],
        )
        write_gzipped_tsv(
            f"/{accession}/functional-annotation/kegg/{accession}_ko_summary.tsv.gz",
            ko_summary_summary_tsv_rows_per_accession[accession],
        )
        write_gzipped_tsv(
            f"/{accession}/functional-annotation/pfam/{accession}_pfam_summary.tsv.gz",
            pfam_summary_summary_tsv_rows_per_accession[accession],
        )

        # pathways-and-systems #
        write_gzipped_tsv(
            f"/{accession}/pathways-and-systems/antismash/{accession}_antismash_summary.tsv.gz",
            antismash_summary_tsv_rows_per_accession[accession],
        )
        write_gzipped_tsv(
            f"/{accession}/pathways-and-systems/kegg/{accession}_kegg_modules_summary.tsv.gz",
            kegg_modules_summary_tsv_rows_per_accession[accession],
        )
        write_gzipped_tsv(
            f"/{accession}/pathways-and-systems/sanntis/{accession}_sanntis_summary.tsv.gz",
            sanntis_modules_summary_tsv_rows_per_accession[accession],
        )
        write_gzipped_tsv(
            f"/{accession}/taxonomy/{accession}.krona.txt.gz",
            taxonomy_tsv_rows_per_accession[accession],
        )

    input_csv = "/input.csv"
    fs.create_file(
        input_csv,
        contents="\n".join(f"{accession},success" for accession in accessions),
    )

    runner = CliRunner()
    result = runner.invoke(
        summarise_analyses,
        ["--assemblies", input_csv, "--study_dir", "/", "-p", "test"],
    )
    assert result.exit_code == 0

    # Verify that the output files exist and have the expected MD5 checksums
    expected_files_md5 = {
        "/test_go_summary.tsv": "3783ef924b0f7d598888334025de1834",
        "/test_goslim_summary.tsv": "3783ef924b0f7d598888334025de1834",
        "/test_interpro_summary.tsv": "ccdb11f77c1631ec3dad99226e791cbb",
        "/test_ko_summary.tsv": "9a0e4ebc8d3c096738ab9be5658fbc0d",
        "/test_kegg_modules_summary.tsv": "c66966de0c38f74ed83919c9fe0dd1b4",
        "/test_pfam_summary.tsv": "314ec298af0a15fbc764b731a120e54d",
        "/test_sanntis_summary.tsv": "97b79aee421aaee39bf6409333218e03",
        "/test_antismash_summary.tsv": "c3105b172ab8c795bb0593e1dc5abe98",
        "/test_taxonomy_summary.tsv": "6260f20e4237d81c00d9787e02653734",
    }

    for file_path, expected_md5 in expected_files_md5.items():
        assert os.path.exists(file_path), f"Output file {file_path} does not exist"
        actual_md5 = calculate_md5(file_path)
        assert (
            actual_md5 == expected_md5
        ), f"MD5 checksum mismatch for {file_path}. Expected: {expected_md5}, Actual: {actual_md5}"
