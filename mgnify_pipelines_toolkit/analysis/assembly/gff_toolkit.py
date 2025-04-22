#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import argparse

from mgnify_pipelines_toolkit.analysis.assembly.gff_annotation_utils import (
    get_ncrnas,
    get_trnas,
    load_annotations,
    load_crispr,
)
from mgnify_pipelines_toolkit.analysis.assembly.gff_file_utils import (
    write_results_to_file,
    print_pseudogene_report,
)


def main():

    (
        gff,
        ipr_file,
        eggnog_file,
        sanntis_file,
        crispr_file,
        amr_file,
        antismash_file,
        gecco_file,
        dbcan_file,
        dbcan_cazys_file,
        defense_finder_file,
        pseudofinder_file,
        rfam_file,
        trnascan_file,
        outfile,
        pseudogene_report_file,
    ) = parse_args()

    # load annotations and add them to existing CDS
    # here header contains leading GFF lines starting with "#",
    # main_gff_extended is a dictionary that contains GFF lines with added in additional annotations
    # fasta is the fasta portion of the original GFF file
    # pseudogene_report_dict is the information on detected pseudogene which can be optionally printed
    # to a separate output file
    header, main_gff_extended, fasta, pseudogene_report_dict = load_annotations(
        gff,
        eggnog_file,
        ipr_file,
        sanntis_file,
        amr_file,
        antismash_file,
        gecco_file,
        dbcan_file,
        dbcan_cazys_file,
        defense_finder_file,
        pseudofinder_file,
    )
    ncrnas = {}
    if rfam_file:
        ncrnas = get_ncrnas(rfam_file)
    trnas = {}
    if trnascan_file:
        trnas = get_trnas(trnascan_file)
    crispr_annotations = {}
    if crispr_file:
        crispr_annotations = load_crispr(crispr_file)

    write_results_to_file(
        outfile, header, main_gff_extended, fasta, ncrnas, trnas, crispr_annotations
    )
    if pseudogene_report_file:
        print_pseudogene_report(pseudogene_report_dict, pseudogene_report_file)


def parse_args():
    parser = argparse.ArgumentParser(
        description="The script extends a user-provided base GFF annotation file by incorporating "
        "information extracted from the user-provided outputs of supplementary annotation tools.",
    )
    parser.add_argument(
        "-g",
        dest="gff_input",
        required=True,
        help="GFF input file containing the base annotation",
    )
    parser.add_argument(
        "-i",
        dest="ips",
        help="InterProScan annotation results (TSV)",
        required=False,
    )
    parser.add_argument(
        "-e",
        dest="eggnog",
        help="EggNOG mapper annotation results (TSV)",
        required=False,
    )
    parser.add_argument(
        "-s",
        dest="sanntis",
        help="SanntiS results",
        required=False,
    )
    parser.add_argument(
        "-c",
        dest="crispr",
        help="CRISPRCasFinder results for the cluster rep (pre-filtered high quality GFF)",
        required=False,
    )
    parser.add_argument(
        "-a",
        dest="amr",
        help="The TSV file produced by AMRFinderPlus",
        required=False,
    )
    parser.add_argument(
        "--antismash",
        help="The GFF file produced by AntiSMASH post-processing script",
        required=False,
    )
    parser.add_argument(
        "--gecco",
        help="The GFF file produced by GECCO",
        required=False,
    )
    parser.add_argument(
        "--dbcan",
        help="The GFF file produced by dbCAN post-processing script that uses cluster annotations",
        required=False,
    )
    parser.add_argument(
        "--dbcan-cazys",
        help="The GFF file produced by dbCAN-CAZYs post-processing script",
        required=False,
    )
    parser.add_argument(
        "--defense-finder",
        help="The GFF file produced by Defense Finder post-processing script",
        required=False,
    )
    parser.add_argument(
        "--pseudofinder",
        help="The GFF file produced by the Pseudofinder post-processing script",
        required=False,
    )
    parser.add_argument("-r", dest="rfam", help="Rfam results", required=False)
    parser.add_argument(
        "-t", dest="trnascan", help="tRNAScan-SE results", required=False
    )
    parser.add_argument("-o", dest="outfile", help="Outfile name", required=True)
    parser.add_argument(
        "--pseudogene-report", help="Pseudogene report filename", required=False
    )

    args = parser.parse_args()
    return (
        args.gff_input,
        args.ips,
        args.eggnog,
        args.sanntis,
        args.crispr,
        args.amr,
        args.antismash,
        args.gecco,
        args.dbcan,
        args.dbcan_cazys,
        args.defense_finder,
        args.pseudofinder,
        args.rfam,
        args.trnascan,
        args.outfile,
        args.pseudogene_report,
    )


if __name__ == "__main__":
    main()
