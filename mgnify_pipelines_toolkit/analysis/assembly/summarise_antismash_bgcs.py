#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
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
import fileinput
import logging
import pandas as pd


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s"
)

ANTISMASH_VERSION = "7.1.x"

f"""
Script parses antismash GFF output and adds descriptions from pre-parsed glossary https://docs.antismash.secondarymetabolites.org/glossary/.
Glossary was taken from version {ANTISMASH_VERSION} and commit dbeeb0e https://github.com/antismash/documentation/blob/master/docs/glossary.md
"""

DESCRIPTIONS = {
    "2dos": "2-deoxy-streptamine aminoglycoside",
    "acyl_amino_acids": "N-acyl amino acid",
    "amglyccycl": "Aminoglycoside/aminocyclitol",
    "aminocoumarin": "Aminocoumarin",
    "aminopolycarboxylic-acid": "Aminopolycarboxylic acid metallophores (doi:10.1039/C8MT00009C)",
    "archaeal-ripp": "Archaeal RiPPs (doi:10.1021/jacs.2c00521 supplemental)",
    "arylpolyene": "Aryl polyene",
    "atropopeptide": "Atropopeptide RiPPs, e.g. scabrirubin and tryptorubin",
    "azoxy-crosslink": "axoxy compounds formed by carboxilic cross-link",
    "azoxy-dimer": "axoxy compounds formed by dimerisation",
    "benzoxazole": "Benzoxazoles",
    "betalactone": "Beta-lactone containing protease inhibitor",
    "blactam": "&beta;-lactam",
    "bottromycin": "Bottromycin",
    "butyrolactone": "Butyrolactone",
    "cdps": "tRNA-dependent cyclodipeptide synthases",
    "crocagin": "Crocagin-like",
    "cyanobactin": "Cyanobactins like patellamide (AY986476)",
    "cyclic-lactone-autoinducer": "agrD-like cyclic lactone autoinducer peptides (AF001782)",
    "cytokinin": "Adenine-type cytokinins, e.g. fusatin and trans-zeatin",
    "darobactin": "Darobactin-like compounds",
    "deazapurine": "Deazapurine",
    "ectoine": "Ectoine",
    "epipeptide": "D-amino-acid containing RiPPs such as yydF (D78193)",
    "fungal_cdps": "Fungal cyclodipeptide synthases",
    "fungal-ripp": "Fungal RiPP with POP or UstH peptidase types and a modification",
    "furan": "Furan",
    "glycocin": "Glycocin",
    "guanidinotides": "Pheganomycin-style protein ligase-containing cluster",
    "hgle-ks": "Heterocyst glycolipid synthase-like PKS",
    "hr-t2pks": "Highly reducing type II PKS like ishigamide and skyllamycin",
    "hserlactone": "Homoserine lactone",
    "hydrogen-cyanide": "Hydrogen cyanide (AF208523, doi:10.1128/jb.182.24.6940-6949.20)",
    "hydroxy-tropolone": "7-hydroxytropolone-like cluster",
    "indole": "Indole",
    "isocyanide": "Isocyanides (doi:10.1093/nar/gkad573)",
    "nrp with isocyanide": "Isocyanides (doi:0.1128/mBio.00785-18)",
    "ladderane": "Ladderane",
    "lanthipeptide class i": "Class I lanthipeptides like nisin",
    "lanthipeptide class ii": "Class II lanthipeptides like mutacin II (U40620)",
    "lanthipeptide class iii": "Class III lanthipeptides like labyrinthopeptin (FN178622)",
    "lanthipeptide class iv": "Class IV lanthipeptides like venezuelin (HQ328852)",
    "lanthipeptide class v": "Glycosylated lanthipeptide/linaridin hybrids like MT210103",
    "lassopeptide": "Lasso peptide",
    "leupeptin": "leupeptin-like compounds",
    "linaridin": "Linear arid peptide such as cypemycin (HQ148718) and salinipeptin (MG788286)",
    "lincosamides": "NRPS-adjacent biosynthesis of lincosamides",
    "lipolanthine": "Lanthipeptide class containing N-terminal fatty acids such as MG673929",
    "melanin": "Melanin",
    "methanobactin": "Copper-chelating/transporting peptides (doi:10.1126/science.aap9437)",
    "microviridin": "Microviridin",
    "mycosporine": "Molecules containing mycosporine-like amino acid",
    "naggn": "N-acetylglutaminylglutamine amide",
    "napaa": "Non-alpha poly-amino acids like e-Polylysin",
    "ni-siderophore": "NRPS-independent, IucA/IucC-like siderophores (*siderophore* prior to 7.0)",
    "nitropropanoic-acid": "3-Nitropropanoic acid (neurotoxin)",
    "nrps": "Non-ribosomal peptide synthetase",
    "nrp-metallophore": "Non-ribosomal peptide metallophores",
    "nucleoside": "Nucleoside",
    "oligosaccharide": "Oligosaccharide",
    "opine-like-metallophore": "Opine-like zincophores like staphylopine (doi:10.1128/mSystems.00554-20)",
    "other": "Cluster containing a secondary metabolite-related protein that does not fit into any other category",
    "pbde": "Polybrominated diphenyl ether",
    "phenazine": "Phenazine",
    "phosphoglycolipid": "Phosphoglycolipid",
    "phosphonate": "Phosphonate",
    "polyhalogenated-pyrrole": "Polyhalogenated pyrrole",
    "polyyne": "Polyyne",
    "ppys-ks": "PPY-like pyrone",
    "prodigiosin": "Serratia-type non-traditional PKS prodigiosin biosynthesis pathway",
    "proteusin": "Proteusin",
    "pufa": "Polyunsaturated fatty acid",
    "pyrrolidine": "Pyrrolidines like described in BGC0001510",
    "ranthipeptide": "Cys-rich peptides (aka. SCIFF: six Cys in fourty-five) like in CP001581:3481278-3502939",
    "ras-ripp": "Streptide-like thioether-bond RiPPs",
    "rcdps": "Fungal Arginine-containing cyclic dipeptides",
    "redox-cofactor": "Redox-cofactors such as PQQ (NC_021985:1458906-1494876)",
    "resorcinol": "Resorcinol",
    "sactipeptide": "Sactipeptide",
    "spliceotide": "RiPPs containing plpX type spliceases (NZ_KB235920:17899-42115)",
    "t1pks": "Type I PKS (Polyketide synthase)",
    "t2pks": "Type II PKS",
    "t3pks": "Type III PKS",
    "terpene": "Terpene",
    "thioamitides": "Thioamitide RiPPs as found in JOBF01000011",
    "thioamide-nrp": "Thioamide-containing non-ribosomal peptide",
    "transat-pks": "Trans-AT PKS",
    "triceptide": "Triceptides",
    "tropodithietic-acid": "Tropodithietic acid",
    "fungal-ripp-like": "Fungal RiPP-likes",
    "nrps-like": "NRPS-like fragment",
    "phosphonate-like": "Phosphonate-like (prior to 7.0 this was the phosphonate rule)",
    "pks-like": "Other types of PKS",
    "ripp-like": "Other unspecified ribosomally synthesised and post-translationally modified peptide product (RiPP)",
    "rre-containing": "RRE-element containing cluster",
    "terpene-precursor": "Compound likely used as a terpene precursor",
    "transat-pks-like": "Trans-AT PKS fragment, with trans-AT domain not found",
    "fatty_acid": "Fatty acid (loose strictness, likely from primary metabolism)",
    "halogenated": "Halogenase-containing cluster, potentially generating a halogenated product",
    "lysine": "Fungal lysine primary metabolism",
    "saccharide": "Saccharide (loose strictness, likely from primary metabolism)",
    "lap": "Linear azol(in)e-containing peptides",
    "mycosporine-like": "Molecules containing mycosporine-like amino acid",
    "thiopeptide": "Thiopeptide",
    "siderophore": "Siderophore",
    "bacteriocin": "Bacteriocin or other unspecified ribosomally synthesised and  post-translationally modified peptide product (RiPP)",
    "fused": "Pheganomycin-style protein ligase-containing cluster",
    "head_to_tail": "Head-to-tail cyclised RiPP (subtilosin-like)",
    "lanthidin": "Glycosylated lanthipeptide/linaridin hybrids like MT210103",
    "lanthipeptide": "Lanthipeptides",
    "tfua-related": "TfuA-related RiPPs",
    "otherks": "Other types of PKS",
    "microcin": "Microcin",
    "cf_saccharide": "Possible saccharide",
    "cf_fatty_acid": "Possible fatty acid",
    "cf_putative": "Putative cluster of unknown type identified  with the ClusterFinder algorithm",
}


def parse_args():
    description = (
        "antiSMASH output summary generator. "
        "Script takes regions from GFF and counts its appearance in annotation. "
        "Output columns contain label, descriptions and count. "
        f"Descriptions were taken from pre-parsed glossary provided on antiSMASH website. "
        f"Current script supports antiSMASH results for version {ANTISMASH_VERSION} and older."
    )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--antismash-gff", help="antiSMASH GFF", required=True)
    parser.add_argument(
        "-o", "--output", help="Antisamsh summary TSV output file.", required=True
    )
    parser.add_argument(
        "-a",
        "--antismash-version",
        help="antiSMASH version that was used to generate GFF",
        required=False,
        default=ANTISMASH_VERSION,
    )
    args = parser.parse_args()
    if args.antismash_version > ANTISMASH_VERSION:
        logging.error(
            "Provided version of antiSMASH is bigger than supported. "
            "Please, make sure you have updated descriptions dictionary. Exit."
        )
        exit(1)
    return args.antismash_gff, args.output


def main():
    input_gff, output_filename = parse_args()
    dict_list = []
    with fileinput.hook_compressed(input_gff, "r") as file_in:
        # TODO: to be merged with the GFF toolkit
        for line in file_in:
            if line.startswith("#"):
                continue
            info = line.strip().split("\t")[8].split(";")
            entry_dict = {}
            for pair in info:
                key, value = pair.split(
                    "=", 1
                )  # Ensure split only occurs at the first '=' occurrence
                entry_dict[key] = value
            dict_list.append(entry_dict)

        # Convert to DataFrame
        df = pd.DataFrame(dict_list)
        df = df[df["product"].notna()]
        df_grouped = (
            df.groupby(["product"]).size().reset_index(name="count")
        ).sort_values(by="count", ascending=False)

        df_grouped = df_grouped.rename(
            columns={
                "product": "label",
            }
        )
        df_grouped["description"] = df_grouped["label"].apply(
            lambda x: ",".join(
                [
                    DESCRIPTIONS.get(cls.strip().lower(), cls.strip())
                    for cls in x.split(",")
                ]
            )
        )
        df_grouped = df_grouped[["label", "description", "count"]]
        df_grouped.to_csv(output_filename, sep="\t", index=False)


if __name__ == "__main__":
    main()
