##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 15:19:48 2025

@author: jlu8


range of scripts to manipulate mgnify data into vector files for creating syncoms with the mimic toolkit - https://github.com/thh32/MiMiC2.git

- generate vector from assembly analysis (pfam)
- generate vector from mgnify



"""


import os
import argparse
import sys
import warnings
import pandas as pd  # type: ignore
import time
from io import BytesIO
import requests
import gzip
import urllib.request
import tarfile
import shutil
from pathlib import Path


# Suppress FutureWarnings
warnings.simplefilter(action="ignore", category=FutureWarning)


# read pfam clan data file into memory
def read_pfamClanInfo2df(pfam_clan_file):

    pfam_clan_df = pd.read_csv(
        pfam_clan_file,
        sep="\t",
        names=["Protein_family", "class", "name", "target.name", "Description"],
    )

    # rename Protein_family to PfamID
    pfam_clan_df.rename(columns={"Protein_family": "PfamID"}, inplace=True)

    return pfam_clan_df


# get mgya_ids from mgys using api
def extract_mgya_from_mgys(mgys_id):

    max_retries = 3  # max retires for url accessing
    retry_delay = 5

    mgya_ids = []  # empty list for storing ids

    # api_Endpoint for getting all analyses information belonging to mgys_id.
    url = f"https://www.ebi.ac.uk/metagenomics/api/v1/studies/{mgys_id}/analyses"

    for attempt in range(max_retries):
        try:
            response = requests.get(url)

            if response.status_code == 200:
                api_endpoint = response.json()

                # get all analyses ids that are from analysis of assemblies
                mgya_ids = [
                    data["id"]
                    for data in api_endpoint["data"]
                    if data["attributes"]["experiment-type"] == "assembly"
                ]

                # get all analyses ids that are from analysis of assemblies
                mgya_ids = [
                    data["id"]
                    for data in api_endpoint["data"]
                    if data["attributes"]["experiment-type"] == "assembly"
                ]

                break
            else:
                pass
                # print(f"Attempt {attempt + 1} failed with status code: {response.status_code}")
        except requests.exceptions.RequestException as e:
            print(f"Attempt {attempt + 1} failed due to error: {e}")

        if attempt < max_retries - 1:
            print(f"Retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)

        else:
            print(
                "All retries failed. Please check the MGYS study accesion or try again later."
            )

    return mgya_ids


# download pfam results for a given MGYA accession
def download_pfam_from_mgya(mgya_id):

    max_retries = 3  # max retires for url accessing
    # retry_delay = 5

    assembly_id = " "
    sample_id = " "

    # api endpoint for downloads
    url = f"https://www.ebi.ac.uk/metagenomics/api/v1/analyses/{mgya_id}/downloads"

    mgya_pfam_df = pd.DataFrame(columns=["protein", "PfamID", "pfam_description"])

    mgya_pfams = []

    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            if response.status_code == 200:
                for downloadable_file in response.json()["data"]:

                    if (
                        "pfam"
                        in downloadable_file["attributes"]["description"][
                            "label"
                        ].lower()
                    ):

                        # url of pfam file
                        url = downloadable_file["links"]["self"]
                        mgya_pfam_df = pd.read_csv(
                            url,
                            sep=",",
                            names=["protein", "PfamID", "pfam_description"],
                        )

                        mgya_pfams = list(mgya_pfam_df["PfamID"])

                        break
            else:
                print(
                    f"Connection attempt failed with status code: {response.status_code}"
                )
        except requests.exceptions.RequestException as e:
            print(f"Connection attempt {attempt + 1} failed due to error: {e}")

        try:
            response = requests.get(
                f"https://www.ebi.ac.uk/metagenomics/api/v1/analyses/{mgya_id}"
            )
            if response.status_code == 200:
                mgya_metadata = response.json()["data"]["relationships"]
                sample_id = mgya_metadata["sample"]["data"]["id"]
                assembly_id = mgya_metadata["assembly"]["data"]["id"]
            else:
                print(
                    f"Attempt {attempt + 1} failed with status code: {response.status_code}"
                )

        except requests.exceptions.RequestException as e:
            print(
                f"Attempt {attempt + 1} failed due to error: {e}. Please try again later"
            )

    return sample_id, assembly_id, mgya_pfams


# for a given mgys_id, download list of pfams for each mgya and get pfams for each mgya into a pfam dict
def download_pfams_from_mgys(mgys_id):

    samples_pfam_dict = {}

    mgya_list = extract_mgya_from_mgys(mgys_id)

    # go through list and download pfams
    for mgya_id in mgya_list:
        print(f"extracting pfams for {mgya_id}")
        sample_id, assembly_id, pfams_list = download_pfam_from_mgya(mgya_id)

        time.sleep(2)

        samples_pfam_dict[mgya_id] = pfams_list

    # print (samples_pfam_dict)

    return samples_pfam_dict


def download_pfams_from_mgya_list(mgya_list):

    samples_pfam_dict = {}

    # go through list and download pfams
    for mgya_id in mgya_list:
        print(f"extracting pfams for {mgya_id}")
        sample_id, assembly_id, pfams_list = download_pfam_from_mgya(mgya_id)

        time.sleep(2)

        samples_pfam_dict[mgya_id] = pfams_list

    return samples_pfam_dict


def download_pfam_clan(version="32.0"):

    # only proceed if pfam version is in list
    pfam_versions = [
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

    pfam_clan_df = []

    if version in pfam_versions:
        try:
            # download pfam clans file of default version 32
            response = requests.get(
                f"https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{version}/Pfam-A.clans.tsv.gz"
            )
            if response.status_code == 200:
                with gzip.open(BytesIO(response.content), "rt") as file:
                    pfam_clan_df = pd.read_csv(
                        file,
                        sep="\t",
                        names=[
                            "Protein_family",
                            "class",
                            "name",
                            "target.name",
                            "Description",
                        ],
                    )

                    pfam_clan_df.rename(
                        columns={"Protein_family": "PfamID"}, inplace=True
                    )

        except requests.exceptions.RequestException as e:
            print(
                f"Pfam clan download failed due to error: {e}. Please try again later"
            )

    else:
        print(
            "This version of pfam is not compatible with the mgnify2mimic_inputs pipeline."
        )
        print(
            "Compatible versions: 32.0, 33.0, 33.1, 34.0, 35.0, 36.0, 37.0, 37.1, 37.2"
        )
        sys.exit(1)

    return pfam_clan_df


def generate_pfam_metagenome_vector(pfam_clan_df, samples_pfam_dict, pfam_vector_file):

    # get list of pfamIDs
    pfam_clan_ids = list(pfam_clan_df["PfamID"])

    # assign pfam dict
    pfam_vector_df = pfam_clan_df[["PfamID"]]

    title = ["PfamID"]
    title.extend(list(samples_pfam_dict.keys()))
    title = "\t".join(title) + "\n"

    print("Generating pfam vector for metagenome")
    try:
        with open(pfam_vector_file, "w") as rfile:
            rfile.write(title)
            for pfam in pfam_clan_ids:
                rline = [
                    "1" if pfam in samples_pfam_dict[mgya] else "0"
                    for mgya in samples_pfam_dict
                ]
                rline.insert(0, pfam)

                # print (rline)

                rline = "\t".join(rline) + "\n"

                rfile.write(rline)

    except FileNotFoundError:
        print(
            f"metagenome file path: {pfam_vector_file} is not correct. Please check check filepath and retry command."
        )

    return pfam_vector_df


# for a given biome, download the functional profiles from the associated catalogue and parse the pfams into a vector
def download_genome_functions(biome, download_directory):

    biomes_dict = {
        "chicken-gut": "v1.0.1",
        # "cow-rumen": "v1.0.1",
        "honeybee-gut": "v1.0.1",
        # "human-gut":"v2.0.2",
        # "human-oral":"v1.0.1",
        "human-vaginal": "v1.0",
        "marine": "v2.0",
        "mouse-gut": "v1.0",
        "non-model-fish-gut": "v2.0",
        # "pig-gut":"v1.0",
        "sheep-rumen": "v1.0",
        "zebrafish-fecal": "v1.0",
    }

    biomes = ", ".join(biomes_dict.keys())

    # if biome is not in biomes_dict (i.e. pangenomes functions doesn't exist at the moment for the specified catalogue)
    if biome not in biomes_dict:
        print(
            f"*** specified biome: no pfam files exist for the {biome} catalogue. Please use one from the following: {biomes}. ***"
        )
    else:
        catalogue_version = biomes_dict[biome]

    ftp_url = (
        f"https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/{biome}/"
        f"{catalogue_version}/pangenome_functions/functional_profiles.tar.gz"
    )

    # downloaded functional folder path
    download_path = os.path.join(download_directory, "functional_profiles.tar.gz")

    # directory of unzipped functional profiles
    directory = "functional_profiles_DB"
    functional_profiles_directory = os.path.join(download_directory, directory)

    try:
        with urllib.request.urlopen(ftp_url) as response:
            if response.status == 200:
                print(
                    f"*** pfams for {biome} catalogue is accessible. Proceeding with download... ***"
                )
                urllib.request.urlretrieve(ftp_url, download_path)

                # extract tar.gz file to the same location and remove
                print(f"Extracting {download_path}...")
                with tarfile.open(download_path, "r:gz") as tar:
                    tar.extractall(
                        os.path.dirname(download_path)
                    )  # Extracts to the same location
                print(
                    f"*** Extraction complete. Files are in '{os.path.dirname(download_path)}'. ***"
                )

                os.remove(download_path)
            else:
                print(f"URL returned status code: {response.status}.")
    except urllib.error.HTTPError as e:
        print(f"HTTP error: {e.code} - {e.reason}")
    except urllib.error.URLError as e:
        print(f"URL error: {e.reason}")
    except Exception as e:
        print(f"Unexpected error: {e}")

    return functional_profiles_directory


# extact pfams from cluster files
def extract_pfams_from_cluster(cluster_file):
    pfam_accessions = []

    genome_df = pd.read_csv(cluster_file, sep="\t", header=1)

    for line in list(genome_df["pfam"].unique()):

        if "," in line:
            line = line.split(",")
            pfam_accessions.extend(line)
        elif line != "-":
            pfam_accessions.append(line)

    return pfam_accessions


# extract pfams from each genome in the functional_profile_folder (biome)
# scripts to generate pfam vector from metagenome
def extract_pfams2dict(pfam_clan_df, functional_profiles_folder):

    functional_profiles = os.listdir(functional_profiles_folder)

    pfam_dict = {pfam: [] for pfam in list(pfam_clan_df["PfamID"].unique())}

    for file in functional_profiles:

        genome_id = file.split("_")[0]
        cluster_file = os.path.join(functional_profiles_folder, file)

        pfam_accessions = extract_pfams_from_cluster(cluster_file)

        # add genome_id to dict
        for pfam in pfam_accessions:
            if pfam in pfam_dict:
                pfam_dict[pfam].append(genome_id)

    # remove keys with empty lists
    pfam_dict = {k: v for k, v in pfam_dict.items() if v != []}

    return pfam_dict


# geenrate pfam vector from folder of catalogue pfam files
def genome_pfam2vector(
    pfam_clan_df, genome_func_profile_path, result_genome_matrix_file
):

    # # get list of pfams (pfam-clans) and list of proteins
    pfam_list = list(pfam_clan_df["PfamID"])

    functional_files = os.listdir(genome_func_profile_path)

    genomes_list = [file.rstrip("_clstr.tsv") for file in functional_files]

    title = ["PfamID"]
    title.extend(genomes_list)
    title = "\t".join(title) + "\n"

    # reads pfams from each genome into dict
    genome_pfam_dict = extract_pfams2dict(pfam_clan_df, genome_func_profile_path)

    print("*** Generating pfam vector for genomes ***")

    with open(result_genome_matrix_file, "w") as rfile:
        rfile.write(title)
        for pfam in pfam_list:
            if pfam in genome_pfam_dict:
                rline = [
                    "1" if genome in genome_pfam_dict[pfam] else "0"
                    for genome in genomes_list
                ]
                rline.insert(0, pfam)
            else:
                rline = ["0" for genome in genomes_list]
                rline.insert(0, pfam)

            rline = "\t".join(rline) + "\n"

            rfile.write(rline)

    print("*** Removing intermediate files ***")
    # remove directory
    shutil.rmtree(genome_func_profile_path)

    return functional_files


# harmonizes the pfam vectors from metagenome sample and genome to one another so they are usable in MiMiC2
def vector_harmonisation(
    sample_vector_file, genome_vector_file, result_directory, pfam_clan_version="32.0"
):

    # download pfam clan info. default-v32
    pfam_clan_df = download_pfam_clan(pfam_clan_version)

    # read genome vector to df
    metagenome_vector_df = pd.read_csv(
        sample_vector_file, sep="\t", header=0
    )  # .set_index('PfamID')

    # read genome_vector to df
    genome_vector_df = pd.read_csv(
        genome_vector_file, sep="\t", header=0
    )  # .set_index('PfamID')

    # merge with old vector_df and fill in na with '0' and remove pfams not found in pfam_clan_df
    new_metagenome_vector_df = pfam_clan_df[["PfamID"]].merge(
        metagenome_vector_df, how="left", on="PfamID"
    )
    new_genome_vector_df = pfam_clan_df[["PfamID"]].merge(
        genome_vector_df, how="left", on="PfamID"
    )

    # file path for vector files
    new_metagenome_vector_file = os.path.join(
        result_directory, "sample_metagenome_pfam_vector_harmonised.tsv"
    )
    new_genome_vector_file = os.path.join(
        result_directory, "genome_pfam_vector_harmonised.tsv"
    )
    pfam_clan_vector_file = os.path.join(
        result_directory, f"Pfam-A.clans-{pfam_clan_version}.tsv"
    )

    # write to tsv file
    new_metagenome_vector_df.to_csv(
        new_metagenome_vector_file, sep="\t", na_rep="0", index=False
    )
    new_genome_vector_df.to_csv(
        new_genome_vector_file, sep="\t", na_rep="0", index=False
    )

    # write pfam_clan to tsv file
    pfam_clan_df.to_csv(pfam_clan_vector_file, sep="\t", index=False)

    return new_metagenome_vector_df, new_genome_vector_df, pfam_clan_df


# check if the input filepath exists
def file_path_check(vector_file):
    # check if sample vector file exists
    try:
        with open(vector_file, "w"):
            file_exists = True

    except OSError as e:
        print(
            f"metagenome file path: {vector_file} is not correct. Please check check filepath and retry command."
        )
        print(e)
        file_exists = False

    return file_exists


def main():
    parser = argparse.ArgumentParser(
        prog="mgnify2mimic inputs",
        usage="mgnify2mimic_inputs --mgnify_study MGYGXXXXXX"
        "--metagenome_vector_prefix mgnify_metagenome_vector.tsv",
        description="The mgnify2mimic_inputs script is designed to support the conversion of MGnify assemblies analysis"
        " and genomes outputs into a format that is usable by the MiMiC. The script downloads all needed files directly "
        "from MGnify's backend and transforms it into binary matrix (presence|absence (1|0) of pfams vector files which "
        "can be used directly in MiMiC2 to estimate the microbial genomes that functionally represent a specific biom "
        "(synthetic community).",
        epilog="The generated vector files can be used to generate synthetic communities using the MiMiC2 program"
        " - https://github.com/thh32/MiMiC2.git",
    )
    parser.add_argument(
        "--mgnify_study",
        "-ms",
        type=str,
        help="MGnify study accession (e.g. MGYS00006156)",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--mgnify_analysis_list",
        "-ma",
        type=str,
        help="provide comma seperated list of MGnify assembly analysis IDs (e.g. MYA0001,MGYA0002,MGYA0003) "
        "or tsv file with line separated list of MGnify assembly analysis IDs",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--metagenome_vector_file",
        "-mfile",
        type=str,
        help="path for sample metagenome vector file. DEFAULT=mgnify_metagenome_vector.tsv",
        required=False,
    )
    parser.add_argument(
        "--genome_vector_file",
        "-gfile",
        type=str,
        help="path for genome vector file.",
        required=False,
    )
    parser.add_argument(
        "--biome",
        "-b",
        type=str,
        help="biome of selected mgnify genome catalogue (e.g. mouse-gut). Used for recreating genome vector file.",
        required=False,
    )
    parser.add_argument(
        "--result_directory",
        "-d",
        type=str,
        help="result directory to store harmonized vector files",
        default=os.getcwd(),
    )
    parser.add_argument(
        "--pfam_clan_version",
        "-p",
        type=str,
        help="pfam version to harmonise metagenome and genomes to. Default pfam 32.0",
        required=False,
        default="32.0",
    )
    parser.add_argument(
        "--harmonize",
        help="",
        required=False,
        default=False,
        action=argparse.BooleanOptionalAction,
    )

    if len(sys.argv) == 1:
        parser.print_help()
        # print (sys.argv)
        sys.exit(1)

    args = parser.parse_args()

    # assign variable
    mgnify_study_id = args.mgnify_study
    mgnify_analysis_list_obj = args.mgnify_analysis_list
    result_directory = args.result_directory
    pfam_clan_version = args.pfam_clan_version

    # assign prefixes for files
    # if args.metagenome_vector_prefix is not None:
    #     sample_pfam_vector_file = args.metagenome_vector_prefix+"_metagenome_pfam_vector.tsv"

    # if args.genome_vector_prefix is not None:
    #     genome_pfam_vector_file = args.genome_vector_prefix+"_genome_pfam_vector.tsv"

    # option 1: vectorise sample metagenomes from mgys_id
    if mgnify_study_id is not None:
        sample_pfam_vector_file = (
            "metagenome_pfam_vector.tsv"
            if args.metagenome_vector_file is None
            else args.metagenome_vector_file
        )

        if file_path_check(sample_pfam_vector_file) is True:
            pfam_clan_df = download_pfam_clan()
            samples_pfam_dict = download_pfams_from_mgys(mgnify_study_id)
            generate_pfam_metagenome_vector(
                pfam_clan_df, samples_pfam_dict, sample_pfam_vector_file
            )

    # option 2: vectorise sample metagenomes from individual mgya_id
    elif mgnify_analysis_list_obj is not None and mgnify_study_id is None:
        if "," in mgnify_analysis_list_obj:
            mgnify_analysis_list = mgnify_analysis_list_obj.split(",")
        elif Path(mgnify_analysis_list_obj).exists():
            try:
                with open(mgnify_analysis_list_obj, "r") as file:
                    mgnify_analysis_list = [
                        line.rstrip() for line in file if line.startswith("MGYA")
                    ]
            except OSError as e:
                print(e)
        else:
            print(
                "mgnify analysis list was incorrectly provided. Please provide with a list of MGnify analysis "
                "accessions (MGYA) of analysed assmeblies. "
            )
            sys.exit(1)

        sample_pfam_vector_file = (
            "metagenome_pfam_vector.tsv"
            if args.metagenome_vector_file is None
            else args.metagenome_vector_file
        )

        if file_path_check(sample_pfam_vector_file) is True:
            pfam_clan_df = download_pfam_clan()
            samples_pfam_dict = download_pfams_from_mgya_list(mgnify_analysis_list)
            generate_pfam_metagenome_vector(
                pfam_clan_df, samples_pfam_dict, sample_pfam_vector_file
            )

    # option 3 harmonize sample vector and genome vector to a given pfam version
    elif args.harmonize is True:
        sample_pfam_vector_file = args.metagenome_vector_file
        genome_pfam_vector_file = args.genome_vector_file

        if (
            Path(sample_pfam_vector_file).exists()
            and Path(genome_pfam_vector_file).exists()
        ):
            new_metagenome_vector_df, new_genome_vector_df, pfam_clan_df = (
                vector_harmonisation(
                    sample_pfam_vector_file,
                    genome_pfam_vector_file,
                    result_directory,
                    pfam_clan_version,
                )
            )
        else:
            print(
                "Please input both the sample pfam vector file (e.g. metagenome_pfam_vector.tsv) and genome "
                "vector file (genome_vector.tsv) for harmonization. The program will harmonize the vectors to "
                "pfam v32.0 by default."
            )

    # option 4: for a given biome, generate pfam vector from mgnify genomes
    elif args.biome is not None:
        biome = args.biome.lower()
        genome_pfam_vector_file = (
            args.genome_vector_file
            if args.genome_vector_file is not None
            else "metagenome_pfam_vector.tsv"
        )

        # proceed if vector genome_pfam_vector_file is accessible
        if (
            file_path_check(genome_pfam_vector_file) is True
            and Path(result_directory).exists() is True
        ):
            pfam_clan_df = download_pfam_clan()
            genome_func_profile_path = download_genome_functions(
                biome, result_directory
            )
            genome_pfam2vector(
                pfam_clan_df, genome_func_profile_path, genome_pfam_vector_file
            )
        else:
            print(
                f"Result directory {result_directory} doesn't exist. Please check the path and retry command."
            )

    print("*** Pipeline complete ***")


if __name__ == "__main__":
    main()
