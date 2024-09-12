from collections import defaultdict
import os
import requests

import pandas as pd
from tqdm import tqdm

url = "https://www.ebi.ac.uk/metagenomics/api/v1"
headers = {"Accept": "application/json"}


def date_parser(date):
    year, month, day = date.split("-")
    return (year, month, day)


def get_metadata_from_run_acc(run_acc):

    res_run = requests.get(f"{url}/runs/{run_acc}", headers=headers)
    sample_acc = res_run.json()["data"]["relationships"]["sample"]["data"]["id"]

    res_samples = requests.get(f"{url}/samples/{sample_acc}/metadata", headers=headers)

    res_samples_studies = requests.get(
        f"{url}/samples/{sample_acc}/studies", headers=headers
    )
    study_acc = res_samples_studies.json()["data"][0]["attributes"][
        "secondary-accession"
    ]

    fields_of_interest = {
        "geographic location (longitude)": "decimalLongitude",
        "geographic location (latitude)": "decimalLatitude",
        "collection date": "collection_date",
        "depth": "Depth",
    }

    res_dict = defaultdict(list)
    res_dict["SampleID"].append(sample_acc)
    res_dict["StudyID"].append(study_acc)
    res_dict["RunID"].append(run_acc)

    present_fields = []
    for val in res_samples.json()["data"]:
        field = val["attributes"]["key"]

        if field in fields_of_interest:
            present_fields.append(field)
            field_val = val["attributes"]["value"]

            if field != "collection date":
                res_dict[fields_of_interest[field]].append(field_val)
            else:
                year, month, day = date_parser(field_val)
                res_dict["Year"].append(year)
                res_dict["Month"].append(month)
                res_dict["Day"].append(day)

    missing_fields = [
        field for field in fields_of_interest if field not in present_fields
    ]

    if len(missing_fields) > 0:
        for field in missing_fields:
            if field != "collection_date":
                res_dict[field].append("NA")
            else:
                res_dict["Year"].append("NA")
                res_dict["Month"].append("NA")
                res_dict["Day"].append("NA")

    res_df = pd.DataFrame.from_dict(res_dict)

    return res_df


def get_all_runs_from_studies(data_path, study_list):

    all_runs = []

    for study in study_list:
        runs_path = f"{data_path}/{study}"
        all_paths = os.listdir(runs_path)

        all_runs.extend([file for file in all_paths if "ERR" in file or "SRR" in file])

    return all_runs


def get_all_asv_data_from_runs(data_path, study_list):

    asv_dict = {}

    no_asv_runs = []
    for study in study_list:

        all_runs = []

        runs_path = f"{data_path}/{study}"
        all_paths = os.listdir(runs_path)
        all_runs.extend([file for file in all_paths if "ERR" in file or "SRR" in file])

        for run in all_runs:
            if os.path.isdir(f"{runs_path}/{run}/asv/"):
                tax_file_path = f"{runs_path}/{run}/asv/{run}_DADA2-SILVA_asv_taxa.tsv"
                count_file_path = f"{runs_path}/{run}/asv/16S-V4-V5/{run}_16S.V4-V5_asv_read_counts.tsv"

                run_tax_df = pd.read_csv(tax_file_path, sep="\t")
                count_df = pd.read_csv(count_file_path, sep="\t")

                merged_df = count_df.merge(run_tax_df, left_on="asv", right_on="ASV")
                merged_df.pop("ASV")

                run_col = [run] * len(merged_df)
                merged_df["RunID"] = run_col
                asv_dict[run] = merged_df
            else:
                no_asv_runs.append(run)

    print(len(no_asv_runs))
    return asv_dict


def get_all_metadata_from_runs(runs):

    run_metadata_dict = defaultdict(dict)

    for run in tqdm(runs):
        res_df = get_metadata_from_run_acc(run)
        run_metadata_dict[run] = res_df

    return run_metadata_dict


def cleanup_taxa(df):

    df.pop("Kingdom")
    cleaned_df = df.rename(columns={"Superkingdom": "Kingdom", "asv": "ASVID"})

    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

    for rank in ranks:
        cleaned_df[rank] = cleaned_df[rank].apply(
            lambda x: x.split("__")[1] if pd.notnull(x) else "NA"
        )

    cleaned_df = cleaned_df[
        [
            "ASVID",
            "StudyID",
            "SampleID",
            "RunID",
            "decimalLongitude",
            "Year",
            "Month",
            "Day",
            "decimalLatitude",
            "depth",
            "Kingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
        ]
    ]

    return cleaned_df


def main():

    data_path = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/atlanteco_v6/analysis"
    # data_path = "/hps/software/users/rdf/metagenomics/service-team/users/chrisata/mgnify-pipelines-toolkit/metadata_asv_test_data"
    study_list = ["ERP107119", "ERP129129", "SRP279131"]

    asv_dict = get_all_asv_data_from_runs(data_path, study_list)
    all_runs = get_all_runs_from_studies(data_path, study_list)
    run_metadata_dict = get_all_metadata_from_runs(all_runs)

    all_merged_df = []

    for run in all_runs:
        if run in asv_dict.keys() and run in run_metadata_dict.keys():
            run_asv_data = asv_dict[run]
            run_metadata = run_metadata_dict[run]
            run_merged_result = run_metadata.merge(run_asv_data, on="RunID")
            all_merged_df.append(run_merged_result)

    final_df = pd.concat(all_merged_df, ignore_index=True)

    final_df = cleanup_taxa(final_df)

    final_df.to_csv("all_merged.csv", index=False, na_rep="NA")


# TODO NEED TO FIX THE NON-NA EMPTYS (',,')


if __name__ == "__main__":
    main()
