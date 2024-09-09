from collections import defaultdict
import requests

import pandas as pd


def date_parser(date):
    year, month, day = date.split("-")
    return (year, month, day)


def main():

    url = "https://www.ebi.ac.uk/metagenomics/api/v1"
    headers = {"Accept": "application/json"}
    run_acc = "ERR2356286"
    # run_acc = "ERR11968495"

    res_run = requests.get(f"{url}/runs/{run_acc}", headers=headers)
    sample_acc = res_run.json()["data"]["relationships"]["sample"]["data"]["id"]

    res_samples = requests.get(f"{url}/samples/{sample_acc}/metadata", headers=headers)

    fields_of_interest = {
        "geographic location (longitude)": "decimalLongitude",
        "geographic location (latitude)": "decimalLatitude",
        "collection date": "collection_date",
        "depth": "Depth",
    }

    res_dict = defaultdict(list)
    res_dict["SampleID"].append(sample_acc)
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
    print(res_df)


if __name__ == "__main__":
    main()
