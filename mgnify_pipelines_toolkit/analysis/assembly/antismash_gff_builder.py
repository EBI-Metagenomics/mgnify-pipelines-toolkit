
import argparse
from collections import defaultdict
import json

import pandas as pd

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="Input json from antiSMASH")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output")

    args = parser.parse_args()

    JSON = args.input
    OUTPUT = args.output

    return JSON, OUTPUT

def main():
    
    JSON, OUTPUT = parse_args()

    res_dict = defaultdict(list)

    res = ''
    with open(JSON, 'r') as json_data:
        res = json.load(json_data)

    antismash_ver = res['version']

    for record in res['records']:
        record_id = record['id']

        for row in record['features']:
            if row['type'] == 'CDS':

                if 'gene_kind' not in row['qualifiers'].keys():
                    continue

                res_dict['contig'].append(record_id)
                res_dict['version'].append(f'antiSMASH:{antismash_ver}')
                res_dict['type'].append('CDS')

                row_start = row['location'].split(':')[0][1:]
                res_dict['start'].append(row_start)
                
                row_end = row['location'].split(':')[1].split(']')[0]
                res_dict['end'].append(row_end)

                res_dict['score'].append(".")

                row_strand = row['location'].split('(')[1][0] # + or -
                res_dict['strand'].append(row_strand)

                res_dict['phase'].append(".")

                note_str = ""

                row_id = row['qualifiers']['ID'][0]
                note_str = f"ID={row_id};"
                
                row_type = row['qualifiers']['gene_kind'][0]
                note_str += f"as_type={row_type};"

                row_functions = ','.join(row['qualifiers']['gene_functions'])
                note_str += f"gene_functions={row_functions}"
                
                res_dict['attributes'].append(note_str)

    res_df = pd.DataFrame.from_dict(res_dict)
    # print(res_df)

    res_df.to_csv("ERZ21659799_test.gff", header=False, index=False, sep='\t')


    # print(res['records'][4604]['features'][1])
    # print(res['records'][4604]['features'][1].keys())
    # print(res['records'][25]['modules']['antismash.detection.genefunctions'])
    # print(res['records'][25]['modules'].keys())

if __name__ == "__main__":
    main()