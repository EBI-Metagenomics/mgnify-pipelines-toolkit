
import argparse
from collections import defaultdict
import json

import pandas as pd

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="Input json from antiSMASH")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output")

    args = parser.parse_args()

    JSON = args.input
    SAMPLE = args.sample
    OUTPUT = args.output

    return JSON, SAMPLE, OUTPUT

def main():
    
    JSON, SAMPLE, OUTPUT = parse_args()

    res_dict = defaultdict(list)

    res = ''
    with open(JSON, 'r') as json_data:
        res = json.load(json_data)

    antismash_ver = res['version']

    note_dict = defaultdict(str)

    for record in res['records']:
        record_id = record['id']

        for row in record['features']:
            if row['type'] == 'CDS':

                if 'antismash.detection.genefunctions' not in record['modules'].keys():
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
                note_str = f"ID={row_id}"
                
                if 'gene_kind' in row['qualifiers'].keys():
                    row_type = row['qualifiers']['gene_kind'][0]
                    note_str += f";as_type={row_type}"

                    # row_functions = ','.join(row['qualifiers']['gene_functions'])
                    # note_str += f"gene_functions={row_functions};"

                else:
                    note_str += f";as_type=other"
                
                note_dict[row_id] += note_str

        if 'antismash.detection.hmm_detection' in record['modules'].keys():
            temp_row = record['modules']['antismash.detection.hmm_detection']['rule_results']['cds_by_protocluster']
            if len(temp_row) > 0:
                for row in temp_row[0][1]:
                    if 'cds_name' in row.keys():
                        row_id = row['cds_name']
                        row_as_clusters = ','.join(list(row['definition_domains'].keys()))
                        note_dict[row_id] += f';as_gene_clusters={row_as_clusters}'

        if 'antismash.detection.genefunctions' in record['modules'].keys():
            for tool in record['modules']['antismash.detection.genefunctions']['tools']:
                note_str = ""
                if tool['tool'] == 'smcogs' and len(tool['best_hits']) > 0:
                    for best_hit in tool['best_hits']:
                        hit_id = tool['best_hits'][best_hit]['hit_id'].split(':')[0]
                        hit_desc = tool['best_hits'][best_hit]['hit_id'].split(':')[1].replace(' ', '_')
                        score = tool['best_hits'][best_hit]['bitscore']
                        eval = tool['best_hits'][best_hit]['evalue']

                        if "as_notes" in note_str:
                            note_str += f'smCOG%3A%20{hit_id}%A{hit_desc}%20%28Score%3A%20{score}%3B%20E-value%3A%20{eval}%29%3B'
                        else:
                            note_str += f';as_notes=smCOG%3A%20{hit_id}%A{hit_desc}%20%28Score%3A%20{score}%3B%20E-value%3A%20{eval}%29%3B'

                        note_dict[best_hit] += note_str
                        break




    note_col = note_dict.values()
    res_dict['attributes'] = note_col

    res_df = pd.DataFrame.from_dict(res_dict)
    res_df.to_csv(f"{SAMPLE}_test.gff", header=False, index=False, sep='\t')

    # print(res['records'][25]['modules']['antismash.modules.clusterblast']['knowncluster']['results'][0]['ranking'][0])
    # print(res['records'][25]['modules']['antismash.detection.hmm_detection']['rule_results']['cds_by_protocluster'][0][1][0]['cds_name'])
    # print(res['records'][25]['modules']['antismash.detection.hmm_detection']['rule_results']['cds_by_protocluster'][0][1][0]['definition_domains'].keys())
    # print(res['records'][25]['modules']['antismash.detection.hmm_detection']['rule_results']['cds_by_protocluster'][0][1][0].keys())

if __name__ == "__main__":
    main()