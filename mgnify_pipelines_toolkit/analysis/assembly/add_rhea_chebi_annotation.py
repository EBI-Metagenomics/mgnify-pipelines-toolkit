import argparse
import gzip
import logging
import re
from pathlib import Path
import sys

import pandas as pd
import requests


def main(input: Path, output: Path, rhea2chebi: Path, up2rhea: Path):
    logging.info("Step 0/4: Checking Rhea-CHEBI mapping file...")
    if not rhea2chebi:
        logging.info("Rhea-CHEBI mapping not provided. Starting download...")
        download_path = Path("./data/")
        try:
            rhea2chebi = download_and_convert_to_tsv(download_path)
        except:
            logging.info("Failed to download file. Exiting...")
            sys.exit(1)
        finally:
            logging.info(f"File downloaded successfully and saved to {rhea2chebi.resolve()}")

    logging.info(f"Step 1/4: Reading input file {input.resolve()}")
    diamond_df = pd.read_csv(input, sep='\t', usecols=['uniref90_ID', 'contig_name'])
    diamond_df[['protein_id', 'contig']] = diamond_df['contig_name'].str.split('-', n=1, expand=True)
    diamond_df['uniref90_rep'] = diamond_df['uniref90_ID'].str.split('_').str[1]

    logging.info(f"Step 2/4: Adding RHEA IDs based on provided file {up2rhea.resolve()}")
    up2rhea_df = pd.read_csv(up2rhea, sep='\t', usecols=['Entry','Rhea ID'])
    up2rhea_df.columns = ['unirefKB_id', 'rhea_id']
    diamond_df = diamond_df.merge(up2rhea_df, left_on='uniref90_rep', right_on="unirefKB_id", how='left')
    diamond_df['rhea_id'] = diamond_df['rhea_id'].str.split()
    diamond_df = diamond_df.explode('rhea_id')

    logging.info(f"Step 3/4: Adding CHEBI reactions based on provided file {rhea2chebi.resolve()}")
    rhea2chebi_df = pd.read_csv(rhea2chebi, sep='\t')
    rhea2chebi_df.columns = ['rhea_id', 'definition', 'chebi_reaction', 'enzyme_id']
    diamond_df = diamond_df.merge(rhea2chebi_df[['rhea_id', 'chebi_reaction']], on='rhea_id', how='left')
    diamond_df = diamond_df.drop(columns=['contig_name', 'unirefKB_id', 'uniref90_rep'])

    logging.info(f"Step 4/4: Saving output table to {output.resolve()}")
    diamond_df = diamond_df[['contig', 'protein_id', 'uniref90_ID', 'rhea_id', 'chebi_reaction']]
    diamond_df.to_csv(output, sep="\t", index=False)

    logging.info("Processed successfully. Exiting.")


def download_and_convert_to_tsv(download_path: Path)  -> Path:
    url = "https://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz"
    
    # Ensure the download directory exists
    download_path.mkdir(parents=True, exist_ok=True)
    
    # Define paths for the created files
    gz_file_path = download_path / "rhea-reactions.txt.gz"
    tsv_file_path = download_path / "rhea-reactions.tsv"
    
    # Download the .gz file
    response = requests.get(url)
    response.raise_for_status()
    with gz_file_path.open('wb') as file:
        file.write(response.content)
    
    # Process the .gz file and convert it to TSV
    with gzip.open(gz_file_path, 'rt') as file_in, open(tsv_file_path, 'w') as file_out:
        input_text = file_in.read()
        pattern = re.compile(r'ENTRY\s+(.*?)\nDEFINITION\s+(.*?)\nEQUATION\s+(.*?)\n(?:ENZYME\s+(.*?))?\n*', re.DOTALL)
        matches = pattern.findall(input_text)
        file_out.write('ENTRY\tDEFINITION\tEQUATION\tENZYME\n')
        for match in matches:
            entry, definition, equation, enzyme = match
            enzyme = enzyme.strip() if enzyme else ''
            file_out.write(f'{entry}\t{definition}\t{equation}\t{enzyme}\n')

    # Clean up the downloaded .gz file
    gz_file_path.unlink()
    return tsv_file_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser("""
                                     Use diamond output file to create a table with Rhea and CHEBI 
                                     reaction annotation for every protein.
                                     """)
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        type=Path, 
        help="DIAMOND results file"
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        type=Path, 
        help="Output TSV file with columns: contig_id, protein_id, UniRef90 cluster, rhea_ids, CHEBI reaction participants"
    )
    parser.add_argument(
        "--rhea2chebi",
        default=None,
        type=Path,
        help="File that maps rhea_idss to CHEBI. Will be downloaded if not provided",
    ) 
    parser.add_argument(
        "--up2rhea",
        required=True,
        type=Path,
        help="File that maps UniProt IDs to Rhea. Must contain at least 2 columns 'Entry' and 'Rhea ID'",
    ) 
    args = parser.parse_args()
    main(args.input, args.output, args.rhea2chebi, args.up2rhea)
