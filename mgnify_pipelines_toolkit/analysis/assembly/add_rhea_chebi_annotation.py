import argparse
import gzip
import logging
import re
from pathlib import Path
import sys

import pandas as pd
import requests


def main(input, output, rhea2chebi, up2rhea):
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
    diamond_df = pd.read_csv(input, sep='\t', usecols=[0, 1], header=0)
    diamond_df.columns = ['uniref90_ID', 'contig_name']
    diamond_df['Uniref90_rep'] = diamond_df['uniref90_ID'].str.split('_').str[1]

    logging.info(f"Step 2/4: Adding RHEA IDs based on provided file {up2rhea.resolve()}")
    up2rhea_df = pd.read_csv(up2rhea, sep='\t', usecols=[0, 2])
    diamond_df = diamond_df.merge(up2rhea_df, left_on='Uniref90_rep', right_on="Entry", how='left')
    diamond_df['Rhea ID'] = diamond_df['Rhea ID'].str.split()
    diamond_df = diamond_df.explode('Rhea ID')

    logging.info(f"Step 3/4: Adding CHEBI reactions based on provided file {rhea2chebi.resolve()}")
    rhea2chebi_df = pd.read_csv(rhea2chebi, sep='\t')
    diamond_df = diamond_df.merge(rhea2chebi_df[['ENTRY', 'EQUATION']], left_on='Rhea ID', right_on='ENTRY', how='left')
    diamond_df = diamond_df.drop(columns=['ENTRY', 'Entry'])

    logging.info(f"Step 4/4: Saving output table to {output.resolve()}")
    diamond_df.to_csv(output, sep="\t", index=False)

    logging.info("Processed successfully. Exiting.")


def download_and_convert_to_tsv(download_path: Path):
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
        help="Output TSV file with columns: contig_id, protein_id, UniRef90 cluster, Rhea ID, CHEBI reaction participants"
    )
    parser.add_argument(
        "--rhea2chebi",
        default=None,
        type=Path,
        help="File that maps RHEA IDs to CHEBI. Will be downloaded if not provided",
    ) 
    parser.add_argument(
        "--up2rhea",
        required=True,
        type=Path,
        help="File that maps UniProt IDs to Rhea. Must contain at least 2 columns 'Entry' and 'Rhea ID'",
    ) 
    args = parser.parse_args()
    main(args.input, args.output, args.rhea2chebi, args.up2rhea)
