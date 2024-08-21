import argparse
from pathlib import Path
import requests
import gzip
import pandas as pd
import re


def main(input, output, rhea_reactions, up2rhea):
    if not rhea_reactions:
        download_path = Path("./data/")
        rhea_reactions = download_and_convert_to_tsv(download_path)

    diamond_df = pd.read_csv(input, sep='\t', usecols=[0, 1], header=0)
    diamond_df.columns = ['uniref90_ID', 'contig_name']

    # Step 2: Create a new column 'Uniref90_rep' by splitting 'uniref90_ID'
    diamond_df['Uniref90_rep'] = diamond_df['uniref90_ID'].str.split('_').str[1]

    # Step 3: Map 'Rhea_reactions' to 'diamond_df' based on 'Uniref90_rep'
    up2rhea_df = pd.read_csv(up2rhea, sep='\t', usecols=[0, 2])
    diamond_df = diamond_df.merge(up2rhea_df, left_on='Uniref90_rep', right_on="Entry", how='left')
    diamond_df['Rhea ID'] = diamond_df['Rhea ID'].str.split()
    diamond_df = diamond_df.explode('Rhea ID')

    # Step 4: Add CHEBI reaction 
    rhea2chebi_df = pd.read_csv(rhea_reactions, sep='\t')
    diamond_df = diamond_df.merge(rhea2chebi_df[['ENTRY', 'EQUATION']], left_on='Rhea ID', right_on='ENTRY', how='left')
    diamond_df = diamond_df.drop(columns=['ENTRY', 'Entry'])
    
    diamond_df.to_csv(output, sep="\t", index=False)


def download_and_convert_to_tsv(download_path: Path):
    # Define the URL and the download path
    url = "https://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz"
    
    # Ensure the download directory exists
    download_path.mkdir(parents=True, exist_ok=True)
    
    # Define paths for the files
    gz_file_path = download_path / "rhea-reactions.txt.gz"
    tsv_file_path = download_path / "rhea-reactions.tsv"
    
    # Download the .gz file
    response = requests.get(url)
    response.raise_for_status()
    
    # Save the file
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
            # Clean up enzyme field if it's empty
            enzyme = enzyme.strip() if enzyme else ''
            file_out.write(f'{entry}\t{definition}\t{equation}\t{enzyme}\n')

    # Clean up the downloaded .gz file
    gz_file_path.unlink()
    return tsv_file_path


if __name__ == "__main__":

    parser = argparse.ArgumentParser("")
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
        help=""
    )
    parser.add_argument(
        "--rhea_reactions",
        default=None,
        type=Path,
        help="File that maps RHEA id to CHEBI",
    ) 
    parser.add_argument(
        "--up2rhea",
        required=True,
        type=Path,
        help="File that maps UniProt id to Rhea",
    ) 
    args = parser.parse_args()
    main(args.input, args.output, args.rhea_reactions, args.up2rhea)