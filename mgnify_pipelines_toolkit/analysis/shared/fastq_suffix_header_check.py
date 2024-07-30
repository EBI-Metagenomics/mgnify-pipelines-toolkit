
import argparse
from collections import defaultdict
import gzip
import json
import logging

logging.basicConfig(level=logging.DEBUG)

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fwd", required=True, type=str, help="Input forward read headers file (PE) OR SE read file")
    parser.add_argument("-r", "--rev", required=False, type=str, help="Input reverse read headers file (PE)")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output")

    args = parser.parse_args()

    FWD = args.fwd
    REV = args.rev
    SAMPLE = args.sample
    OUTPUT = args.output

    return FWD, REV, SAMPLE, OUTPUT

def choose_open_func(file_path):

    open_func = open

    if file_path[-2:] == "gz":
        open_func = gzip.open

    return open_func

def main():
    
    FWD, REV, SAMPLE, OUTPUT = parse_args()

    files_to_parse = []

    if "_1" in FWD:
        if REV == None:
            logging.error("No reverse file given, yet given forward file has the \"_1\" suffix implying it's paired-end. " +
                           "Either supply the reverse file, or supply a single-end file.")
            exit(1)
        elif "_2" not in REV:
            logging.error("The expected suffix \"_2\" for a supplied reverse file is missing. Please verify your inputs.")
            exit(1)
        else:
            files_to_parse = [FWD, REV]

    else:
        files_to_parse = [FWD]

    open_func = choose_open_func(FWD) # Choose between gzip.open() and open() by checking the file extension
    reads_with_err = defaultdict(list)

    for file in files_to_parse:

        header_str = ""

        if "_1" in file:
            header_str = "/1"
        elif "_2" in file:
            header_str = "/2"
        else:
            header_str = "/1" # SE files still have "/1" in the headers
        
        for counter, line in enumerate(open_func(file)):

            if counter % 4 == 0: # Only do stuff every four lines to hit the header
                line = line.decode("ascii").strip()
                curr_read_strand = line[-2:]

                if curr_read_strand != header_str:
                    reads_with_err[file].append(line)
                    reads_with_err["total"].append(1)

    if len(reads_with_err) != 0:

        num_of_reads_with_err = len(reads_with_err["total"])
        reads_with_err["total"] = num_of_reads_with_err

        logging.error(f"Found {num_of_reads_with_err} reads with header strands that don't match file suffix. See log file at {OUTPUT}/{SAMPLE}_suffix_header_err.json")
        
        with open(f"{OUTPUT}/{SAMPLE}_suffix_header_err.json", 'w') as fw: # Writes JSON file containing the headers of reads with errors
            json.dump(reads_with_err, fw)

    else:
        with open(f"{OUTPUT}/{SAMPLE}_suffix_header_err.json", 'w') as fw: # Creates an empty file if there are no errors
            print("No errors.")

if __name__ == "__main__":
    main()

