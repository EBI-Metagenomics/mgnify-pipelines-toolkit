import argparse
import os

from mgnify_pipelines_toolkit.analysis.assembly.go_utils import parse_ips_file

description = "Go slim pipeline."
parser = argparse.ArgumentParser(description=description)
parser.add_argument(
    "-i", "--ips_input", help="InterProScan result file.", required=True
)
parser.add_argument("-o", "--output", help="GO summary output file.", required=True)
args = parser.parse_args()

IPS_INPUT = args.ips_input
OUTPUT = args.output


def create_gaf_file(gaf_input_file_path, go_id_set):
    """
    :param gaf_input_file_path:
    :param go2proteinDict:
    :return: nothing
    """
    with open(gaf_input_file_path, "w") as file:
        file.write("!gaf-version: 2.1\n")
        file.write("!Project_name: EBI Metagenomics\n")
        file.write("!URL: http://www.ebi.ac.uk/metagenomics\n")
        file.write("!Contact Email: metagenomics-help@ebi.ac.uk\n")
        for go_id in go_id_set:
            gaf_file_entry_line_str = "EMG\t{0}\t{1}\t\t{2}\tPMID:12069591\tIEA\t\t{3}\t\t\tprotein\ttaxon:1310605\t{4}\t{5}\t\t".format(
                go_id, "GO", go_id, "P", "20160528", "InterPro"
            )
            file.write("" + gaf_file_entry_line_str + "\n")


def main():

    if not os.stat(IPS_INPUT).st_size == 0:
        # Parse InterProScan result file; map protein accessions and GO terms
        print("Parsing InterProScan input: " + IPS_INPUT)
        go2protein_count_dict = parse_ips_file(IPS_INPUT)
        print("Finished parsing.")

        gaf_output_path = f"{OUTPUT}_ips_annotations.gaf"
        print("Creating GAF file: " + gaf_output_path)
        go_id_set = go2protein_count_dict.keys()
        create_gaf_file(gaf_output_path, go_id_set)
        print("Finished GAF file generation.")


if __name__ == "__main__":
    main()
