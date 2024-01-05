
import argparse
import sys

from Bio import Seq, SeqIO

def parse_args(argv=None):

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="Path to finalised primer list fasta file")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path")
    args = parser.parse_args(argv)
  
    _INPUT = args.input
    _SAMPLE = args.sample
    _OUTPUT = args.output

    return _INPUT, _SAMPLE, _OUTPUT

def main(argv=None):
    
    _INPUT, _SAMPLE, _OUTPUT = parse_args(argv)

    primers_dict = SeqIO.to_dict(SeqIO.parse(_INPUT, "fasta"))
    
    for primer_key in primers_dict.keys():

        primer = primers_dict[primer_key]
        primer_name = primer.name

        if "R" in primer_name:
            primers_dict[primer_key].seq = primer.seq.reverse_complement()

    SeqIO.write(primers_dict.values(), f"{_OUTPUT}/{_SAMPLE}_rev_comp_se_primers.fasta", "fasta")


if __name__ == "__main__":
    main()