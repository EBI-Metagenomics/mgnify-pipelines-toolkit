import argparse
from itertools import product
from pathlib import Path

from Bio import SeqIO


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_primers",
        required=True,
        type=str,
        help="Input primers to generate permutations for due to IUPAC ambiguous codes",
    )
    parser.add_argument("-p", "--prefix", required=True, type=str, help="Output prefix")

    args = parser.parse_args()

    input_path = args.input_primers
    prefix = args.prefix

    return input_path, prefix


def permute_seq(seq):

    ambiguous_bases_dict = {
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "T", "G"],
    }

    seq_template = []

    for base in seq:
        if base in ["A", "C", "T", "G"]:
            seq_template.append(base)
        else:
            seq_template.append(ambiguous_bases_dict[base])

    seq_permutations = []
    for combo in product(*seq_template):
        seq_permutations.append("".join(combo))

    return seq_permutations


def make_primer_permutations(primers_dict, prefix):

    with open(f"{prefix}_permuted_primers.fasta", "w") as fw:
        for primer_name, seq in primers_dict.items():

            primer_seq = seq.seq
            fw.write(f">{primer_name}\n{primer_seq}\n")

            if primer_name == "F_auto" or primer_name[-1] == "F":
                strand = "F"
            elif primer_name == "R_auto" or primer_name[-1] == "R":
                strand = "R"

            seq_permutations = permute_seq(primer_seq)

            for counter, permuted_seq in enumerate(seq_permutations, 1):
                variant_name = f"{primer_name}_variant_{counter}_{strand}"
                fw.write(f">{variant_name}\n{permuted_seq}\n")


def main():

    input_path, prefix = parse_args()
    primers_dict = SeqIO.to_dict(SeqIO.parse(Path(input_path), "fasta"))
    make_primer_permutations(primers_dict, prefix)


if __name__ == "__main__":
    main()
