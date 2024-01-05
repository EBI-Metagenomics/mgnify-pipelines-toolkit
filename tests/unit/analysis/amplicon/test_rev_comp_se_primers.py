
import hashlib

import pytest

from mgnify_pipelines_toolkit.analysis.amplicon.rev_comp_se_primers import parse_args
from mgnify_pipelines_toolkit.analysis.amplicon.rev_comp_se_primers import main

def test_right_arguments():
    args = [
        "-i",
        "input",
        "-s",
        "sample",
        "-o",
        "output"
    ]
    parse_args(args)

    assert True

def test_wrong_arguments():
    with pytest.raises(SystemExit):
        parse_args(["--chromosome", "X"])
    
def test_argparse():
    args = [
        "-i",
        "input",
        "-s",
        "sample",
        "-o",
        "output"
    ]

    assert parse_args(args) == ("input", "sample", "output")

def test_main_output():

    args = [
        "-i",
        "./tests/fixtures/sequences/f_and_r_primers.fasta",
        "-s",
        "sample",
        "-o",
        "./"
    ]

    main(args)
    out_md5 = hashlib.md5(open('./sample_rev_comp_se_primers.fasta','rb').read()).hexdigest()
    assert out_md5 == "aa3443aec8e07139360a77c35f69e326"