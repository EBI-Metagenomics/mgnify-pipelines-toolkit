import logging
import pandas as pd
from mgnify_pipelines_toolkit.analysis.amplicon.mapseq_to_asv_table import parse_label, parse_mapseq, process_blank_tax_ends, mapseq_to_asv_table
from mgnify_pipelines_toolkit.constants.tax_ranks import (
    SHORT_SILVA_TAX_RANKS,
    SILVA_TAX_RANKS,
)


def test_parse_label():
    """
    Test that parse_label returns the correct short and long ranks for a given database label.
    """
    short_ranks, long_ranks = parse_label("DADA2-SILVA")
    assert short_ranks == [f"{rank}__" for rank in SHORT_SILVA_TAX_RANKS]
    assert long_ranks == SILVA_TAX_RANKS


def test_parse_mapseq():
    """
    Test that parse_mapseq correctly converts MAPseq taxonomy strings into a structured DataFrame.
    """
    short_ranks = ["k__", "p__", "c__"]
    long_ranks = ["Kingdom", "Phylum", "Class"]
    mseq_data = {"query": ["ASV1", "ASV2"], "tax": ["k__Bacteria;p__Firmicutes", "k__Archaea"]}
    mseq_df = pd.DataFrame(mseq_data)

    res_df = parse_mapseq(mseq_df, short_ranks, long_ranks)

    assert len(res_df) == 2
    assert list(res_df.columns) == ["ASV", "Kingdom", "Phylum", "Class"]

    assert res_df.iloc[0]["ASV"] == "ASV1"
    assert res_df.iloc[0]["Kingdom"] == "k__Bacteria"
    assert res_df.iloc[0]["Phylum"] == "p__Firmicutes"
    assert res_df.iloc[0]["Class"] == "c__"

    assert res_df.iloc[1]["ASV"] == "ASV2"
    assert res_df.iloc[1]["Kingdom"] == "k__Archaea"
    assert res_df.iloc[1]["Phylum"] == "p__"
    assert res_df.iloc[1]["Class"] == "c__"


def test_process_blank_tax_ends():
    """
    Test that process_blank_tax_ends correctly replaces trailing empty ranks with 'NA'.
    """
    ranks = ["k__", "p__", "c__"]
    df = pd.DataFrame(
        {
            "ASV": ["ASV1", "ASV2", "ASV3"],
            "Kingdom": ["k__Bacteria", "k__Archaea", "k__"],
            "Phylum": ["p__Firmicutes", "p__", "p__"],
            "Class": ["c__", "c__", "c__"],
        }
    )

    processed_df = process_blank_tax_ends(df, ranks)

    # ASV1: Kingdom=Bacteria, Phylum=Firmicutes, Class=c__ (empty) -> Class should be NA
    assert processed_df.iloc[0]["Class"] == "NA"
    assert processed_df.iloc[0]["Phylum"] == "p__Firmicutes"

    # ASV2: Kingdom=Archaea, Phylum=p__ (empty), Class=c__ (empty) -> Phylum and Class should be NA
    assert processed_df.iloc[1]["Class"] == "NA"
    assert processed_df.iloc[1]["Phylum"] == "NA"
    assert processed_df.iloc[1]["Kingdom"] == "k__Archaea"

    # ASV3: All empty -> Phylum and Class should be NA, Kingdom should stay k__
    assert processed_df.iloc[2]["Class"] == "NA"
    assert processed_df.iloc[2]["Phylum"] == "NA"
    assert processed_df.iloc[2]["Kingdom"] == "k__"


def test_mapseq_to_asv_table_empty_input(tmp_path, monkeypatch, caplog):
    """
    Test that mapseq_to_asv_table handles an empty input file correctly.
    Similar to tests in test_dwc_summary_generator.py.
    """
    caplog.set_level(logging.INFO)

    input_file = tmp_path / "input.mseq"
    input_file.touch()

    sample = "test_sample"
    label = "DADA2-SILVA"

    # Need to change directory because the function currently hardcodes "./" output
    monkeypatch.chdir(tmp_path)

    mapseq_to_asv_table(input_file, label, sample)

    output_file = tmp_path / f"{sample}_{label}_asv_taxa.tsv"
    assert output_file.exists()
    assert output_file.stat().st_size == 0
    assert f"Run {sample}'s mapseq file {input_file} exists but is empty. Skipping." in caplog.text
