# dbCAN examples from pipeline results 6.0 and 6.1

These are verbatim, single-gene extracts for `ERZ19210871_100019_1` from the
`ERZ19210871` result folders labelled `6.0` and `6.1`. The labels identify the
result folders, not dbCAN versions: the runs used dbCAN 5.1.2 and 5.2.9,
respectively. Both runs used mgnify-pipelines-toolkit 1.4.9 for the final GFF.

Each folder contains:

- `overview.tsv`: the header and matching row from
  `functional-annotation/dbcan/ERZ19210871_dbcan_overview.tsv.gz`.
- `sub_hmm.tsv`: the header and matching HMM hit from
  `functional-annotation/dbcan/ERZ19210871_dbcan_sub_hmm.tsv.gz`.
- `annotation_summary.gff`: the GFF version directive and matching CDS from
  `annotation-summary/ERZ19210871_annotation_summary.gff.gz`.

The overview EC field changes from
`3.2.1.-:8|3.2.1.113:1|3.2.1.114:1|3.2.1.24:1` in 6.0 to
`3.2.1.24:1;3.2.1.113:1;3.2.1.114:1;3.2.1.-:8` in 6.1.
The four ECs are the same; their order and separators differ. The 6.1 overview
also adds a `Substrate` column.

The recorded 6.0 GFF preserves all four ECs. The recorded 6.1 GFF contains only
`eC_number=3.2.1.24`: the old parser split on `|` and then took the text before
the first `:`, losing the other ECs in the semicolon-separated group.
These GFF files document the original outputs, including that bug; they are
not expected outputs for the corrected parser.

The regression test uses each recorded CDS for coordinates and its HMM hit for
substrates, then checks that the corrected parser exports all four ECs in
overview order, without evidence counts, while retaining the family and substrate.
