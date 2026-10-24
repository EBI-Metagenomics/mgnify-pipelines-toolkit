# dbCAN fixtures

All fixtures are organised by dbCAN version. The files directly inside `5.1.2/`
are the original dbCAN 5.1.2 fixtures used by the substrate parsing and checksum
tests.

The `assembly-analysis-pipeline` subfolders contain verbatim, single-gene extracts for `ERZ19210871_100019_1` from the
`ERZ19210871` results generated with assembly-analysis pipeline versions 6.0.5
and 6.1.x, using dbCAN 5.1.2 and 5.2.9, respectively.

The examples are organised by dbCAN version and source pipeline:

```text
5.1.2/assembly-analysis-pipeline/
5.2.9/assembly-analysis-pipeline/
```

Both runs used mgnify-pipelines-toolkit 1.4.9 for the final GFF.

Each `assembly-analysis-pipeline` folder contains:

- `overview.tsv`: the header and matching row from
  `functional-annotation/dbcan/ERZ19210871_dbcan_overview.tsv.gz`.
- `sub_hmm.tsv`: the header and matching HMM hit from
  `functional-annotation/dbcan/ERZ19210871_dbcan_sub_hmm.tsv.gz`.
- `annotation_summary.gff`: the GFF version directive and matching CDS from
  `annotation-summary/ERZ19210871_annotation_summary.gff.gz`.

The overview EC field changes from
`3.2.1.-:8|3.2.1.113:1|3.2.1.114:1|3.2.1.24:1` in dbCAN 5.1.2 to
`3.2.1.24:1;3.2.1.113:1;3.2.1.114:1;3.2.1.-:8` in dbCAN 5.2.9.
The four ECs are the same; their order and separators differ. The dbCAN 5.2.9 overview
also adds a `Substrate` column.

The recorded assembly-analysis pipeline 6.0.5 GFF preserves all four ECs.
The recorded assembly-analysis pipeline 6.1.x GFF contains only
`eC_number=3.2.1.24`: the old parser split on `|` and then took the text before
the first `:`, losing the other ECs in the semicolon-separated group.
These GFF files document the original outputs, including that bug; they are
not expected outputs for the corrected parser.

The regression test uses each recorded CDS for coordinates and its HMM hit for
substrates, then checks that the corrected parser exports all four ECs in
overview order, without evidence counts, while retaining the family and substrate.
