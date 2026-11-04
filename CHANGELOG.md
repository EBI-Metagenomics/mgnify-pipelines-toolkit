# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.5.6] - 2026-11-04

### Added

- Support for run_dbCAN 5.2.9 in the GFF generation script

### Changed

- Reject dbCAN overview rows with unexpected gene IDs instead of silently skipping them. This script is only compatible with ERZ and MGYG accessions ATM, we will revisit this in the future.

### Fixed

- Parse reordered dbCAN overview columns and mixed EC separators in dbCAN 5.2.9, preserving EC order and repeated accessions.
- Filter missing substrate tokens and preserve empty TSV fields when parsing dbCAN HMM results.
- Report missing overview columns and incomplete rows.

## [1.5.5] - 2026-10-19

### Added

- Support compressed inputs in the CGC [#163](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/163)
- add bam as valid option for reads [#164](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/164)

## [1.5.4] - 2026-07-13

### Changed

- EMG-9948 Handle empty amplicon result files [#162](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/162)

## [1.5.3] - 2026-07-10

### Added

- Feature/automated docker container [#158](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/158)
- Support compressed inputs on go_utils and summarise_goslims [#160](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/160)

### Changed

- Update container to version 1.5.1 [#159](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/159)

## [1.5.2] - 2026-06-30

### Changed

- Improve robustness of amplicon summaries generating scripts [#157](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/157)

## [1.5.1] - 2026-05-19

### Added

- Add reads submission support to webin cli wrapper [#155](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/155)

## [1.5.0] - 2026-05-05

### Fixed

- Fix/dwcr empty files support [#154](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/154)

## [1.4.24] - 2026-04-23

### Added

- Add missing `__init__.py` in the shared proteins module [#153](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/153)

### Fixed

- Fix contig_id parsing in RHEA script [#152](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/152)

## [1.4.23] - 2026-04-16

### Added

- Add unclassified column to asv count table [#149](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/149)

## [1.4.22] - 2026-04-15

### Added

- support more than one manifest in webin_cli_handler [#148](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/148)

## [1.4.21] - 2026-03-20

### Added

- Add asvs_only category for amplicon output summary table. [#146](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/146)

### Fixed

- bugfix: hardcoded defaults for JAVA_HEAP_SIZE_INITIAL and JAVA_HEAP_SIZE_MAX caused OUT_OF_MEM [#147](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/147)

## [1.4.20] - 2026-03-18

### Added

- Add antismash_sideload_schemas to dependencies [#145](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/145)

### Changed

- Refactor error handling and increase verbosity in webin-cli-wrapper [#144](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/144)

### Fixed

- Debug/antismash gff builder [#143](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/143)

## [1.4.19] - 2026-03-12

### Fixed

- Debug/json validation [#142](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/142)

## [1.4.18] - 2026-03-06

### Added

- bgc_mapper script [#138](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/138)

### Changed

- Refactor FASTA file handling and improve readability [#140](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/140)
- Modules / packages layout for GFF and BGC code [#141](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/141)

## [1.4.17] - 2026-02-24

### Added

- gbk generator script [#135](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/135)
- Zenodo downloader script [#137](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/137)
- feat(combined_gene_caller_merge): Support provision of version string of the tools to annotate GFF source column with tool+version [#136](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/136)

### Changed

- Have the webin cli version in the GH action instead of vars [#139](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/139)

## [1.4.16] - 2026-01-28

### Fixed

- Fix/dwcr skip unclassified line [#134](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/134)

## [1.4.15] - 2026-01-23

### Changed

- NCBITaxa - with domain instead of superkingdom - asa study summaries [#133](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/133)

### Fixed

- Fix/dwcr generator dblabels [#132](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/132)

## [1.4.14] - 2026-01-21

### Changed

- Use d__ instead of sk__ in for contigs taxa (krona) [#130](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/130)
- For the assembly analysis pipeline - contigs taxonomic assignation we settled with `domain` and `d__` instead of `sk__` (superkingdom). This keeps our results aligned with the NCBI taxonomy (for more details - https://ncbiinsights.ncbi.nlm.nih.gov/2024/06/04/changes-ncbi-taxonomy-classifications/)

## [1.4.13] - 2026-01-20

### Fixed

- Bugfix: replace superkingdom to domain [#129](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/129)

## [1.4.12] - 2026-01-13

### Added

- Add sorting of output FASTA files in CGC to mitigate pyhmmer bug [#127](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/127)
- Add taxids to DwC-R generator, and generate ITSoneDB + UNITE + SILVA-LSU closed reference summaries [#123](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/123)

### Changed

- Pin versions of deps and bump version 1.4.9 -> 1.4.10 [#122](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/122)

### Fixed

- Fix failure on empty mapseq output file [#126](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/126)

## [1.4.11] - 2026-01-07

### Added

- Add error handling for empty mapseq file [#124](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/124)

### Changed

- DbCAN output changed the substarates separator from , to ; [#120](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/120)

### Fixed

- Bugfix/empty mapseq to asv file [#125](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/125)

## [1.4.10] - 2025-12-28

### Changed

- Pin versions of deps and bump version 1.4.9 -> 1.4.10 [#122](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/122)

## [1.4.9] - 2025-12-23

### Changed

- DbCAN output changed the substarates separator from , to ; [#120](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/120)

## [1.4.8] - 2025-11-28

### Changed

- For the assembly analysis summary, antismash needs to be optional [#116](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/116)

### Fixed

- Fix attribute formatting for CDS features [#114](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/114)
- Fix to deal with empty inputs [#115](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/115)

## [1.4.6] - 2025-11-25

### Fixed

- Fix on amrintegrator script formatting for CDS features [#114](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/114)
- Fix on amrintegrator script to deal with empty inputs [#115](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/115)

## [1.4.5] - 2025-11-13

### Changed

- Allow any sequence ID in classify_var_regions.py [#111](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/111)
- Use better logging in AMP scripts, includes some small refactors [#109](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/109)

### Fixed

- Bugfix/primer validation single end [#112](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/112)

## [1.4.4] - 2025-11-12

### Changed

- Improve logging of primer validation, improve validation logic [#108](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/108)

## [1.4.3] - 2025-11-11

### Added

- webin_cli handler for assembly_uploader and genome_uploader [#100](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/100)

### Changed

- Merge dev to main [#107](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/107)

## [1.4.2] - 2025-11-11

### Changed

- Amplicon and raw-reads study summaries are named *_study_summary.tsv (now assemblies too) [#105](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/105)

### Fixed

- Bugfix/handle empty antismash [#104](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/104)

## [1.4.1] - 2025-11-06

### Added

- Feature/amr integrator [#101](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/101)

### Fixed

- Fix lack of output in make_asv_count_table.py in rare case [#102](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/102)

## [1.3.0] - 2025-10-09

### Changed

- Big schemas refactor pydantic+pandera => pandera only [#98](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/98)
- Panderas only Schemas for validation of pipeline results. [#99](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/99)

## [1.2.11] - 2025-09-29

### Fixed

- Prevent non-empty tsv output on empty dictionary [#97](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/97)

## [1.2.10] - 2025-09-23

### Added

- Added dada2_stats_fail to amplicon results schema [#96](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/96)

## [1.2.9] - 2025-09-17

### Changed

- AntiSMASH GFF builder fix to work with Bakta-based outputs [#93](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/93)

## [1.2.8] - 2025-08-19

### Changed

- Updated raw-reads summary generation/merging scripts [#91](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/91)
- Dwcsummary generator changes [#90](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/90)

## [1.2.7] - 2025-07-11

### Added

- added raw-reads-analysis-pipeline study summary generator [#77](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/77)

### Changed

- change column names and filenames for raw-reads-analysis-pipeline [#88](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/88) (and [#89](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/89))

## [1.2.6] - 2025-07-07

### Changed

- loosen "useful sequences" threshold [#86](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/86)

## [1.2.5] - 2025-07-02

### Fixed

- fix bug related to not parsing properly the original primer name if it appears in the cmsearch file [#85](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/85)

## [1.2.4] - 2025-07-02

### Added

- add new script that will permute primers because of ambiguous IUPAC codes [#83](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/83)

## [1.2.3] - 2025-07-01

Release notes were not provided.

## [1.2.2] - 2025-06-30

### Fixed

- Prevent failure if, among SILVA's or PR's taxonomic annotations, none match detected asvs [#80](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/80)

## [1.2.1] - 2025-06-17

### Removed

- remove anything related to scripts that now makeup PIMENTO [#78](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/78)

## [1.2.0] - 2025-05-16

### Changed

- antiSMASH v8 support for ASA study summary [#75](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/75)

## [1.1.2] - 2025-05-08

### Changed

- Refactor dbcan_process_cazys based on new dbcan version output [#72](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/72)

## [1.1.1] - 2025-05-02

### Changed

- Study summary generator for assembly analysis

## [1.1.0] - 2025-04-29

### Changed

- make another fix to classify_var_regions [#66](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/66)

## [1.0.9] - 2025-04-28

### Fixed

- Fix another classify_var_region bug [#65](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/65)

## [1.0.8] - 2025-04-24

### Fixed

- Fix classify_var_regions bug [#64](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/64)

## [1.0.7] - 2025-04-23

### Fixed

- fix hook_compressed in process_dbcan scripts [#63](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/63)

## [1.0.6] - 2025-04-22

### Fixed

- Fix argparse issues with gff-toolkit [#62](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/62)

## [1.0.5] - 2025-04-17

### Changed

- BGCs summary for assembly-analysis-pipeline [#58](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/58)
- Gff toolkit [#34](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/34)

### Fixed

- Hotfix/rename columns bgcs [#60](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/60)

## [1.0.4] - 2025-03-27

### Changed

- get_subunits fasta header ID slash fix [#53](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/53)
- Dependency upgrades [#54](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/54)

### Fixed

- Fix parsing of taxdump files [#55](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/55)

## [1.0.3] - 2025-03-25

### Added

- New script added to produce KRONA txt input file from CAT results [#47](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/47)

### Changed

- Tweak the go sumarries - as sorted TSVs [#51](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/51)
- ncrna genomes pipeline scripts + tests using executables [#50](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/50)

## [1.0.2] - 2025-03-21

### Added

- Add some functionality needed for general RNA detection process [#48](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/48)

## [1.0.1] - 2025-03-15

### Changed

- Study summary generator implementation improvements [#45](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/45)

### Fixed

- Fix the combined gene caller merge script [#46](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/46)

## [1.0.0] - 2025-03-10

### Added

- Feature/gff output hard [#42](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/42)
- Add DarwinCore-ready generating summary script [#41](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/41)

### Changed

- Push dev into main for version 0.2.2 [#37](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/37)
- update upload-artifact to v4 [#39](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/39)
- Dev [#40](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/40)
- cgc_merge is updated to generate GFF output and parse annotations from GFF input [#35](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/35)
- Dev [#44](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/44)

### Fixed

- Ported missing * -> X replacement in the CGC merge script [#43](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/43)

## [0.2.2-latest] - 2025-03-04

### Changed

- Push dev into main for version 0.2.2 [#37](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/37)
- update upload-artifact to v4 [#39](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/39)
- Dev [#40](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/40)

## [0.2.1] - 2025-01-28

### Changed

- PR [#31](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/31)
- PR [#32](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/32)

## [0.2.0] - 2025-01-24

### Added

- add combined-gene-caller-merge [#29](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/29)
- Add goslim scripts [#28](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/28)

## [0.1.9] - 2025-01-10

### Added

- Add new script to produce table with Rhea IDs and CHEBI reactions from DIAMOND output [#20](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/20)
- New script: study_summary_generator + Pydantic/Pandera models/schemas + Click [#27](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/27)

### Changed

- Merging Sonya's Rhea script [#25](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/25)

## [0.1.8] - 2024-10-23

### Fixed

- Fixed ASV discrepancy bug

## [0.1.7] - 2024-10-08

### Fixed

- Fixing some bugs in get_subunits and primer_val_classification

## [0.1.6] - 2024-09-02

### Added

- New features, all coming from [#19](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/19)

### Changed

- antismash_gff_builder unit test
- get_mpt_version script for getting current version of the toolkit

### Removed

- remove exit codes from suffix_header_checks - this would cause it to fail in a Nextflow context

## [0.1.5] - 2024-08-16

### Changed

- [#8](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/8)
- [#9](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/9)
- [#14](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/pull/14)

## [0.1.4] - 2024-07-12

### Changed

- Changes to make_asv_count_tables.py to generate ASV read count files

## [0.1.3] - 2024-05-24

### Changed

- This release contains the final changes to the toolkit to support the initial release of the V6 Amplicon pipeline (PR here: https://github.com/EBI-Metagenomics/asv-gen/pull/7).

## [0.1.2] - 2024-04-26

### Added

- Added support for MAPseq assignments of ASVs

## [0.1.1] - 2024-02-22

### Changed

- Made optimisations for fetch_mcp, reorganised constants, fixed read_count function + more unit tests

## [0.1.0] - 2024-02-19

### Added

- Add standard_primer_matching alias I had forgotten

## [0.0.9] - 2024-02-19

Release notes were not provided.

## [0.0.8] - 2024-02-16

### Changed

- Made some final changes for the V6 amplicon pipeline, should now be ready to be pushed to Bioconda for the first time

## [0.0.7] - 2024-02-06

### Added

- Added unit tests and changed the automation workflows a bit to add an independent pytest step, removed unnecessary bash binaries

## [0.0.6-alpha] - 2024-01-04 (pre-release)

Release notes were not provided.

## [0.0.6] - 2024-01-04 (pre-release)

### Changed

- Initial release

[Unreleased]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/compare/v1.5.5...HEAD
[1.5.5]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.5.5
[1.5.4]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.5.4
[1.5.3]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.5.3
[1.5.2]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.5.2
[1.5.1]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.5.1
[1.5.0]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.5.0
[1.4.24]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.24
[1.4.23]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.23
[1.4.22]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.22
[1.4.21]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.21
[1.4.20]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.20
[1.4.19]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.19
[1.4.18]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.18
[1.4.17]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.17
[1.4.16]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.16
[1.4.15]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.15
[1.4.14]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.14
[1.4.13]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.13
[1.4.12]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.12
[1.4.11]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.11
[1.4.10]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.10
[1.4.9]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.9
[1.4.8]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.8
[1.4.6]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.6
[1.4.5]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.5
[1.4.4]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.4
[1.4.3]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.3
[1.4.2]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.2
[1.4.1]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.4.1
[1.3.0]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.3.0
[1.2.11]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.11
[1.2.10]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.10
[1.2.9]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.9
[1.2.8]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.8
[1.2.7]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.7
[1.2.6]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.6
[1.2.5]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.5
[1.2.4]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.4
[1.2.3]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.3
[1.2.2]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.2
[1.2.1]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.1
[1.2.0]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.2.0
[1.1.2]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.1.2
[1.1.1]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.1.1
[1.1.0]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.1.0
[1.0.9]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.0.9
[1.0.8]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.0.8
[1.0.7]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.0.7
[1.0.6]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.0.6
[1.0.5]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.0.5
[1.0.4]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.0.4
[1.0.3]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.0.3
[1.0.2]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.0.2
[1.0.1]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.0.1
[1.0.0]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v1.0.0
[0.2.2-latest]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.2.2-latest
[0.2.1]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.2.1
[0.2.0]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.2.0
[0.1.9]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.1.9
[0.1.8]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.1.8
[0.1.7]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.1.7
[0.1.6]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v.0.1.6
[0.1.5]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.1.5
[0.1.4]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.1.4
[0.1.3]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.1.3
[0.1.2]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.1.2
[0.1.1]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.1.1
[0.1.0]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.1.0
[0.0.9]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.0.9
[0.0.8]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v.0.0.8
[0.0.7]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v.0.0.7
[0.0.6-alpha]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.0.6-alpha
[0.0.6]: https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit/releases/tag/v0.0.6
