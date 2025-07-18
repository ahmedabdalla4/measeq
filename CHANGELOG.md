# phac-nml/MeaSeq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.3.0 - [2025-07-18]

### `Added`

- CI tests and Linting
- Providence to final `Measeq_Report.html` file
- Approximate genome position annotation to final report figures
- nf-validation to work with IRIDA Next
- nf-iridanext to work with IRIDA Next
- nf-prov to allow some more providence options
- Hashes to novel DSId calls so that in the same run new DSIds will match
  - And they will match in other locations as well

### `Adjusted`

- Output directories for Nanopore process in the `modules.config` file
- Updated `artic` version from `1.6.2` --> `1.7.4` for nanopore pipeline
- Ambiguous position handling for illumina data
  - Specifically for rare postions where there was a low-supported INDEL along with a SNP
- Negative control default string to add in 'en'
- Samtools depth `meta1` to `meta` as it was breaking the IRIDA-Next plugin
- Minimum nextflow version required to `24.10.0`

### `Removed`

- nf-schema and associated workflows

## v0.2.1 - [2025-06-12]

### `Added`

- New plots to the summary table for the final report
  - Summary sequencing depth plots
  - No longer stand-alone as it can get a bit large, a smaller report will be made later
- Allow primer schemes to have different direction extensions in the names
  - `_L`, `_R`, `_FORWARD`, `_REVERSE`, `_F`
  - Added in proper errors for if the primer file was not formatted correctly
- Readme updates for primers

### `Bugfixes`

- `genome_completeness` fixed to `genome_completeness_percent` in the negative control checking

### `Removed`

- N plots of the final report to help lower size and speed up opening

## v0.2.0 - [2025-06-06]

### `Added`

- All of the pipeline has been rewritten in Nextflow
- Illumina paired-end sequencing workflow added
  - Freebayes for variant calling over ivar variants/consensus previously
- Nanopore (initial) workflow added
  - clair3
- DSId assignment added when using `--dsid fasta` parameter
  - Based on full sequence match
- Summary outputs added
  - Amplicon summary report
  - The current Rmarkdown report needs to be fixed for the new outputs

### `Deprecated`

- Current MeaSeq script that utilized viralrecon depreciated to make the whole pipeline nextflow

## v0.1.0 - [2025-04-08]

### `Added`

- MeaSeq pipeline created and initial code added
