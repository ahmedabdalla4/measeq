# phac-nml/MeaSeq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.3.0] - 2025-07-18

### `Added`

- CI tests and Linting [PR #6](https://github.com/phac-nml/measeq/pull/6)
- Providence to final `Measeq_Report.html` file [PR #5](https://github.com/phac-nml/measeq/pull/5)
- Approximate genome position annotation to final report figures [PR #5](https://github.com/phac-nml/measeq/pull/5)
- `nf-validation` to work with IRIDA Next [PR #6](https://github.com/phac-nml/measeq/pull/6)
- `nf-iridanext` to work with IRIDA Next [PR #6](https://github.com/phac-nml/measeq/pull/6)
- `nf-prov` to allow some more providence options [PR #5](https://github.com/phac-nml/measeq/pull/5)
- Hashes to novel DSId calls so that in the same run new DSIds will match [PR #5](https://github.com/phac-nml/measeq/pull/6)
  - And they will match in other locations as well
- `min_indel_threshold` parameter to set the minimum indel threshold required to call an indel [PR #6](https://github.com/phac-nml/measeq/pull/6)

### `Adjusted`

- Output directories for Nanopore process in the `modules.config` file [PR #6](https://github.com/phac-nml/measeq/pull/6)
- Updated `artic` version from `1.6.2` --> `1.7.4` for nanopore pipeline [PR #6](https://github.com/phac-nml/measeq/pull/6)
- Ambiguous position handling for illumina data [PR #5](https://github.com/phac-nml/measeq/pull/5)
  - Specifically for rare postions where there was a low-supported INDEL along with a SNP
- Negative control default string to add in 'en' [PR #5](https://github.com/phac-nml/measeq/pull/5)
- Samtools depth `meta1` to `meta` as it was breaking the IRIDA-Next plugin [PR #6](https://github.com/phac-nml/measeq/pull/6)
- Minimum nextflow version required to `24.10.0` [PR #5](https://github.com/phac-nml/measeq/pull/5)
- Samplesheet input added a `sample_name` column to work with IRIDA Next [PR #6](https://github.com/phac-nml/measeq/pull/6)

### `Removed`

- nf-schema and associated workflows [PR #6](https://github.com/phac-nml/measeq/pull/6)

## [v0.2.1-dev] - 2025-06-12

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

## [v0.2.0-dev] - 2025-06-06

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

## v0.1.0 - 2025-04-08

### `Added`

- MeaSeq pipeline created and initial code added

[v0.3.0]: https://github.com/phac-nml/measeq/releases/tag/0.3.0
[v0.2.1-dev]: https://github.com/phac-nml/measeq/releases/tag/0.2.1-dev
[v0.2.0-dev]: https://github.com/phac-nml/measeq/releases/tag/0.2.0-dev
