# phac-nml/MeaSeq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.4.3] - 2025-11-06

Adjusting alignment filtering for both nanopore and illumina data to remove supplementary reads and secondary reads. The supplementary reads were rarely adding in artifacts to the final consensus sequence including snps and indels that only they contained. This lowers read counts slightly

Also added the updates to the clair3 vcf filter found in Artic to add in the SNP allele frequency cutoff along with the high frameshift indel cutoff. There was also the addition of a RefCall filter to prevent those from being reported in the final report table (Note that they don't affect the consensus result at all, just the report visualization)

Overall, the update should fix some rare artifacts to improve the accuracy of the final results

### `Added`

- Samtools sort nf-core module to replace samtools index module in illumina pipeline [PR #21](https://github.com/phac-nml/measeq/pull/21)
  - As `sort` was removed from the process `BWAMEM2_MEM` itself to use view and filter specific alignments alignments
  - Output remains the same as `<sample>.sorted.bam` with the alignment filtering having been done

### `Adjusted`

- `BWAMEM2_MEM` process to run samtools view with `-f 3 -F 2048` to just keep mapped and properly paired alignments along with removing supplementary [PR #21](https://github.com/phac-nml/measeq/pull/21)
- `MINIMAP2_ALIGN` process added in `-F 2308` to remove unmapped, not primary alignment, and supplementary alignment [PR #21](https://github.com/phac-nml/measeq/pull/21)
- Tests adjusted to conform with the overall lower mapped reads with the change in both variant calling subworkflows [PR #21](https://github.com/phac-nml/measeq/pull/21)
  - There was a single variant change in nanopore data test
    - Lower frequency variant at 0.47 AF
  - A couple more Ns overall at the start and end of genome for both
  - Slightly lower depth
- `nanoq` removed the max read length from `modules.json`
- `cs_vcf_filter.py` adjusted to match the source more with adding in min frameshift quality (30 instead of 50) and min allele freq [PR #21](https://github.com/phac-nml/measeq/pull/21)
  - Added a REFCALL filter so that refcalls were just removed as well

## [v0.4.2] - 2025-10-23

Exposing more parameters to allow users more options in adjusting illumina data analyses with parameters

### `Added`

- New parameters to control illumina variant calling [PR #18](https://github.com/phac-nml/measeq/pull/18)
  - `--ivar_trim_min_read_length`
  - `--min_alt_fraction_freeabyes`
  - `--min_variant_qual_freebayes`
- Minimum and Maximum thresholds for parameters in the `nextflow_schema` [PR #18](https://github.com/phac-nml/measeq/pull/18)
- Example DSID fasta file [PR #19](https://github.com/phac-nml/measeq/pull/19)

### `Fixes`

- Fixed the final report provenance typo [PR #18](https://github.com/phac-nml/measeq/pull/18), [Issue #17](https://github.com/phac-nml/measeq/issues/17)

## [v0.4.1] - 2025-09-09

Small addition of Picard MarkDuplicates workflow along with some new tests

### `Added`

- nf-core Picard MarkDuplicates workflow as an optional parameter/workflow to use for illumina data [PR #15](https://github.com/phac-nml/measeq/pull/15)
  - Along with this, added the bam stats samtools workflow to run even when the picard workflow isn't to keep outputs the same

## [v0.4.0] - 2025-09-03

### `Added`

- 3 columns to final Excel and CSV file [PR #12](https://github.com/phac-nml/measeq/pull/12)
  - N450_completeness
  - N450_mean_depth
  - N450_status
- Fail tracking for all samples and runs where every sample fails [PR #12](https://github.com/phac-nml/measeq/pull/12)

### `Adjusted`

- Adjusted the DSId final report page to have a list of all individual sample calls and the summary data shown [PR #12](https://github.com/phac-nml/measeq/pull/12)
- Added the N450 status to the initial report summary page [PR #12](https://github.com/phac-nml/measeq/pull/12)
- Always run the DSId check even with no database fasta file [PR #12](https://github.com/phac-nml/measeq/pull/12)
  - Everything will be labeled as `novel-hash` but it will group them up
- Updated CI tests [PR #12](https://github.com/phac-nml/measeq/pull/12)

## [v0.3.2] - 2025-08-01

### `Added`

- DSId tab to final HTML report for DSId summary information if a DSId file was given as input [PR #10](https://github.com/phac-nml/measeq/pull/10)
- New `overall.xlsx` final QC file based on adding in a few more columns [PR #10](https://github.com/phac-nml/measeq/pull/10)
  - `genome_fasta` and `N450_fasta` that contain fasta formatted sequence data
    - This broke the CSV file version parsing so that remains unchanged to get metadata to IRIDANext

### `Adjusted`

- Fixed a bug with MeaSeq Report summary table not correctly linking to certain samples [PR #10](https://github.com/phac-nml/measeq/pull/10)
- Removed `versions.yml` from being created during final report as versions already were reported [PR #10](https://github.com/phac-nml/measeq/pull/10)
- Adjusted code for final report row formatting to match throughout Rmd files [PR #10](https://github.com/phac-nml/measeq/pull/10)
- DSId assignment changes [PR #10](https://github.com/phac-nml/measeq/pull/10)
  - Hash adjusted to 7 characters
  - Semi-Complete removed to just be Incomplete
- New container for the MAKE_FINAL_QC_CSV step to remove artic container [PR #10](https://github.com/phac-nml/measeq/pull/10)
  - Adds in openpyxl

## [v0.3.1] - 2025-07-29

### `Adjusted`

- Fixed a bug where the irida json file wasn't being populated with metadata [PR #8](https://github.com/phac-nml/measeq/pull/8)
- Adjusted Illumina nf-test to use the `sample_name` field of the samplesheet to test that it works [PR #8](https://github.com/phac-nml/measeq/pull/8)

### `Removed`

- Normalized median read depth plot for now as it overlaps the full depth one too much [PR #8](https://github.com/phac-nml/measeq/pull/8)

## [v0.3.0] - 2025-07-25

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

[v0.4.3]: https://github.com/phac-nml/measeq/releases/tag/0.4.3
[v0.4.2]: https://github.com/phac-nml/measeq/releases/tag/0.4.2
[v0.4.1]: https://github.com/phac-nml/measeq/releases/tag/0.4.1
[v0.4.0]: https://github.com/phac-nml/measeq/releases/tag/0.4.0
[v0.3.2]: https://github.com/phac-nml/measeq/releases/tag/0.3.2
[v0.3.1]: https://github.com/phac-nml/measeq/releases/tag/0.3.1
[v0.3.0]: https://github.com/phac-nml/measeq/releases/tag/0.3.0
[v0.2.1-dev]: https://github.com/phac-nml/measeq/releases/tag/0.2.1-dev
[v0.2.0-dev]: https://github.com/phac-nml/measeq/releases/tag/0.2.0-dev
