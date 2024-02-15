# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.1.1]
### Changed
- Reduced the memory requested by some processes to avoid failing in WSL (since there is slightly less memory available in WSL than specified in `.wslconfig`).

## [v1.1.0]
### Changed
- Some formatting changes to github issue template.

### Added
- Produce a MMI index file.
- `--reference_mmi_file` option to use a pre-generated MMI index file as reference.

### Removed
- The limit of `20` for the `--threads` parameter.

## [v1.0.3]
### Fixed
- Fix regression in depth plots that concatenated the curves of the different samples, rather than displaying them as a multi-line plot

### Changed
- BAM tags from uBAM inputs are now carried over to the resulting BAM files.

## [v1.0.2]
### Fixed
- Regression that caused failing on compressed references.
- The `alignReads` process requesting too little memory in some cases.

### Changed
- Reduced the minimum memory requirement from 16 to 12 GB.

## [v1.0.1]
### Fixed
- The workflow failing due to commas in reference sequence names.

### Changed
- How samples, reference files, and reference sequence names are listed in the summary section at the beginning of the report.

## [v1.0.0]
### Added
- Memory requirements for each process.

### Changed
- Reworked docs to follow new layout.

## [v0.6.3]
### Fixed
- Mangled depth plots when there are multiple reference sequences.
- Report generation failing when there is only a single read or a small number of reads with near-identical mean quality for a sample or reference file.

## [v0.6.2]
### Fixed
- Report generation failing when a sample name begins with the name of another sample (e.g. 'sample_A' and 'sample_A_2').

### Removed
- Default local executor CPU and RAM limits.

### Changed
- Names of barcoded directories in the sample sheet now need to be of format `barcodeXY`.

## [v0.6.1]
### Fixed
- Workflow failing when using a large number of reference sequences.

## [v0.6.0]
### Removed
- `--ubam` option. `--bam` can now be used for both BAM and uBAM files. The workflow will determine if files are aligned or not (and align them against the provided reference in that case).

## [v0.5.3]
### Fixed
- Read length histogram only displaying a small number of bins when there are a few outlier reads a lot longer than the other reads.
- configure-jbrowse breaking on unescaped spaces

### Changed
- x-axis limits for accuracy, mean read quality, and read alignment coverage histograms to be more dynamic.
### Fixed
- Workflow will no longer crash when running with `--bam` on an input directory containing more than one `.bam` file.

## [v0.5.2]
### Changed
- Removed no longer used `--concat_fastq` parameter.

## [v0.5.1]
### Changed
- Updated GitHub issue templates to force capture of more information.
- Example command to use demo data.

### Fixed
- Tooltips in depth plots not showing.

## [v0.5.0]
### Changed
- Bumped minimum required Nextflow version to 22.10.8.
- Enum choices are enumerated in the `--help` output.
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice.

## [v0.4.1]
### Fixed
- Workflow aborting on `fastcat_or_mv` process.

## [v0.4.0]
### Changed
- Replaced `--threads` option with `--mapping_threads` and `--sorting_threads`, which control the number of threads used during the alignment process.
    - `--mapping_threads` controls the number of threads used by `minimap2`.
    - `--sorting_threads` controls the number of threads used to sort the aligned reads.
    - The total number of threads used by the alignment process is the sum of the two values.
    - Other processes use a hard-coded number of threads ranging between 1 and 3.

### Added
- Parameters `--minimap_args` and `--minimap_preset` to expose additional `minimap2` options to the user.
    - For RNA data sets, `--minimap_preset` can be set to `'rna'` to automatically configure the workflow accordingly (`'dna'` is the default preset).
    - Advanced users can provide `--minimap_args` to pass additional overriding arguments to `minimap2`

## [v0.3.6]
### Added
- Configuration for running demo data in AWS

## [v0.3.5]
### Fixed
- Bug crashing the report when running on AWS without a `--counts` file.

## [v0.3.4]
### Changed
- Now uses ONT Public License.
- Report now uses dropdown menus instead of tabs.

### Fixed
- Missing `seqkit` in `getVersions`process.

## [v0.3.3]
### Removed
- `-y` flag from `minimap2` command

## [v0.3.2]
### Changed
- format to 'directory-path' for parameters fastq, bam, ubam, references

## [v0.3.1]
### Fixed
- missing header for 'Useful links' in docs
- description about references in schema (now only mentions an input directory)

## [v0.3.0]
### Changed
- uses bamstats instead of mapula
- uses ezcharts for report

### Removed
- legacy option 'demultiplex'

## [v0.2.4]
### Fixed
- sample_sheet format in schema to expect a file

## [v0.2.3]
### Changed
- Updated description in manifest

## [v0.2.2]
### Changed
- Harmonized line plot colours in report.
- Expanded explanation for coverage plots.

## [v0.2.1]
### Changed
- Changed plot layout and margins to avoid overflowing plots

## [v0.2.0]
### Added
- Workflow will now output a JBrowse2 `jbrowse.json` configuration

### Changed
- Output combined reference file to `out_dir`
- `-profile conda` is no longer supported, users should use `-profile standard` (Docker) or `-profile singularity` instead
- Removed option for specifying report suffix
- Restructured workflow parameter schema

## [v0.1.9]
### Added
- Input params and handling for bam and ubam formats

### Updated
- Bumped base container to v0.2.0

## [v0.1.8]
### Changed
- Fastqingress metadata map

### Fixed
- Set out_dir option type to ensure output is written to correct directory on Windows.

### Added
- Argument Parser for fastqingress.

## [v0.1.7]
### Fixed
- Coloring with less than 3 samples

## [v0.1.6]
### Fixed
- run id and barcode output correctly

## [v0.1.5]
### Added
- concat_fastq boolean parameter

### Changed
- Better help text on cli

## [v0.1.4]
### Fixed
- Mosdepth 0 step

### Added
- Depth coverage steps parameters

## [v0.1.3]
### Fixed
- Cumulative coverage plotting incorrect numbers

## [v0.1.2]
### Added
- Cumulative coverage plot

## [v0.1.1]
### Changed.
- reference can be either a directory or single file.
- output one merged CSV vs one for each barcode.
- speed up a few steps including mosdepth and report creation.

## [v0.1.0]
### Fixed.
- run_id in mapula output json.
- Only accept certain format files as references.
- reduce storage required for workspace.

### Added.
- Handling for no alignments.
- Integration with EPI2ME Labs notebook environment.

## [v0.0.9]
### Added
- Error message if no references in directory provided.
- Singularity profile.
- Ping telemetry file.

### Fixed
- Calculate depth coverage graph steps based on length of reference.

### Changed
- Sample name to sample id

## [v0.0.8]
### Added
- Option to add suffix to HTML report name.
- Unmapped QC statistics
- Depth coverage graph per reference

### Changed
- Help message now uses JSON schema
- Updated fastqingress

## [v0.0.7]
### Fixed
- Correct conda profile environment file path

## [v0.0.6]
### Fixed
- Remove erroneous --prefix messages
- Increase default batch_size to 1000
- Increase default max local executor cpus to 8

## [v0.0.5]
### Changed
- Retag of v0.0.4, updated sample reports

## [v0.0.4]
### Changed
- Make prefix optional

## [v0.0.3]
### Added
- Barcode awarenesss support with --demultiplex flag (requires guppy_barcoder to be installed)
- Output naming via new required --prefix argument

## [v0.0.2]
### Changed
- Standardised report name.
- Make docker executor default.

## [v0.0.1]
* Initial release

### Added
- Basic running of alignment workflow and reporting
