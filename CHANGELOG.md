# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Fixed.
- run_id in mapula output json.
- Only accept certain format files as references
- reduce storage required for workspace

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

Initial release

### Added
- Basic running of alignment workflow and reporting
