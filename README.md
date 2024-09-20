# Alignment workflow

Align Nanopore reads and visualize mapping statistics.



## Introduction

This workflow provides an easy way to align Oxford Nanopore reads and gather mapping
stats either locally for small amounts of data or at scale in a distributed
environment such as a cluster or the cloud.

> This workflow contains minimal functionality that is duplicated in many of our more specialised workflows.
> Please consider using one of these alternative workflows before using this one: you very likely do not need
> to use this workflow.

In brief, it will perform the following:
* Combine all reference files in the directory passed to `--references`.
* Align input reads (passed as FASTQ or unaligned BAM files) against the reference (Note that BAM files with aligned reads can be used as well; these will skip the alignment step and only stats and the report will be produced).
* Create alignment stats.
* Calculate depth of coverage along the reference sequences (this step can be skipped if requested).
* Create an HTML report to illustrate the results.




## Compute requirements

Recommended requirements:

+ CPUs = 12
+ Memory = 32GB

Minimum requirements:

+ CPUs = 6
+ Memory = 12GB

Approximate run time: 0.5-5 minutes per sample (depending on number of reads, length of reference, and available compute).

ARM processor support: True




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop)
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-alignment --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-alignment
```

A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-alignment/wf-alignment-demo.tar.gz
tar -xzvf wf-alignment-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-alignment \
	--fastq 'wf-alignment-demo/fastq' \
	--references 'wf-alignment-demo/references' \
	-profile standard
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).




## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts either FASTQ or BAM files as input.

The FASTQ or BAM input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ or BAM file; (ii) the path to a top-level directory containing FASTQ or BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ or BAM files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| references | string | Path to a directory containing FASTA reference files. | Accepted file extensions are '.fasta', '.fna', '.ffn', '.faa', '.frn', '.fa', '.txt', '.fa.gz', '.fna.gz', '.frn.gz', '.ffn.gz', '.fasta.gz'. In addition, an MMI index file can be provided to make the workflow run faster using the option `--reference_mmi_file`. |  |
| reference_mmi_file | string | Path to an MMI index file to be used as reference. | Accepted file extension is '.mmi'. The references parameter is still required if this is provided. Note that some minimap2 alignment options are set by the reference MMI and cannot be overridden. |  |
| counts | string | Path to a CSV file containing expected counts as a control. | The expected counts CSV file must contain columns named 'reference' and 'expected_counts' in order to be valid. the 'reference' column should contain names matching the names of reference sequences within the fasta files provided using --references. |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |
| prefix | string | Optional prefix attached to each of the output filenames. | Output filename format will be `<prefix>-filename.ext`. |  |
| per_read_stats | boolean | Generate Bamstats per-read stats. | With this option, the workflow will produce detailed per-read alignment stats emitted as gzipped TSV file. As these files can get quite large, it is recommended to only request them when necessary. | False |


### Advanced options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| depth_coverage | boolean | Calculate depth coverage statistics and include them in the report. | This step can be a computational bottleneck. Set this to false if your reference sequences are >50mb to speed things up. | True |
| minimap_preset | string | Pre-defined parameter sets for `minimap2`, covering most common use cases. | Available parameter sets are: 'dna' (`-ax map-ont`), 'rna' (`-ax splice -uf`). | dna |
| minimap_args | string | String of command line arguments to be passed on to `minimap2`. | This overrides the options defined by `--minimap_preset` and allows for running the alignment step in a more customized way. |  |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Number of CPU threads to use for the alignment step. | The alignment process will run with this many threads (note that the memory used by minimap2 scales with the number of threads). The total CPU resources used by the workflow are constrained by the Nextflow executor configuration. | 4 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | wf-alignment-report.html | Report for all samples | aggregated |
| Combined references | combined-refs.fasta | FASTA file containing all input references. | aggregated |
| Combined references index | combined-refs.fasta.fai | Index file for combined references FASTA. | aggregated |
| Combined references MMI index | combined-refs.mmi | Minimap2 index file for combined references FASTA. | aggregated |
| Per-read alignment stats | {{ alias }}.readstats.tsv.gz | Bamstats per-read output TSV file (compressed with gzip). | per-sample |
| Per-reference alignment stats | {{ alias }}.flagstat.tsv | Bamstats flagstat output TSV file. | per-sample |
| Alignment accuracy histogram | {{ alias }}-histograms/accuracy.hist | Bamstats alignment accuracy histogram TSV file. | per-sample |
| Alignment coverage histogram | {{ alias }}-histograms/coverage.hist | Bamstats alignment coverage histogram TSV file. | per-sample |
| Read length histogram (mapped) | {{ alias }}-histograms/length.hist | Bamstats read length histogram TSV file (for mapped reads). | per-sample |
| Read length histogram (unmapped) | {{ alias }}-histograms/length.unmap.hist | Bamstats read length histogram TSV file (for unmapped reads). | per-sample |
| Read quality histogram (mapped) | {{ alias }}-histograms/quality.hist | Bamstats read quality histogram TSV file (for mapped reads). | per-sample |
| Read quality histogram (unmapped) | {{ alias }}-histograms/quality.unmap.hist | Bamstats read quality histogram TSV file (for unmapped reads). | per-sample |
| Alignments BAM file | {{ alias }}.sorted.aligned.bam | BAM file with alignments of filtered input reads against the combined references. | per-sample |
| Alignments index file | {{ alias }}.sorted.aligned.bam.bai | Index for alignments BAM file. | per-sample |
| IGV config JSON file | igv.json | JSON file with IGV config options to be used by the EPI2ME Desktop Application. | aggregated |




## Pipeline overview

### 1. Combine reference files

All reference files in the directory passed to `--references` are concatenated.

### 2. Align reads

Input reads are aligned against the combined reference with [Minimap2](https://github.com/lh3/minimap2). If BAM files are used as input (with `--bam`), only reads in files without a reference in the SAM header are aligned. For other BAM files this step is skipped.

### 3. Create alignment stats

[Bamstats](https://github.com/epi2me-labs/fastcat#bamstats) is used to create per-read and per-reference alignment stats from the BAM files.

### 4. Calculate depth of coverage

Depth of coverage along the reference sequences is determined with [Mosdepth](https://github.com/brentp/mosdepth) (using 200 windows per reference sequence). To speed up the workflow, this step can be skipped by adding `--depth-coverage false`.




## Troubleshooting

+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ Please see [here](https://labs.epi2me.io/trouble-shooting/) for how to resolve some common Nextflow issues and [here](https://labs.epi2me.io/how-to-exits/) for how to interpret command exit codes.




## FAQ's

*I cannot select a single reference file in the EPI2ME desktop app.* - When running the workflow via the desktop app, you need to provide a directory with reference files. If you only have a single file, you can create a directory to place your reference file inside and select this with the reference input option.

*How are the values in the `acc` column (and other metrics) in the per-read output stats calculated?* -
For details on the per-read stats output files, please refer to the [fastcat/bamstats documentation](https://github.com/epi2me-labs/fastcat#output-format).

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-alignment/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).




## Related blog posts

- [How to align your data](https://labs.epi2me.io/how-to-align/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.




