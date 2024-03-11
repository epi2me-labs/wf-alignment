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
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |


