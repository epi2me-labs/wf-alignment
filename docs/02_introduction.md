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
