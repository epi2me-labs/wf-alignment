# wf-alignment






## Introduction

This workflow provides an easy way to align Oxford Nanopore reads and gather mapping
stats either locally for small amounts of data or at scale in a distributed
environment such as a cluster.




## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/)
to provide isolation of the required software. Both methods are automated
out-of-the-box provided either docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-alignment --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* A sorted, indexed BAM file containing alignments.
* A CSV containing various mapping stats.
* An HTML report with visualisations of the mapping stats.




## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://docs.sylabs.io/guides/latest/user-guide/)
* [minimap2](https://github.com/lh3/minimap2)