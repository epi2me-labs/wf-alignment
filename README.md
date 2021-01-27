# Mapping Workflow

This repository contains a nextflow workflow and associated Docker
container build for performing mapping (aka alignment) with minimap2
from basecalls (.fastq) and a reference file (.fasta). The outputs are
a sorted BAM file containing the alignments made, as well as a JSON
file that contains basic statistics which can easily be accessed
programmatically to produce visualisations. In addition, a minimal
report is pre-constructed using this data, in the form of an HTML file
which you can open in the browser. This workflow, using the power of
nextflow, should be able to scale horizontally across compute, e.g.
making use of a cluster.

## Roadmap

> The pipeline is functional, however the report generation feature
> will be significantly improved in an upcoming update. In addition, 
> more command lime options will be added in due course, in addition
> to improved documentation.

## Quickstart

```bash
# build the container
CONTAINER_TAG=mapping
docker build -t ${CONTAINER_TAG} -f Dockerfile  .

# run the pipeline with the test data
nextflow run workflow.nf \
    -w mapping_docker/workspace 
    -profile withdocker
    --fastq test_data/fastq_data/
    --reference test_data/reference.fasta
    --threads 4 --out_dir mapping
```

The output of the pipeline will be found in `./mapping` for the above
example.


## Useful links

* [minimap2](https://github.com/lh3/minimap2)