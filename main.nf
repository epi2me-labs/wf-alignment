#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage(){
    log.info """
wf-alignment

Usage:
    nextflow run epi2melabs/wf-alignment [options]

Script Options:
    --fastq        DIR     Path to directory containing FASTQ files (required)
    --references   DIR     Path to directory containing FASTA reference files (required)
    --counts       FILE    Path to a CSV file containing expected counts (optional)
    --demultiplex  BOOL    Provide this flag to enable demultiplexing of the data (optional)
    --out_dir      DIR     Path for output (default: $params.out_dir)
    --threads      INT     Number of threads per process for alignment and sorting steps (4)
    --batch        INT     Determines how many fastq to split into each parallel job (100)
    --prefix       STR     The prefix attached to each of the output filenames (optional)
    --help

Notes:
    The expected counts CSV file must contain columns named 'reference' 
    and 'expected_counts' in order to be valid. the 'reference' column
    should contain names matching the names of reference sequences within
    the fasta files provided using --references.
"""
}

def displayParamError(msg) {
    helpMessage()
    println("")
    println(msg)
    exit 1
}


if (params.help) {
    helpMessage()
    exit 1
}
if (!params.fastq || !params.references) {
    displayParamError("Error: `--fastq` and `--references` are required")
}


process demultiplexReads {
    cpus params.threads
    input:
        file "reads_*.fastq"
    output:
        path "guppy_barcoder/{**,.}/*.fastq", glob: true
    """
    guppy_barcoder -t $task.cpus -i . --save_path guppy_barcoder
    """
}


process combineReferences {
    label "wfalignment"
    cpus 1
    input:
        file "reference_*_.fasta"
    output:
        file "combined.fasta"
    """
    cat reference_*_.fasta > "combined.fasta"
    """
}


process alignReads {
    label "wfalignment"
    cpus params.threads
    input:
        file "reads_*.fastq"
        file reference_files
        file combined
        file counts
    output:
        path "sorted.aligned.bam", emit: sorted
        path "mapula.json", emit: json
    script:
        def counts_arg = counts.name != 'NO_COUNTS' ? "-c ${counts}" : ""
    """
    fastq_header_to_SAM_tags.py *.fastq \
    | minimap2 -y -t $task.cpus -ax map-ont $combined - \
    | mapula count $counts_arg -r $reference_files -s fasta barcode run_id -f json -p \
    | samtools sort -o sorted.aligned.bam -
    """
}


process mergeBAM {
    label "wfalignment"
    cpus params.threads
    input:
        file "sorted.aligned._*_.bam"
    output:
        file "merged.sorted.aligned.bam"
    """
    samtools merge -@ $task.cpus merged.sorted.aligned.bam sorted.aligned._*_.bam
    """
}


process indexBAM {
    label "wfalignment"
    cpus 1
    input:
        file merged_BAM
    output:
        file "merged.sorted.aligned.bam.bai"
    """
    samtools index $merged_BAM
    """
}


process gatherStats {
    label "wfalignment"
    cpus 1
    input:
        file "mapula_*_.json"
        file counts
    output:
        path "merged.mapula.csv", emit: merged_mapula_csv
        path "merged.mapula.json", emit: merged_mapula_json
    script:
        def counts_arg = counts.name != 'NO_COUNTS' ? "-c ${counts}" : ""
    """
    mapula merge mapula_*_.json $counts_arg -f all -n merged.mapula
    """
}


process plotStats {
    label "wfalignment"
    cpus 1
    input:
        file merged_mapula_json
        file counts
    output:
        file "report.html"
    script:
        def counts_arg = counts.name != 'NO_COUNTS' ? "-c ${counts}" : ""
    """
    aplanat mapula -n report $merged_mapula_json $counts_arg
    """
}


// workflow module
workflow pipeline {
    take:
        fastq
        references
        counts
    main:
        // Get fastq files from dir path
        fastq_files = channel
            .fromPath("${fastq}{**,.}/*.fastq", glob: true)
            .buffer( size: params.batch, remainder: true )

        // Demux if enabled
        if ( params.demultiplex )
            fastq_files = demultiplexReads(fastq_files)

        // Get reference fasta files from dir path
        references_files = channel
            .fromPath("${references}/*", glob: true)

        // Cat the references together for alignment
        combined = combineReferences(references_files.collect())

        // Align the reads to produce bams and stats
        aligned = alignReads(
            fastq_files, references_files.collect(), combined, counts)
        merged = mergeBAM(aligned.sorted.collect())
        indexed = indexBAM(merged)

        // Merge the stats together to get a csv out
        stats = gatherStats(alignReads.out.json.collect(), counts)
        report = plotStats(stats.merged_mapula_json, counts)
    emit:
        merged = merged
        indexed = indexed
        merged_mapula_csv = stats.merged_mapula_csv
        merged_mapula_json = stats.merged_mapula_json
        report = report
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*", saveAs: { 
        f -> params.prefix ? "${params.prefix}-${f}" : "${f}" }
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


// entrypoint workflow
workflow {
    // Acquire fastq directory
    fastq = file(params.fastq, type: "dir", checkIfExists: true)
    // Acquire reference files
    references = file(params.references, type: "dir", checkIfExists: true)
    counts = file(params.counts, checkIfExists: params.counts == 'NO_COUNTS' ? false : true)
    // Run pipeline
    results = pipeline(fastq, references, counts)
    output(results.merged.concat( 
        results.indexed, results.merged_mapula_csv, 
        results.merged_mapula_json, results.report ))
}
