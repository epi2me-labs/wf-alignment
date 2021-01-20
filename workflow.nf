#!/usr/bin/env nextflow

params.help = ""
params.fastq = ""
params.reference = ""
params.threads = 3
params.output = "./output"

if(params.help) {
    log.info ''
    log.info 'Workflow: Alignment'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run workflow.nf [options]'
    log.info ''
    log.info 'Script Options: '
    log.info '    --fastq		DIR	Path to FASTQ files'
    log.info '    --reference   FILE	Path to reference genome'
    log.info '    --threads		INT	Number of threads to use'
    log.info '    --output		DIR	Directory to write output files'
    log.info '    --help		BOOL	Display help message'
    log.info ''

    return
}

Channel
    .fromPath("$params.fastq/*")
    .buffer( size: 300, remainder: true )
    .set { reads }

process alignAndCountReads {
    publishDir "${params.output}", pattern: "*.json"
    input:
    tuple fastqs from reads
    output:
    file "sorted.aligned.bam" into sorted_BAM
    file "mapping-stats.json" into stats_JSON
    """
    minimap2 -y -t $task.cpus -ax map-ont $params.reference $fastqs \
    | gather_mapping_stats -j mapping-stats.json \
    | samtools sort -o sorted.aligned.bam -
    """
}

process bamMerge {
    publishDir "${params.output}"
    input:
    tuple bams from sorted_BAM
    output:
    file "merged.sorted.aligned.bam" into merged_BAM
    """
    samtools merge merged.sorted.aligned.bam $bams
    """
}

process indexBAM {
    publishDir "${params.output}"
    input:
    file bam from merged_BAM
    output:
    file "merged.sorted.aligned.bam.bai" into indexed_BAM
    """
    samtools index $bam
    """
}

process createAlignmentReport {
    publishDir "${params.output}"
    input:
    file stats from stats_JSON
    output:
    file "VisualiseMappingStats_Report.html" into report_HTML
    file "VisualiseMappingStats_Export.json" into export_JSON
    """
    visualise_mapping_stats -j $stats -o .
    """
}