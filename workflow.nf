#!/usr/bin/env nextflow

params.help = ""
params.fastq = ""
params.reference = ""
params.threads = 3
params.out_dir = "./output"

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
    log.info '    --out_dir		DIR	Directory to write output files'
    log.info '    --help		BOOL	Display help message'
    log.info ''

    return
}

Channel
    .fromPath("$params.fastq/*")
    .buffer( size: 300, remainder: true )
    .set { reads }

process alignAndCountReads {
    input:
    file "*_reads.fastq" from reads
    output:
    file "sorted.aligned.bam" into sorted_BAM
    file "mapping-stats.json" into stats_JSON
    """
    minimap2 -y -t $task.cpus -ax map-ont $params.reference *_reads.fastq \
    | gather_mapping_stats -j mapping-stats.json \
    | samtools sort -o sorted.aligned.bam -
    """
}

process bamMerge {
    publishDir "${params.out_dir}", mode: 'copy'
    input:
    file 'sorted.aligned_*_.bam' from sorted_BAM.collect()
    output:
    file "merged.sorted.aligned.bam" into merged_BAM
    """
    samtools merge merged.sorted.aligned.bam sorted.aligned_*_.bam
    """
}

process indexBAM {
    publishDir "${params.out_dir}", mode: 'copy'
    input:
    file bam from merged_BAM
    output:
    file "merged.sorted.aligned.bam.bai" into indexed_BAM
    """
    samtools index $bam
    """
}

process combineJSON {
    publishDir "${params.out_dir}", mode: 'copy'
    input:
    file 'mapping-stats_*_.json' from stats_JSON.collect()
    output:
    file "merged.mapping-stats.json" into merged_JSON
    """
    combine_mapping_stats -j mapping-stats_*_.json
    """
}

process createAlignmentReport {
    publishDir "${params.out_dir}", mode: 'copy'
    input:
    file stats from merged_JSON
    output:
    file "VisualiseMappingStats_Report.html" into report_HTML
    file "VisualiseMappingStats_Export.json" into export_JSON
    """
    visualise_mapping_stats -j $stats -o .
    """
}