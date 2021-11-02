#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 


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


process fastcatUncompress {
    label "wfalignment"
    cpus params.threads
    input:
        tuple file(directory), val(sample_name) 
    output:
        path "*.fastq", emit: fastq
        env SAMPLE_NAME, emit: sample_name
    """
    fastcat -H -s ${sample_name} -r ${sample_name}.stats -x ${directory}  >> ${sample_name}.fastq
    SAMPLE_NAME="${sample_name}"
    """
}


process alignReads {
    label "wfalignment"
    cpus params.threads
    input:
        file fastq
        file reference_files
        file combined
        file counts
    output:
        path "*.sorted.aligned.bam", emit: sorted
        path "*.mapula.json", emit: json
        path "*.unmapped.stats", emit: unmapped_stats 
    script:
        def counts_arg = counts.name != 'NO_COUNTS' ? "-c ${counts}" : ""
        def sampleName = fastq.simpleName

    """
    minimap2 -y -t $task.cpus -ax map-ont $combined $fastq \
    | mapula count $counts_arg -r $reference_files -s fasta run_id barcode -f json -p \
    | samtools sort -@ $task.cpus -o ${sampleName}.sorted.aligned.bam - 
    mv mapula.json ${sampleName}.mapula.json
    bamtools split -in ${sampleName}.sorted.aligned.bam -mapped
    (bedtools bamtofastq -i *UNMAPPED.bam -fq unmapped.fq && fastcat -s unmapped.fq -r ${sampleName}.unmapped.stats -x unmapped.fq >> uncompressed.fastq) \
    || touch ${sampleName}.unmapped.stats
    
    """
}

   
process mergeBAM {
    label "wfalignment"
    cpus params.threads
    input:
        file sorted_aligned_bam
    output:
        file "*.merged.sorted.aligned.bam"
    script:
        def sampleName = sorted_aligned_bam.simpleName
    """
    samtools merge -@ $task.cpus ${sampleName}.merged.sorted.aligned.bam $sorted_aligned_bam
    
    """
}




process readDepth {
    label "wfalignment"
    cpus 1
    input:
        file merged_BAM
    output:
        path "*.bed.gz"
        file "*merged.sorted.aligned.bam.bai"
    script:
        def sampleName = merged_BAM.simpleName
    """
    samtools index $merged_BAM
    mosdepth -n --fast-mode --by 5 ${sampleName}.merged $merged_BAM
    """
}


process gatherStats {
    label "wfalignment"
    cpus 1
    input:
        file mapula_json
        file counts
    output:
        path "*merged.mapula.csv", emit: merged_mapula_csv
        path "*merged.mapula.json", emit: merged_mapula_json
    script:
        def counts_arg = counts.name != 'NO_COUNTS' ? "-c ${counts}" : ""
        def sampleName = mapula_json.simpleName
    """
    mapula merge $mapula_json $counts_arg -f all -n ${sampleName}.merged.mapula
    """
}


process getVersions {
    label "wfalignment"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    aplanat --version | sed 's/ /,/' >> versions.txt
    python -c "import pysam; print(pysam.__version__)" | sed 's/^/spoa,/'  >> versions.txt
    """
}


process getParams {
    label "wfalignment"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process plotStats {
    label "wfalignment"
    cpus params.threads
    input:
        path "merged_mapula_json/*"
        file counts
        file "bed_file/*"
        file reference_files
        path "unmapped_stats/*"
        path "versions/*"
        path "params.json"
        file sample_names

    output:
        path "*.html", optional: true, emit: report
    script:
        def counts_arg = counts.name != 'NO_COUNTS' ? "-c ${counts}" : ""
        def report_name = "wf-alignment-" + params.report_name
    """
    report.py merged_mapula_json/* $counts_arg --report_name '$report_name' \
    --bedfile bed_file/* \
    --references $reference_files \
    --unmapped_stats unmapped_stats/* \
    --params params.json \
    --versions versions \
    --sample_names $sample_names

    """
}


// workflow module
workflow pipeline {
    take:
        fastq
        references
        counts
    main:
        
        //uncompress aplotnd combine fastq's if multiple files
        uncompressed = fastcatUncompress(fastq)

        // Get reference fasta files from dir path
        references_files = channel
            .fromPath("${references}/*", glob: true)

        // Cat the references together for alignment
        combined = combineReferences(references_files.collect())

        // Align the reads to produce bams and stats
        aligned = alignReads(
           uncompressed.fastq, references_files.collect(), combined, counts)
        merged = mergeBAM(aligned.sorted)
        depth = readDepth(merged)
        workflow_params = getParams()
        software_versions = getVersions()


        // Merge the stats together to get a csv out
        stats = gatherStats(alignReads.out.json, counts)
        sample_names = uncompressed.sample_name.collectFile(name: 'sample_names.csv', newLine: true)

        report = plotStats(stats.merged_mapula_json.collect(), counts,
        depth[0].collect(), references_files.collect(),
        aligned.unmapped_stats.collect(),  software_versions, workflow_params,
        sample_names)
    emit:
        merged = merged
        indexed = depth[1]
        merged_mapula_csv = stats.merged_mapula_csv
        merged_mapula_json = stats.merged_mapula_json
        report = report.report
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
WorkflowMain.initialise(workflow, params, log)
workflow {
    // Acquire fastq directory
    fastq = fastq_ingress(
        params.fastq, params.out_dir, params.samples, params.sanitize_fastq)
    // Acquire reference files
    references = file(params.references, type: "dir", checkIfExists: true)
    counts = file(params.counts, checkIfExists: params.counts == 'NO_COUNTS' ? false : true)
    // Run pipeline
    results = pipeline(fastq, references, counts)
    output(results.merged.concat( 
        results.indexed, results.merged_mapula_csv, 
        results.merged_mapula_json, results.report ))
}
