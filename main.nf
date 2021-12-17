#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 
include { start_ping; end_ping } from './lib/ping'



def nameIt(ch) {
            return ch.map { it -> return tuple("$it".split(/\./)[-2], it) }
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


process fastcatUncompress {
    label "wfalignment"
    cpus params.threads
    input:
        tuple path(directory), val(sample_id), val(type)
    output:
        path "*.fastq.gz", emit: fastq
        env SAMPLE_ID, emit: sample_id
    """
    fastcat -H -s ${sample_id} -r temp.stats -x ${directory} > temp.fastq 
    fastq_header_to_SAM_tags.py temp.fastq -o ${sample_id}.reads.fastq
    gzip ${sample_id}.reads.fastq
    SAMPLE_ID="${sample_id}"
    rm -rf *temp*
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
    (bedtools bamtofastq -i *UNMAPPED.bam -fq temp.unmapped.fq \
    && fastcat -s temp.unmapped.fq -r ${sampleName}.unmapped.stats -x temp.unmapped.fq >> temp.uncompressed.fastq) \
    || touch ${sampleName}.unmapped.stats
    rm -rf *temp*
    rm -rf *UNMAPPED*
    """
}

   
process mergeBAM {
    label "wfalignment"
    cpus params.threads
    input:
        file sorted_aligned_bam
    output:
        path "*.merged.sorted.aligned.bam"
    script:
        def sampleName = sorted_aligned_bam.simpleName
    """
    samtools merge -@ $task.cpus ${sampleName}.merged.sorted.aligned.bam $sorted_aligned_bam
    """
}


process indexBam {
    label "wfalignment"
    cpus 1
    input:
        file merged_BAM
    output:
        file "*merged.sorted.aligned.bam.bai"
    script:
        def sampleName = merged_BAM.simpleName
    """
    samtools index $merged_BAM
    """
}


process splitByReference {
    label "wfalignment"
    cpus params.threads
    input:
        file merged_bam
    output:
       path "*.sorted.aligned.*.bam", emit: per_ref_bam
    script:
        def sampleName = merged_bam.simpleName
    """
    samtools index $merged_bam
    bamtools split -in $merged_bam -reference -refPrefix "$sampleName".
    """
}


process refLengths {
    label "wfalignment"
    cpus params.threads
    input:
        path reference
        each path(indexed_bam)
    output:
        path "lengths.csv"
    """
    seqkit fx2tab --length --name --only-id $reference > lengths.txt  
    echo 'name,lengths' > lengths.csv
    tr -s '[:blank:]' ',' <lengths.txt >> lengths.csv
    """
}


process readDepthPerRef {
    depth_threads = {params.threads >= 4  ? 4 : 1}
    label "wfalignment"
    cpus depth_threads
    input:
        tuple val(ref_name), val(ref_len), file(ref_bam)
    output:
        path "*.bed"
    script:
        def sampleName = "$ref_bam".split(/\./)[4]
        int test = "$ref_len".toInteger()
        def steps = Math.round(Math.floor(test/200))
        if (steps == 0){
            steps = 1;
        }
    
    """
    samtools index $ref_bam
    mosdepth -n --fast-mode --by $steps -t $task.cpus ${sampleName}.${ref_name}.temp $ref_bam
    zgrep '$ref_name' ${sampleName}.${ref_name}.temp.regions.bed.gz > ${sampleName}.${ref_name}.bed
    rm -rf *temp*
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
        file reference_files
        path "unmapped_stats/*"
        path "versions/*"
        path "params.json"
        file sample_ids
        path "depth_beds/*"

    output:
        path "*.html", optional: true, emit: report
    script:
        def counts_arg = counts.name != 'NO_COUNTS' ? "-c ${counts}" : ""
        def report_name = "wf-alignment-" + params.report_name
    """
    report.py merged_mapula_json/* $counts_arg --report_name '$report_name' \
    --references $reference_files \
    --unmapped_stats unmapped_stats/* \
    --params params.json \
    --versions versions \
    --sample_names $sample_ids
    """
}


// workflow module
workflow pipeline {
    take:
        fastq
        references
        counts
    main:
        
        //uncompress and combine fastq's if multiple files
        uncompressed = fastcatUncompress(fastq)

        // Cat the references together for alignment
        combined = combineReferences(references.collect())

        // Align the reads to produce bams and stats
        aligned = alignReads(
           uncompressed.fastq, references.collect(), combined, counts)   
        merged = mergeBAM(aligned.sorted)
        indexed_bam = indexBam(merged)

        // Split bam by reference
        split_by_ref = splitByReference(merged)

        // Find reference lengths
        ref_length = refLengths(combined, indexed_bam)

        // Output tuples containing ref_id, length
        length_ch = ref_length[0].splitCsv(header:true)
                    .map{ row-> tuple(row.name, row.lengths) }
        // Output tuples containing ref_id, bam_file
        named_bams = nameIt(split_by_ref.per_ref_bam.flatten())
        // Join tuples to make ref_id, length, bam_file
        ref_len_bam = length_ch.join(named_bams)

        // Find read_depth per reference/bam file
        depth_per_ref = readDepthPerRef(ref_len_bam)

        // get params & versions
        workflow_params = getParams()
        software_versions = getVersions()


        // Merge the stats together to get a csv out
        stats = gatherStats(alignReads.out.json, counts)
        sample_ids = uncompressed.sample_id.collectFile(name: 'sample_ids.csv', newLine: true)

        report = plotStats(stats.merged_mapula_json.collect(), counts,
        references.collect(),
        aligned.unmapped_stats.collect(),  software_versions, workflow_params,
        sample_ids, depth_per_ref.collect())
    emit:
        merged = merged
        indexed = indexed_bam
        merged_mapula_csv = stats.merged_mapula_csv
        merged_mapula_json = stats.merged_mapula_json
        report = report.report
        telemetry = workflow_params
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
    // Start ping
    start_ping()
    // Acquire fastq directory
    fastq = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet, params.sanitize_fastq)
    // Acquire reference files
    references = file(params.references, type: "dir", checkIfExists:true);
    reference_files = []
    extensions = ["fasta", "fna", "ffn", "faa", "frn", "fa", "txt", "fa.gz", "fna.gz", "frn.gz", "ffn.gz", "fasta.gz"]
    for (ext in extensions) {
        reference_files += file(references.resolve("*.${ext}"), type: 'file', maxdepth: 1)
    }
    if (reference_files.size() == 0) {
            println('Error: No references found in the directory provided.')
            exit 1 
    }       
    else {
        reference_files = channel
            .fromPath(reference_files)
    }
    counts = file(params.counts, checkIfExists: params.counts == 'NO_COUNTS' ? false : true)
    // Run pipeline
    results = pipeline(fastq, reference_files, counts)
    output(results.merged.concat( 
        results.indexed, results.merged_mapula_csv, 
        results.merged_mapula_json, results.report ))
    // End ping
    end_ping(pipeline.out.telemetry)
}



