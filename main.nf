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
        tuple path(directory), val(meta)
    output:
        path "*.fastq.gz", emit: fastq
        env SAMPLE_ID, emit: sample_id
    """
    fastcat -H -s ${meta.sample_id} -r temp.stats -x ${directory} > ${meta.sample_id}.reads.fastq
    gzip ${meta.sample_id}.reads.fastq
    SAMPLE_ID="${meta.sample_id}"
    """
}


process nameFastq {
    label "wfalignment"
    cpus params.threads
    input:
        tuple path(directory), val(sample_id), val(type)
    output:
        path "*.fastq.gz", emit: fastq
        env SAMPLE_ID, emit: sample_id
    """
    SAMPLE_ID="${sample_id}"
    rm -rf *temp*
    cp $directory/* .
    for f in * ; do mv -- "\$f" "${sample_id}.\$f" ; done
    (gzip ${sample_id}.*.f*q) || echo "already gzipped"
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


process refLengths {
    label "wfalignment"
    cpus 1
    input:
        path reference
    output:
        path "lengths.csv"
    """
    seqkit fx2tab --length --name --only-id $reference > lengths.txt
    echo 'name,lengths' > lengths.csv
    tr -s '[:blank:]' ',' <lengths.txt >> lengths.csv
    """
}

process addStepsColumn {
    label "wfalignment"
    cpus 1
    input:
        path "lengths.csv"
    output:
        path "lengths_with_steps.csv"
    '''
    #!/usr/bin/env python
    import pandas as pd
    all = pd.read_csv('lengths.csv')
    all["step"] = all["lengths"]//200
    all = all.replace(0, 1)
    all.to_csv('lengths_with_steps.csv', index=False, header=False)
    '''

}


process readDepthPerRef {
    depth_threads = {params.threads >= 4  ? 4 : params.threads}
    label "wfalignment"
    cpus depth_threads
    input:
        file ref_len
        file alignment
    output:
        path "*bed*"
    script:
       sampleName = alignment.simpleName
       def part = "${alignment}".split(/\./)[1]
    """
    samtools index $alignment
    while IFS=, read -r name lengths steps; do
        mosdepth -n --fast-mode --by "\$steps" --chrom "\$name" -t $task.cpus ${sampleName}."\$name".temp $alignment
    done < $ref_len
    cat *.regions.bed.gz > ${sampleName}.all_regions.bed.gz
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
    samtools --version | (head -n 1 && exit 0) | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    aplanat --version | sed 's/ /,/' >> versions.txt
    python -c "import pysam; print(pysam.__version__)" | sed 's/^/pysam,/'  >> versions.txt
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


process mergeCSV {
    label "wfalignment"
    cpus 1
    input:
        path files
    output:
        path 'final_merged.csv'
    """
    awk '(NR == 1) || (FNR > 1)' $files > final_merged.csv
    """
}


process getRefNames {
    label "wfalignment"
    maxForks 1
    cpus 1
    input:
        each path(reference)
    output:
        path "*.txt"
    """
    seqkit fx2tab --name --only-id $reference > ${reference.baseName}.txt
    """
}


// workflow module
workflow pipeline {
    take:
        fastq
        references
        counts
    main:

        // Ready the optional file
        OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

        //uncompress and combine fastq's if multiple files
        if (params.concat_fastq){
            prepared_fastq = fastcatUncompress(fastq)
        }
        else{
            prepared_fastq = nameFastq(fastq)
        }

        // Cat the references together for alignment
        combined = combineReferences(references.collect())

        // Align the reads to produce bams and stats
        aligned = alignReads(
           prepared_fastq.fastq,
           references.collect(),
           combined,
           counts)
        merged_bams = mergeBAM(aligned.sorted)
        indexed_bam = indexBam(merged_bams)

        // Find reference lengths
        ref_length = refLengths(combined)

        // add steps column for use with mosdepth
        ref_steps = addStepsColumn(ref_length[0])

        // Find read_depth per reference/bam file
        depth_per_ref = file(OPTIONAL, type: "file")
        if (params.depth_coverage){
            depth_per_ref = readDepthPerRef(ref_steps, aligned.sorted)
        }

        // get params & versions
        workflow_params = getParams()
        software_versions = getVersions()

        //get references names faster that pysam
        ref = getRefNames(references.flatten())

        // Merge the stats together to get a csv out
        stats = gatherStats(alignReads.out.json, counts)
        sample_ids = prepared_fastq.sample_id.collectFile(
            name: 'sample_ids.csv', newLine: true)
        merged_csv = mergeCSV(stats.merged_mapula_csv.collect())

        report = plotStats(
            stats.merged_mapula_json.collect(), counts,
            ref.collect(),
            aligned.unmapped_stats.collect(),
            software_versions, workflow_params,
            sample_ids, depth_per_ref.collect())
    emit:
        merged = merged_bams
        indexed = indexed_bam
        merged_mapula_csv = merged_csv
        merged_mapula_json = stats.merged_mapula_json
        report = report.report
        telemetry = workflow_params
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    label "wfalignment"
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*", saveAs: {
        f -> params.prefix ? "${params.prefix}-${f}" : "${f}" }
    input:
        path fname
    output:
        path fname
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
    fastq = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "sanitize": params.sanitize_fastq,
        "output":params.out_dir])
    extensions = ["fasta", "fna", "ffn", "faa", "frn", "fa", "txt", "fa.gz", "fna.gz", "frn.gz", "ffn.gz", "fasta.gz"]

    // Acquire reference files
    input = file(params.references);
    if (input.isDirectory()) {
        reference_files = []
        references = file(params.references, type: "dir", checkIfExists:true);
        for (ext in extensions) {
            reference_files += file(references.resolve("*.${ext}"), type: 'file', maxdepth: 1)
        }
    }
    else {
        reference_files = file(params.references, type: "file", checkIfExists:true);
    }


    if (reference_files.size() == 0) {
            println('Error: No references found in the directory provided.')
            exit 1
    }
    else {
          reference_files = channel.fromPath(reference_files)
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
