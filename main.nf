#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress; } from "./lib/ingress"
include { process_references } from "./subworkflows/process_references"


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")
MINIMAP_ARGS_PRESETS = [
    "dna": "-ax map-ont",
    "rna": "-ax splice -uf"
]

process alignReads {
    label "wfalignment"
    cpus params.mapping_threads + params.sorting_threads
    input:
        tuple val(meta), path(input)
        path combined_refs
        val is_xam
        val minimap_args
    output:
        tuple val(meta), path(bam_name)
    script:
        def sample_name = meta["alias"]
        bam_name = "${sample_name}.sorted.aligned.bam"
    """
    ${is_xam ? "samtools fastq -T '*' $input" : "cat $input"} \
    | minimap2 -t $params.mapping_threads $minimap_args $combined_refs - \
    | samtools sort -@ ${params.sorting_threads - 1} -o $bam_name -
    """
}

process indexBam {
    label "wfalignment"
    cpus 1
    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path(bam), path("*.bai")
    script:
    """
    samtools index $bam
    """
}

process bamstats {
    label "wfalignment"
    cpus 2
    input:
        tuple val(meta), path(bam), path(index)
    output:
        path "*.readstats.tsv", emit: read_stats
        path "*.flagstat.tsv", emit: flagstat
    script:
        def sample_name = meta["alias"]
    """
    bamstats $bam -s $sample_name -u -f ${sample_name}.flagstat.tsv -t $task.cpus \
    > ${sample_name}.readstats.tsv
    """
}

process addStepsColumn {
    // TODO: we don't need 200 windows for very short references; find heuristics for
    // determining window length / number for such cases
    label "wfalignment"
    cpus 1
    input: path "lengths.csv"
    output: path "lengths_with_steps.csv"
    """
    #!/usr/bin/env python
    import pandas as pd
    all = pd.read_csv('lengths.csv')
    all["step"] = all["lengths"]//200
    all = all.replace(0, 1)
    all.to_csv('lengths_with_steps.csv', index=False, header=False)
    """
}

process readDepthPerRef {
    // TODO: check if parallelisation with `xargs` or `parallel` is more efficient
    label "wfalignment"
    cpus 3
    input:
        tuple val(meta), path(alignment), path(index)
        path ref_len
    output:
        path outfname
    script:
        def sample_name = meta["alias"]
        outfname = "${sample_name}.all_regions.bed.gz"
    """
    while IFS=, read -r name lengths steps; do
        mosdepth -n --fast-mode --by "\$steps" --chrom "\$name" -t $task.cpus \
            ${sample_name}."\$name".temp $alignment \
        || echo "No alignments for "\$name""
    done < $ref_len

    # use `find` instead of a glob as this cannot fail with `Argument list too long`
    find *.regions.bed.gz -exec cat {} + > $outfname
    # remove all the temp files
    rm -rf ${sample_name}.*.temp*
    """
}

process makeReport {
    label "wfalignment"
    cpus 1
    input:
        path "readstats/*"
        path "flagstat/*"
        path "refnames/*"
        path depths, stageAs: "depths/*"
        path counts
        path versions
        path params
    output:
        path "*.html"
    script:
    String depth_args = "--depths_dir depths"
    // we need to check against `.baseName` here because Nextflow includes the staging
    // directory in the `.name` of a `TaskPath`
    if (!(depths instanceof List) && depths.baseName == OPTIONAL_FILE.name) {
        depth_args = ""
    }
    String counts_args = (counts.name == OPTIONAL_FILE.name) ? "" : "--counts $counts"
    """
    workflow-glue report \
        --name wf-alignment \
        --stats_dir readstats \
        --flagstat_dir flagstat \
        --refnames_dir refnames \
        --versions $versions \
        --params $params \
        $depth_args \
        $counts_args
    """
}


process getVersions {
    label "wfalignment"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python --version | tr -s ' ' ',' | tr '[:upper:]' '[:lower:]' > versions.txt
    seqkit version | sed 's/ /,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | (head -n 1 && exit 0) | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    ezcharts --version | sed 's/ /,/' >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    bgzip --version | head -n1 | sed -E 's/(.*) /\\1,/' >> versions.txt
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


// workflow module
workflow pipeline {
    take:
        sample_data
        refs
        counts
        depth_coverage
    main:
        // get params & versions
        workflow_params = getParams()
        software_versions = getVersions()

        // handle references
        refs = process_references(params.references)

        sample_data = sample_data
        | map { meta, path, stats -> [meta, path] }

        String minimap_args

        if (params.bam) {
            ch_branched = sample_data.branch { meta, bam ->
                to_align: meta["is_unaligned"]
                aligned: true
            }
            ch_to_align = ch_branched.to_align
            // `xam_ingress` sorts the BAMs, so we don't have to
            bam = ch_branched.aligned
        } else {
            // FASTQ input
            ch_to_align = sample_data
            bam = Channel.empty()
        }

        // run minimap
        if (! MINIMAP_ARGS_PRESETS.containsKey(params.minimap_preset)) {
            error "'--minimap_preset' needs to be one of " +
                "${MINIMAP_ARGS_PRESETS.keySet()}."
        }
        minimap_args = params.minimap_args ?: \
            MINIMAP_ARGS_PRESETS[params.minimap_preset]
        bam = bam
        | mix(
            alignReads(ch_to_align, refs.combined, params.bam as boolean, minimap_args)
        )
        | indexBam

        // get stats
        stats = bamstats(bam)

        // determine read_depth per reference / bam file if requested
        depth_per_ref = Channel.of(OPTIONAL_FILE)
        if (depth_coverage) {
            // add step column to ref lengths
            ref_lengths_with_steps = addStepsColumn(refs.lengths_combined)
            depth_per_ref = readDepthPerRef(bam, ref_lengths_with_steps)
        }

        report = makeReport(
            stats.read_stats.collect(),
            stats.flagstat.collect(),
            refs.names_per_ref_file.collect(),
            depth_per_ref.collect(),
            counts,
            software_versions,
            workflow_params,
        )
    emit:
        alignments = bam.map { it[1] }
        indices = bam.map{ it[2] }
        per_read_stats = stats.read_stats
        per_file_stats = stats.flagstat
        report
        params_json = workflow_params
        software_versions
        combined_ref = refs.combined
        combined_ref_index = refs.combined_index
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


process configure_jbrowse {
    label "wfalignment"
    input:
        path(alignments)
        path(indexes)
        path(reference)
        path(ref_idx)
    output:
        path("jbrowse.json")
    script:
    ArrayList alignment_args = []
    int i = 0;
    for(a in alignments) {
        // don't be fooled into iterating over bam.size() here
        // when the cardinality is 1, bam.size() returns the filesize of the bam!
        this_bam = a
        this_bai = indexes[i]
        alignment_args << "--alignment '${params.out_dir}/${this_bam.name}' '${params.out_dir}/${this_bai.name}'"
        i++;
    }
    String alignment_args_str = alignment_args.join(' ')
    """
    workflow-glue configure_jbrowse \
        --reference '${reference}' '${params.out_dir}/${reference.name}' '${params.out_dir}/${ref_idx.name}' \
        ${alignment_args_str} > jbrowse.json
    """
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

    Map ingress_args = [
        "sample": params.sample,
        "sample_sheet": params.sample_sheet,
        "analyse_unclassified": params.analyse_unclassified,
        "stats": false,
    ]

    // get input data
    if (params.fastq) {
        sample_data = fastq_ingress(ingress_args + ["input": params.fastq])
    } else {
        sample_data = xam_ingress(
            ingress_args + ["input": params.bam, "keep_unaligned": true]
        )
    }

    counts = file(params.counts ?: OPTIONAL_FILE, checkIfExists: true)

    // Run pipeline
    results = pipeline(
        sample_data, params.references, counts, params.depth_coverage
    )

    // create jbrowse file
    jb2_conf = configure_jbrowse(
        results.alignments.collect(),
        results.indices.collect(),
        results.combined_ref,
        results.combined_ref_index
    )
    output(jb2_conf.concat(results))
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
