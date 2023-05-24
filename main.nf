#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { ingress } from "./lib/ingress"
include { process_references } from "./subworkflows/process_references"


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")
MINIMAP_ARGS_PRESETS = [
    "dna": "-ax map-ont",
    "rna": "-ax splice -uf"
]

process alignReads {
    label "wfalignment"
    cpus params.threads
    input:
        tuple val(meta), path(input)
        path combined_refs
        val input_type
        val minimap_args
    output:
        tuple val(meta), path(bam_name)
    script:
        def sample_name = meta["alias"]
        bam_name = "${sample_name}.sorted.aligned.bam"
        if (!(input_type in ["fastq", "ubam"])) {
            error "`input_type` must be either 'fastq' or 'ubam'."
        }
    """
    ${(input_type == "fastq") ? "cat $input" : "samtools fastq -T '*' $input"} \
    | minimap2 -t $task.cpus $minimap_args $combined_refs - \
    | samtools sort -@ $task.cpus -o $bam_name -
    """
}

process sortInputBam {
    label "wfalignment"
    cpus params.threads
    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path(sorted_bam_name)
    script:
        def sample_name = meta["alias"]
        sorted_bam_name = "${sample_name}.sorted.aligned.bam"
    """
    samtools sort -@ $task.cpus -o $sorted_bam_name $bam
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
    cpus params.threads
    input:
        tuple val(meta), path(bam), path(index)
    output:
        path "*.readstats.tsv", emit: read_stats
        path "*.flagstat.tsv", emit: flagstat
    script:
        def sample_name = meta["alias"]
    """
    bamstats $bam -s $sample_name -u -f ${sample_name}.flagstat.tsv \
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
    depth_threads = {params.threads >= 4  ? 4 : params.threads}
    label "wfalignment"
    cpus depth_threads
    input:
        tuple val(meta), path(alignment)
        path ref_len
    output:
        path outfname
    script:
        def sample_name = meta["alias"]
        outfname = "${sample_name}.all_regions.bed.gz"
    """
    samtools index $alignment
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
    cpus params.threads
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
        input_type
        refs
        counts
        depth_coverage
    main:
        // get params & versions
        workflow_params = getParams()
        software_versions = getVersions()

        // handle references
        refs = process_references(params.references)

        def bam
        String minimap_args
        if (input_type in ["fastq", "ubam"]) {
            // run minimap
            if (! MINIMAP_ARGS_PRESETS.containsKey(params.minimap_preset)) {
                error "'--minimap_preset' needs to be one of " +
                    "${MINIMAP_ARGS_PRESETS.keySet()}."
            }
            minimap_args = params.minimap_args ?: \
                MINIMAP_ARGS_PRESETS[params.minimap_preset]
            bam = alignReads(sample_data, refs.combined, input_type, minimap_args)
        } else {
            // input is bam
            bam = sortInputBam(sample_data)
        }

        // index the bam
        bam = indexBam(bam)

        // get stats
        stats = bamstats(bam)

        // determine read_depth per reference / bam file if requested
        depth_per_ref = Channel.of(OPTIONAL_FILE)
        if (depth_coverage) {
            // add step column to ref lengths
            ref_lengths_with_steps = addStepsColumn(refs.lengths_combined)
            depth_per_ref = readDepthPerRef(
                bam.map { it[0..1] }, ref_lengths_with_steps
            )
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
        alignment_args << "--alignment ${params.out_dir}/${this_bam.name} ${params.out_dir}/${this_bai.name}"
        i++;
    }
    String alignment_args_str = alignment_args.join(' ')
    """
    workflow-glue configure_jbrowse \
        --reference ${reference} ${params.out_dir}/${reference.name} ${params.out_dir}/${ref_idx.name} \
        ${alignment_args_str} > jbrowse.json
    """
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

    // macke sure that one of `--fastq`, `--bam`, or `--ubam` was given
    def input_type = ['fastq', 'bam', 'ubam'].findAll { params[it] }
    if (input_type.size() != 1) {
        error "Only provide one of '--fastq', '--bam' or '--ubam'."
    }

    input_type = input_type[0]

    // get input data
    sample_data = ingress([
        "input": params[input_type],
        "input_type": input_type,
        "sample": params.sample,
        "sample_sheet": params.sample_sheet,
        "analyse_unclassified": params.analyse_unclassified,
        "fastcat_stats": false
    ])

    counts = file(params.counts ?: OPTIONAL_FILE, checkIfExists: true)

    // Run pipeline
    results = pipeline(
        sample_data, input_type, params.references, counts, params.depth_coverage
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

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }
    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
