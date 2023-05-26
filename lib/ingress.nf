import java.nio.file.NoSuchFileException

import ArgumentParser

EXTENSIONS = [
    fastq: ["fastq", "fastq.gz", "fq", "fq.gz"],
    bam: ["bam"],
    ubam: ["ubam"]
]

/**
 * Take a map of input arguments, find valid inputs, and return a channel
 * with elements of `[metamap, path-to-input-file]`. The file might be a concatenation
 * of multiple files if the original input was a top-level directory or a directory with
 * sub-directories which in turn contain input files. Can find either FASTQ, BAM, or
 * uBAM files.
 *
 * @param arguments: map with arguments containing
 *  - "input": path to either: (i) single input file, (ii) top-level directory
 *     containing input files, (iii) directory containing sub-directories which contain
 *     input files
 *  - "input_type": string of either "fastq", "bam", or "ubam".
 *  - "sample": string to name a single sample
 *  - "sample_sheet": path to CSV sample sheet
 *  - "analyse_unclassified": boolean whether to keep unclassified reads
 * @return Channel of `[Map(alias, barcode, type, ...), Path]`. The first element is a
 *  map with metadata and the second is the path to the (potentially concatenated) file.
 */
def ingress(Map arguments) {
    // check arguments
    Map margs = parse_arguments(arguments)
    if (!EXTENSIONS.containsKey(margs["input_type"])) {
        error "Input type needs to be one of ${EXTENSIONS.keySet()}"
    }
    ArrayList extensions = EXTENSIONS[margs["input_type"]]
    // Get a channel with valid input data [metamap, input_path]. It will be filled by
    // the data of the three input types (single file or dir with fastq or subdirs with
    // fastq).
    def ch_input = get_valid_inputs(margs, extensions).map {
        // we could have directories with a single valid file in the input --> "unwrap"
        // these files so that we only have single files or directories containing
        // multiple files in the channel
        meta, path ->
        if (path.isDirectory()) {
            List fq_files = get_fq_files_in_dir(path, extensions)
            if (fq_files.size() == 1) {
                path = fq_files[0]
            }
        }
        [meta, path]
    }
    if (margs.input_type == "fastq") {
        return fastcat_or_mv(ch_input)
    } else {
        // BAM or uBAM
        return concatBam_or_mv(ch_input, extensions)
    }
}

/**
 * Process to either (i) run `fastcat` on a directory with FASTQ files or (ii) rename or
 * compress a single FASTQ file.
 *
 * @param: tuple of metadata and path to input directory / file
 * @return: tuple of metadata and path to concatenated / renamed / gzipped file
 */
process fastcat_or_mv {
    label "wfalignment"
    cpus 3
    input:
        tuple val(meta), path(input)
    output:
        tuple val(meta), path(out)
    script:
        out = "reads.fastq.gz"
        int bgzip_threads = task.cpus - 1
        if (input.isFile()) {
            // if the file is already gzipped, only rename; otherwise compress
            if (input.name.endsWith(".gz")) {
                // we need to take into account that the file could already be named
                // "reads.fastq.gz" in which case `mv` would fail
                """
                [ "$input" == "$out" ] || mv $input $out
                """
            } else {
                """
                cat $input | bgzip -@ $bgzip_threads > $out
                """
            }
        } else if (input.isDirectory()) {
            """
            fastcat $input | bgzip -@ $bgzip_threads > $out
            """
        } else {
            error "$input is neither file nor directory."
        }
}

/**
 * Process to either (i) concatenate BAM or uBAM files in a directory with
 * `samtools cat` or (ii) rename a single BAM / uBAM file.
 *
 * @param: tuple of metadata and path to input directory / file
 * @param: List of allowed extensions (e.g. ['bam'])
 * @return: tuple of metadata and path to concatenated / renamed file
 */
process concatBam_or_mv {
    label "wfalignment"
    cpus 1
    input:
        tuple val(meta), path(input)
        val(extensions)
    output:
        tuple val(meta), path(out)
    script:
        out = "reads.bam"
        exts = "{" + extensions.join(",") + "}"
        if (input.isFile()) {
            """
            mv $input $out
            """
        } else if (input.isDirectory()) {
            """
            samtools cat -o $out ${input}/*.$exts
            """
        } else {
            error "$input is neither file nor directory."
        }
}


/**
 * Parse input arguments for `fastq_ingress`.
 *
 * @param arguments: map with input arguments (see `fastq_ingress` for details)
 * @return: map of parsed arguments
 */
Map parse_arguments(Map arguments) {
    def parser = new ArgumentParser(
        args:["input", "input_type"],
        kwargs:["sample": null,
                "sample_sheet": null,
                "analyse_unclassified": false,
                "fastcat_stats": false,
                "watch_path": false],
        name: "fastq_ingress")
    return parser.parse_args(arguments)
}


/**
 * Find valid inputs based on the input type.
 *
 * @param margs: parsed arguments (see `fastq_ingress` for details)
 * @return: channel of `[metamap, input-path]`; `input-path` can be the path to
 *  a single FASTQ file or to a directory containing FASTQ files
 */
def get_valid_inputs(Map margs, ArrayList extensions){
    log.info "Checking fastq input."
    Path input
    try {
        input = file(margs.input, checkIfExists: true)
    } catch (NoSuchFileException e) {
        error "Input path $margs.input does not exist."
    }
    boolean dir_has_subdirs = false
    boolean dir_has_fastq_files = false
    // declare resulting input channel
    def ch_input
    // handle case of `input` being a single file
    if (input.isFile()) {
        // the `fastcat` process can deal with directories or single file inputs
        ch_input = Channel.of(
            [create_metamap([alias: margs["sample"] ?: input.simpleName]), input])
    } else if (input.isDirectory()) {
        // input is a directory --> we accept two cases: (i) a top-level directory with
        // fastq files and no sub-directories or (ii) a directory with one layer of
        // sub-directories containing fastq files
        dir_has_fastq_files = get_fq_files_in_dir(input, extensions) as boolean
        // find potential subdirectories (this list can be empty)
        ArrayList sub_dirs = file(input.resolve("*"), type: "dir")
        dir_has_subdirs = sub_dirs as boolean
        // deal with first case (top-lvl dir with fastq files and no sub-directories)
        if (dir_has_fastq_files){
            if (dir_has_subdirs) {
                error "Input directory '$input' cannot contain " +
                    "${margs["input_type"]} files and sub-directories."
            }
            ch_input = Channel.of(
                [create_metamap([alias: margs["sample"] ?: input.name]), input])
        } else {
            // deal with the second case (sub-directories with fastq data) --> first
            // check whether we actually found sub-directories
            if (!dir_has_subdirs) {
                error "Input directory '$input' must contain either " +
                    "${margs["input_type"]} files or sub-directories."
            }
            // make sure that there are no sub-sub-directories with FASTQ files and that
            // the sub-directories actually contain fastq files)
            if (sub_dirs.any {
                def subsubdirs = file(it.resolve("*"), type: "dir")
                subsubdirs.any { get_fq_files_in_dir(it, extensions) }
            }) {
                error "Input directory '$input' cannot contain more than one level " +
                    "of sub-directories with ${margs["input_type"]} files."
            }
            ArrayList sub_dirs_with_fastq_files = sub_dirs.findAll {
                get_fq_files_in_dir(it, extensions) as boolean
            }
            if (!sub_dirs_with_fastq_files) {
                error "No ${margs["input_type"]} files were found in the " +
                    "sub-directories of '$input'."
            }
            // remove directories called 'unclassified' unless otherwise specified
            if (!margs.analyse_unclassified) {
                sub_dirs_with_fastq_files = sub_dirs_with_fastq_files.findAll {
                    it.name != "unclassified"
                }
            }
            // filter based on sample sheet in case one was provided
            if (margs.sample_sheet) {
                // get channel of entries in the sample sheet
                def ch_sample_sheet = get_sample_sheet(file(margs.sample_sheet))
                // get the intersection of both channels
                def ch_intersection = Channel.fromPath(sub_dirs_with_fastq_files).map {
                    [it.name, it]
                }.join(ch_sample_sheet.map{[it.barcode, it]}, remainder: false)
                // TODO: we should let the user know if a sample present in the sample
                // sheet didn't have a matching sub-directory. We could throw an error
                // or print a warning, but ideally we would return an extra channel from
                // `fastq_ingress` with a summary about what has been found (in terms of
                // sample sheet and sub-dirs) --> we'll do this later

                // put metadata into map and return
                ch_input = ch_intersection.map {barcode, path, sample_sheet_entry -> [
                    create_metamap(sample_sheet_entry), path
                ]}
            } else {
                ch_input = Channel.fromPath(sub_dirs_with_fastq_files).map {
                    [create_metamap([alias: it.name]), it]
                }
            }
        }
    } else {
        error "Input $input appears to be neither a file nor a directory."
    }
    // a sample sheet only makes sense in the case of a directory with
    // sub-directories
    if (margs.sample_sheet && !dir_has_subdirs) {
        error "Sample sheet was provided, but input does not contain " +
            "sub-directories with ${margs["input_type"]} files."
    }
    return ch_input
}


/**
 * Create a map that contains at least these keys: `[alias, barcode, type]`.
 * `alias` is required, `barcode` and `type` are filled with default values if
 * missing. Additional entries are allowed.
 *
 * @param kwargs: map with input parameters; must contain `alias`
 * @return: map(alias, barcode, type, ...)
 */
Map create_metamap(Map arguments) {
    def parser = new ArgumentParser(
        args: ["alias"],
        kwargs: [
            "barcode": null,
            "type": "test_sample",
        ],
        name: "create_metamap",
    )
    return parser.parse_known_args(arguments)
}


/**
 * Get the fastq files in the directory (non-recursive).
 *
 * @param dir: path to the target directory
 * @return: list of found fastq files
 */
ArrayList get_fq_files_in_dir(Path dir, ArrayList extensions) {
    return extensions.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
}


/**
 * Check the sample sheet and return a channel with its rows if it is valid.
 *
 * @param sample_sheet: path to the sample sheet CSV
 * @return: channel of maps (with values in sample sheet header as keys)
 */
def get_sample_sheet(Path sample_sheet) {
    // If `validate_sample_sheet` does not return an error message, we can assume that
    // the sample sheet is valid and parse it. However, because of Nextflow's
    // asynchronous magic, we might emit values from `.splitCSV()` before the
    // error-checking closure finishes. This is no big deal, but undesired nonetheless
    // as the error message might be overwritten by the traces of new nextflow processes
    // in STDOUT. Thus, we use the somewhat clunky construct with `concat` and `last`
    // below. This lets the CSV channel only start to emit once the error checking is
    // done.
    ch_err = validate_sample_sheet(sample_sheet).map {
        // check if there was an error message
        if (it) error "Invalid sample sheet: $it"
        it
    }
    // concat the channel holding the path to the sample sheet to `ch_err` and call
    // `.last()` to make sure that the error-checking closure above executes before
    // emitting values from the CSV
    return ch_err.concat(Channel.fromPath(sample_sheet)).last().splitCsv(
        header: true, quote: '"'
    )
}


/**
 * Python script for validating a sample sheet. The script will write messages
 * to STDOUT if the sample sheet is invalid. In case there are no issues, no
 * message is emitted.
 *
 * @param: path to sample sheet CSV
 * @return: string (optional)
 */
process validate_sample_sheet {
    label "wfalignment"
    input: path csv
    output: stdout
    """
    workflow-glue check_sample_sheet $csv
    """
}
