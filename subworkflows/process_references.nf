process combine {
    label "wfalignment"
    cpus 1
    memory "2 GB"
    input: path "reference*.fasta"
    output:
        path outfname, emit: fasta
        path "*fai", emit: index
    script:
    outfname = "combined_refs.fasta"
    """
    find -name 'reference*.fasta' -exec zcat -f {} + > $outfname

    # make sure all sequence IDs are unique
    if [ "\$(grep "^>" $outfname | sort | uniq -d)" ]; then
        echo "Sequence IDs in the reference files must be unique." 1>&2
        exit 1
    fi
    samtools faidx $outfname
    """
}

process fx2tab {
    label "wfalignment"
    cpus 1
    memory { reference.size() > 1e9 ? "16 GB" : "2 GB" }
    input:
        path reference
    output:
        path "*.names.txt", emit: names
        path "*.lengths.tsv", emit: lengths
    script:
    """
    seqkit fx2tab --length --name --only-id $reference > fx2tab.out
    cut -f1 fx2tab.out > ${reference}.names.txt
    echo -e 'name\\tlengths' > ${reference}.lengths.tsv
    cat fx2tab.out >> ${reference}.lengths.tsv
    """
}

workflow process_references {
    List extensions = [
        "fasta", "fna", "ffn", "faa", "frn", "fa", "txt",
        "fa.gz", "fna.gz", "frn.gz", "ffn.gz", "fasta.gz"
    ]
    take:
        input
    main:
        // get the reference files
        Path input = file(input, checkIfExists: true)
        List ref_files
        if (input.isDirectory()) {
            // we got a directory with one or multiple references
            ref_files = extensions.collect {
                file(input.resolve("*.$it"), type: "file")
            }.flatten()
            if (ref_files.size() == 0) {
                    error "No references found in ${input}."
            }
        }
        else {
            // the reference is a single file --> make sure it has an expected extension
            if (!extensions.any { input.name.endsWith(it) }) {
                error "Reference file $input does not have an " +
                    "accepted extension ($extensions)."
            }
            ref_files = [input]
        }
        fx2tab(Channel.of(ref_files).flatten())
        combine(ref_files)
    emit:
        combined = combine.out.fasta
        combined_index = combine.out.index
        names_per_ref_file = fx2tab.out.names
        lengths_per_ref_file = fx2tab.out.lengths
        lengths_combined = fx2tab.out.lengths.collectFile(
            name: "combined_lengths.tsv", keepHeader: true
            // we need to call `.first()` to get a value channel (`.collectFile()`
            // always returns a queue channel, even when it only produces a single file)
        ).first()
}
