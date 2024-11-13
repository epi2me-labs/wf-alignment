Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | wf-alignment-report.html | Report for all samples | aggregated |
| Combined references | combined_refs.fasta | FASTA file containing all input references. | aggregated |
| Combined references index | combined_refs.fasta.fai | Index file for combined references FASTA. | aggregated |
| Combined references MMI index | combined_refs.mmi | Minimap2 index file for combined references FASTA. | aggregated |
| Per-read alignment stats | {{ alias }}.readstats.tsv.gz | Bamstats per-read output TSV file (compressed with gzip). | per-sample |
| Per-reference alignment stats | {{ alias }}.flagstat.tsv | Bamstats flagstat output TSV file. | per-sample |
| Alignment accuracy histogram | {{ alias }}-histograms/accuracy.hist | Bamstats alignment accuracy histogram TSV file. | per-sample |
| Alignment coverage histogram | {{ alias }}-histograms/coverage.hist | Bamstats alignment coverage histogram TSV file. | per-sample |
| Read length histogram (mapped) | {{ alias }}-histograms/length.hist | Bamstats read length histogram TSV file (for mapped reads). | per-sample |
| Read length histogram (unmapped) | {{ alias }}-histograms/length.unmap.hist | Bamstats read length histogram TSV file (for unmapped reads). | per-sample |
| Read quality histogram (mapped) | {{ alias }}-histograms/quality.hist | Bamstats read quality histogram TSV file (for mapped reads). | per-sample |
| Read quality histogram (unmapped) | {{ alias }}-histograms/quality.unmap.hist | Bamstats read quality histogram TSV file (for unmapped reads). | per-sample |
| Alignments BAM file | {{ alias }}.sorted.aligned.bam | BAM file with alignments of filtered input reads against the combined references. | per-sample |
| Alignments index file | {{ alias }}.sorted.aligned.bam.bai | Index for alignments BAM file. | per-sample |
| IGV config JSON file | igv.json | JSON file with IGV config options to be used by the EPI2ME Desktop Application. | aggregated |
