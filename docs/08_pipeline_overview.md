### 1. Combine reference files

All reference files in the directory passed to `--references` are concatenated.

### 2. Align reads

Input reads are aligned against the combined reference with [Minimap2](https://github.com/lh3/minimap2). If BAM files are used as input (with `--bam`), only reads in files without a reference in the SAM header are aligned. For other BAM files this step is skipped.

### 3. Create alignment stats

[Bamstats](https://github.com/epi2me-labs/fastcat#bamstats) is used to create per-read and per-reference alignment stats from the BAM files.

### 4. Calculate depth of coverage

Depth of coverage along the reference sequences is determined with [Mosdepth](https://github.com/brentp/mosdepth) (using 200 windows per reference sequence). To speed up the workflow, this step can be skipped by adding `--depth-coverage false`.
