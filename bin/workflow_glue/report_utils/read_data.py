"""Read data for report."""

import os

from ezcharts.components import fastcat
import pandas as pd


def length_hist(data_dir):
    """Load bamstats length histogram data."""
    lengths_mapped = fastcat.load_histogram(data_dir, "length")
    try:
        lengths_unmapped = fastcat.load_histogram(data_dir, "length.unmap")
    except FileNotFoundError:
        lengths_unmapped = None
    return fastcat.sum_hists((lengths_mapped, lengths_unmapped))


def flagstat(sample_dir):
    """Load bamstats flagstat."""
    return fastcat.load_bamstats_flagstat(sample_dir / "bamstats.flagstat.tsv")


def depths(data_dir):
    """Parse mosdepth results file if there is one."""
    depths_file = data_dir / "depth.all_regions.bed.gz"
    if not depths_file.exists():
        return None
    depths = pd.read_csv(
        depths_file,
        sep="\t",
        header=None,
        names=["ref", "start", "end", "depth"],
        dtype={"ref": str, "start": int, "end": int, "depth": float},
    )
    return depths


def refnames(refnames_dir):
    """Read files mapping reference sequence IDs to reference file names.

    :param refnames_dir: directory containing files with ref. names
    :return: `dict` mapping reference sequence IDs to reference files
    """
    refname2reffile = {"*": "unmapped"}
    for ref_name_file in os.listdir(refnames_dir):
        with open(f"{refnames_dir}/{ref_name_file}", "r") as f:
            for ref_name in f:
                refname2reffile[ref_name.strip()] = os.path.basename(
                    ref_name_file.split(".names")[0]
                )
    return refname2reffile


def counts(counts_file):
    """Read expected counts data.

    :param counts_file: CSV file with expected counts (needs columns `["reference",
        "expected_count"]`)
    :raises ValueError: throw error if one of the required columns is missing
    :return: `pd.Series` with the expected counts
    """
    counts = pd.read_csv(
        counts_file,
        dtype={"Reference": str, "expected_count": float})
    counts.columns = [col.lower() for col in counts.columns]
    if not ("reference" in counts.columns and "expected_count" in counts.columns):
        raise ValueError(
            (
                "Counts CSV must (at least) contain the columns 'reference' and "
                "'expected_count' (capitalisation is ignored)."
            )
        )
    counts["reference"] = counts["reference"].astype(str)
    return counts.set_index("reference")["expected_count"].squeeze()
