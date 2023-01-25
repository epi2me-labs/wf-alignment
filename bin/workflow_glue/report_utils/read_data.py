"""Read data for report."""
import os
import sys

import pandas as pd
from pandas.api import types as pd_types

CATEGORICAL = pd_types.CategoricalDtype(ordered=True)


def bamstats(stats_dir):
    """Read bamstats per-read stats.

    :param stats_dir: directory with bamstats per-read stats output files
    :return: `pd.DataFrame` with bamstats per-read stats data
    """
    # these files can be quite large; so only keep relevant columns and store the string
    # columns as categorical
    relevant_stats_cols_dtypes = {
        "name": str,
        "sample_name": CATEGORICAL,
        "ref": CATEGORICAL,
        "coverage": float,
        "ref_coverage": float,
        "read_length": int,
        "mean_quality": float,
        "acc": float,
    }
    input_files = os.listdir(stats_dir)
    dfs = []
    for fname in input_files:
        df = pd.read_csv(
            f"{stats_dir}/{fname}",
            sep="\t",
            index_col=0,
            usecols=relevant_stats_cols_dtypes,
            dtype=relevant_stats_cols_dtypes,
        )
        if df.empty:
            sys.stderr.write(
                f"Per-read stats file '{fname}' is empty."
                "Were there reads in the corresponding input file?"
            )
            continue
        dfs.append(df)
    # check that we actually got any reads
    if not dfs:
        raise ValueError(
            "No reads were found in any of the per-read stats files "
            f"{input_files}. Could it be that all original input files were empty?"
        )
    # pd.concat will revert the `dtype` of categorical columns back to `object` if the
    # concatenated columns don't contain the same set of items --> we use
    # `union_categoricals` to avoid that
    cat_cols = [k for k, v in relevant_stats_cols_dtypes.items() if v == CATEGORICAL]
    for col in cat_cols:
        uc = pd_types.union_categoricals([df[col] for df in dfs], sort_categories=True)
        for df in dfs:
            df[col] = pd.Categorical(df[col], categories=uc.categories, ordered=True)
    return pd.concat(dfs)


def flagstat(flagstat_dir):
    """Read bamstats per-file stats.

    :param stats_dir: directory with bamstats per-file stats output files
    :return: `pd.DataFrame` with bamstats per-file stats data
    """
    input_files = os.listdir(flagstat_dir)
    dfs = []
    for fname in input_files:
        df = pd.read_csv(f"{flagstat_dir}/{fname}", sep="\t")
        if df.empty:
            sys.stderr.write(
                f"Per-file stats file '{fname}' is empty."
                "Were there reads in the corresponding input file?"
            )
            continue
        dfs.append(df)
    if not dfs:
        raise ValueError(
            "No entries were found in any of the per-file stats files "
            f"{input_files}. Could it be that all original input files were empty?"
        )
    return pd.concat(dfs)


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


def depths(depths_dir, sample_names):
    """Read depth data.

    :param depths_dir: directory with `mosdepth` output files
    :param sample_names: list of sample names for which we found reads
    :return: `pd.DataFrame` with columns `["ref", "start", "end", "depth",
        "sample_name"]`
    """
    dfs = []
    for sample_name in sample_names:
        (fname,) = [x for x in os.listdir(depths_dir) if x.startswith(sample_name)]
        df = pd.read_csv(
            f"{depths_dir}/{fname}",
            sep="\t",
            header=None,
            names=["ref", "start", "end", "depth"],
        )
        dfs.append(df.eval(f'sample_name = "{sample_name}"'))
    return pd.concat(dfs).astype({"sample_name": CATEGORICAL, "ref": CATEGORICAL})


def counts(counts_file):
    """Read expected counts data.

    :param counts_file: CSV file with expected counts (needs columns `["reference",
        "expected_count"]`)
    :raises ValueError: throw error if one of the required columns is missing
    :return: `pd.Series` with the expected counts
    """
    counts = pd.read_csv(counts_file)
    counts.columns = [col.lower() for col in counts.columns]
    if not ("reference" in counts.columns and "expected_count" in counts.columns):
        raise ValueError(
            (
                "Counts CSV must (at least) contain the columns 'reference' and "
                "'expected_count' (capitalisation is ignored)."
            )
        )
    return counts.set_index("reference")["expected_count"].squeeze()
