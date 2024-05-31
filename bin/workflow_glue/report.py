#!/usr/bin/env python
"""Report using ezcharts."""

from pathlib import Path

from ezcharts.components.reports import labs
import pandas as pd

from .report_utils import read_data, sanitizer, sections  # noqa: ABS101

from .util import get_named_logger, wf_parser  # noqa: ABS101

sanitizer = sanitizer.Sanitizer()


def main(args):
    """Run entry point."""
    logger = get_named_logger("report")

    per_sample_dirs = sorted(args.data.glob("*"))

    flagstat_df = pd.concat(
        (
            read_data.flagstat(sample_dir).assign(sample_name=sample_dir.name)
            for sample_dir in per_sample_dirs
        )
    ).astype({"sample_name": "category"})

    # read the ref names (i.e. get a dict mapping ref names to the ref file)
    refname2reffile = read_data.refnames(args.refnames_dir)
    ref_files = sorted(set(refname2reffile.values()) - set(["unmapped"]))
    ref_seqs = sorted(set(refname2reffile.keys()) - set(["*"]))
    # read depth info if available
    try:
        depth_df = pd.concat(
            (read_data.depths(d).assign(sample_name=d.name) for d in per_sample_dirs)
        ).astype({"sample_name": "category"})
    except AttributeError:
        depth_df = None
    # read counts if available
    counts = read_data.counts(args.counts) if args.counts is not None else None

    # add a column with the respective ref. files to the stats dataframes
    for df in (flagstat_df, depth_df):
        if df is None:
            continue
        try:
            df["ref_file"] = (
                df["ref"].apply(lambda ref: refname2reffile[ref]).astype("category")
            )
        except KeyError as e:
            (missing_ref,) = e.args
            raise ValueError(
                f"Reference '{missing_ref}' not found in the provided "
                f"reference files {ref_files}."
            )

    # create the report
    report = labs.LabsReport(
        "wf-alignment report",
        "wf-alignment",
        args.params,
        args.versions,
        workflow_version=args.wf_version,
    )
    # add sections
    sections.summary(report, per_sample_dirs, ref_files, ref_seqs, flagstat_df)
    sections.seqsum(report, per_sample_dirs)
    if depth_df is not None:
        sections.depths(report, depth_df)
    if counts is not None:
        sections.counts(report, flagstat_df, counts, sanitizer)
    # "de-sanitize" the report and write to the output file
    report_fname = "wf-alignment-report.html"
    with open(report_fname, "w") as report_file:
        report_file.write(sanitizer.desanitise_report(str(report)))

    logger.info(f"Written sanitized report to '{report_fname}'.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "--data",
        type=Path,
        help="directory with per-sample data (with a sub-directory for each sample)",
    )
    parser.add_argument(
        "--refnames_dir",
        help="directory with files containing reference names",
    )
    parser.add_argument(
        "--counts",
        required=False,
        help=(
            "CSV file with expected counts "
            "(columns: Reference, expected_count, expected_length)"
        ),
    )
    parser.add_argument(
        "--params",
        default=None,
        help="CSV file with workflow parameters",
    )
    parser.add_argument(
        "--versions",
        help="CSV file with software versions",
    )
    parser.add_argument(
        "--wf-version",
        default="unknown version",
        help="Workflow version",
    )
    return parser
