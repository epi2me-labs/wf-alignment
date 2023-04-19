#!/usr/bin/env python
"""Report using ezcharts."""

from ezcharts.components.reports import labs

from .report_utils import read_data, sanitizer, sections  # noqa: ABS101

from .util import get_named_logger, wf_parser  # noqa: ABS101

sanitizer = sanitizer.Sanitizer()


def main(args):
    """Run entry point."""
    logger = get_named_logger("report")

    # read input stats data
    stats_df = read_data.bamstats(args.stats_dir)
    sample_names = stats_df["sample_name"].cat.categories
    flagstat_df = read_data.flagstat(args.flagstat_dir)
    # read the ref names (i.e. get a dict mapping ref names to the ref file)
    refname2reffile = read_data.refnames(args.refnames_dir)
    ref_files = sorted(set(refname2reffile.values()) - set(['unmapped']))
    ref_seqs = sorted(set(refname2reffile.keys()) - set(['*']))
    # read depth info if available
    depth_df = None
    if args.depths_dir is not None:
        sample_names = stats_df["sample_name"].cat.categories
        depth_df = read_data.depths(args.depths_dir, sample_names)
    # read counts if available
    counts = read_data.counts(args.counts) if args.counts is not None else None

    # add a column with the respective ref. files to the stats dataframes
    for df in (stats_df, flagstat_df, depth_df):
        if df is None:
            continue
        try:
            df["ref_file"] = (
                df["ref"]
                .apply(lambda ref: refname2reffile[ref])
                .astype(read_data.CATEGORICAL)
            )
        except KeyError as e:
            (missing_ref,) = e.args
            raise ValueError(
                f"Reference '{missing_ref}' not found in the provided "
                f"reference files {ref_files}."
            )

    # extract the mapped reads and some other metrics used in the report sections
    stats_df_mapped = stats_df.query('ref != "*"')

    # create the report
    report = labs.LabsReport(
        f"{args.name} report",
        args.name,
        args.params,
        args.versions,
    )
    # add sections
    sections.summary(report, sample_names, ref_files, ref_seqs, stats_df, flagstat_df)
    sections.quality(report, stats_df, sanitizer)
    sections.accuracy(report, stats_df_mapped, sanitizer)
    sections.read_coverage(report, stats_df_mapped, sanitizer)
    if depth_df is not None:
        sections.depths(report, depth_df)
    if counts is not None:
        sections.counts(report, flagstat_df, counts, sanitizer)
    # "de-sanitize" the report and write to the output file
    report_fname = f"{args.name}-report.html"
    with open(report_fname, "w") as report_file:
        report_file.write(sanitizer.desanitise_report(str(report)))

    logger.info(f"Written sanitized report to '{report_fname}'.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "--name",
        help="report name",
    )
    parser.add_argument(
        "--stats_dir",
        help="directory with `bamstats` per-read stats",
    )
    parser.add_argument(
        "--flagstat_dir",
        help="directory with `bamstats` per-file stats",
    )
    parser.add_argument(
        "--refnames_dir",
        help="directory with files containing reference names",
    )
    parser.add_argument(
        "--depths_dir",
        required=False,
        help="directory with depth files",
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
        # TODO: implement
        "--revision",
        default="unknown",
        help="git branch/tag of the executed workflow",
    )
    parser.add_argument(
        # TODO: implement
        "--commit",
        default="unknown",
        help="git commit of the executed workflow",
    )
    return parser
