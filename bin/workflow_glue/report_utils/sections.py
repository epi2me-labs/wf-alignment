"""Sections for ezcharts report."""
import numpy as np
import pandas as pd
from scipy import stats

import dominate.tags as dom_tags  # noqa: I100,I202
import dominate.util as dom_util

import ezcharts as ezc  # noqa: I202
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets import Grid, Progress, Tabs

THEME = "epi2melabs"


def histogram_with_mean_and_median(
    series,
    title=None,
    x_axis_name=None,
    y_axis_name=None,
    bins=100,
    round_digits=1,
):
    """Create ezcharts histogram showing the mean and median underneath the plot title.

    :param series: `pd.Series` with data to plot
    :param title: plot title, defaults to None
    :param x_axis_name: x axis label, defaults to None
    :param y_axis_name: y axis label, defaults to None
    :param bins: number of bins, defaults to 100
    :param round_digits: number of decimals to round the mean and median values to,
        defaults to 1
    :raises ValueError: Raise error if `series` is not a `pd.Series`
    :return: the histogram (`ezcharts.plots.Plot`)
    """
    if not isinstance(series, pd.Series):
        raise ValueError("`series` must be `pd.Series`.")

    plt = ezc.histplot(data=series, bins=bins)
    plt.title = dict(
        text=title,
        subtext=(
            f"Mean: {series.mean().round(round_digits)}. "
            f"Median: {series.median().round(round_digits)}"
        ),
    )
    if x_axis_name is not None:
        plt.xAxis.name = x_axis_name
    if y_axis_name is not None:
        plt.yAxis.name = y_axis_name
    return plt


def dropdown_with_histograms(
    data,
    groupby_column,
    data_column,
    plot_title_func,
    x_axis_name,
    y_axis_name,
    sanitizer,
):
    """Create multiple tabs with a histogram in each.

    Take a `pd.DataFrame`, perform a groupby on the `groupby_column` and then
    generate a tab with a histogram of the values in `data_column` for each group.

    :param data: `pd.DataFrame` containing at least `groupby_column` and `data_column`
    :param groupby_column: column to group by
    :param data_column: column of data to plot in the histogram
    :param plot_title_func: function taking the name of the group and returning the plot
        title
    :param x_axis_name: x axis label
    :param y_axis_name: y axis label
    :param sanitizer: `sanitizer.Sanitizer` for sanitizing strings passed to ezcharts
    """
    tabs = Tabs()
    with tabs.add_dropdown_menu():
        for grp_name, df in data.groupby(groupby_column, observed=True):
            grp_name = sanitizer(grp_name)
            plot_title = plot_title_func(grp_name)
            with tabs.add_dropdown_tab(grp_name):
                plt = histogram_with_mean_and_median(
                    df[data_column].dropna(),
                    title=plot_title,
                    x_axis_name=x_axis_name,
                    y_axis_name=y_axis_name,
                )
                EZChart(plt, theme=THEME)


def sub_heading(string):
    """Create sub heading in report.

    :param string: string to show in heading
    """
    with dom_tags.div():
        dom_tags.h5(string, {"class": "mb-0 pb-3 pt-3"})


def plural_s(string, n):
    """Add 's' to a string if `n > 1`.

    :param string: input string
    :param n: number to check if greater than 1
    :return: string with extra 's' to indicate plural in case `n` was greater than 1
    """
    return string + ("s" if n > 1 else "")


def get_summary_table(
    stats_df, flagstat_df, n_ref_seqs, n_reads_total, n_bases_total, secondary=False
):
    """Create table with summary statistics.

    :param stats_df: `pd.DataFrame` with bamstats per-read stats
    :param flagstat_df: `pd.DataFrame` with bamstats per-file stats
    :param n_ref_seqs: total number of reference sequences
    :param n_reads_total: total number of reads
    :param n_bases_total: total number of sequenced bases
    :param secondary: _description_, defaults to False
    """
    # get metrics
    n_detected_ref_seqs = (
        flagstat_df.query('ref != "*"').groupby("ref", observed=True)["total"].sum() > 0
    ).sum()
    n_reads = stats_df.shape[0]
    n_unmapped = flagstat_df["unmapped"].sum()
    n_bases = stats_df["read_length"].sum()
    # get percentages
    perc_detected_refs = n_detected_ref_seqs / n_ref_seqs * 100
    perc_reads = n_reads / n_reads_total * 100
    perc_unmapped = n_unmapped / n_reads * 100
    perc_bases = n_bases / n_bases_total * 100
    # get alignments per ref file (the table will have one row per ref file)
    alignments_per_ref_file = (
        flagstat_df.groupby("ref_file")["primary"].sum().drop("unmapped")
    ).sort_index()

    # function for progress bar cell showing percentage values
    def percentage_table_cell(value):
        dom_tags.td(
            Progress(
                value_min=0,
                value_max=100,
                value_now=round(value, 1),
                bar_cls="bg-secondary" if secondary else "",
                height=30,
            )
        )

    # create table
    with dom_tags.table(cls="table align-middle"):
        # table header
        with dom_tags.thead():
            dom_tags.th("Metric")
            dom_tags.th("Value")
            dom_tags.th("Percentage")
        # row for ref seqs
        with dom_tags.tr():
            dom_tags.td("Detected reference sequences")
            dom_tags.td(f"{n_detected_ref_seqs:,}")
            percentage_table_cell(perc_detected_refs)
        # row for reads
        with dom_tags.tr():
            dom_tags.td("Reads")
            dom_tags.td(f"{n_reads:,}")
            percentage_table_cell(perc_reads)
        # one row per ref with aligned reads
        for ref_file, primary_alignments in alignments_per_ref_file.items():
            perc_alignments = primary_alignments / n_reads * 100
            with dom_tags.tr():
                dom_tags.td(f"Reads aligned to '{ref_file}'")
                dom_tags.td(f"{primary_alignments:,}")
                percentage_table_cell(perc_alignments)
        # unmapped reads
        with dom_tags.tr():
            dom_tags.td("Unmapped reads")
            dom_tags.td(f"{n_unmapped:,}")
            percentage_table_cell(perc_unmapped)
        # bases
        with dom_tags.tr():
            dom_tags.td("Bases")
            dom_tags.td(f"{n_bases:,}")
            percentage_table_cell(perc_bases)


def summary(report, sample_names, ref_files, ref_seqs, stats_df, flagstat_df):
    """Create summary section.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param sample_names: collection of sample names
    :param ref_files: collection of reference file names
    :param ref_seqs: collection of reference sequence names
    :param stats_df: `pd.DataFrame` with bamstats per-read stats
    :param flagstat_df: `pd.DataFrame` with bamstats per-file stats
    """
    with report.add_section("Summary", "Summary"):
        # basic stats + intro
        dom_tags.p(
            """
            This report contains visualisations of statistics that can help in
            understanding the results from the wf-alignment workflow. Each section
            contains different plots or tables, and in general the results are broken
            down by sample or the reference file to which alignments were made. You can
            quickly jump to an individual section with the links in the header bar.
            """
        )
        ref_seqs_str = ", ".join(ref_seqs[:7]) + (", ..." if len(ref_seqs) > 7 else "")
        dom_util.raw(
            f"""
            <b>{len(sample_names)} {plural_s("sample", len(sample_names))}:</b><br>
            {', '.join(sample_names)}<br><br>
            <b>{len(ref_files)} {plural_s("reference file", len(ref_files))}:</b><br>
            {', '.join(ref_files)}<br><br>
            <b>{len(ref_seqs)} {plural_s("reference sequence", len(ref_seqs))}:</b><br>
            {ref_seqs_str}<br><br>
            """
        )
        n_reads_total = stats_df.shape[0]
        n_bases_total = stats_df["read_length"].sum()
        tabs = Tabs()
        # one tab for summary stats of all samples combined; one tab with dropdown per
        # sample (if we only got a single sample, the second tab won't be a dropdown
        # menu but rather a regular tab)
        with tabs.add_tab("total"):
            # summary table for all samples
            get_summary_table(
                stats_df,
                flagstat_df,
                len(ref_seqs),
                n_reads_total,
                n_bases_total,
                secondary=False,
            )
        with tabs.add_dropdown_menu(title="per sample"):
            for sample_name in sample_names:
                with tabs.add_dropdown_tab(sample_name):
                    # summary table for an individual sample
                    get_summary_table(
                        stats_df.query(f"sample_name == '{sample_name}'"),
                        flagstat_df.query(f"sample_name == '{sample_name}'"),
                        len(ref_seqs),
                        n_reads_total,
                        n_bases_total,
                        secondary=True,
                    )

        # fastcat / bamstats summary
        sub_heading("Reads Summary")
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for sample_name, sample_df in stats_df.groupby(
                "sample_name", observed=True
            ):
                with tabs.add_dropdown_tab(sample_name):
                    with Grid():
                        # read length plot
                        max_read_length_to_show = (
                            (sample_df["read_length"] / 1000)
                            .quantile(0.99, interpolation="lower")
                            .round()
                        )
                        plt = fastcat.read_length_plot(sample_df)
                        plt.xAxis.max = max_read_length_to_show
                        EZChart(plt, theme=THEME)
                        # base yield plot
                        plt = fastcat.base_yield_plot(sample_df)
                        # the base yield plot has kb as x-axis unit
                        plt.xAxis.max = max_read_length_to_show
                        plt.tooltip = None
                        EZChart(plt, theme=THEME)


def quality(report, stats_df, sanitizer):
    """Create read quality section.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param stats_df: `pd.DataFrame` with bamstats per-read stats
    :param sanitizer: `sanitizer.Sanitizer` for sanitizing strings passed to ezcharts
    """
    with report.add_section("Mean read quality", "Quality"):
        with Grid():
            # quality per sample
            with dom_tags.div():
                sub_heading("Per sample:")
                dropdown_with_histograms(
                    data=stats_df,
                    groupby_column="sample_name",
                    data_column="mean_quality",
                    plot_title_func=lambda _: None,
                    x_axis_name="Quality",
                    y_axis_name="Number of reads",
                    sanitizer=sanitizer,
                )
            # quality per ref file
            with dom_tags.div():
                sub_heading("Per ref.file:")
                dropdown_with_histograms(
                    data=stats_df,
                    groupby_column="ref_file",
                    data_column="mean_quality",
                    plot_title_func=lambda _: None,
                    x_axis_name="Quality",
                    y_axis_name="Number of reads",
                    sanitizer=sanitizer,
                )


def accuracy(report, stats_df_mapped, sanitizer):
    """Create mapping accuracy section.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param stats_df_mapped: `pd.DataFrame` with bamstats per-read stats (aligned reads
        only)
    :param sanitizer: `sanitizer.Sanitizer` for sanitizing strings passed to ezcharts
    """
    with report.add_section("Alignment accuracy", "Accuracy"):
        with Grid():
            # accuracy per sample
            with dom_tags.div():
                sub_heading("Per sample:")
                dropdown_with_histograms(
                    data=stats_df_mapped,
                    groupby_column="sample_name",
                    data_column="acc",
                    plot_title_func=lambda _: None,
                    x_axis_name="Accuracy [%]",
                    y_axis_name="Number of reads",
                    sanitizer=sanitizer,
                )
            # accuracy per ref file
            with dom_tags.div():
                sub_heading("Per ref.file:")
                dropdown_with_histograms(
                    data=stats_df_mapped,
                    groupby_column="ref_file",
                    data_column="acc",
                    plot_title_func=lambda _: None,
                    x_axis_name="Accuracy [%]",
                    y_axis_name="Number of reads",
                    sanitizer=sanitizer,
                )


def read_coverage(report, stats_df_mapped, sanitizer):
    """Create read coverage section.

    This section contains histograms showing the portion of individual reads covered by
    alignments.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param stats_df_mapped: `pd.DataFrame` with bamstats per-read stats (aligned reads
        only)
    :param sanitizer: `sanitizer.Sanitizer` for sanitizing strings passed to ezcharts
    """
    with report.add_section("Read coverage", "Read cov."):
        dom_tags.p(
            """
            These histograms show how much of an individual aligned read was covered by
            the alignment.
            """
        )
        with Grid():
            with dom_tags.div():
                # coverage per sample
                sub_heading("Per sample:")
                dropdown_with_histograms(
                    data=stats_df_mapped,
                    groupby_column="sample_name",
                    data_column="coverage",
                    plot_title_func=lambda _: None,
                    x_axis_name="Coverage [%]",
                    y_axis_name="Number of reads",
                    sanitizer=sanitizer,
                )
            with dom_tags.div():
                sub_heading("Per ref.file:")
                # coverage per ref (not per ref file)
                dropdown_with_histograms(
                    data=stats_df_mapped,
                    groupby_column="ref_file",
                    data_column="coverage",
                    plot_title_func=lambda _: None,
                    x_axis_name="Coverage [%]",
                    y_axis_name="Number of reads",
                    sanitizer=sanitizer,
                )


def ref_coverage(report, stats_df_mapped, sanitizer):
    """Create ref coverage section.

    This section contains histograms showing the portion of reference
    sequences covered by alignments of individual reads.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param stats_df_mapped: `pd.DataFrame` with bamstats per-read stats (aligned reads
        only)
    :param sanitizer: `sanitizer.Sanitizer` for sanitizing strings passed to ezcharts
    """
    with report.add_section("Reference coverage", "Ref. cov."):
        dom_tags.p(
            """
            These histograms show how much of the reference was covered by an individual
            alignment (i.e. by a single read mapped to it). These values will of course
            vary greatly for references with different lengths.
            """
        )
        with Grid():
            for ref_file, df in stats_df_mapped.groupby("ref_file", observed=True):
                with dom_tags.div():
                    sub_heading(f"Reference file '{ref_file}':")
                    dropdown_with_histograms(
                        data=df,
                        groupby_column="ref",
                        data_column="ref_coverage",
                        plot_title_func=lambda _: None,
                        x_axis_name="Coverage [%]",
                        y_axis_name="Number of reads",
                        sanitizer=sanitizer,
                    )


def get_relative_cumulative_depths(depth_df):
    """Calculate percentage of cumulative depth.

    Returns a `pd.DataFrame` in long format with the columns `['ref_percentage',
    'depth', 'sample_name']`. The values in `ref_percentage` indicate how much of the
    reference sequence were sequenced to at least this depth.

    :param depth_df: `pd.DataFrame` with depth data
    :return: `pd.DataFrame` with columns `['ref_percentage', 'depth', 'sample_name']`
    """
    cumul_depth_dfs = []
    for sample_name, df in depth_df.groupby("sample_name", observed=True):
        df = df.eval("step = end - start")
        df = df.sort_values("depth", ascending=False)
        df = df.eval("ref_percentage = step.cumsum()")[["ref_percentage", "depth"]]
        df["ref_percentage"] /= df["ref_percentage"].iloc[-1]
        df["ref_percentage"] *= 100
        # set ref_percentage at min/mox depth to 0 and 100 to make the plots look nicer
        df["ref_percentage"].iloc[0] = 0
        df.iloc[-1] = [100, 0]
        df["sample_name"] = sample_name
        cumul_depth_dfs.append(df)
    return pd.concat(cumul_depth_dfs)


def depths(report, depth_df):
    """Create depth section.

    This section contains a plot with depth coverage vs. genome position on the left and
    relative cumulative coverage (i.e. percentage of genome covered to at least a
    certain depth vs. depth) on the right.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param depth_df: `pd.DataFrame` with depth data.
    """
    with report.add_section("Depth of coverage", "Depth"):
        dom_tags.p(
            """
            This section illustrates the depth of coverage of the reference genomes. The
            left plot shows coverage vs. genomic position (note that the coordinates on
            the x-axis are the positions along the concatenated reference including all
            reference sequences in the respective reference file). The right plot shows
            the cumulative fraction of the reference that was covered to at least a
            certain depth.
            """
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for ref_file, df_ref_file in depth_df.groupby("ref_file"):
                # prepare data for depth vs coordinate plot
                df_depth_vs_coords = (
                    df_ref_file.eval("mean_pos = (start + end) / 2")
                    .eval("step = end - start")
                    .reset_index()
                )
                ref_lengths = df_depth_vs_coords.groupby("ref", observed=True)[
                    "end"
                ].last()
                total_ref_starts = ref_lengths.cumsum().shift(1, fill_value=0)
                df_depth_vs_coords["total_mean_pos"] = df_depth_vs_coords.groupby(
                    "ref", observed=True, group_keys=False
                )["mean_pos"].apply(lambda s: s + total_ref_starts[s.name])
                # prepare data for cumulative depth plot
                df_cumul_depth = get_relative_cumulative_depths(df_ref_file)
                with tabs.add_dropdown_tab(ref_file):
                    # depth vs genomic coordinate plot on the left and cumulative depth
                    # plot on the right
                    with Grid():
                        # TODO: maybe smooth this line a bit
                        plt = ezc.lineplot(
                            data=df_depth_vs_coords.round(2),
                            x="total_mean_pos",
                            y="depth",
                            hue="sample_name",
                        )
                        plt.title = {"text": "Coverage along reference"}
                        plt.xAxis.name = "Position along reference"
                        plt.yAxis.name = "Sequencing depth"
                        for s in plt.series:
                            s.showSymbol = False
                        plt.tooltip = {"trigger": "item"}
                        EZChart(plt, theme=THEME)
                        # now the cumulative depth plot
                        plt = ezc.lineplot(
                            data=df_cumul_depth.round(2),
                            x="depth",
                            y="ref_percentage",
                            hue="sample_name",
                        )
                        plt.xAxis.scale = True
                        plt.title = {"text": "Cumulative coverage"}
                        plt.xAxis.name = "Sequencing depth"
                        plt.yAxis.name = "Percentage of reference"
                        for s in plt.series:
                            s.showSymbol = False
                        plt.tooltip = {"trigger": "item"}
                        EZChart(plt, theme=THEME)


def counts(report, flagstat_df, counts, sanitizer):
    """Create counts section.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param flagstat_df: `pd.DataFrame` with bamstats per-file stats
    :param counts: `pd.Series` of expected counts (with ref. sequence IDs as index)
    :param sanitizer: `sanitizer.Sanitizer` for sanitizing strings passed to ezcharts
    """
    exp_obs_counts_df = (
        flagstat_df.query('ref != "*"')[["ref", "sample_name", "ref_file", "primary"]]
        .rename(columns={"primary": "obs"})
        .set_index("ref")
    )
    exp_obs_counts_df["exp"] = counts
    exp_obs_counts_df.dropna(inplace=True)

    with report.add_section("Read count control", "Control"):
        dom_tags.p(
            r"""
            When a file with expected read counts was provided, this plot shows the
            observed vs. expected counts for each sample / reference sequence
            combination.
            """
        )
        no_plot_generated = True
        # have one group of tabs per ref file
        for ref_file, ref_file_df in exp_obs_counts_df.groupby(
            "ref_file", observed=True
        ):
            sub_heading(f"Reference file '{ref_file}':")
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for sample_name, sample_df in ref_file_df.groupby(
                    "sample_name", observed=True
                ):
                    n_detected_refs = (sample_df["obs"] > 0).sum()
                    log_counts_df = np.log10(sample_df[["exp", "obs"]] + 1)
                    spear_r, spear_p = stats.spearmanr(
                        log_counts_df["exp"], log_counts_df["obs"]
                    )
                    pears_r, pears_p = stats.pearsonr(
                        log_counts_df["exp"], log_counts_df["obs"]
                    )
                    with tabs.add_dropdown_tab(sample_name):
                        plt = ezc.scatterplot(
                            data=log_counts_df.round(2), x="exp", y="obs"
                        )
                        plt.title = {
                            "text": "Expected vs. observed counts",
                            "subtext": sanitizer(
                                f"Spearman's: {spear_r:.2f} (p={spear_p:.1e}); "
                                f"Pearson's: {pears_r:.2f} (p={pears_p:.1e}); "
                                "detected refs: "
                                f"{n_detected_refs} / {sample_df.shape[0]}"
                            ),
                        }
                        plt.xAxis.name = "log10 of expected counts"
                        plt.yAxis.name = "log10 of observed counts"
                        plt.xAxis.scale = True
                        plt.yAxis.scale = True

                        # add ref IDs as tooltips
                        plt.dataset[0].source = (
                            log_counts_df.reset_index()
                            .eval("hue = None")[["exp", "obs", "hue", "ref"]]
                            .values
                        )
                        plt.dataset[0].dimensions = ["x", "y", "hue", "tooltip"]
                        plt.series[0].encode["tooltip"] = "tooltip"
                        plt.tooltip = {"trigger": "item"}

                        EZChart(plt, theme=THEME)
                        no_plot_generated = False
        if no_plot_generated:
            dom_tags.p(
                """
                No observations were found for this sample / reference combination.
                """
            )
