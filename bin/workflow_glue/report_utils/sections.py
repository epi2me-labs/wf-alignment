"""Sections for ezcharts report."""

from bokeh.models import HoverTool
from bokeh.models import Title
import numpy as np
import pandas as pd
from scipy import stats
import dominate.tags as dom_tags  # noqa: I100,I202
import dominate.util as dom_util

import ezcharts as ezc  # noqa: I202
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets import Grid, Progress, Tabs

from . import read_data  # noqa: ABS101

THEME = "epi2melabs"


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
    n_reads,
    n_bases,
    flagstat_df,
    n_ref_seqs,
    n_reads_total,
    n_bases_total,
    secondary=False,
):
    """Create table with summary statistics.

    :param n_reads: number of reads
    :param n_bases: number of sequenced bases
    :param flagstat_df: `pd.DataFrame` with bamstats per-file stats
    :param n_ref_seqs: total number of reference sequences
    :param n_reads_total: total number of reads
    :param n_bases_total: total number of sequenced bases
    :param secondary: set `bar_cls` of progress bars to 'bg-secondary'
    """
    # get metrics
    n_detected_ref_seqs = (
        flagstat_df.query('ref != "*"').groupby("ref", observed=True)["total"].sum() > 0
    ).sum()
    n_unmapped = flagstat_df["unmapped"].sum()
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


def summary(report, per_sample_dirs, ref_files, ref_seqs, flagstat_df):
    """Create introduction section.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param per_sample_dirs: collection of directories with per-sample result files (i.e.
        one directory per sample)
    :param ref_files: collection of reference file names
    :param ref_seqs: collection of reference sequence names
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
        sample_names = [d.name for d in per_sample_dirs]
        # we don't want to list all sample names / ref files / ref seq IDs as there
        # might be hundreds; we thus only show the first 7 and `, ...`
        sample_names_str = "&emsp;".join(sample_names[:7]) + (
            "&emsp;..." if len(sample_names) > 7 else ""
        )
        ref_files_str = "&emsp;".join(ref_files[:7]) + (
            "&emsp;..." if len(ref_files) > 7 else ""
        )
        ref_seqs_str = "&emsp;".join(ref_seqs[:7]) + (
            "&emsp;..." if len(ref_seqs) > 7 else ""
        )
        dom_util.raw(
            f"""
            <b>{len(sample_names)} {plural_s("sample", len(sample_names))}:</b><br>
            {sample_names_str}<br><br>
            <b>{len(ref_files)} {plural_s("reference file", len(ref_files))}:</b><br>
            {ref_files_str}<br><br>
            <b>{len(ref_seqs)} {plural_s("reference sequence", len(ref_seqs))}:</b><br>
            {ref_seqs_str}<br><br>
            """
        )

        # read the length histograms to get the number of reads and bases for each
        # sample (and the totals as well)
        length_hist_dict = {
            sample_dir.name: read_data.length_hist(sample_dir)
            for sample_dir in per_sample_dirs
        }
        n_reads_dict = {
            sample_name: length_hist_dict[sample_name]["count"].sum()
            for sample_name in sample_names
        }
        n_bases_dict = {
            sample_name: (length_hist_dict[sample_name].eval("start * count")).sum()
            for sample_name in sample_names
        }

        n_reads_total = sum(hist["count"].sum() for hist in length_hist_dict.values())
        n_bases_total = sum(
            (hist["start"] * hist["count"]).sum() for hist in length_hist_dict.values()
        )
        tabs = Tabs()
        # one tab for summary stats of all samples combined; one tab with dropdown per
        # sample (if we only got a single sample, the second tab won't be a dropdown
        # menu but rather a regular tab)
        with tabs.add_tab("total"):
            # summary table for all samples
            get_summary_table(
                n_reads_total,
                n_bases_total,
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
                        n_reads_dict[sample_name],
                        n_bases_dict[sample_name],
                        flagstat_df.query(f"sample_name == '{sample_name}'"),
                        len(ref_seqs),
                        n_reads_total,
                        n_bases_total,
                        secondary=True,
                    )


def seqsum(report, per_sample_dirs):
    """Get ezcharts fastcat.SeqSummary.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param per_sample_dirs: collection of directories with per-sample result files (i.e.
        one directory per sample)
    """
    with report.add_section("Read and alignment summary statistics", "Read stats"):
        dom_tags.p(
            """
            Read quality and length in addition to alignment accuracy and coverage are
            illustrated in the plots below.
            """
        )
        fastcat.SeqSummary(
            tuple(per_sample_dirs), sample_names=tuple(d.name for d in per_sample_dirs)
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
                df_depth_vs_coords = df_ref_file.eval("step = end - start")
                # for each sample, get the mean position of each depth window (while
                # concatenating all the sequences in that ref file)
                df_depth_vs_coords = df_depth_vs_coords.groupby("sample_name").apply(
                    lambda df: df.eval("total_mean_pos=step.cumsum() - step / 2")
                )
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
                            title="Coverage along reference",
                            marker=False,
                        )
                        plt._fig.xaxis.axis_label = "Position along reference"
                        plt._fig.yaxis.axis_label = "Sequencing depth"
                        EZChart(plt, theme=THEME)
                        # now the cumulative depth plot
                        plt = ezc.lineplot(
                            data=df_cumul_depth.round(2),
                            x="depth",
                            y="ref_percentage",
                            hue="sample_name",
                            title="Cumulative coverage",
                            marker=False,
                        )
                        # plt.xAxis.scale = True
                        plt._fig.xaxis.axis_label = "Sequencing depth"
                        plt._fig.yaxis.axis_label = "Percentage of reference"
                        EZChart(plt, theme=THEME)


def counts(report, flagstat_df, counts):
    """Create counts section.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param flagstat_df: `pd.DataFrame` with bamstats per-file stats
    :param counts: `pd.Series` of expected counts (with ref. sequence IDs as index)
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
                    log_counts_df = np.log10(sample_df[["exp", "obs"]] + 1).round(2)
                    spear_r, spear_p = stats.spearmanr(
                        log_counts_df["exp"], log_counts_df["obs"]
                    )
                    pears_r, pears_p = stats.pearsonr(
                        log_counts_df["exp"], log_counts_df["obs"]
                    )
                    log_counts_df = log_counts_df.reset_index()
                    with tabs.add_dropdown_tab(sample_name):
                        plt = ezc.scatterplot(
                            data=log_counts_df,
                            x="exp",
                            y="obs",
                            hue="ref",
                            title="Expected vs. observed counts",
                        )
                        plt._fig.add_layout(Title(text=(
                                f"Spearman's: {spear_r:.2f} (p={spear_p:.1e}); "
                                f"Pearson's: {pears_r:.2f} (p={pears_p:.1e}); "
                                "detected refs: "
                                f"{n_detected_refs} / {sample_df.shape[0]}"
                            ), text_font_style="italic"), 'above')
                        plt._fig.xaxis.axis_label = "log10 of expected counts"
                        plt._fig.yaxis.axis_label = "log10 of observed counts"

                        # Hide legend because could be many entries.
                        plt._fig.legend.visible = False
                        # IDs as tooltips
                        # Add ref to the plot to be able to access later in hover
                        # There is one glyd per ref, which is the hue
                        # Replace with the corresponding info based on the index value
                        for gly in plt._fig.renderers:
                            # add more properties for hover later
                            gly.data_source.data['ref'] = log_counts_df.loc[
                                list(gly.data_source.to_df().index)]['ref']

                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [
                            ("exp", "$x"),
                            ("obs", "$y"),
                            ("ref", "@ref"),
                        ]
                        EZChart(plt, theme=THEME)
                        no_plot_generated = False
        if no_plot_generated:
            dom_tags.p(
                """
                No observations were found for this sample / reference combination.
                """
            )
