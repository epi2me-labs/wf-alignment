#!/usr/bin/env python
"""Report component for displaying information from mapula's json output."""

import argparse
import json
import math
import os
import sys

from aplanat import hist, points
from aplanat import report
from aplanat.components import depthcoverage, fastcat
from aplanat.components import simple as scomponents
from aplanat.report import HTMLSection
from aplanat.util import Colors
from bokeh.layouts import gridplot, layout
from bokeh.models import ColumnDataSource, CustomJS, Select
from bokeh.models import Div, HoverTool, Panel, Slope, Tabs, Title
from bokeh.palettes import Category20c
from bokeh.plotting import figure
from bokeh.transform import cumsum
import markdown
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


class PlotMappingStats(HTMLSection):
    """Build an aplanat dashboard from mapula's json output."""

    def __init__(
        self,
        json,
        counts,
        references,
        unmapped,
        sample_name,
        depth_file,
    ):
        """Load the json file and output the dashboard."""
        super().__init__()
        self.json = json
        self.md = markdown.Markdown()
        self.counts = counts
        self.references = references
        self.unmapped = unmapped
        self.sample_name = sample_name
        self.depth_file = depth_file
        self.data = self.load_data(self.json)
        self.report = self.build_report(
            self.counts,
            self.references,
            self.unmapped,
            self.sample_name,
            self.depth_file,
            **self.data
        )

    @staticmethod
    def load_data(path):
        """Load_data."""
        try:
            with open(path) as data:
                try:
                    return json.load(data)
                except json.decoder.JSONDecodeError:
                    sys.stdout.write('Error loading data from {}'.format(
                        path))
                    raise
        except IOError:
            sys.stdout.write("Path {} cannot be opened".format(path))
            raise

    def build_report(
        self,
        counts,
        references,
        unmapped,
        sample_name,
        depth_file,
        **data
    ):
        """Build_report."""
        tabs = []

        # Plot summary
        summary_tab = self.build_summary_tab(data)
        tabs.append(summary_tab)

        # Plot accuracy
        acc_tab = self.build_accuracy_tab(data)
        tabs.append(acc_tab)

        # Plot quality
        qual_tab = self.build_quality_tab(data, unmapped)
        tabs.append(qual_tab)

        # Plot lengths
        len_tab = self.build_length_tab(data)
        tabs.append(len_tab)

        # Plot coverage
        cov_tab = self.build_coverage_tab(data)
        tabs.append(cov_tab)

        # Plot coverage
        depth_tab = self.depth_graph(references, depth_file)
        tabs.append(depth_tab)

        # Plot whole genome coverage
        cov_tab = self.build_cumulative_coverage_tab(references, depth_file)
        tabs.append(cov_tab)

        # Plot control
        if counts:
            con_tab = self.build_control_tab(data, counts)
            tabs.append(con_tab)

        # Tabula rasa
        tabs = list(filter(None, tabs))
        panel = Tabs(tabs=tabs)
        self.plot(panel)
        return panel

    def depth_graph(self, references, depth_file):
        """Create depth vs position graph."""
        # create list of graphs
        graphs = []
        # create dictionary of {ref: depths}
        if str(depth_file) == 'None':
            return None
        depth_file.columns = ['ref', 'start', 'end', 'depth']
        all_ref = dict(tuple(depth_file.groupby(['ref'])))
        ref_keys = list(all_ref.keys())
        for name, isoforms in references.items():
            refs_present = list(set(isoforms).intersection(ref_keys))
            graph_dic = {}
            for entry in refs_present:
                v = all_ref[entry]
                graph_dic[entry] = {
                    'x': v['start'],
                    'y': v['depth']}
            if len(graph_dic) > 0:
                drop_down = list(graph_dic.keys())
                first_plot = graph_dic[drop_down[0]]
                plot = figure()
                source = ColumnDataSource(data=first_plot)
                select = Select(
                    title="Select reference for " + name + ':',
                    value=drop_down[0], options=drop_down)
                plot = figure(
                    plot_height=350, y_range=None, x_range=None,
                    x_axis_label='Position along reference',
                    y_axis_label='Sequencing depth / Bases',
                    title=str(name))
                plot.step(x='x', y='y', source=source)
                plot.xaxis.formatter.use_scientific = False
                select.js_on_change('value', CustomJS(
                    args=dict(
                            source=source,
                            source_graph=graph_dic,
                            select=select), code="""
                            const new_data = Object.assign({}, source.data)
                            source.data = source_graph[select.value]
                            """))
                drop_down_plot = layout([select, plot])
                graphs.append(drop_down_plot)
        alignment_grid = gridplot(
            graphs, ncols=2,
            sizing_mode="stretch_width")
        text = (
            "This tab contains visualisations of "
            "the depth of coverage of alignments across the reference."
            " If there are no Depth graphs for"
            " a particular reference"
            " it is likely to mean no alignments were found."
        )
        plots = [
            [self.get_description(text)],
            [alignment_grid]]
        main = layout(plots, sizing_mode="scale_width")
        return Panel(child=main, title="Depth")

    def build_summary_tab(self, data):
        """Build_summary_tab."""
        alignment_plots = [
            self.plot_alignment_summary(value)
            for value in data.values()
            if value['primary_count']
        ]

        if not alignment_plots:
            return None

        alignment_grid = gridplot(
            alignment_plots, ncols=2,
            sizing_mode="stretch_width")

        text1 = (
            "This report contains visualisations of "
            "statistics that can help in understanding "
            "the results from the wf-alignment workflow. "
            "Each tab above contains a different type of "
            "plot, and in general the results are broken "
            "down by the reference genome to which "
            "alignments were made. On this tab, the donut "
            "plot displays the proportion of mapped to "
            "unmapped reads (hover over for counts) and "
            "the bar charts give a breakdown of the different "
            "alignment types per-reference."
        )
        plots = [
            [self.get_description(text1)],
            [self.plot_base_pairs(data)],
            [alignment_grid]
        ]

        main = layout(plots, sizing_mode="scale_width")
        return Panel(child=main, title="Summary")

    def build_accuracy_tab(self, data):
        """Build_accuracy_tab."""
        accuracy_plots = [
            self.plot_accuracy_distribution(value)
            for value in data.values()
            if sum(value['alignment_accuracies'])
        ]

        if not accuracy_plots:
            return None

        accuracy_grid = gridplot(
            accuracy_plots, ncols=2,
            sizing_mode="stretch_width")

        text = (
            "This tab contains visualisations of "
            "percentage alignment accuracy broken "
            "down by aligned reference. Naturally, "
            "no accuracy measurements are available "
            "for unmapped reads."
        )
        plots = [
            [self.get_description(text)],
            [accuracy_grid]
        ]

        main = layout(plots, sizing_mode="scale_width")
        return Panel(child=main, title="Accuracy")

    def build_coverage_tab(self, data):
        """Build_coverage_tab."""
        coverage_plots = [
            self.plot_coverage_distribution(value)
            for value in data.values()
            if sum(value['alignment_coverages'])
        ]

        if not coverage_plots:
            return None

        coverage_grid = gridplot(
            coverage_plots, ncols=2,
            sizing_mode="stretch_width")

        text = (
            "This tab contains visualisations of "
            "reference-coverage (i.e. the proportion of "
            "the reference that a given alignment covers "
            "ex. soft-clipping). The cov80 gives the number "
            "of alignments which cover more than 80 percent "
            "of the reference."
        )
        plots = [
            [self.get_description(text)],
            [coverage_grid]
        ]

        main = layout(plots, sizing_mode="scale_width")
        return Panel(child=main, title="Coverage")

    def build_cumulative_coverage_tab(self, references, depth_df):
        """Build cumulative coverage tab."""
        if str(depth_df) == 'None':
            return None
        depth_df = depth_df.copy()
        depth_df.columns = ['ref', 'start', 'end', 'depth']
        # Map Genome names onto the depth frame
        genome_name_map = {}
        for genome, contigs in references.items():
            genome_name_map.update({contig: genome for contig in contigs})
        depth_df['genome'] = depth_df.ref.map(genome_name_map)

        plots = []
        for genome, df in depth_df.groupby('genome'):
            plots.append(
                depthcoverage.cumulative_depth_from_bed(df, title=genome))

        coverage_grid = gridplot(
            plots, ncols=2,
            sizing_mode="stretch_width")

        text = (
            "Plot(s) showing cumulative read depth across reference genome."
        )
        plots = [
            [self.get_description(text)],
            [coverage_grid]
        ]

        main = layout(plots, sizing_mode="scale_width")
        return Panel(child=main, title="Cumulative genome coverage")

    def build_length_tab(self, data):
        """Build_length_tab."""
        length_plots = [
            self.plot_read_length_distribution(value)
            for value in data.values()
            if sum(value['read_lengths'])
        ]

        if not length_plots:
            return None

        length_grid = gridplot(
            length_plots, ncols=2,
            sizing_mode="stretch_width")

        text = (
            "This tab contains visualisations of "
            "read length. The n50 is defined as the "
            "length N for which 50 percent of all "
            "bases in the sequence are in a sequence "
            "of length L < N."
        )
        plots = [
            [self.get_description(text)],
            [length_grid]
        ]

        main = layout(plots, sizing_mode="scale_width")
        return Panel(child=main, title="Length")

    def build_quality_tab(self, data, unmapped):
        """Build_quality_tab."""
        quality_plots = [
            self.plot_qscore_distribution(value)
            for value in data.values()
            if sum(value['aligned_qualities'])
        ]
        fastcat_unmapped = fastcat.read_quality_plot(unmapped)
        fastcat_unmapped.title.text = 'Unmapped reads'
        quality_plots.append(fastcat_unmapped)

        if not quality_plots:
            return None

        quality_grid = gridplot(
            quality_plots, ncols=2,
            sizing_mode="stretch_width")

        text = (
            "This tab contains visualisations of "
            "avg read quality. Each data point "
            "represents the mean per-base Phred quality "
            "score for a given read. Quality scores are "
            "related to the underlying error probabilities."
        )
        plots = [
            [self.get_description(text)],
            [quality_grid]
        ]

        main = layout(plots, sizing_mode="scale_width")
        return Panel(child=main, title="Quality")

    def build_control_tab(self, data, counts):
        """Build_control_tab."""
        counts_df = pd.read_csv(counts)
        counts_df.columns = map(str.lower, counts_df.columns)

        if not all(
            col in counts_df.columns
            for col in ['reference', 'expected_count']
        ):
            sys.stdout.write(
                "[Error]: Supplied counts file does not "
                "contain the required columns, "
                "'reference,expected_count'"
            )
            sys.exit(1)

        counts_df = counts_df.set_index('reference')
        required = ['spearmans_rho', 'observed_references']
        text = (
            "This tab becomes available if you "
            "have supplied an expected counts csv to the "
            "wf-alignment workflow. Providing expected counts "
            "permits us to calculate the correlation between "
            "the number of observations for each of a given set "
            "of reference sequences against their expected counts. "
        )
        plots = [
            [self.get_description(text)],
        ]

        for key, value in data.items():
            if not all(value.get(req) for req in required):
                continue

            observed_expected_corr = (
                self.plot_observed_vs_expected_correlations(
                    value, counts_df))
            detected_vs_undetected = (
                self.plot_detected_vs_undetected_references(
                    value, counts_df))
            observations_vs_expected_molarity = (
                self.plot_observations_vs_expected_molarity(
                    value, counts_df))

            plots.append([Div(text=f'<h3>{key}</h3>')])
            plots.append([observed_expected_corr, detected_vs_undetected])
            plots.append([observations_vs_expected_molarity])

        main = layout(plots, sizing_mode="scale_width")
        return Panel(child=main, title="Control")

    def plot_base_pairs(self, data):
        """Plot_base_pairs."""
        base_pairs = {
            self.abbreviate_name(k, data): v['base_pairs']
            for k, v in data.items()
        }

        data = pd.Series(base_pairs).reset_index(
            name='value').rename(columns={'index': 'group'})
        data['angle'] = data['value']/data['value'].sum() * 2 * math.pi
        if Category20c.get(len(base_pairs)):
            data['color'] = Category20c[len(base_pairs)]
        else:
            data['color'] = Colors.cerulean
        data['percentage'] = (data['value']/data['value'].sum())
        data['megabases'] = (data['value'] / 1000000)

        plot = figure(
            plot_height=200,
            title="Base pairs breakdown",
            tools="hover",
            tooltips="@group: @value",
            x_range=(-0.5, 1.0),
            toolbar_location=None
        )

        plot.annular_wedge(
            x=0,
            y=1,
            inner_radius=0.1,
            outer_radius=0.2,
            start_angle=cumsum('angle', include_zero=True),
            end_angle=cumsum('angle'),
            line_color="white",
            fill_color='color',
            legend_field='group',
            source=data
        )

        plot.axis.axis_label = None
        plot.axis.visible = False
        plot.grid.grid_line_color = None
        plot.min_border_left = 0
        plot.legend.visible = True

        hover = plot.select({'type': HoverTool})
        hover.tooltips = [
            ('group', '@group'),
            ('percentage', '@percentage{0.00%}'),
            ('count (mb)', '@megabases{0.00}')
        ]

        total = "Total: {:.2f}mb".format(data['megabases'].sum())
        self.add_plot_title(plot, {}, total)
        self.style_plot(plot)
        return plot

    def plot_alignment_summary(self, data):
        """Plot_alignment_summary."""
        labels = ['Primary', 'Secondary', 'Supplementary']
        counts = [
            data['primary_count'],
            data['secondary_count'],
            data['supplementary_count']
        ]

        if not sum(counts):
            return

        plot = figure(
            x_range=labels,
            plot_height=350,
            output_backend="webgl",
            toolbar_location=None,
            x_axis_label='Outcome',
            y_axis_label='Number of alignments',
            tools=""
        )
        plot.vbar(
            x=labels, top=counts,
            width=0.9, fill_color=Colors.cerulean)
        plot.xgrid.grid_line_color = None
        plot.y_range.start = 0

        total = f"Total observations: {data['observations']}"
        self.style_plot(plot)
        self.add_plot_title(plot, data, total)
        return plot

    def plot_accuracy_distribution(self, data):
        """Plot_accuracy_distribution."""
        counts = []
        for index, item in enumerate(data['alignment_accuracies']):
            for i in range(item):
                counts.append(index/10)

        plot = hist.histogram(
            [counts],
            bins=100,
            height=300,
            width=400,
            xlim=(0, 100),
            colors=[Colors.cerulean],
            x_axis_label='Accuracy %',
            y_axis_label='Number of reads'
        )

        median = f"Median: {data['median_accuracy']}%"
        self.style_plot(plot)
        self.add_plot_title(plot, data, median)
        return plot

    def plot_coverage_distribution(self, data):
        """Plot_coverage_distribution."""
        counts = []
        for index, item in enumerate(data['alignment_coverages']):
            for i in range(item):
                counts.append(index)

        plot = hist.histogram(
            [counts],
            bins=101,
            height=300,
            width=400,
            xlim=(0, 100),
            colors=[Colors.cerulean],
            x_axis_label='Coverage %',
            y_axis_label='Number of reads',
        )

        cov80_perc = round(data['cov80_percent'], 2)
        cov80 = (
            f"Cov80: {data['cov80_count']} ({cov80_perc}%)"
        )
        self.style_plot(plot)
        self.add_plot_title(plot, data, cov80)
        return plot

    def plot_qscore_distribution(self, data):
        """Plot_qscore_distribution."""
        counts = []
        for index, item in enumerate(data['aligned_qualities']):
            for i in range(item):
                counts.append(index/10)

        plot = hist.histogram(
            [counts],
            bins=600,
            height=300,
            width=400,
            xlim=(0, 30),
            colors=[Colors.cerulean],
            x_axis_label='Avg Quality',
            y_axis_label='Number of reads',
        )

        median = f"Median: {int(data['median_quality'])}"
        self.style_plot(plot)
        self.add_plot_title(plot, data, median)
        return plot

    def plot_read_length_distribution(self, data):
        """Plot_read_length_distribution."""
        counts = []
        max_length = 0
        for index, item in enumerate(data['read_lengths']):
            for i in range(item):
                counts.append(index * 50)
                length = index * 50
                if length > max_length:
                    max_length = length

        plot = hist.histogram(
            [counts],
            bins=1000,
            height=300,
            width=400,
            xlim=(0, max_length + 200),
            colors=[Colors.cerulean],
            x_axis_label='Length',
            y_axis_label='Number of reads',
        )

        n50 = f"Read n50: {int(data['read_n50'])}"
        self.style_plot(plot)
        self.add_plot_title(plot, data, n50)
        return plot

    def plot_observed_vs_expected_correlations(self, value, counts):
        """Plot_observed_vs_expected_correlations."""
        names = []
        obs = []
        exps = []

        countdict = counts.to_dict()['expected_count']
        for name, exp in countdict.items():
            names.append(name)
            exps.append(exp)
            obs.append(value['observed_references'].get(name, 0))

        data = pd.DataFrame(
            zip(names, obs, exps),
            columns=['Name', 'Observed', 'Expected']
        )

        data['log_obs'] = [math.log(y+1, 10) for y in data['Observed']]
        data['log_exp'] = [math.log(y+1, 10) for y in data['Expected']]

        model = LinearRegression().fit(
            np.array(data['log_exp'].values).reshape(-1, 1),
            data['log_obs']
        )

        regression_line = Slope(
            gradient=model.coef_[0],
            y_intercept=model.intercept_,
            line_color=Colors.cerulean
        )

        plot = points.points(
            [data['log_exp'].tolist()],
            [data['log_obs'].tolist()],
            height=350,
            tools="",
            toolbar_location=None,
            output_backend="webgl",
            colors=[Colors.cerulean],
            x_axis_label='log10(Expected Count)',
            y_axis_label='log10(Observed Count)',
            title='Expected vs Observed'
        )

        corrs = (
            "Spearmans: {:.2f}, (p = {:.2f})  "
            "Pearsons: {:.2f}, (p = {:.2f})".format(
                value['spearmans_rho'],
                value['spearmans_rho_pval'],
                value['pearson'],
                value['pearson_pval']
            ))
        plot.add_layout(regression_line)
        self.add_plot_title(plot, {}, corrs)
        self.style_plot(plot)
        return plot

    def plot_observations_vs_expected_molarity(self, value, counts):
        """Plot_observations_vs_expected_molarity."""
        names = []
        obs = []
        exps = []

        countdict = counts.to_dict()['expected_count']
        for name, exp in countdict.items():
            names.append(name)
            exps.append(exp)
            obs.append(value['observed_references'].get(name, 0))

        data = pd.DataFrame(
            zip(names, obs, exps),
            columns=['Name', 'Observed', 'Expected']
        )

        data['log_obs'] = [math.log(y+1, 10) for y in data['Observed']]
        data['log_exp'] = [math.log(y+1, 10) for y in data['Expected']]

        data.sort_values('Expected', inplace=True)

        plot = figure(
            x_range=data['Name'].tolist(),
            plot_height=200,
            toolbar_location=None,
            tools="",
            x_axis_label='Reference',
            y_axis_label='log10(Observed Count)',
            output_backend="webgl",
            title='Observations ordered by increasing expected count'
        )
        plot.vbar(
            x=data['Name'].tolist(),
            top=data['log_obs'].values, width=0.9,
            fill_color=Colors.cerulean
        )
        plot.xgrid.grid_line_color = None
        plot.y_range.start = 0

        plot.xaxis.major_label_orientation = 3.14/2
        self.style_plot(plot)
        return plot

    def plot_detected_vs_undetected_references(self, value, counts_df):
        """plot_detected_vs_undetected_references."""
        names = ['Detected', 'Not Detected']
        detected = value['observed_reference_count']
        undetected = len(counts_df)
        data = [detected, undetected]

        plot = figure(
            x_range=names,
            toolbar_location=None,
            tools="",
            plot_height=350,
            x_axis_label='Outcome',
            y_axis_label='Number of alignments',
            output_backend="webgl",
            title='Detected vs Undetected references'
        )
        plot.vbar(
            x=names, top=data, width=0.9,
            fill_color=Colors.cerulean
        )
        plot.xgrid.grid_line_color = None
        plot.y_range.start = 0

        detected = f"Detected: {detected} / {undetected}"
        self.add_plot_title(plot, {}, detected)
        self.style_plot(plot)
        return plot

    def get_description(self, desc):
        """Get Description."""
        styles = [
            "display:block;",
            "width:100%;",
            "padding:25px 0 0 0;",
            "font-size: 16px;",
            "margin-bottom: 0;"
        ]
        return Div(
            text=f'<p style="{"".join(styles)}" class="lead">{desc}</p>'
        )

    def style_plot(self, plot):
        """Style plots."""
        plot.margin = (10, 10, 40, 10)
        plot.background_fill_alpha = 0

    def abbreviate_name(self, name, data):
        """Abbreviate the name."""
        item = data[name]
        fasta = item['fasta']
        barcode = item['barcode']
        run_id = item['run_id']

        if len(fasta) > 20:
            fasta = fasta[0:20] + '...'

        if len(run_id) > 20:
            run_id = run_id[0:20] + '...'

        return f"{fasta} / {barcode} / {run_id}"

    def add_plot_title(self, plot, data, *extra_lines):
        """Annotate plot with titles."""
        for line in extra_lines:
            plot.add_layout(
                Title(
                    text=line,
                    text_line_height=0.5,
                    text_font_style='italic',
                    text_color='grey'
                ),
                'above',
            )
        for identkey in ['source', 'fasta', 'reference', 'read_group', 'barcode']:  # noqa
            if not data.get(identkey):
                continue
            plot.add_layout(
                Title(
                    text=str(identkey)+':'+data[identkey],
                    text_line_height=0.5,
                ),
                'above'
            )


def get_references(references):
    """Grab multiple  name references out of reference files."""
    reference_sets = {}
    for fname in references:
        name = os.path.splitext(fname)[0]
        with open(fname) as f:
            ref_iso = [line.rstrip('\n') for line in f]
        reference_sets[name] = ref_iso
    return reference_sets


def gather_sample_files(sample_names):
    """Check files exist for the report."""
    sample_names = pd.read_csv(sample_names, header=None).iloc[:, 0].tolist()
    sample_names.sort()
    sample_files = {}
    for sample_name in sample_names:
        mapula_json = os.path.join(
            'merged_mapula_json',
            sample_name + '.merged.mapula.json')
        unmapped_stats = os.path.join(
            'unmapped_stats', sample_name + '.unmapped.stats')
        depth_beds = os.path.join(
            'depth_beds', sample_name + '.all_regions.bed.gz')
        expected_files = {
            'Mapula file': mapula_json,
            'Unmapped file': unmapped_stats,
            'depth beds': depth_beds}
        final_files = {
            'Mapula file': mapula_json,
            'Unmapped file': unmapped_stats,
            'depth beds': depth_beds}
        for name, file in expected_files.items():
            if os.path.exists(file):
                pass
            else:
                final_files[name] = 'None'
                sys.stdout.write(
                    'Missing {0} required for report for: {1}'.format(
                        name, sample_name))
        sample_files[sample_name] = final_files
    return sample_files


def main():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser(
        'Visualise mapula output',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        help="Input .json file to plot graphs from.",
        nargs='+',
        dest="json")
    parser.add_argument(
        '-c',
        help=(
            "Expected counts CSV. "
            "Required columns: reference,expected_count."),
        dest="counts",
        default=None,
        required=False)
    parser.add_argument(
        "--report_name",
        help="Report file name.",
        dest="name",
        required=False,
        default="report")
    parser.add_argument(
        "--references", nargs='+',
        help="reference files")
    parser.add_argument(
        "--params", default=None,
        help="A csv containing the parameter key/values")
    parser.add_argument(
        "--versions",
        help="directory contained CSVs containing name,version.")
    parser.add_argument(
        "--unmapped_stats", nargs='+', default=None,
        help="Fastcat stats for unmapped")
    parser.add_argument(
        "--revision", default='unknown',
        help="revision number")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit number")
    parser.add_argument(
        "--sample_names", nargs='+',
        help="sample names or barcodes")
    args = parser.parse_args()
    report_doc = report.WFReport(
        "Workflow Alignment Report",
        "wf-alignment",
        revision=args.revision,
        commit=args.commit)
    section = report_doc.add_section()
    section.markdown(
        "Results generated through the wf-alignment "
        "nextflow workflow provided by Oxford Nanopore "
        "Technologies")
    section.markdown('### Alignment statistics')
    section.markdown(
        "Statistics are provided for each barcode "
        "or sample name if a sample sheet was provided.")
    section.markdown(
        "Alignment was done using the following reference files:")
    all_references = get_references(args.references)
    for ref in args.references[::-1]:
        section.markdown('- ' + str(ref))
    sample_files = gather_sample_files(args.sample_names[0])
    for name, values in sample_files.items():
        if values['Unmapped file'] == "None":
            unmapped = pd.DataFrame(
                data=None, index=None, columns=[
                    'read_id', 'filename',
                    'sample_name', 'read_length',
                    'mean_quality', 'channel',
                    'read_number', 'start_time'],
                dtype=None, copy=None)
        else:
            unmapped = pd.read_csv(values['Unmapped file'], sep='\t')
        if str(values['depth beds']) != 'None':
            depth_sample = pd.read_csv(
                values['depth beds'], sep='\t', header=None)
        else:
            depth_sample = 'None'
        stats_panel = PlotMappingStats(
                json=values['Mapula file'],
                counts=args.counts,
                references=all_references,
                unmapped=unmapped,
                sample_name=name,
                depth_file=depth_sample)
        report_doc.add_section().markdown('##' + name)
        report_doc.add_section(section=stats_panel)
    # Versions and params
    report_doc.add_section(
        section=scomponents.version_table(args.versions))
    report_doc.add_section(
        section=scomponents.params_table(args.params))
    report_doc.write(args.name + '.html')


if __name__ == "__main__":
    main()
