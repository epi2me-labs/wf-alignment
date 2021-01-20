import os
import sys
import json
import argparse
from math import pi
import pandas as pd
from typing import List
from aplanat import hist
from bokeh.embed import json_item
from bokeh.plotting import figure
from bokeh.transform import cumsum
from aplanat.report import HTMLReport
from bokeh.palettes import Category20c


class VisualiseMappingStats(object):

    EXPORT = 'VisualiseMappingStats_Export.json'
    REPORT = 'VisualiseMappingStats_Report.html'

    def __init__(
        self,
        title: str,
        outdir: str,
        json_path: str,
    ) -> None:
        self.title = title
        self.outdir = outdir
        self.json_path = json_path
        self.data = self.load_data(self.json_path)

        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        self.summary_plot = self.plot_alignment_summary()
        self.accuracy_plot = self.plot_accuracy_distribution()
        self.qscore_plot = self.plot_qscore_distribution()
        self.length_plot = self.plot_read_length_distribution()
        self.depth_frame = self.get_read_depth_dataframe()

        self.build_report(
            title=title,
            outdir=outdir,
            plots=[
                self.summary_plot, 
                self.accuracy_plot, 
                self.qscore_plot, 
                self.length_plot
            ],
            tables=[
                self.depth_frame
            ]
        )

        self.export_data(
            outdir=outdir,
            summary_plot=self.summary_plot,
            accuracy_plot=self.accuracy_plot,
            qscore_plot=self.qscore_plot,
            length_plot=self.length_plot
        )

    @staticmethod
    def load_data(
        path: str,
    ) -> dict:
        try:
            with open(path) as data:
                try:
                    return json.load(data)
                except json.decoder.JSONDecodeError:
                    print('Error loading data from {}'.format(path))
                    raise
        except IOError:
            print("Path {} cannot be opened".format(path))
            raise

    def export_data(
        self, 
        outdir: str,
        **plots: List[figure]
    ):
        path = os.path.join(outdir, self.EXPORT)

        export_data = {}
        for k, v in plots.items():
            export_data[k] = json_item(v)

        with open(path, 'w') as out:
            out.write(json.dumps(export_data))

    def build_report(
        self, 
        title: str,
        outdir: str,
        plots: List[figure],
        tables: List[figure]
    ):
        path = os.path.join(outdir, self.REPORT)
        report = HTMLReport(title=title)

        for plot in plots:
            report.plot(plot)

        for table in tables:
            report.table(table)

        report.render()
        report.write(path)
        
    def plot_alignment_summary(self):
        x = {
            'Primary alignments': self.data['primary_count'],
            'Secondary alignments': self.data['secondary_count'],
            'Supplementary alignments': self.data['supplementary_count'],
            'Unmapped': self.data['unmapped_count'],
        }

        data = pd.Series(x).reset_index(name='value').rename(
            columns={'index':'type'}
        )
        data['angle'] = data['value']/data['value'].sum() * 2*pi
        data['color'] = Category20c[len(x)]

        plot = figure(
            plot_height=350, 
            title="Alignments by type", 
            toolbar_location=None,
            tools="hover", 
            tooltips="@type: @value", 
            x_range=(-0.5, 1.0)
        )

        plot.wedge(
            x=0, 
            y=1, 
            radius=0.4,
            start_angle=cumsum('angle', include_zero=True),
            end_angle=cumsum('angle'),
            line_color="white", 
            fill_color='color',
            legend_field='type',
            source=data
        )

        plot.axis.axis_label=None
        plot.axis.visible=False
        plot.grid.grid_line_color = None
        return plot

    def plot_accuracy_distribution(self):
        counts = []
        for index, item in enumerate(self.data['alignment_accuracies']):
            for i in range(item):
                counts.append(index/10)
        
        plot = hist.histogram(
            [counts],
            bins=100,
            height=300,
            x_axis_label='Accuracy %',
            y_axis_label='Count',
            title='Distribution of average alignment accuracy'
        )
        return plot

    def plot_qscore_distribution(self):
        counts = []
        for index, item in enumerate(self.data['read_qualities']):
            for i in range(item):
                counts.append(index/10)
        
        plot = hist.histogram(
            [counts],
            bins=600,
            height=300,
            x_axis_label='Accuracy %',
            y_axis_label='Count',
            title='Distribution of mean qscore per read'
        )
        return plot

    def plot_read_length_distribution(self):
        counts = []
        for index, item in enumerate(self.data['read_lengths']):
            for i in range(item):
                counts.append(index * 50)
        
        plot = hist.histogram(
            [counts],
            bins=1000,
            height=300,
            x_axis_label='Length',
            y_axis_label='Count',
            title='Distribution of read lengths'
        )
        return plot

    def get_read_depth_dataframe(self):
        df = pd.DataFrame.from_dict(
            self.data['contigs'], 
            orient='index', 
            columns=['length', 'base_pairs']
        )
        return df


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Visualise Mapping Stats")

    parser.add_argument(
        '-t',
        '--title',
        required=False,
        default='Visualise Mapping Stats -> Report'
    )

    parser.add_argument(
        '-o',
        '--outdir',
        required=False,
        default='./output/'
    )

    parser.add_argument(
        '-j',
        '--json',
        required=False,
        default='mapping_stats.json'
    )

    args = parser.parse_args(sys.argv[1:])
    return args


def run_main():
    args = parse_arguments()
    VisualiseMappingStats(
        title=args.title,
        outdir=args.outdir,
        json_path=args.json,
    )


if __name__ == '__main__':
    run_main()
