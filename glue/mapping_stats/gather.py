import sys
import pysam
import argparse
from typing import TypedDict, List
from mapping_stats.common import MappingStats


# Qscore lookup table
LOOKUP = [pow(10, -.1 * q) for q in range(100)]


class GatherMappingStats(MappingStats):

    def __init__(
        self,
        sam_path: str,
        out_path: str,
        json_path: str,
    ) -> None:
        """
        This class accepts a SAM or BAM files and
        iterates over each of the records within,
        calculating cumulative stats as it goes.
        After all records have been processed, it
        dumps the data to a JSON file.
        """
        super(GatherMappingStats, self).__init__()

        self.json_path = json_path
        self.out_path = out_path
        self.sam_path = sam_path

        self.alignments = pysam.AlignmentFile(
            sam_path, "r"
        )

        outfile = pysam.AlignmentFile(
            out_path, 
            "w", 
            template=self.alignments
        )

        self.data = self.load_data(json_path)
        self.update_data(outfile)
        self.write_data(json_path, self.data)

    def update_data(
        self,
        outfile: pysam.AlignmentFile = None
    ) -> None:
        """
        Carries out the heavy lifting of scanning
        each record within the input SAM or BAM and
        updates the internal json object with various
        metrics. In addition, if provided, will write
        each record to an output AlignmentFile,
        permitting records to be piped out to another
        location.
        """
        for aln in self.alignments.fetch(until_eof=True):
            if outfile:
                outfile.write(aln)

            if aln.is_supplementary:
                self.supplementary_count += 1
                continue

            if aln.is_secondary:
                self.secondary_count += 1
                continue

            # General stats for primary and unmapped reads
            self.total_base_pairs += aln.query_length
            self.read_count += 1

            self.update_read_length_dist(aln)
            self.update_read_n50()

            self.update_read_quality_dist(aln)
            self.update_median_quality()

            if aln.is_unmapped:
                self.unmapped_count += 1
                continue

            self.add_or_update_contig(aln, self.alignments)
            self.update_alignment_accuracy_dist(aln)
            self.update_median_accuracy()

            self.primary_count += 1


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Gather Mapping Stats")

    parser.add_argument(
        '-s',
        '--SAM',
        default=sys.stdin,
        help="Input sam/bam file. (default: stdin)"
    )

    parser.add_argument(
        "-o", 
        "--OUT",
        default='-',
        help="Output sam file. (default: stdout)"
    )

    parser.add_argument(
        '-j',
        '--JSON',
        required=False,
        default='mapping-stats.json'
    )

    args = parser.parse_args(sys.argv[1:])
    return args


def run_main():
    args = parse_arguments()
    GatherMappingStats(
        sam_path=args.SAM,
        out_path=args.OUT,
        json_path=args.JSON,
    )


if __name__ == '__main__':
    run_main()
