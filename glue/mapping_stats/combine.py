import os
import sys
import argparse
from operator import add
from typing import List
from mapping_stats.common import MappingStats


class CombineMappingStats(MappingStats):

    def __init__(
        self,
        out_path: str,
        json_paths: List[str],
    ) -> None:
        super(CombineMappingStats, self).__init__()

        self.out_path = out_path
        self.json_paths = [
            os.path.abspath(p) 
            for p in json_paths
        ]

        self.data = self.DEFAULTS
        self.combine_data()
        self.write_data(out_path, self.data)

    def combine_data(
        self,
    ):
        for path in self.json_paths:
            data = self._load_data(path)

            self.alignment_count += data[
                self.ALN_COUNT
            ]
            self.read_count += data[
                self.READ_COUNT
            ]
            self.total_base_pairs += data[
                self.TOTAL_BASE_PAIRS
            ]
            self.unmapped_count += data[
                self.UNMAPPED_COUNT
            ]
            self.primary_count += data[
                self.PRIMARY_COUNT
            ]
            self.secondary_count += data[
                self.SECONDARY_COUNT
            ]
            self.supplementary_count += data[
                self.SUPPLEMENTARY_COUNT
            ]
            self.data[self.ALN_ACCS] = list(
                map(
                    add, 
                    self.alignment_accuracies, 
                    data[self.ALN_ACCS]
                ) 
            )
            self.data[self.READ_QUALS] = list(
                map(
                    add, 
                    self.read_qualities, 
                    data[self.READ_QUALS]
                ) 
            )
            self.data[self.READ_LENS] = list(
                map(
                    add, 
                    self.read_lengths, 
                    data[self.READ_LENS]
                ) 
            )

            for k, v in data[self.CONTIGS].items():
                self._add_or_update_contig(
                    name=k,
                    base_pairs=v['base_pairs'],
                    length=v['length']
                )

        self.update_read_n50()
        self.update_median_quality()
        self.update_median_accuracy()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Combine Mapping Stats")

    parser.add_argument(
        "-o", 
        "--OUT",
        default="merged.mapping-stats.json",
        help=(
            "Output merged json file. "
            "(default: merged.mapping-stats.json)"
        )
    )

    parser.add_argument(
        '-j',
        '--JSON',
        nargs='*',
        required=False,
        default=['mapping-stats.json']
    )

    args = parser.parse_args(sys.argv[1:])
    return args


def run_main():
    args = parse_arguments()
    CombineMappingStats(
        out_path=args.OUT,
        json_paths=args.JSON,
    )


if __name__ == '__main__':
    run_main()
