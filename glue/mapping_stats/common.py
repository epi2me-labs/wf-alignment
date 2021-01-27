import os
import math
import json
import pysam
from collections.abc import Iterable
from typing import Union, TypedDict, List


# Qscore lookup table
LOOKUP = [pow(10, -.1 * q) for q in range(100)]


class MappingStats(object):
    ALN_COUNT = 'alignment_count'
    READ_COUNT = 'read_count'
    TOTAL_BASE_PAIRS = 'total_base_pairs'

    UNMAPPED_COUNT = 'unmapped_count'
    PRIMARY_COUNT = 'primary_count'
    SECONDARY_COUNT = 'secondary_count'
    SUPPLEMENTARY_COUNT = 'supplementary_count'

    ALN_ACCS = 'alignment_accuracies'
    MED_ACC = 'median_accuracy'
    READ_QUALS = 'read_qualities'
    MED_QUAL = 'median_quality'
    READ_LENS = 'read_lengths'
    READ_N50 = 'read_n50'

    CONTIGS = 'contigs'

    class Contig(TypedDict):
        length: int
        base_pairs: int

    DEFAULTS = {
        # General stats
        ALN_COUNT: 0,
        READ_COUNT: 0,
        TOTAL_BASE_PAIRS: 0,
        # Mapping stats
        UNMAPPED_COUNT: 0,
        PRIMARY_COUNT: 0,
        SECONDARY_COUNT: 0,
        SUPPLEMENTARY_COUNT: 0,
        # Distributions
        ALN_ACCS: [0] * 1001,
        MED_ACC: 0,
        READ_QUALS: [0] * 600,
        MED_QUAL: 0,
        READ_LENS: [0] * 1000,
        READ_N50: 0,
        # Other
        CONTIGS: {}
    }

    def __init__(self):
        self.data = {}

    def _load_data(
        self,
        path: str,
    ):
        with open(path) as data:
            try:
                return json.load(data)
            except json.decoder.JSONDecodeError:
                print('Error loading data file {}.'.format(
                    path
                ))
                raise

    def load_data(
        self,
        path: str,
        default: dict = DEFAULTS
    ) -> dict:
        """
        Attempts to find and open a file located
        at self.json_path. If it cannot be found,
        for instance because it does not exist, 
        this method will return a default dict of
        keys and starting values.
        """
        if path is not None and os.path.exists(path):
            return self._load_data(path)
        return default

    def write_data(
        self, 
        path: str,
        data: dict,
    ) -> None:
        """
        Writes self.data out to a file located
        at self.json_path.
        """
        with open(path, 'w') as out:
            json.dump(data, out)

    #
    # Alignment count
    #
    @property
    def alignment_count(self):
        return self.data[self.ALN_COUNT]

    @alignment_count.setter
    def alignment_count(self, value):
        self.data[self.ALN_COUNT] = value

    #
    # Read count
    #
    @property
    def read_count(self):
        return self.data[self.READ_COUNT]

    @read_count.setter
    def read_count(self, value):
        self.data[self.READ_COUNT] = value

    #
    # Base pairs
    #
    @property
    def total_base_pairs(self):
        return self.data[self.TOTAL_BASE_PAIRS]

    @total_base_pairs.setter
    def total_base_pairs(self, value):
        self.data[self.TOTAL_BASE_PAIRS] = value

    #
    # Unmapped reads
    #
    @property
    def unmapped_count(self):
        return self.data[self.UNMAPPED_COUNT]

    @unmapped_count.setter
    def unmapped_count(self, value):
        self.data[self.UNMAPPED_COUNT] = value

    #
    # Primary alignments
    #
    @property
    def primary_count(self):
        return self.data[self.PRIMARY_COUNT]

    @primary_count.setter
    def primary_count(self, value):
        self.data[self.PRIMARY_COUNT] = value

    #
    # Secondary alignments
    #
    @property
    def secondary_count(self):
        return self.data[self.SECONDARY_COUNT]

    @secondary_count.setter
    def secondary_count(self, value):
        self.data[self.SECONDARY_COUNT] = value

    #
    # Supplementary alignments
    #
    @property
    def supplementary_count(self):
        return self.data[self.SUPPLEMENTARY_COUNT]

    @supplementary_count.setter
    def supplementary_count(self, value):
        self.data[self.SUPPLEMENTARY_COUNT] = value

    #
    # Read qualities
    #
    @property
    def read_qualities(self):
        return self.data[self.READ_QUALS]

    @read_qualities.setter
    def read_qualities_setter(self, value):
        print(
            "Warning: do not set qualities "
            "distribution directly. (But if "
            "you really want to, use "
            "obj.data['read_qualities'] "
            "= value)."
        )
        return

    def _update_read_quality_dist(
        self, 
        quality: Union[float, int, None]
    ) -> None:
        self.read_qualities[int(quality / 0.1)] += 1
    
    def update_read_quality_dist(
        self, 
        aln: pysam.AlignedSegment
    ) -> None:
        quality = self.get_alignment_mean_qscore(
            aln.query_qualities
        )
        self._update_read_quality_dist(quality)

    #
    # Median quality
    #
    @property
    def median_quality(self):
        return self.data[self.MED_QUAL]

    @median_quality.setter
    def median_quality(self, value):
        self.data[self.MED_QUAL] = value

    def update_median_quality(self):
        self.median_quality = self.get_median_from_frequency_dist(
            self.read_qualities, 0.1
        )

    #
    # Alignment accuracies
    #
    @property
    def alignment_accuracies(self):
        return self.data[self.ALN_ACCS]

    @alignment_accuracies.setter
    def alignment_accuracies_setter(self, value):
        print(
            "Warning: do not set accuracies "
            "distribution directly. (But if "
            "you really want to, use "
            "obj.data['alignment_accuracies'] "
            "= value)."
        )
        return

    def _update_alignment_accuracy_dist(
        self, 
        accuracy: Union[float, int]
    ) -> None:
        self.alignment_accuracies[int(accuracy / 0.1)] += 1

    def update_alignment_accuracy_dist(
        self, 
        aln: pysam.AlignedSegment
    ) -> None:
        accuracy = self.get_alignment_accuracy(aln) or 0
        self._update_alignment_accuracy_dist(accuracy)

    #
    # Median accuracy
    #
    @property
    def median_accuracy(self):
        return self.data[self.MED_ACC]

    @median_accuracy.setter
    def median_accuracy(self, value):
        self.data[self.MED_ACC] = value

    def update_median_accuracy(self):
        self.median_accuracy = self.get_median_from_frequency_dist(
            self.alignment_accuracies, 0.1
        )

    #
    # Read lengths
    #
    @property
    def read_lengths(self):
        return self.data[self.READ_LENS]

    @read_lengths.setter
    def read_lengths(self, value):
        print(
            "Warning: do not set read lengths "
            "distribution directly. (But if "
            "you really want to, use "
            "obj.data['read_lengths'] = value)."
        )
        return

    def _update_read_length_dist(
        self, 
        length: int
    ):
        self.read_lengths[int(length / 50)] += 1

    def update_read_length_dist(
        self, 
        aln: pysam.AlignedSegment
    ) -> None:
        self._update_read_length_dist(aln.query_length)

    #
    # Read N50
    #
    @property
    def read_n50(self):
        return self.data[self.READ_N50]

    @read_n50.setter
    def read_n50(self, value):
        self.data[self.READ_N50] = value

    def update_read_n50(self):
        self.read_n50 = self.get_n50_from_frequency_dist(
            self.read_lengths, 50, self.total_base_pairs
        )

    #
    # Contigs
    #
    @property
    def contigs(self):
        return self.data[self.CONTIGS]

    @contigs.setter
    def contigs(self, value):
        print(
            "Warning: do not set contigs "
            "distribution directly. (But if "
            "you really want to, use "
            "obj.data['contigs'] = value)."
        )
        return

    def _add_or_update_contig(
        self,
        name: str,
        base_pairs: int,
        length: int
    ):
        """
        Either gets or creates an entry in self.contigs 
        for the contig matching name. If found, adds the 
        base pairs to the total count.
        """
        entry = self.contigs.get(name)
    
        if entry:
            self.contigs[
                name
            ]['base_pairs'] += base_pairs
            return

        self.contigs[name] = self.Contig(
            length=length,
            base_pairs=base_pairs
        )

    def add_or_update_contig(
        self,
        aln: pysam.AlignedSegment,
        alignments: pysam.AlignmentFile
    ) -> None:
        """
        Given an AlignedSegment, gets the reference
        it was aligned to, and either gets or creates
        an entry in self.contigs for that reference.
        If found, adds the query's base pairs to the
        total count.
        """
        name = aln.reference_name
        self._add_or_update_contig(
            name=name,
            base_pairs=aln.query_alignment_length,
            length=alignments.get_reference_length(name)
        )

    #
    # Utilities
    #
    @staticmethod
    def get_median_from_frequency_dist(
        arr: Iterable,
        width: Union[int, float]
    ):
        """
        Returns the median value from an array whose 
        positions represent frequency counts of values 
        at that index in a given range.
        """
        arr_sum: int = sum(arr)
        half_way_pos = arr_sum / 2
        is_odd = bool(arr_sum % 2)

        lower = None

        accumulator = 0
        for idx, count in enumerate(arr):
            accumulator += count

            if accumulator > half_way_pos:
                if lower is not None:
                    avg = ((lower * width) + (idx * width)) / 2
                    return float(format(avg, '.2f'))

                return float(format(idx * width, '.2f'))

            if accumulator == half_way_pos:
                if is_odd:
                    return float(format(idx * width, '.2f'))

                if lower is None:
                    lower = idx
                    continue

    @staticmethod
    def get_alignment_accuracy(
        alignment: pysam.AlignedSegment
    ):
        """
        Returns the percentage accuracy of a given aligned 
        segment as a float.
        """
        if alignment.is_unmapped:
            return None

        el_count = [0] * 10
        for el in alignment.cigartuples:
            el_count[el[0]] += el[1]

        # Number of ambiguous bases
        nn = 0
        if alignment.has_tag("nn"):
            nn = alignment.get_tag("nn")

        # NM = #mismatches + #I + #D + #ambiguous_bases
        nm = alignment.get_tag("NM")
        dels = el_count[2]
        ins = el_count[1]

        mismatches = nm - dels - ins - nn
        matches = el_count[0] - mismatches

        accuracy = float(matches) / (matches + nm) * 100

        return accuracy

    @staticmethod
    def get_n50_from_frequency_dist(
        arr: Iterable,
        width: int,
        total: int
    ) -> Union[float, int]:
        """
        Calculates the N50 from an array whose positions represent
        frequency counts of values at that index in a given range.

        N50 is defined as:

        Length N for which 50% of all bases in the sequences
        are in a sequence of length L < N

        Returns the approximate 'size' of the value at which
        the N50 is achieved.
        """
        n50 = 0
        cumulative_value = 0
        half_total = total / 2

        for bin_number, bin_count in enumerate(arr):
            bin_value = (bin_number * width) + (width / 2)

            cumulative_value += bin_value * bin_count
            if cumulative_value >= half_total:
                n50 = bin_value
                break

        return n50

    @staticmethod
    def get_alignment_mean_qscore(
        scores: List[int]
    ) -> Union[float, None]:
        """
        Returns the phred score corresponding to the mean of
        the probabilities associated with the phred scores 
        provided.
        """
        if scores is None:
            return None

        if not scores:
            return 0.0

        sum_prob = 0.0
        for val in scores:
            sum_prob += LOOKUP[val]

        mean_prob = sum_prob / len(scores)

        return -10.0 * math.log10(mean_prob)