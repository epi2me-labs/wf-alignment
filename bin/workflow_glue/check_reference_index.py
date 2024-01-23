#!/usr/bin/env python
"""Parse references from MMI index file."""


import struct
from .util import get_named_logger, wf_parser  # noqa: ABS101

MM_IDX_MAGIC = b"MMI\2"


def parse_mmi(mmi_file):
    """Get references from the MMI file."""
    # Minimap2 Index Description: .mmi -> binary file.
    # Created in:
    # https://github.com/lh3/minimap2/blob/fc2e1607d73ea4870e6f396697c79051aff23eed/index.c
    # see this Github Issue: https://github.com/lh3/minimap2/issues/820
    # It contains 4 parts:
    # 1) Magic string: ("MMI\x02", MM_IDX_MAGIC in the source) -> first four bytes.
    # 2) Followed by 5 integers of constants: #L273
    #        mi = mm_idx_init(int w, int k, int bucket_bits, int B, int flag);
    #        - w-> minimizer width w
    #        - k-> window length k,
    #        - bucket_bits (b)-> hardcoded value 14: #L403, 421
    #        - B-> number of sequences ?? -> L418
    #        - flags: is_hpc, name...
    # 3) REFERENCES IDS: The reference information is stored: #L427C2-L427C2
    #        - name: the name in ascii
    #        - offset: sum_len
    #        - len: the length as integer.
    #        - is_alt: ALT contigs
    #    This is repeated for each sequence.
    # 4) SEQUENCES
    sequences = dict()
    # open file in read binary mode
    with open(mmi_file, "rb") as fd:
        # type I: convert C values to python
        # https://docs.python.org/3/library/struct.html#format-characters
        magic = fd.read(4)
        if magic != MM_IDX_MAGIC:
            raise ValueError(f"{mmi_file} does not appear to be a minimap2 index.")
        # w, k, bucket bits, n_seqs, flags
        _, _, _, n_seq, _ = struct.unpack("5I", fd.read(20))
        # contents is then (len name, name, len seq)
        for _ in range(n_seq):
            name_length = struct.unpack("B", fd.read(1))[0]
            name = fd.read(name_length).decode("ascii")
            length = struct.unpack("I", fd.read(4))[0]
            sequences[name] = length
    return sequences


def get_fai_references(fai_file):
    """Get references from the FASTA fai file."""
    with open(fai_file, 'r') as fasta_ids:
        # check references
        return {i.split('\t')[0] for i in fasta_ids.readlines()}


def main(args):
    """Run entry point."""
    logger = get_named_logger("check_reference_index")
    fasta_names = get_fai_references(args.fasta_fai)
    sequences = set(parse_mmi(args.mmi_file).keys())
    if len(sequences) != len(fasta_names):
        raise Exception(
            f"Number of sequences in the MMI file ({len(sequences)}) doesn't match "
            f"the number of sequences provided references ({len(fasta_names)}).")
    if not sequences - fasta_names:
        logger.info("All the MMI references are in the FASTA reference")
    else:
        raise Exception(
            "The next references found in the MMI file are not in the FASTA file: "
            f"{sequences - fasta_names}")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "--mmi_file",
        help="MMI index file"
    )
    parser.add_argument(
        "--fasta_fai",
        help="Fai index file."
    )
    return parser
