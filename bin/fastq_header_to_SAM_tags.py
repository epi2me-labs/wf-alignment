#!/usr/bin/env python
"""Re-write fastq headers to be SAM tag compatible."""
import os
import re
import sys
import glob
import argparse
import pysam


def errprint(*args, **kwargs):
    """Prints to stderr."""
    sys.stderr.write(*args, **kwargs)
    sys.stderr.write('\n')


def parse_barcode(read_comment):
    """Parse barcode from FASTQ read comment."""
    if not read_comment:
        return "none"

    match = re.search(
        'barcode=([0-9]+)',
        read_comment
    )

    if match:
        return "barcode{:02d}".format(
            int(match.group(1))
        )

    return "unclassified"


def parse_runid(read_comment):
    """Parses the runid from a fastq comment."""
    if not read_comment:
        return "none"

    match = re.search(
        'runid=([0-9A-Za-z]+)',
        read_comment
    )

    if match:
        return match.group(1)

    return "none"


def rewrite_fastq_header(path, out_file):
    """Rewrite header for each read in a fastq file."""
    with pysam.FastxFile(path) as fastq_file:

        # Iterate over fastq reads
        for entry in fastq_file:

            # Pass the barcoding info from the read comment
            comment = ""
            if entry.comment:
                barcode = parse_barcode(entry.comment)
                runid = parse_runid(entry.comment)

                if barcode:
                    comment = "BC:Z:" + barcode

                if runid:
                    runid = "RD:Z:" + runid
                    comment = '\t'.join([comment, runid]) if comment else runid

            # Set the barcoding info on the read and write it
            entry.comment = comment
            out_file.write(str(entry) + "\n")


def get_file_names(path, endswith=".fastq"):
    """Discover files at the given path."""
    if path == "/dev/stdin":
        yield path
        return

    if path.endswith(endswith):
        yield path
        return

    if os.path.isdir(path):
        for file in glob.glob(os.path.join(path, f"*{endswith}")):
            yield file
        return

    errprint("Could not find {}".format(path))


def find_and_rewrite_fastq_headers(fastq_paths, out_file):
    """Discover fastq files and rewrites their headers."""
    with open(out_file, mode='w') as open_out_file:
        for path in fastq_paths:
            for filename in get_file_names(path):
                rewrite_fastq_header(filename, open_out_file)


def parse_args(argv):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Find and re-write fastq headers",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-o",
        "--output-file",
        dest="OUT_FILE",
        type=str,
        default="/dev/stdout",
        help="Output file. (default: stdout)"
    )

    parser.add_argument(
        "FASTQ",
        nargs="*",
        type=str,
        help="FASTQ files or folders containing FASTQ files",
        default=['/dev/stdin']
    )

    args = parser.parse_args(argv)

    return args


def main(argv=sys.argv[1:]):
    """Run main function."""
    args = parse_args(argv)
    find_and_rewrite_fastq_headers(args.FASTQ, args.OUT_FILE)


if __name__ == '__main__':
    main()
