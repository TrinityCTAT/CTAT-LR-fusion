#!/usr/bin/env python3

import sys
import logging
import argparse
import pysam


def main():

    parser = argparse.ArgumentParser(
        description="filter mm2 bam for chimeric alignments",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_bam",
        type=str,
        required=True,
        help="input bam filename",
    )

    parser.add_argument(
        "--output_bam",
        type=str,
        required=True,
        help="output bam filename",
    )

    args = parser.parse_args()

    input_bam_filename = args.input_bam
    output_bam_filename = args.output_bam

    bamreader = pysam.AlignmentFile(input_bam_filename, "rb")

    if (not "SO" in bamreader.header.as_dict()["HD"]) or (
        bamreader.header.as_dict()["HD"]["SO"] != "unsorted"
    ):
        raise RuntimeError(
            "Error, file: {} must be coordinate sorted".format(input_bam_filename)
        )

    bamwriter = pysam.AlignmentFile(output_bam_filename, "wb", template=bamreader)

    process_bam(bamreader, bamwriter)

    sys.exit(0)


def process_bam(bam_reader, bam_writer):

    prev_read_name = ""
    reads = list()

    for read in bam_reader:
        read_name = read.query_name
        if read_name != prev_read_name:
            evaluate_chimeric_read_candidates(reads, bam_writer)
            reads.clear()

        prev_read_name = read_name
        reads.append(read)

    # get last one.
    evaluate_chimeric_read_candidates(reads, bam_writer)

    return


def evaluate_chimeric_read_candidates(reads, bam_writer):

    if len(reads) < 2:
        return

    # ensure there's a supplementary alignment
    has_supplementary_alignment = False
    for read in reads:
        if read.is_supplementary:
            has_supplementary_alignment = True
            break

    if not has_supplementary_alignment:
        return

    # write candidates
    for read in reads:
        bam_writer.write(read)

    return


if __name__ == "__main__":
    main()
