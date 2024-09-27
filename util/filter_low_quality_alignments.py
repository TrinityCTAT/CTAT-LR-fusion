#!/usr/bin/env python

import sys, os, re
import pysam
import argparse
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="__add_descr__",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_namesorted_bam",
        type=str,
        required=True,
        help="bam file, reads sorted by name NOT COORDINATE",
    )
    parser.add_argument(
        "--output_bam", type=str, required=True, help="output bam filename"
    )

    parser.add_argument(
        "--min_quality",
        "-Q",
        dest="min_quality",
        type=int,
        required=True,
        help="minimum quality value for primary alignments",
    )

    parser.add_argument(
        "--debug", required=False, action="store_true", default=False, help="debug mode"
    )

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    args = parser.parse_args()

    input_bam_file = args.input_namesorted_bam
    output_bam_filename = args.output_bam

    min_quality = args.min_quality

    logger.info("-parsing {}".format(input_bam_file))
    samreader = pysam.AlignmentFile(input_bam_file, "rb")

    if (not "SO" in samreader.header.as_dict()["HD"]) or samreader.header.as_dict()[
        "HD"
    ]["SO"] != "queryname":
        raise RuntimeError(
            "Error, file: {} must be sorted by queryname (samtools sort -n)  not coordinate".format(
                input_bam_filename
            )
        )

    bamwriter = pysam.AlignmentFile(output_bam_filename, "wb", template=samreader)

    logger.info("-writing to {}".format(output_bam_filename))

    prev_query_name = ""

    reads_for_single_query = list()

    total_reads = 0
    reads_written = 0

    for read in samreader:
        total_reads += 1
        query_name = read.query_name
        if query_name != prev_query_name:
            # start a new one.
            if len(reads_for_single_query) > 0:
                reads_written += process_alignments(
                    reads_for_single_query, bamwriter, min_quality
                )
                reads_for_single_query.clear()

        reads_for_single_query.append(read)
        prev_query_name = query_name

    if len(reads_for_single_query) > 0:
        reads_written += process_alignments(
            reads_for_single_query, bamwriter, min_quality
        )

    frac_reads_written = reads_written / total_reads

    logger.info(
        "done. {}/{} = {:.1f}% reads passed quality filtering".format(
            reads_written, total_reads, frac_reads_written * 100
        )
    )


def process_alignments(reads_list, bam_writer, min_quality):

    assert len(reads_list) > 0, "Error, no reads in reads_list"

    query_name = reads_list[0].query_name

    # should have at least one primary
    primary_aligns = list()
    secondary_aligns = list()

    for read in reads_list:
        if not read.is_secondary:
            primary_aligns.append(read)
        else:
            secondary_aligns.append(read)

    if len(primary_aligns) == 0:
        logger.warn("no primary alignment identified for {}".format(query_name))
        return 0

    if len(primary_aligns) > 1:
        # take the maximum mapq score
        primary_aligns = sorted(
            primary_aligns, key=lambda x: x.mapping_quality, reverse=True
        )

    primary_mapping_quality = primary_aligns[0].mapping_quality

    if primary_mapping_quality < min_quality:
        # not reporting records.
        return 0

    else:
        # report the alignments.
        for read in reads_list:
            bam_writer.write(read)

    return len(reads_list)


if __name__ == "__main__":
    main()
