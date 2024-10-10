#!/usr/bin/env python3

import sys, os, re
import pysam
import argparse
import csv
import logging

csv.field_size_limit(sys.maxsize)

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
        "--input_bam", type=str, required=True, help="LRF-phase-1.mm2.sorted.bam"
    )

    parser.add_argument(
        "--candidates_with_reads",
        type=str,
        required=True,
        help="chimeric_read_candidates.FI_listing.with_reads",
    )

    parser.add_argument(
        "--candidates_filter_fail",
        type=str,
        required=True,
        help="chimeric_read_candidates.FI_listing.wAnnot.annot_filter.fail",
    )

    parser.add_argument(
        "--output_bam",
        required=True,
        type=str,
        help="output bam file containined the reads alignments of interest",
    )

    args = parser.parse_args()

    input_bam_file = args.input_bam
    candidates_with_reads_filename = args.candidates_with_reads
    candidates_filter_fail_filename = args.candidates_filter_fail
    output_bam_filename = args.output_bam

    failed_fusions = get_failed_fusions_neighbors_overlap(
        candidates_filter_fail_filename
    )

    reads_want = get_fusion_read_names(candidates_with_reads_filename, failed_fusions)

    extract_wanted_reads_into_bam(input_bam_file, reads_want, output_bam_filename)

    sys.exit(0)


def extract_wanted_reads_into_bam(input_bam_file, reads_want, output_bam_filename):

    bamreader = pysam.AlignmentFile(input_bam_file, "rb")

    bamwriter = pysam.AlignmentFile(output_bam_filename, "wb", template=bamreader)

    logger.info("-writing reads to {}".format(output_bam_filename))

    read_counter = 0

    for read in bamreader:
        if read.query_name in reads_want:
            bamwriter.write(read)
            read_counter += 1

    logger.info("-wrote {} reads to bam".format(read_counter))

    bamreader.close()
    bamwriter.close()

    return


def get_fusion_read_names(candidates_with_reads_filename, failed_fusions):

    logger.info(
        "-getting reads for failed fusions from: {}".format(
            candidates_with_reads_filename
        )
    )

    reads_want = set()

    with open(candidates_with_reads_filename) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["#FusionName"] in failed_fusions:
                reads = row["reads"].split(",")
                for read in reads:
                    reads_want.add(read)

    return reads_want


def get_failed_fusions_neighbors_overlap(candidates_filter_fail_filename):

    logger.info(
        "-getting neighbor overlap fusions from: {}".format(
            candidates_filter_fail_filename
        )
    )

    failed_fusions = set()

    with open(candidates_filter_fail_filename) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            m = re.search("NEIGHBORS_OVERLAP", row["Annot_Fail_Reason"])
            if m is not None:
                fusion_name = row["#FusionName"]
                failed_fusions.add(fusion_name)

    return failed_fusions


if __name__ == "__main__":
    main()
