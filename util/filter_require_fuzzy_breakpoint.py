#!/usr/bin/env python3

import sys, os, re
import intervaltree as itree
import gzip
from collections import defaultdict
import csv
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

UTILDIR = os.path.dirname(os.path.abspath(__file__))

REF_EXONS_FILE = os.path.sep.join(
    [
        UTILDIR,
        "../../resources/GRCh38.ref_annot_exon_coords.tsv.gz",
    ]
)
FUZZY = 5


def main():

    usage = "usage: {} preds.collected\n\n".format(sys.argv[0])

    if len(sys.argv) < 2:
        exit(usage)

    preds_file = sys.argv[1]

    logger.info("-building itrees for exon and fuzzy breakpoints")
    chrom_itrees = build_fuzzy_breakpoint_itree(REF_EXONS_FILE)

    preds_reader = csv.DictReader(open(preds_file, "rt"), delimiter="\t")
    preds_writer = csv.DictWriter(
        sys.stdout,
        fieldnames=preds_reader.fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )

    logger.info("-filtering breakpoints")
    preds_writer.writeheader()

    for row in preds_reader:
        breakpoint = row["breakpoint"]
        break_lend, break_rend = breakpoint.split("--")
        chrom_lend, coord_lend = break_lend.split(":")
        chrom_rend, coord_rend = break_rend.split(":")

        coord_lend = int(coord_lend)
        coord_rend = int(coord_rend)

        if (len(chrom_itrees[chrom_lend][coord_lend : coord_lend + 1]) > 0) and (
            len(chrom_itrees[chrom_rend][coord_rend : coord_rend + 1]) > 0
        ):
            # within fuzzy break distance of exon boundaries
            preds_writer.writerow(row)

    logger.info("-done")

    sys.exit(0)


def build_fuzzy_breakpoint_itree(ref_exons_file):

    chr_itrees = defaultdict(lambda: itree.IntervalTree())

    with gzip.open(ref_exons_file, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            chrom, lend, rend = line.split("\t")
            lend = int(lend)
            rend = int(rend)
            chr_itrees[chrom][(lend - FUZZY - 1) : (lend + FUZZY + 1)] = True
            chr_itrees[chrom][(rend - FUZZY - 1) : (rend + FUZZY + 1)] = True

    return chr_itrees


if __name__ == "__main__":
    main()
