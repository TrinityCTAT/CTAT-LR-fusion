#!/usr/bin/env python3

import sys, os, re
import intervaltree as itree
import gzip
from collections import defaultdict
import csv
import logging
import argparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


FUZZY = 5


def main():

    parser = argparse.ArgumentParser(
        description="Filter fusions for those having breakpoints within +/- distance from refernece annotations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--fusions", type=str, required=True, help="fusion predictions tsv file"
    )
    parser.add_argument(
        "--genome_lib_dir",
        type=str,
        default=os.environ.get("CTAT_GENOME_LIB", None),
        required=False,
        help="path to CTAT genome lib dir",
    )
    parser.add_argument(
        "--fuzzy_dist",
        type=int,
        required=False,
        default=FUZZY,
        help="distance allowed +/- reference exon boundaries",
    )

    args = parser.parse_args()

    genome_lib_dir = args.genome_lib_dir
    preds_file = args.fusions
    fuzzy_dist = args.fuzzy_dist

    if genome_lib_dir is None:
        raise RuntimeError("must specify --genome_lib_dir")

    # build exon trees
    logger.info("-building itrees for exon and fuzzy breakpoints")
    ref_annot_gtf_exons = os.path.join(genome_lib_dir, "ref_annot.gtf.mini.sortu")
    chrom_itrees = build_fuzzy_breakpoint_itree(ref_annot_gtf_exons, fuzzy_dist)

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

        break_lend = row["LeftBreakpoint"]
        break_rend = row["RightBreakpoint"]

        chrom_lend, coord_lend, strand_lend = break_lend.split(":")
        chrom_rend, coord_rend, strand_rend = break_rend.split(":")

        coord_lend = int(coord_lend)
        coord_rend = int(coord_rend)

        if (len(chrom_itrees[chrom_lend][coord_lend : coord_lend + 1]) > 0) and (
            len(chrom_itrees[chrom_rend][coord_rend : coord_rend + 1]) > 0
        ):
            # within fuzzy break distance of exon boundaries
            preds_writer.writerow(row)

    logger.info("-done")

    sys.exit(0)


def build_fuzzy_breakpoint_itree(ref_exons_file, fuzzy_dist):

    chr_itrees = defaultdict(lambda: itree.IntervalTree())

    with open(ref_exons_file, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            chrom, lend, rend = vals[0], vals[3], vals[4]
            lend = int(lend)
            rend = int(rend)
            chr_itrees[chrom][(lend - fuzzy_dist - 1) : (lend + fuzzy_dist + 1)] = True
            chr_itrees[chrom][(rend - fuzzy_dist - 1) : (rend + fuzzy_dist + 1)] = True

    return chr_itrees


if __name__ == "__main__":
    main()
