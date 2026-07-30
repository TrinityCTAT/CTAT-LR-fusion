#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
import pandas as pd
import csv
from decimal import Decimal, ROUND_CEILING


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


def main():
    
    parser = argparse.ArgumentParser(description="filtering fusion calls based on read counts", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--fusions_input", type=str, required=True, help="fusions input file")
    parser.add_argument("--filtered_fusions_output", type=str, required=True, help="name for filtered fusions output file")
    
    parser.add_argument("--min_num_LR", default=1, type=int, help="min number of long reads with canonical splice support")
    parser.add_argument("--min_LR_novel_junction_support", type=int, default=2, help="min number of long reads with non-canonical splice support")
    parser.add_argument("--min_J", type=int, default=1, help="min number of Illumina junction reads with canonical splice breakpoints")
    parser.add_argument("--min_sumJS", type=int, default=1, help="min number of Illumina reads supporting junction and spanning frags summed")
    parser.add_argument("--min_novel_junction_support", type=int, default=3, help="min number of junction reads with non-canonical splice support")
    parser.add_argument("--min_FFPM", type=float, default=0.1, help="min FFPM value for long or short reads.  If short reads >= min_FFPM and long reads < min_FFPM, still reported")
    parser.add_argument("--num_LR_total", type=int, required=True, help="total number of long reads used to compute fusion-pair FFPM")
    parser.add_argument("--min_frac_dom_iso", type=float, default=0.05, help="minimum fraction of the dominant long-read isoform count included in fusion-pair FFPM")

    args = parser.parse_args()
    if args.num_LR_total <= 0:
        parser.error("--num_LR_total must be a positive integer")

    fusions_input_filename = args.fusions_input
    fusions_output_filename = args.filtered_fusions_output

    min_num_LR = args.min_num_LR
    min_LR_novel_junction_support = args.min_LR_novel_junction_support
    min_J = args.min_J
    min_sumJS = args.min_sumJS
    min_novel_junction_support = args.min_novel_junction_support
    min_FFPM = args.min_FFPM
    num_LR_total = args.num_LR_total
    min_frac_dom_iso = args.min_frac_dom_iso

    data = pd.read_csv(fusions_input_filename, sep="\t", quotechar='"')

    LR_ok = (
        (
            (data.SpliceType == "ONLY_REF_SPLICE")
            & (data.num_LR >= min_num_LR)
        )
        | (data.num_LR >= min_LR_novel_junction_support)
    ).fillna(False)

    LR_supported_counts = data.num_LR.where(LR_ok)
    max_pair_LR = LR_supported_counts.groupby(
        data["#FusionName"], sort=False, dropna=False
    ).transform("max")
    if min_frac_dom_iso == 0:
        LR_contributor = LR_ok
    else:
        LR_contributor = (
            LR_ok & (data.num_LR >= min_frac_dom_iso * max_pair_LR)
        ).fillna(False)

    pair_LR_count = (
        data.num_LR.where(LR_contributor, 0)
        .groupby(data["#FusionName"], sort=False, dropna=False)
        .transform("sum")
    )
    min_pair_LR_reads = int(
        (
            Decimal(str(min_FFPM))
            * Decimal(num_LR_total)
            / Decimal("1000000")
        ).to_integral_value(rounding=ROUND_CEILING)
    )
    pair_LR_abundant = (pair_LR_count >= min_pair_LR_reads).fillna(False)

    if 'JunctionReadCount' in data.columns:
        row_SR_abundant = (data.FFPM >= min_FFPM).fillna(False)
        SR_ok = (
            (
                (
                    (data.JunctionReadCount >= min_J)
                    & (data.SpliceType == "ONLY_REF_SPLICE")
                )
                | (data.JunctionReadCount >= min_novel_junction_support)
            )
            & (data.JunctionReadCount + data.SpanningFragCount >= min_sumJS)
        ).fillna(False)
        keep = (
            (LR_ok & (pair_LR_abundant | row_SR_abundant))
            | (SR_ok & row_SR_abundant)
        )
    else:
        keep = LR_ok & pair_LR_abundant

    data_filtered = data[keep]
    
    

    data_filtered.to_csv(fusions_output_filename, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    sys.exit(0)

    


if __name__=='__main__':
    main()

