#!/usr/bin/env python3

import sys, os, re
import pandas as pd
import logging
import argparse


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description="merges the mm2 fusion and FI fusion reports", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--mm2_fusions", type=str, required=True, help="minimap2-fusion report")
    parser.add_argument("--FI_fusions", type=str, required=True, help="FI-fusion report")
    parser.add_argument("--output_file", type=str, required=True, help="output filename for merged data table")

    args = parser.parse_args()

    mm2_fusions_filename = args.mm2_fusions
    FI_fusions_filename = args.FI_fusions
    output_filename = args.output_file

    logger.info("-parsing {}".format(mm2_fusions_filename))
    mm2_df = pd.read_csv(mm2_fusions_filename, sep="\t")

    logger.info("-parsing {}".format(FI_fusions_filename))
    FI_df = pd.read_csv(FI_fusions_filename, sep="\t")


    logger.info("-merging data frames.")
    merged_df = pd.merge(mm2_df, FI_df,
                         on=['#FusionName', 'LeftLocalBreakpoint', 'RightLocalBreakpoint', 'LeftBreakpoint', 'RightBreakpoint', 'SpliceType'],
                         how='outer')
    

    logger.info("-writing output: {}".format(output_filename))
    merged_df.to_csv(output_filename, sep="\t", index=False, na_rep="NA")

    sys.exit(0)

if __name__=='__main__':
    main()
    
