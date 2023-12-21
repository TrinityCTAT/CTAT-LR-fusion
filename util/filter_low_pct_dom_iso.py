#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
import pandas as pd
import csv


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


def main():
    
    parser = argparse.ArgumentParser(description="filtering fusion calls based on read counts", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--fusions_input", type=str, required=True, help="fusions input file")
    parser.add_argument("--filtered_fusions_output", type=str, required=True, help="name for filtered fusions output file")
    
    parser.add_argument("--min_frac_dom_iso", type=float, default=0.05, help="min fraction expression of dominant fusion isoform")
    
    args = parser.parse_args()

    fusions_input_filename = args.fusions_input
    fusions_output_filename = args.filtered_fusions_output
    min_frac_dom_iso = args.min_frac_dom_iso

    data = pd.read_csv(fusions_input_filename, sep="\t", quotechar='"')

    if data.shape[0] == 0:
        logger.info("no fusion entries to filter.\n")
        data.to_csv(fusions_output_filename, sep="\t", index=False, quoting=csv.QUOTE_NONE)
        sys.exit(0)
    

    def filter_frac_dom_iso (group_df):
        group_df['max_LR_FFPM'] = group_df['LR_FFPM'].max()
        group_df['frac_dom_iso'] = group_df['LR_FFPM'] / group_df['max_LR_FFPM']
        group_df['above_frac_dom_iso'] =  group_df['frac_dom_iso'] >= min_frac_dom_iso

        return(group_df)
    
    data = data.groupby('#FusionName', group_keys=False).apply(filter_frac_dom_iso).reset_index(drop=True)

    filtered_out_fusions = data[ ~ data['above_frac_dom_iso' ] ]

    retained_fusions = data[ data['above_frac_dom_iso' ] ]  

    # sort descending by long read support
    
    filtered_out_fusions = filtered_out_fusions.sort_values(by=['num_LR'], ascending=False)

    retained_fusions = retained_fusions.sort_values(by=['num_LR'], ascending=False)

    # retain any with Illumina read support
    if "FFPM" in filtered_out_fusions.columns.tolist():
        recovered_FI_fusions = filtered_out_fusions[filtered_out_fusions["FFPM"] > 0]
        num_recovered_fusions = recovered_FI_fusions.shape[0]
        logger.info("-recovering Illumina supported fusions: {}".format(num_recovered_fusions))
        if num_recovered_fusions > 0:
            retained_fusions = pd.concat([retained_fusions, recovered_FI_fusions])
            filtered_out_fusions = filtered_out_fusions[~ (filtered_out_fusions["FFPM"] > 0) ] 

            
    filtered_out_fusions.to_csv(fusions_output_filename + ".removed_below_min_frac_dom_iso", sep="\t", index=False, quoting=csv.QUOTE_NONE) 
    retained_fusions.to_csv(fusions_output_filename, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    logger.info("-filter_low_pct_dom_iso.py removed low dom iso frac fusions: " + str(filtered_out_fusions.shape[0]))
    logger.info("-filter_low_pct_dom_iso.py RETAINED above min dom iso frac fusions: " + str(retained_fusions.shape[0]))
    
    
    sys.exit(0)

    


if __name__=='__main__':
    main()

