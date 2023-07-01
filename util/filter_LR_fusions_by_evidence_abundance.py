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
    
    parser.add_argument("--min_num_LR", default=1, type=int, help="min number of long reads with canonical splice support")
    parser.add_argument("--min_LR_novel_junction_support", type=int, default=2, help="min number of long reads with non-canonical splice support")
    parser.add_argument("--min_J", type=int, default=1, help="min number of Illumina junction reads with canonical splice breakpoints")
    parser.add_argument("--min_sumJS", type=int, default=1, help="min number of Illumina reads supporting junction and spanning frags summed")
    parser.add_argument("--min_novel_junction_support", type=int, default=3, help="min number of junction reads with non-canonical splice support")
    parser.add_argument("--min_FFPM", type=float, default=0.1, help="min FFPM value for short reads")

    args = parser.parse_args()

    fusions_input_filename = args.fusions_input
    fusions_output_filename = args.filtered_fusions_output

    min_num_LR = args.min_num_LR
    min_LR_novel_junction_support = args.min_LR_novel_junction_support
    min_J = args.min_J
    min_sumJS = args.min_sumJS
    min_novel_junction_support = args.min_novel_junction_support
    min_FFPM = args.min_FFPM

    data = pd.read_csv(fusions_input_filename, sep="\t", quotechar='"')


    if 'JunctionReadCount' in data.columns:
        # filter based on long or short read results:
        data_filtered = data[
            (    # long read criteria

                 ( (data.SpliceType == "ONLY_REF_SPLICE") & (data.num_LR >= min_num_LR) )
                  |
                 (data.num_LR >= min_LR_novel_junction_support)
            )
                |
            (    # short read criteria
              
                (
                    ( (data.JunctionReadCount >= min_J) & (data.SpliceType == "ONLY_REF_SPLICE"))
                            |
                            (data.JunctionReadCount >= min_novel_junction_support)
                )
                    &
                (data.JunctionReadCount + data.SpanningFragCount >= min_sumJS)
                    &
                (data.FFPM >= min_FFPM)
            )
            ]                     

    else:
        # filter just based on long reads
        data_filtered = data[
            
              (  (data.SpliceType == "ONLY_REF_SPLICE") & (data.num_LR >= 1) )
                  |
                (data.num_LR >= 3)
            ]
    
    

    data_filtered.to_csv(fusions_output_filename, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    sys.exit(0)

    


if __name__=='__main__':
    main()

