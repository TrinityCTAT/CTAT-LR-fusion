#!/usr/bin/env python3

import sys, os, re
import pysam
from collections import defaultdict
import pandas as pd

def main():

    usage = "usage: {} jaffal.simreads.fq.gz\n\n"

    if len(sys.argv) < 2:
        exit(usage)

    fastq_file = sys.argv[1]

    fusion_read_support = defaultdict(int)

    with pysam.FastxFile(fastq_file) as fh:
        for entry in fh:
            fusion_gene_pair = re.search(r"(\S+)\|.*--(\S+)\|.*", entry.comment)
            fusion_gene_pair = "--".join(fusion_gene_pair.groups())
            
            fusion_read_support[fusion_gene_pair] += 1


    df = pd.DataFrame(fusion_read_support.items(), columns=['fusion_name', 'num_reads'])

    df.sort_values('num_reads', ascending=False, inplace=True)

    df.to_csv(sys.stdout, sep="\t", index=False)

if __name__=='__main__':
    main()
