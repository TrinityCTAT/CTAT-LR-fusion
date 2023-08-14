#!/usr/bin/env python3

import sys, os, re
import pysam
import csv
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def main():

    usage = "\n\n\tusage: {} fusions.full.tsv mas-seq.bam\n\n".format(sys.argv[0])
    if len(sys.argv) < 3:
        exit(usage)

    fusions_filename = sys.argv[1]
    bam_filename = sys.argv[2]

    logger.info(f"-parsing {fusions_filename} for long-read fusion support")
    read_to_fusion = parse_fusion_reads(fusions_filename)

    logger.info(f"-parsing {bam_filename} to capture cell barcodes for fusion reads")
    read_name_to_cell_barcode = extract_cell_barcodes_for_fusion_reads(read_to_fusion, bam_filename)

    # report
    print("\t".join(["read_name", "cell_barcode", "fusion_name"]))
    for fusion_read_name, fusion_name in read_to_fusion.items():
        cell_barcode = read_name_to_cell_barcode.get(fusion_read_name, "NA")
        print("\t".join([fusion_read_name, cell_barcode, fusion_name]))

    sys.exit(0)
    
    

def extract_cell_barcodes_for_fusion_reads(read_to_fusion, bam_filename):

    read_name_to_cell_barcode = dict()

    samreader = pysam.AlignmentFile(bam_filename, "rb", check_sq=False)
    for read in samreader:
        read_name = read.query_name
        if read_name in read_to_fusion:
            fusion_name = read_to_fusion[read_name]
            cell_barcode = read.get_tag("CB", "Z")[0]
            read_name_to_cell_barcode[read_name] = cell_barcode
            logger.info(f"{read_name} -> {cell_barcode} -> {fusion_name}")

    return read_name_to_cell_barcode



def parse_fusion_reads(fusions_filename):

    read_to_fusion = dict()

    with open(fusions_filename, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fusion_name = row['#FusionName']
            num_LR = int(row['num_LR'])
            reads = row['LR_accessions']

            if num_LR == 0:
                continue

            LR_accs = reads.split(",")
            for LR_acc in LR_accs:
                read_to_fusion[LR_acc] = fusion_name

    return read_to_fusion


if __name__=='__main__':
    main()
