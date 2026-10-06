#!/usr/bin/env python3

# Rewrites fastq/fasta records that carry single-cell tags in the header comment, eg.
#     @readname<TAB>CB:Z:GTACAACAGGAGAGTA<TAB>UB:Z:AAGCGAAGAGAG
# (as written by 'samtools fastq -T CB,UB')
# so that the cell barcode and umi are encoded in the read name instead:
#     @GTACAACAGGAGAGTA^AAGCGAAGAGAG^readname
# matching the read naming used by 10x_ubam_to_fastq.py and sc-Kinnex_ubam_to_fastq.py
# Records lacking a CB or UB tag get 'NA' in that position.
# Output is written to stdout.

import sys, os, re
import gzip
import argparse
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="encode CB and UB header tags into the read name as CB^UMI^readname",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--reads", type=str, required=True, help="reads in fastq or fasta format (can be gzipped), or '-' for stdin")
    args = parser.parse_args()

    reads_file = args.reads

    if reads_file == "-":
        fh = sys.stdin
    elif reads_file.endswith(".gz"):
        fh = gzip.open(reads_file, "rt")
    else:
        fh = open(reads_file, "rt")

    first_line = fh.readline()
    if first_line == "":
        logger.info("no reads found in {}".format(reads_file))
        sys.exit(0)
    elif first_line.startswith("@"):
        is_fastq = True
    elif first_line.startswith(">"):
        is_fastq = False
    else:
        raise RuntimeError("Error, not recognizing {} as fastq or fasta format".format(reads_file))

    out = sys.stdout

    num_records = 0
    num_missing_CB = 0
    num_missing_UB = 0

    line = first_line
    while line:

        header = line.rstrip("\n")
        new_header, has_CB, has_UB = encode_header(header)
        num_records += 1
        if not has_CB:
            num_missing_CB += 1
        if not has_UB:
            num_missing_UB += 1

        out.write(new_header + "\n")

        if is_fastq:
            # sequence, '+', and quality lines
            for i in range(3):
                rec_line = fh.readline()
                if not rec_line:
                    raise RuntimeError("Error, truncated fastq record for: {}".format(header))
                out.write(rec_line)
            line = fh.readline()
        else:
            # sequence lines up to the next header
            line = fh.readline()
            while line and not line.startswith(">"):
                out.write(line)
                line = fh.readline()

    fh.close()

    logger.info(
        "encoded CB^UMI into read names for {} records ({} lacking CB, {} lacking UB)".format(
            num_records, num_missing_CB, num_missing_UB
        )
    )

    sys.exit(0)


def encode_header(header):

    prefix = header[0]
    fields = header[1:].split()
    read_name = fields[0]

    cell_barcode = "NA"
    umi = "NA"
    for field in fields[1:]:
        if field.startswith("CB:Z:"):
            cell_barcode = re.sub("-1$", "", field[5:])
        elif field.startswith("UB:Z:"):
            umi = field[5:]

    new_header = prefix + "^".join([cell_barcode, umi, read_name])

    return new_header, cell_barcode != "NA", umi != "NA"


if __name__ == "__main__":
    main()
