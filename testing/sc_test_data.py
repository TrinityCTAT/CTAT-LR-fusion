#!/usr/bin/env python3

# Builds synthetic single-cell inputs from the bulk test data by adding fake CB:Z: and UB:Z: tags
# (derived from a hash of the read name, so a read gets the same tags in the fastq and the bam),
# and checks that the resulting fusion predictions carry CB^UMI^readname read accessions.
#
#   sc_test_data.py make_fastq transcripts.fq.gz transcripts.sc.fq.gz
#   sc_test_data.py make_bam   transcripts.mm2.bam transcripts.sc.mm2.bam
#   sc_test_data.py check      ctat_LR_fusion_outdir.sc_fq/ctat-LR-fusion.fusion_predictions.tsv

import sys, os, re
import gzip
import hashlib
import subprocess
import csv

NUM_CELLS = 4


def get_sc_tags(read_name):
    digest = hashlib.md5(read_name.encode()).digest()
    cell_barcode = to_nucs(hashlib.md5(b"cell" + bytes([digest[0] % NUM_CELLS])).digest(), 16)
    umi = to_nucs(digest[1:], 12)
    return cell_barcode, umi


def to_nucs(byte_vals, length):
    return "".join("ACGT"[b % 4] for b in byte_vals)[:length]


def make_fastq(input_fq, output_fq):
    with gzip.open(input_fq, "rt") as fh, gzip.open(output_fq, "wt") as ofh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                read_name = line[1:].split()[0]
                cell_barcode, umi = get_sc_tags(read_name)
                line = "@{}\tCB:Z:{}\tUB:Z:{}\n".format(read_name, cell_barcode, umi)
            ofh.write(line)


def make_bam(input_bam, output_bam):
    sam_reader = subprocess.Popen(["samtools", "view", "-h", input_bam], stdout=subprocess.PIPE, text=True)
    bam_writer = subprocess.Popen(["samtools", "view", "-b", "-o", output_bam, "-"], stdin=subprocess.PIPE, text=True)
    for line in sam_reader.stdout:
        if not line.startswith("@"):
            line = line.rstrip("\n")
            cell_barcode, umi = get_sc_tags(line.split("\t", 1)[0])
            line = "{}\tCB:Z:{}\tUB:Z:{}\n".format(line, cell_barcode, umi)
        bam_writer.stdin.write(line)
    bam_writer.stdin.close()
    if sam_reader.wait() != 0 or bam_writer.wait() != 0:
        raise RuntimeError("Error, samtools failed converting {}".format(input_bam))


def check(fusion_predictions_file):
    num_fusions = 0
    num_accessions = 0
    with open(fusion_predictions_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            num_fusions += 1
            for accession in row["LR_accessions"].split(","):
                num_accessions += 1
                m = re.match(r"^([ACGT]{16})\^([ACGT]{12})\^(.+)$", accession)
                if not m:
                    raise RuntimeError("Error, LR accession not encoded as CB^UMI^readname: {}".format(accession))
                if (m.group(1), m.group(2)) != get_sc_tags(m.group(3)):
                    raise RuntimeError("Error, LR accession has wrong CB^UMI for its read: {}".format(accession))

    if num_fusions == 0:
        raise RuntimeError("Error, no fusions reported in {}".format(fusion_predictions_file))

    print(
        "OK - {} fusions with {} LR accessions, all encoded as CB^UMI^readname".format(num_fusions, num_accessions)
    )


def main():
    usage = "\n\n\tusage: {} (make_fastq in.fq.gz out.fq.gz | make_bam in.bam out.bam | check fusion_predictions.tsv)\n\n".format(
        sys.argv[0]
    )
    if len(sys.argv) < 3:
        exit(usage)

    action = sys.argv[1]
    if action == "make_fastq" and len(sys.argv) == 4:
        make_fastq(sys.argv[2], sys.argv[3])
    elif action == "make_bam" and len(sys.argv) == 4:
        make_bam(sys.argv[2], sys.argv[3])
    elif action == "check" and len(sys.argv) == 3:
        check(sys.argv[2])
    else:
        exit(usage)


if __name__ == "__main__":
    main()
