#!/usr/bin/env python3

# Reads SAM (with header) from stdin, sorted or grouped by read name (eg. 'samtools sort -N'),
# and rewrites each read name as CB^UMI^readname using the CB:Z: and UB:Z: tags,
# matching the read naming from encode_cb_umi_in_read_names.py.
# All records for a read get the same name: tags from the primary record are used
# (as 'samtools fasta' extracts the primary record), falling back to any record of that read carrying them.
# Missing values are set to 'NA'.
# Output SAM is written to stdout.

import sys, os, re
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def main():

    out = sys.stdout

    num_reads = 0
    num_missing_CB = 0
    num_missing_UB = 0

    read_group = []
    prev_read_name = None

    for line in sys.stdin:
        if line.startswith("@"):
            out.write(line)
            continue

        read_name = line.split("\t", 1)[0]
        if read_name != prev_read_name and read_group:
            has_CB, has_UB = write_renamed_read_group(read_group, out)
            num_reads += 1
            num_missing_CB += not has_CB
            num_missing_UB += not has_UB
            read_group = []

        read_group.append(line.rstrip("\n").split("\t"))
        prev_read_name = read_name

    if read_group:
        has_CB, has_UB = write_renamed_read_group(read_group, out)
        num_reads += 1
        num_missing_CB += not has_CB
        num_missing_UB += not has_UB

    logger.info(
        "encoded CB^UMI into read names for {} reads ({} lacking CB, {} lacking UB)".format(
            num_reads, num_missing_CB, num_missing_UB
        )
    )

    sys.exit(0)


def get_tag_val(sam_fields, tag):
    tag_prefix = tag + ":Z:"
    for field in sam_fields[11:]:
        if field.startswith(tag_prefix):
            return field[len(tag_prefix) :]
    return None


def is_primary(sam_fields):
    flag = int(sam_fields[1])
    return (flag & 0x900) == 0


def get_read_group_tag_val(read_group, tag):

    # prefer primary record
    for sam_fields in read_group:
        if is_primary(sam_fields):
            val = get_tag_val(sam_fields, tag)
            if val is not None:
                return val

    # any record
    for sam_fields in read_group:
        val = get_tag_val(sam_fields, tag)
        if val is not None:
            return val

    return None


def write_renamed_read_group(read_group, out):

    cell_barcode = get_read_group_tag_val(read_group, "CB")
    umi = get_read_group_tag_val(read_group, "UB")

    has_CB = cell_barcode is not None
    has_UB = umi is not None

    cell_barcode = re.sub("-1$", "", cell_barcode) if has_CB else "NA"
    umi = umi if has_UB else "NA"

    for sam_fields in read_group:
        sam_fields[0] = "^".join([cell_barcode, umi, sam_fields[0]])
        out.write("\t".join(sam_fields) + "\n")

    return has_CB, has_UB


if __name__ == "__main__":
    main()
