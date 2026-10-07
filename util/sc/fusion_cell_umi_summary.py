#!/usr/bin/env python3

# Summarizes single-cell representation of ctat-LR-fusion predictions, with umi deduplication.
#
# Requires the LR_accessions in the ctat-LR-fusion.fusion_predictions.tsv to be named as
#     CB^UMI^readname
# as done by ctat-LR-fusion for single-cell inputs (or by the util/sc/*_to_fastq.py and
# encode_cb_umi_in_*read_names.py utilities).
#
# UMIs are deduplicated within each (fusion, cell). Long-read UMIs frequently differ by
# substitutions and by indel-induced shifts (eg. TCAATCCAATCT vs. CAATCCAATCTT), so UMIs are
# grouped by Levenshtein edit distance: UMIs are visited from most to least abundant (by read count),
# and each is assigned to the first already-established umi group whose representative UMI
# is within --max_umi_edit_dist, otherwise it founds a new group. This is a heuristic; use
# --max_umi_edit_dist 0 for exact-match UMI counting.
#
# Outputs:
#   {prefix}.fusion_summary.tsv         one row per fusion: read, umi, and cell counts
#   {prefix}.cell_fusion_counts.tsv     one row per (fusion, cell): read and umi counts
#   {prefix}.read_umi_assignments.tsv   one row per (fusion, LR_accession): cell, umi, and assigned umi group

import sys
import csv
import argparse
import logging
from collections import defaultdict, Counter

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)

csv.field_size_limit(sys.maxsize)

BREAKPOINT_COLS = ["LeftBreakpoint", "RightBreakpoint", "SpliceType"]


def main():

    parser = argparse.ArgumentParser(
        description="summarize cell and umi-deduplicated support per fusion from single-cell ctat-LR-fusion predictions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fusions", type=str, required=True,
                        help="ctat-LR-fusion.fusion_predictions.tsv (full version, including LR_accessions as CB^UMI^readname)")
    parser.add_argument("--output_prefix", type=str, required=True, help="prefix for output files")
    parser.add_argument("--max_umi_edit_dist", type=int, default=2,
                        help="max Levenshtein distance for grouping UMIs within the same cell and fusion (0 = exact match only)")
    parser.add_argument("--by_breakpoint", action="store_true",
                        help="report per fusion breakpoint (LeftBreakpoint, RightBreakpoint, SpliceType) instead of per fusion gene pair")
    args = parser.parse_args()

    if args.max_umi_edit_dist < 0:
        raise ValueError("--max_umi_edit_dist must be >= 0")

    key_cols = ["#FusionName"] + (BREAKPOINT_COLS if args.by_breakpoint else [])

    fusion_info, fusion_to_reads = parse_fusion_reads(args.fusions, key_cols)

    logger.info("-grouping UMIs for {} fusion entries".format(len(fusion_to_reads)))

    summary_rows = []
    cell_rows = []
    read_rows = []

    for key, reads in fusion_to_reads.items():

        # cell -> umi -> [LR_accessions]
        cell_to_umi_reads = defaultdict(lambda: defaultdict(list))
        num_reads_lacking_cb_umi = 0
        for LR_acc, (cell_barcode, umi) in reads.items():
            if cell_barcode == "NA" or umi == "NA":
                num_reads_lacking_cb_umi += 1
                read_rows.append(list(key) + [cell_barcode, umi, "NA", LR_acc])
                continue
            cell_to_umi_reads[cell_barcode][umi].append(LR_acc)

        num_umis_exact = 0
        num_umis = 0
        max_cell_umis = 0

        for cell_barcode in sorted(cell_to_umi_reads):
            umi_to_reads = cell_to_umi_reads[cell_barcode]
            umi_counts = Counter({umi: len(accs) for umi, accs in umi_to_reads.items()})
            umi_to_group = group_umis(umi_counts, args.max_umi_edit_dist)

            cell_num_reads = sum(umi_counts.values())
            cell_num_umis_exact = len(umi_counts)
            cell_num_umis = len(set(umi_to_group.values()))

            num_umis_exact += cell_num_umis_exact
            num_umis += cell_num_umis
            max_cell_umis = max(max_cell_umis, cell_num_umis)

            cell_rows.append(list(key) + [cell_barcode, cell_num_reads, cell_num_umis_exact, cell_num_umis])

            for umi in sorted(umi_to_reads):
                for LR_acc in sorted(umi_to_reads[umi]):
                    read_rows.append(list(key) + [cell_barcode, umi, umi_to_group[umi], LR_acc])

        num_cells = len(cell_to_umi_reads)
        frac_umis_top_cell = "{:.3f}".format(max_cell_umis / num_umis) if num_umis > 0 else "NA"

        left_gene, right_gene, annots = fusion_info[key]
        summary_rows.append(list(key) + [left_gene, right_gene, len(reads), num_reads_lacking_cb_umi,
                                         num_umis_exact, num_umis, num_cells, max_cell_umis,
                                         frac_umis_top_cell, annots])

    # most broadly represented fusions first
    nkey = len(key_cols)
    summary_rows.sort(key=lambda x: (-x[nkey + 6], -x[nkey + 5], -x[nkey + 2], x[:nkey]))
    cell_rows.sort(key=lambda x: (x[:nkey], -x[nkey + 3], x[nkey]))
    read_rows.sort(key=lambda x: x)

    key_header = ["FusionName"] + key_cols[1:]

    write_tsv(args.output_prefix + ".fusion_summary.tsv",
              key_header + ["LeftGene", "RightGene", "num_LR", "num_LR_lacking_CB_or_UMI",
                            "num_UMIs_exact", "num_UMIs", "num_cells", "max_cell_UMIs",
                            "frac_UMIs_top_cell", "annots"],
              summary_rows)

    write_tsv(args.output_prefix + ".cell_fusion_counts.tsv",
              key_header + ["cell_barcode", "num_LR", "num_UMIs_exact", "num_UMIs"],
              cell_rows)

    write_tsv(args.output_prefix + ".read_umi_assignments.tsv",
              key_header + ["cell_barcode", "umi", "umi_group", "LR_accession"],
              read_rows)

    logger.info("-done. {} fusion entries, {} (fusion, cell) pairs, {} UMIs exact, {} UMIs after grouping".format(
        len(summary_rows), len(cell_rows),
        sum(row[nkey + 4] for row in summary_rows),
        sum(row[nkey + 5] for row in summary_rows)))

    sys.exit(0)


def parse_fusion_reads(fusions_filename, key_cols):

    # key -> (LeftGene, RightGene, annots); annots taken from the breakpoint entry with the most reads
    fusion_info = dict()
    fusion_info_num_LR = dict()

    # key -> LR_accession -> (cell_barcode, umi)
    fusion_to_reads = defaultdict(dict)

    logger.info("-parsing {}".format(fusions_filename))

    with open(fusions_filename, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for col in key_cols + ["LeftGene", "RightGene", "LR_accessions"]:
            if col not in reader.fieldnames:
                raise RuntimeError("Error, missing column {} in {}. Use the full (not abridged) fusion_predictions.tsv".format(
                    col, fusions_filename))

        for row in reader:
            key = tuple(row[col] for col in key_cols)
            accessions = [acc for acc in row["LR_accessions"].split(",") if acc not in ("", ".", "NA")]

            num_LR = len(accessions)
            if key not in fusion_info or num_LR > fusion_info_num_LR[key]:
                fusion_info[key] = (row["LeftGene"], row["RightGene"], row.get("annots", "."))
                fusion_info_num_LR[key] = num_LR

            for acc in accessions:
                parts = acc.split("^", 2)
                if len(parts) != 3:
                    raise RuntimeError(
                        "Error, LR accession {} for {} is not formatted as CB^UMI^readname. ".format(acc, key[0])
                        + "Single-cell barcodes and umis must be encoded in the read names.")
                cell_barcode, umi, read_name = parts
                # reads can be repeated across breakpoint entries of the same fusion; count each read once
                fusion_to_reads[key][acc] = (cell_barcode, umi)

    return fusion_info, fusion_to_reads


def group_umis(umi_counts, max_edit_dist):
    """
    umi_counts: Counter of umi -> num reads (all from one cell and one fusion)
    returns dict: umi -> representative umi of its group
    """

    # most abundant first; ties broken lexically for deterministic output
    ordered_umis = sorted(umi_counts, key=lambda u: (-umi_counts[u], u))

    representatives = []
    umi_to_group = dict()

    for umi in ordered_umis:
        assigned = None
        if max_edit_dist > 0:
            for rep in representatives:
                if levenshtein_within(umi, rep, max_edit_dist):
                    assigned = rep
                    break
        if assigned is None:
            representatives.append(umi)
            assigned = umi
        umi_to_group[umi] = assigned

    return umi_to_group


def levenshtein_within(a, b, max_dist):
    """ returns True if edit distance between a and b is <= max_dist """

    if abs(len(a) - len(b)) > max_dist:
        return False

    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        curr = [i]
        for j, cb in enumerate(b, 1):
            curr.append(min(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + (ca != cb)))
        if min(curr) > max_dist:
            return False
        prev = curr

    return prev[-1] <= max_dist


def write_tsv(filename, header, rows):

    with open(filename, "wt") as ofh:
        ofh.write("\t".join(header) + "\n")
        for row in rows:
            ofh.write("\t".join(str(x) for x in row) + "\n")

    logger.info("-wrote {}".format(filename))


if __name__ == "__main__":
    main()
