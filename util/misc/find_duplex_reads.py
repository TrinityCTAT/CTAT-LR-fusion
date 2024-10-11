#!/usr/bin/env python3

import sys, os, re
import pysam
from ssw import AlignmentMgr
import argparse


def main():

    min_read_seq_length = 1000
    trim_from_ends = 100

    parser = argparse.ArgumentParser(
        description="identify reads that align to their reverse complement sequence and may be duplex type",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--fastq", type=str, required=True, help="fastq file")

    parser.add_argument(
        "--trim_ends_length",
        type=int,
        default=trim_from_ends,
        required=False,
        help="first trim this length from each end of the read sequence",
    )

    parser.add_argument(
        "--min_read_length",
        type=int,
        default=min_read_seq_length,
        required=False,
        help="only analyze reads at least this length",
    )

    parser.add_argument(
        "--show_pretty_alignment",
        action="store_true",
        default=False,
        help="show the pretty alignment",
    )

    args = parser.parse_args()

    fastq_filename = args.fastq
    trim_from_ends = args.trim_ends_length
    min_read_seq_length = args.min_read_length
    show_pretty_alignment_flag = args.show_pretty_alignment

    # For alignment method, see: https://libnano.github.io/ssw-py/quickstart.html#alignment
    # and using alignment scoring as in: https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/blast_score.html

    align_mgr = AlignmentMgr(
        match_score=2,
        mismatch_penalty=3,
    )

    fastq_reader = pysam.FastqFile(fastq_filename)

    for read in fastq_reader:

        accession = read.name
        sequence = read.sequence

        seqlen = len(sequence)

        if seqlen < min_read_seq_length:
            continue

        # trim off the ends
        sequence = sequence[trim_from_ends:]
        sequence = sequence[: -1 * trim_from_ends]

        revseq = reverse_complement(sequence)

        align_mgr.set_read(sequence)
        align_mgr.set_reference(revseq)

        alignment = align_mgr.align(gap_open=5, gap_extension=2)
        # print(alignment)

        revseq_start = alignment.reference_start
        revseq_end = alignment.reference_end

        # revcomp back the match coordinates
        revseq_start, revseq_end = seqlen - revseq_end + 1, seqlen - revseq_start + 1

        seq_start = alignment.read_start
        seq_end = alignment.read_end

        match_len = min(seq_end - seq_start, revseq_end - revseq_start)

        frac_read_self_aligns = match_len / seqlen

        print(
            "\t".join(
                [
                    accession,
                    str(seqlen),
                    "{}-{}".format(seq_start, seq_end),
                    "{}-{}".format(revseq_start, revseq_end),
                    str(match_len),
                    "{:.2f}".format(frac_read_self_aligns),
                ]
            )
        )

        if show_pretty_alignment_flag:
            align_mgr.print_result(alignment)

    sys.exit(0)


def reverse_complement(sequence):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[base] for base in sequence[::-1])


if __name__ == "__main__":
    main()
