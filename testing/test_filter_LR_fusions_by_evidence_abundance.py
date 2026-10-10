#!/usr/bin/env python3

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
FILTER_SCRIPT = REPO_ROOT / "util" / "filter_LR_fusions_by_evidence_abundance.py"
LONG_COLUMNS = [
    "marker",
    "#FusionName",
    "num_LR",
    "SpliceType",
    "LR_accessions",
    "LR_FFPM",
]
COMBINED_COLUMNS = LONG_COLUMNS + ["JunctionReadCount", "SpanningFragCount", "FFPM"]


class EvidenceAbundanceFilterTest(unittest.TestCase):

    def run_filter(self, rows, columns=LONG_COLUMNS, *, num_lr_total=300_000_000,
                   min_frac_dom_iso=0.05):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_path = Path(tmpdir) / "output.tsv"
            with input_path.open("w", newline="") as input_fh:
                writer = csv.DictWriter(input_fh, fieldnames=columns, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            result = subprocess.run(
                [
                    "python3",
                    str(FILTER_SCRIPT),
                    "--fusions_input",
                    str(input_path),
                    "--filtered_fusions_output",
                    str(output_path),
                    "--num_LR_total",
                    str(num_lr_total),
                    "--min_frac_dom_iso",
                    str(min_frac_dom_iso),
                    "--min_FFPM",
                    "0.1",
                    "--min_num_LR",
                    "1",
                    "--min_LR_novel_junction_support",
                    "2",
                    "--min_J",
                    "1",
                    "--min_sumJS",
                    "1",
                    "--min_novel_junction_support",
                    "3",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            with output_path.open(newline="") as output_fh:
                reader = csv.DictReader(output_fh, delimiter="\t")
                output_rows = list(reader)
                self.assertEqual(reader.fieldnames, columns)
                return output_rows

    def test_sums_lr_support_across_breakpoint_isoforms_only_within_pair(self):
        rows = [
            self.long_row("pair-a", "PAIR", 18, "read-a"),
            self.long_row("pair-b", "PAIR", 18, "read-b"),
            self.long_row("different-a", "DIFFERENT-A", 18, "read-c"),
            self.long_row("different-b", "DIFFERENT-B", 18, "read-d"),
        ]

        output = self.run_filter(rows)

        self.assertEqual([row["marker"] for row in output], ["pair-a", "pair-b"])

    def test_excludes_invalid_and_subdominant_isoforms_from_pair_total(self):
        rows = [
            self.long_row("invalid-valid", "INVALID", 29, "read-a"),
            self.long_row(
                "invalid-novel",
                "INVALID",
                1,
                "read-b",
                splice_type="INCL_NON_REF_SPLICE",
            ),
            self.long_row("dominant", "DOMINANCE", 29, "read-c"),
            self.long_row("subdominant", "DOMINANCE", 1, "read-d"),
        ]

        self.assertEqual(self.run_filter(rows), [])
        dominance_disabled = self.run_filter(rows[2:], min_frac_dom_iso=0)
        self.assertEqual(
            [row["marker"] for row in dominance_disabled],
            ["dominant", "subdominant"],
        )

    def test_preserves_row_level_short_read_filtering(self):
        rows = [
            self.combined_row("mixed-lr-a", "MIXED", 18, "read-a"),
            self.combined_row("mixed-lr-b", "MIXED", 18, "read-b"),
            self.combined_row(
                "mixed-sr-low", "MIXED", "", "", junction=1, spanning=20, ffpm=0.06
            ),
            self.combined_row(
                "sr-split-a", "SR-SPLIT", "", "", junction=1, spanning=1, ffpm=0.06
            ),
            self.combined_row(
                "sr-split-b", "SR-SPLIT", "", "", junction=1, spanning=1, ffpm=0.06
            ),
            self.combined_row(
                "sr-exact", "SR-EXACT", "", "", junction=1, spanning=1, ffpm=0.1
            ),
            self.combined_row(
                "sr-rescues-lr", "SR-RESCUE", 1, "read-c", junction=0, spanning=0, ffpm=0.2
            ),
        ]

        output = self.run_filter(rows, COMBINED_COLUMNS)

        self.assertEqual(
            [row["marker"] for row in output],
            ["mixed-lr-a", "mixed-lr-b", "sr-exact", "sr-rescues-lr"],
        )

    def test_exact_lr_pair_threshold_and_empty_input(self):
        rows = [
            self.long_row("exact-a", "EXACT", 15, "read-a"),
            self.long_row("exact-b", "EXACT", 15, "read-b"),
            self.long_row("below-a", "BELOW", 14, "read-c"),
            self.long_row("below-b", "BELOW", 15, "read-d"),
        ]

        output = self.run_filter(rows)

        self.assertEqual([row["marker"] for row in output], ["exact-a", "exact-b"])
        self.assertEqual(self.run_filter([]), [])

    @staticmethod
    def long_row(marker, fusion_name, num_lr, accession, *,
                 splice_type="ONLY_REF_SPLICE"):
        lr_ffpm = "" if num_lr == "" else float(num_lr) / 300
        return {
            "marker": marker,
            "#FusionName": fusion_name,
            "num_LR": num_lr,
            "SpliceType": splice_type,
            "LR_accessions": accession,
            "LR_FFPM": lr_ffpm,
        }

    @classmethod
    def combined_row(cls, marker, fusion_name, num_lr, accession, *, junction=0,
                     spanning=0, ffpm=0):
        row = cls.long_row(marker, fusion_name, num_lr, accession)
        row.update(
            {
                "JunctionReadCount": junction,
                "SpanningFragCount": spanning,
                "FFPM": ffpm,
            }
        )
        return row


if __name__ == "__main__":
    unittest.main()
