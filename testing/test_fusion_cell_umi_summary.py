#!/usr/bin/env python3

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SUMMARY_SCRIPT = REPO_ROOT / "util" / "sc" / "fusion_cell_umi_summary.py"
COLUMNS = [
    "#FusionName",
    "num_LR",
    "LeftGene",
    "LeftBreakpoint",
    "RightGene",
    "RightBreakpoint",
    "SpliceType",
    "LR_accessions",
    "annots",
]


class FusionCellUmiSummaryTest(unittest.TestCase):

    def run_summary(self, rows, *extra_args):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            prefix = Path(tmpdir) / "out"
            with input_path.open("w", newline="") as input_fh:
                writer = csv.DictWriter(input_fh, fieldnames=COLUMNS, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            result = subprocess.run(
                [
                    "python3",
                    str(SUMMARY_SCRIPT),
                    "--fusions",
                    str(input_path),
                    "--output_prefix",
                    str(prefix),
                    *extra_args,
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            outputs = dict()
            for suffix in ("fusion_summary", "cell_fusion_counts", "read_umi_assignments"):
                with open(f"{prefix}.{suffix}.tsv", newline="") as fh:
                    outputs[suffix] = list(csv.DictReader(fh, delimiter="\t"))
            return outputs

    def summary_by_fusion(self, outputs):
        return {row["FusionName"]: row for row in outputs["fusion_summary"]}

    def test_umi_dedup_within_cell(self):
        accs = [
            "CELLA^AAAACCCCGGGG^r1",
            "CELLA^AAAACCCCGGGG^r2",  # exact duplicate
            "CELLA^AAACCCCGGGGT^r3",  # shifted by one base: edit distance 2
            "CELLA^AAAACCCCGGGT^r4",  # substitution: edit distance 1
            "CELLA^TTTTGGGGAAAA^r5",  # distinct molecule
        ]
        outputs = self.run_summary([self.row("A--B", accs)])
        summary = self.summary_by_fusion(outputs)["A--B"]

        self.assertEqual(summary["num_LR"], "5")
        self.assertEqual(summary["num_UMIs_exact"], "4")
        self.assertEqual(summary["num_UMIs"], "2")
        self.assertEqual(summary["num_cells"], "1")

        groups = {row["LR_accession"]: row["umi_group"] for row in outputs["read_umi_assignments"]}
        self.assertEqual(groups["CELLA^AAACCCCGGGGT^r3"], "AAAACCCCGGGG")
        self.assertEqual(groups["CELLA^TTTTGGGGAAAA^r5"], "TTTTGGGGAAAA")

    def test_exact_umi_mode(self):
        accs = ["CELLA^AAAACCCCGGGG^r1", "CELLA^AAAACCCCGGGT^r2"]
        outputs = self.run_summary([self.row("A--B", accs)], "--max_umi_edit_dist", "0")
        self.assertEqual(self.summary_by_fusion(outputs)["A--B"]["num_UMIs"], "2")

    def test_same_umi_different_cells_not_merged(self):
        accs = [
            "CELLA^AAAACCCCGGGG^r1",
            "CELLB^AAAACCCCGGGG^r2",
            "CELLC^AAAACCCCGGGG^r3",
            "CELLC^AAAACCCCGGGG^r4",
        ]
        outputs = self.run_summary([self.row("A--B", accs)])
        summary = self.summary_by_fusion(outputs)["A--B"]

        self.assertEqual(summary["num_UMIs"], "3")
        self.assertEqual(summary["num_cells"], "3")
        self.assertEqual(summary["max_cell_UMIs"], "1")

        cell_counts = {row["cell_barcode"]: row for row in outputs["cell_fusion_counts"]}
        self.assertEqual(cell_counts["CELLC"]["num_LR"], "2")
        self.assertEqual(cell_counts["CELLC"]["num_UMIs"], "1")

    def test_missing_cb_or_umi_excluded(self):
        accs = [
            "CELLA^AAAACCCCGGGG^r1",
            "NA^AAAACCCCGGGG^r2",
            "CELLB^NA^r3",
        ]
        outputs = self.run_summary([self.row("A--B", accs)])
        summary = self.summary_by_fusion(outputs)["A--B"]

        self.assertEqual(summary["num_LR"], "3")
        self.assertEqual(summary["num_LR_lacking_CB_or_UMI"], "2")
        self.assertEqual(summary["num_UMIs"], "1")
        self.assertEqual(summary["num_cells"], "1")

    def test_breakpoints_aggregated_unless_requested(self):
        rows = [
            self.row("A--B", ["CELLA^AAAACCCCGGGG^r1", "CELLB^TTTTGGGGAAAA^r2"], left_brkpt="chr1:100:+"),
            self.row("A--B", ["CELLA^AAAACCCCGGGG^r3", "CELLA^AAAACCCCGGGG^r1"], left_brkpt="chr1:200:+"),
            self.row("C--D", ["CELLA^GGGGAAAACCCC^r4"]),
        ]

        outputs = self.run_summary(rows)
        summary = self.summary_by_fusion(outputs)
        # read r1 is listed under both breakpoints, counted once
        self.assertEqual(summary["A--B"]["num_LR"], "3")
        self.assertEqual(summary["A--B"]["num_UMIs"], "2")
        self.assertEqual(summary["A--B"]["num_cells"], "2")
        # most broadly represented fusion first
        self.assertEqual([row["FusionName"] for row in outputs["fusion_summary"]], ["A--B", "C--D"])

        outputs = self.run_summary(rows, "--by_breakpoint")
        brkpt_rows = [row for row in outputs["fusion_summary"] if row["FusionName"] == "A--B"]
        self.assertEqual(len(brkpt_rows), 2)
        self.assertEqual(sorted(row["num_LR"] for row in brkpt_rows), ["2", "2"])

    def test_requires_cb_umi_encoded_accessions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            with input_path.open("w", newline="") as input_fh:
                writer = csv.DictWriter(input_fh, fieldnames=COLUMNS, delimiter="\t")
                writer.writeheader()
                writer.writerow(self.row("A--B", ["plain_read_name"]))

            result = subprocess.run(
                ["python3", str(SUMMARY_SCRIPT), "--fusions", str(input_path),
                 "--output_prefix", str(Path(tmpdir) / "out")],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("CB^UMI^readname", result.stderr)

    @staticmethod
    def row(fusion_name, accessions, *, left_brkpt="chr1:100:+"):
        left_gene, right_gene = fusion_name.split("--")
        return {
            "#FusionName": fusion_name,
            "num_LR": len(accessions),
            "LeftGene": left_gene,
            "LeftBreakpoint": left_brkpt,
            "RightGene": right_gene,
            "RightBreakpoint": "chr2:500:-",
            "SpliceType": "ONLY_REF_SPLICE",
            "LR_accessions": ",".join(accessions),
            "annots": ".",
        }


if __name__ == "__main__":
    unittest.main()
