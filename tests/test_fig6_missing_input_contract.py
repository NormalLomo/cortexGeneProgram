from __future__ import annotations

import csv
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class Fig6MissingInputContractTest(unittest.TestCase):
    def test_public_docs_define_the_missing_object_without_claiming_reproducibility(self) -> None:
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        sources = (ROOT / "DATA_SOURCES.tsv").read_text(encoding="utf-8")
        for token in (
            "decay_per_program_full.csv",
            "panelA_pairs_long.csv",
            "panelA_medians.csv",
            "SHA256",
            "MISSING_INPUT_AND_PRODUCER",
        ):
            self.assertIn(token, readme + sources)
        self.assertIn("Fig. 6a cannot be reproduced", readme)

    def test_fig6_asset_is_blocked_and_the_unverifiable_prep_script_is_absent(self) -> None:
        with (ROOT / "final_asset_map.tsv").open(newline="", encoding="utf-8") as handle:
            rows = {row["asset_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(rows["Fig6"]["source_status"], "BLOCKED_MISSING_INPUT_PRODUCER")
        self.assertFalse((ROOT / "02_main_figures/Fig6/source/prep_panelA_60.py").exists())


if __name__ == "__main__":
    unittest.main()
