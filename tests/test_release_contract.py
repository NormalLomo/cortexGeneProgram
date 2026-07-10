from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ASSET_MAP = ROOT / "final_asset_map.tsv"

EXPECTED_ASSETS = {
    *(f"Fig{number}" for number in range(1, 9)),
    *(f"FigS{number}" for number in range(1, 11)),
    *(f"TableS{number}" for number in range(1, 7)),
    *(f"SuppData{number}" for number in range(1, 7)),
}
FORBIDDEN_FILENAMES = {
    "BUILD_REPORT.md",
    "LICENSE_DECISION_REQUIRED.txt",
    "NUMERIC" + "_20_REVIEW.tsv",
    "AGENT_MANIFEST.md",
}


class ReleaseContractTest(unittest.TestCase):
    def test_final_asset_map_is_complete_and_has_runnable_entrypoints(self) -> None:
        with ASSET_MAP.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["asset_id"] for row in rows}, EXPECTED_ASSETS)
        for row in rows:
            self.assertTrue((ROOT / row["directory"]).is_dir(), row["asset_id"])
            self.assertTrue((ROOT / row["entrypoint"]).is_file(), row["asset_id"])

    def test_final_asset_map_records_real_asset_specific_chains(self) -> None:
        with ASSET_MAP.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            self.assertIn(" -> ", row["source_chain"], row["asset_id"])
            self.assertNotIn("public_asset_entrypoint", row["entrypoint"], row["asset_id"])
            self.assertNotEqual(Path(row["entrypoint"]).name, "run.py", row["asset_id"])

    def test_public_tree_has_no_forbidden_release_files(self) -> None:
        names = {path.name for path in ROOT.rglob("*") if path.is_file()}
        self.assertFalse(FORBIDDEN_FILENAMES & names)

    def test_smoke_uses_the_asset_chain_validator_not_the_removed_wrapper_api(self) -> None:
        smoke = (ROOT / "SMOKE_TESTS.sh").read_text(encoding="utf-8")
        self.assertIn("python run_release.py --validate", smoke)
        self.assertNotIn('run_release.py --asset "$asset" --dry-run', smoke)

    def test_table_s1_builder_writes_only_the_resource_comparison_schema(self) -> None:
        builder = ROOT / "04_tables_and_supplementary_data/TableS1_resource_comparison/source/01_build_resource_comparison.py"
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "TableS1.tsv"
            completed = subprocess.run(
                [sys.executable, str(builder), "--output-tsv", str(output)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 7)
        self.assertEqual(
            list(rows[0]),
            [
                "existing_resource_axis",
                "representative_resources",
                "what_they_enable",
                "remaining_limitation",
                "what_this_study_adds",
                "source_note",
            ],
        )


if __name__ == "__main__":
    unittest.main()
