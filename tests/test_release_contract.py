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
        rows_by_id = {row["asset_id"]: row for row in rows}
        for row in rows:
            self.assertIn(" -> ", row["source_chain"], row["asset_id"])
            self.assertNotIn("public_asset_entrypoint", row["entrypoint"], row["asset_id"])
            self.assertNotEqual(Path(row["entrypoint"]).name, "run.py", row["asset_id"])

        self.assertEqual(rows_by_id["FigS9"]["description"], "Disease sensitivity")
        self.assertEqual(rows_by_id["FigS10"]["description"], "Cognition association")
        self.assertEqual(
            rows_by_id["FigS9"]["source_chain"],
            "fig8_analysis.py -> fig8_program_disease.py -> ed_program_disease_supp.py",
        )
        self.assertEqual(
            rows_by_id["FigS10"]["source_chain"],
            "01_program_cognition.py -> 02_render.py",
        )
        cognition_sources = (
            ROOT / "02_main_figures/Fig8/source/fig7_program_cognition.py",
            ROOT / "03_supplementary_figures/FigS10/source/01_program_cognition.py",
            ROOT / "03_supplementary_figures/FigS10/source/02_render.py",
        )
        for source_path in cognition_sources:
            source = source_path.read_text(encoding="utf-8")
            self.assertIn("Fig. S10", source, source_path)
            self.assertIn("figS10_cognition.pdf", source, source_path)
            self.assertNotIn("Fig. S9", source, source_path)
            self.assertNotIn("figS9_cognition.pdf", source, source_path)
            self.assertNotIn("fig7_program_cognition.pdf", source, source_path)
            self.assertNotIn("fig7_program_cognition.png", source, source_path)

    def test_figs9_chain_rejects_aging_and_cognition_nodes(self) -> None:
        with ASSET_MAP.open(newline="", encoding="utf-8") as handle:
            rows = {row["asset_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        chain = rows["FigS9"]["source_chain"].lower()
        self.assertNotIn("aging", chain)
        self.assertNotIn("cognition", chain)
        source_names = {
            path.name.lower()
            for path in (ROOT / "03_supplementary_figures/FigS9/source").iterdir()
            if path.is_file()
        }
        self.assertFalse(any("aging" in name or "cognition" in name for name in source_names))

    def test_figs10_chain_rejects_disease_and_aging_nodes(self) -> None:
        with ASSET_MAP.open(newline="", encoding="utf-8") as handle:
            rows = {row["asset_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        chain = rows["FigS10"]["source_chain"].lower()
        self.assertNotIn("disease", chain)
        self.assertNotIn("aging", chain)
        source_names = {
            path.name.lower()
            for path in (ROOT / "03_supplementary_figures/FigS10/source").iterdir()
            if path.is_file()
        }
        self.assertFalse(any("disease" in name or "aging" in name for name in source_names))

    def test_public_tree_has_no_forbidden_release_files(self) -> None:
        names = {path.name for path in ROOT.rglob("*") if path.is_file()}
        self.assertFalse(FORBIDDEN_FILENAMES & names)

    def test_smoke_uses_the_asset_chain_validator_not_the_removed_wrapper_api(self) -> None:
        smoke = (ROOT / "SMOKE_TESTS.sh").read_text(encoding="utf-8")
        self.assertIn("python run_release.py --validate", smoke)
        self.assertNotIn('run_release.py --asset "$asset" --dry-run', smoke)
        self.assertIn("mktemp -d", smoke)
        self.assertIn("CORTEX_SMOKE_DISPOSABLE=1", smoke)

    def test_release_has_a_prudent_ignore_policy(self) -> None:
        ignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
        for token in ("__pycache__/", "config/paths.env", "/data/", "/results/", "*.h5ad", "*.rds", "*.sqlite*"):
            self.assertIn(token, ignore)

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
