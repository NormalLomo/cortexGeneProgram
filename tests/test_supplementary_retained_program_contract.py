from __future__ import annotations

import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class SupplementaryRetainedProgramContractTest(unittest.TestCase):
    def test_suppdata4_title_is_native_and_finalizer_only_reembeds(self) -> None:
        source_dir = ROOT / "04_tables_and_supplementary_data/SuppData4_area_program_matrix/source"
        renderer = (source_dir / "02_area_program_matrix_54.R").read_text(encoding="utf-8")
        finalizer = (source_dir / "03_relabel_reembed_450dpi.py").read_text(encoding="utf-8")

        expected_title = (
            "54 retained cortical gene programs (six cohort-technical components excluded), "
            "ordered by cross-region variability (eta-squared)"
        )
        self.assertIn(expected_title, renderer)
        self.assertNotIn("K=60", renderer)
        self.assertNotIn("add_redact_annot", finalizer)
        self.assertNotIn("insert_text", finalizer)
        self.assertNotIn("source-title", finalizer)
        self.assertIn("get_pixmap", finalizer)

    def test_canonical_scripts_encode_the_retained_program_shapes(self) -> None:
        expected_snippets = {
            "03_supplementary_figures/FigS3/source/01_program_batch_check.py": (
                'assert tab["program"].nunique() == 54',
                "PROGRAM_SUMMARY_ZIP",
            ),
            "03_supplementary_figures/FigS5/source/00_build_released_summary_inputs.py": (
                "assert matrix.shape == (22, 54)",
                "mat_program_subclass.tsv",
            ),
            "03_supplementary_figures/FigS5/source/01_cluster_confusion.py": (
                "assert len(MAP) == 54",
                "cell_program_region_subclass.parquet",
                "spatial_bin50_program_score_SCT.parquet",
            ),
            "03_supplementary_figures/FigS6/source/analysis/xregion_m1_expression_auroc.py": (
                "assert len(mapping) == 54",
                "assert len(perprog) == 54",
                "P-1",
            ),
            "03_supplementary_figures/FigS8/source/02_by_region.py": (
                "assert T.shape == (len(Rg), 10, 54, 54)",
                "assert FP.shape[1] == 1431",
            ),
        }
        for relative_path, snippets in expected_snippets.items():
            text = (ROOT / relative_path).read_text(encoding="utf-8")
            for snippet in snippets:
                self.assertIn(snippet, text, relative_path)

    def test_asset_map_traces_each_recomputed_figure_to_its_canonical_script(self) -> None:
        with (ROOT / "final_asset_map.tsv").open(newline="", encoding="utf-8") as handle:
            rows = {row["asset_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

        expected = {
            "FigS3": ("01_program_batch_check.py",),
            "FigS5": (
                "00_build_released_summary_inputs.py",
                "01_cluster_confusion.py",
                "02_render.py",
            ),
            "FigS6": ("analysis/xregion_m1_expression_auroc.py",),
            "FigS8": ("02_by_region.py",),
        }
        for asset_id, tokens in expected.items():
            for token in tokens:
                self.assertIn(token, rows[asset_id]["source_chain"])
            self.assertEqual(rows[asset_id]["source_status"], "SOURCE_CHAIN_DOCUMENTED")


if __name__ == "__main__":
    unittest.main()
