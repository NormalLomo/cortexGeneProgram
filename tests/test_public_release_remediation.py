from __future__ import annotations

import csv
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {9, 18, 19, 35, 52, 57}


def load_workflow_module(name: str):
    path = ROOT / "workflow" / f"{name}.py"
    if not path.is_file():
        raise AssertionError(f"missing workflow module: {path.relative_to(ROOT)}")
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_root_module(name: str):
    path = ROOT / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PublicReleaseRemediationTest(unittest.TestCase):
    def test_canonical_entrypoints_use_a_portable_root_contract(self) -> None:
        entrypoints = (
            "02_main_figures/Fig2/source/07_build_final.sh",
            "02_main_figures/Fig7/source/compose.py",
            "02_main_figures/Fig8/source/06_compose.py",
            "02_main_figures/Fig8/source/07_finalize_labels.py",
        )
        for relative in entrypoints:
            path = ROOT / relative
            self.assertTrue(path.is_file(), relative)
            text = path.read_text(encoding="utf-8")
            self.assertIn("CORTEX_PROGRAM_CANONICAL_ROOT", text, relative)
            self.assertIn("--canonical-root", text, relative)

        for figure_id in ("Fig2", "Fig7", "Fig8"):
            for path in (ROOT / "02_main_figures" / figure_id / "source").rglob("*"):
                if path.is_file() and path.suffix in {".py", ".R", ".sh"}:
                    self.assertNotIn("__PRIVATE_CANONICAL_ROOT__", path.read_text(encoding="utf-8"), path)

        root_contract = load_workflow_module("root_contract")
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(ValueError):
                root_contract.resolve_canonical_root(None)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch.dict(os.environ, {"CORTEX_PROGRAM_CANONICAL_ROOT": str(root)}, clear=True):
                self.assertEqual(root_contract.resolve_canonical_root(None), root.resolve())
            self.assertEqual(root_contract.resolve_canonical_root(str(root)), root.resolve())

    def test_dag_validator_rejects_missing_intermediate_code_nodes(self) -> None:
        release = load_root_module("run_release")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "figures/FigX/source"
            source.mkdir(parents=True)
            for name in ("01_stage.py", "02_transform.R", "03_terminal.py"):
                (source / name).write_text("# fixture\n", encoding="utf-8")
            asset_map = root / "final_asset_map.tsv"
            with asset_map.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("asset_id", "directory", "entrypoint", "source_chain"),
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "asset_id": "FigX",
                        "directory": "figures/FigX",
                        "entrypoint": "figures/FigX/source/03_terminal.py",
                        "source_chain": "source input -> 01_stage.py -> 02_transform.R -> 03_terminal.py",
                    }
                )

            self.assertEqual(release.validate_asset_map(asset_map, root), [])
            (source / "02_transform.R").unlink()
            failures = release.validate_asset_map(asset_map, root)
            self.assertIn("FigX: missing source-chain node: 02_transform.R", failures)

    def test_main_and_supplementary_54_program_schemas_are_equivalent(self) -> None:
        contract = load_workflow_module("program_schema_contract")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            main_path = directory / "main.tsv"
            supplementary_path = directory / "supplementary.tsv"

            main_rows = []
            supplementary_rows = []
            new_id = 1
            for old_id in range(1, 61):
                if old_id in EXCLUDED:
                    main_rows.append({"old_P": f"P{old_id}", "new_P": "EXCLUDED"})
                    continue
                main_rows.append({"old_P": f"P{old_id}", "new_P": f"P{new_id}"})
                supplementary_rows.append({"cnmf_component": old_id, "new_P": f"P{new_id}"})
                new_id += 1

            for path, rows, fields in (
                (main_path, main_rows, ("old_P", "new_P")),
                (supplementary_path, supplementary_rows, ("cnmf_component", "new_P")),
            ):
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                    writer.writeheader()
                    writer.writerows(rows)

            self.assertEqual(
                contract.assert_equivalent_54_program_schemas(main_path, supplementary_path),
                54,
            )

            supplementary_rows[0]["new_P"] = "P54"
            with supplementary_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=("cnmf_component", "new_P"), delimiter="\t")
                writer.writeheader()
                writer.writerows(supplementary_rows)
            with self.assertRaises(ValueError):
                contract.assert_equivalent_54_program_schemas(main_path, supplementary_path)

    def test_release_excludes_internal_history_and_versioned_workers(self) -> None:
        self.assertFalse((ROOT / "AGENT_MANIFEST.md").exists())
        for path in ROOT.rglob("*"):
            if not path.is_file() or path.parts[-2:-1] == ("tests",):
                continue
            self.assertNotIn("hardening_m", path.name, path.relative_to(ROOT))
            self.assertNotRegex(path.name, r"_v\d+\.(?:py|R|sh|tsv)$", path.relative_to(ROOT))

    def test_readme_and_metadata_identify_the_source_release(self) -> None:
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertIn("## Research Background", readme)
        self.assertIn("## Fig. 1A", readme)
        self.assertIn("assets/Fig1A_workflow_overview.png", readme)
        self.assertIn("version `2026.07.12`", readme)
        self.assertIn("10.5281/zenodo.21245200", readme)
        self.assertIn("version: 2026.07.12", citation)
        self.assertIn("date-released: 2026-07-12", citation)
        self.assertEqual(citation.count("family-names:"), 17)

    def test_receipt_claims_match_the_included_evidence(self) -> None:
        status = (ROOT / "PUBLIC_SOURCE_STATUS.md").read_text(encoding="utf-8")
        receipt = (ROOT / "source_receipts/README.md").read_text(encoding="utf-8")
        sources = (ROOT / "DATA_SOURCES.tsv").read_text(encoding="utf-8")
        self.assertIn("provider metadata receipt", status)
        self.assertIn("No recheck log", status)
        self.assertIn("No recheck log", receipt)
        self.assertNotIn("rechecked", status.lower())
        self.assertNotIn("rechecked", receipt.lower())
        self.assertIn("provider metadata receipt", sources)

    def test_release_text_scan_rejects_prohibited_historical_terms(self) -> None:
        scanner = load_workflow_module("release_text_scan")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "producer.py"
            forbidden_term = "reco" + "very"
            source.write_text(f'label = "{forbidden_term}"\n', encoding="utf-8")
            self.assertEqual(scanner.scan_release_text(root), [source])
            source.write_text('label = "enrichment"\n', encoding="utf-8")
            self.assertEqual(scanner.scan_release_text(root), [])

            source.write_text(
                "\n".join(
                    (
                        'note = "current SUBMIT Fig. 8"',
                        'note = "cohort/donor hardening"',
                        'note = "results/hardening/VERDICTS.json"',
                        'note = "bio-MAJ1 fix 2026-06-25"',
                        'note = "FONT UNIFY (W-figfont-unify 2026-06-26)"',
                        'note = "worker P0"',
                    )
                )
                + "\n",
                encoding="utf-8",
            )
            self.assertEqual(scanner.scan_release_text(root), [source])

        self.assertEqual(scanner.scan_release_text(ROOT), [])

    def test_release_text_scan_checks_unmapped_shipped_code(self) -> None:
        scanner = load_workflow_module("release_text_scan")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "mapped.py").write_text('label = "public"\n', encoding="utf-8")
            legacy = root / "unmapped.py"
            legacy.write_text(
                'output = "' + "figure" + "_release/SUBMISSION" + "_final/" + "_final" + "_checks" + '"\n',
                encoding="utf-8",
            )
            (root / "final_asset_map.tsv").write_text(
                "asset_id\tdirectory\tentrypoint\tsource_chain\n"
                "FigX\t.\tmapped.py\tinput -> mapped.py\n",
                encoding="utf-8",
            )

            self.assertEqual(scanner.scan_release_text(root), [legacy])

    def test_retained_program_docstrings_distinguish_raw_and_retained_shapes(self) -> None:
        expected = {
            "03_supplementary_figures/FigS5/source/01_cluster_confusion.py": (
                "raw 60-component inputs",
                "retained 54-program outputs",
            ),
            "03_supplementary_figures/FigS6/source/analysis/xregion_m1_expression_auroc.py": (
                "raw 60-component inputs",
                "retained 54-program analysis",
            ),
            "03_supplementary_figures/FigS8/source/02_by_region.py": (
                "raw 60-component tensor",
                "retained 54 by 54 outputs",
            ),
        }
        for relative, snippets in expected.items():
            text = (ROOT / relative).read_text(encoding="utf-8")
            for snippet in snippets:
                self.assertIn(snippet, text, relative)

    def test_public_main_and_supplementary_schema_contracts_match(self) -> None:
        main = load_root_module("02_main_figures/Fig5/source/program_id_contract")
        expected_exclusions = set(main.EXCLUDED_OLD_IDS)
        self.assertEqual(expected_exclusions, EXCLUDED)
        self.assertEqual(main.EXPECTED_SOURCE_PROGRAMS, 60)
        self.assertEqual(main.EXPECTED_RETAINED_PROGRAMS, 54)

        expected_assertions = (
            "assert len(",
            "list(range(1, 55))",
            "{9, 18, 19, 35, 52, 57}",
        )
        supplementary_sources = (
            "03_supplementary_figures/FigS5/source/01_cluster_confusion.py",
            "03_supplementary_figures/FigS6/source/analysis/xregion_m1_expression_auroc.py",
            "03_supplementary_figures/FigS8/source/02_by_region.py",
        )
        for relative in supplementary_sources:
            source = (ROOT / relative).read_text(encoding="utf-8")
            for assertion in expected_assertions:
                self.assertIn(assertion, source, relative)


if __name__ == "__main__":
    unittest.main()
