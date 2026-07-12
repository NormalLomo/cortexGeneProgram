from __future__ import annotations

import csv
import importlib.util
import tempfile
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "02_main_figures/Fig5/source/program_id_contract.py"
FIG5_R = ROOT / "02_main_figures/Fig5/source/figB_large.R"


def load_contract():
    spec = importlib.util.spec_from_file_location("fig5_program_id_contract", CONTRACT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_map(path: Path, excluded: set[int]) -> None:
    next_id = 1
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("old_P", "new_P"), delimiter="\t")
        writer.writeheader()
        for old_id in range(1, 61):
            if old_id in excluded:
                writer.writerow({"old_P": old_id, "new_P": "EXCLUDED"})
            else:
                writer.writerow({"old_P": old_id, "new_P": next_id})
                next_id += 1


class Fig5ProgramContractTest(unittest.TestCase):
    def test_map_accepts_only_the_canonical_54_program_repertoire(self) -> None:
        module = load_contract()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "program_renumber_map.tsv"
            write_map(path, set(module.EXCLUDED_OLD_IDS))
            contract = module.load_program_contract(path)

        self.assertEqual(len(contract.retained_old_ids), 54)
        self.assertEqual(set(contract.old_to_new.values()), set(range(1, 55)))
        self.assertEqual(set(contract.excluded_old_ids), {9, 18, 19, 35, 52, 57})

    def test_map_rejects_a_different_exclusion_set(self) -> None:
        module = load_contract()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "program_renumber_map.tsv"
            write_map(path, {9, 18, 19, 35, 52, 58})
            with self.assertRaisesRegex(ValueError, "excluded old IDs"):
                module.load_program_contract(path)

    def test_fig5a_source_asserts_a_54_by_54_matrix(self) -> None:
        source = FIG5_R.read_text(encoding="utf-8")
        self.assertIn("stopifnot(identical(dim(M), c(54L, 54L)))", source)
        self.assertNotIn("60x60 heatmap", source)


if __name__ == "__main__":
    unittest.main()
