from __future__ import annotations

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "01_data_processing/00_mouse_snrna_preprocessing/00_fix_subset_mouse.py"


class MousePreprocessingContractTest(unittest.TestCase):
    def test_uses_a_conservative_configurable_dense_block_and_current_successor(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn('parser.add_argument("--block-cells", type=int, default=10_000)', source)
        self.assertIn("01_prepare_mouse_raw_counts.py", source)
        self.assertNotIn("01_prepareData_snRNA.R", source)


if __name__ == "__main__":
    unittest.main()
