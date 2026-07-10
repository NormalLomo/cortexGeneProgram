from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
VALIDATE_HAN = ROOT / "01_data_processing/02_spatial_preprocessing/03_validate_han_sections.py"
CONVERT = ROOT / "01_data_processing/02_spatial_preprocessing/01_convert_h5ad_to_h5seurat.R"
PREPARE_RCTD = ROOT / "01_data_processing/02_spatial_preprocessing/02_prepare_rctd_objects.R"


class SpatialContractTest(unittest.TestCase):
    def test_han_manifest_rejects_an_empty_section_set(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "han_sections.tsv"
            manifest.write_text("section_id\trelative_file\tsha256\n", encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, str(VALIDATE_HAN), "--manifest", str(manifest), "--input-root", str(root)],
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("at least one", completed.stderr)

    def test_spatial_h5ad_to_rctd_conversion_steps_are_present(self) -> None:
        self.assertTrue(CONVERT.is_file())
        self.assertTrue(PREPARE_RCTD.is_file())
        self.assertIn("Convert(", CONVERT.read_text(encoding="utf-8"))
        self.assertIn("saveRDS", PREPARE_RCTD.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
