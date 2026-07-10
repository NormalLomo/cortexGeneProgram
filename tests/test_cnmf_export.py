from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
EXPORTER = ROOT / "01_data_processing/00_human_snrna_discovery/03_export_cnmf_outputs.py"


class CnmfExportContractTest(unittest.TestCase):
    def test_exports_cnmf_171_native_density_threshold_filenames(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run_dir = root / "cnmf_human_k60" / "human_cnmf"
            output_dir = root / "cnmf_human_k60"
            run_dir.mkdir(parents=True)
            suffix = "k_60.dt_0_15"
            for stem in ("gene_spectra_score", "gene_spectra_tpm", "usages"):
                consensus = ".consensus" if stem == "usages" else ""
                (run_dir / f"human_cnmf.{stem}.{suffix}{consensus}.txt").write_text(
                    "feature\t1\nexample\t1.0\n", encoding="utf-8"
                )

            completed = subprocess.run(
                [
                    sys.executable,
                    str(EXPORTER),
                    "--cnmf-run-dir",
                    str(run_dir),
                    "--output-dir",
                    str(output_dir),
                    "--name",
                    "human_cnmf",
                    "--components",
                    "60",
                    "--local-density-threshold",
                    "0.15",
                ],
                check=False,
                capture_output=True,
                text=True,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue((output_dir / "human_k60_gene_spectra_score.tsv").is_file())
            manifest = json.loads((output_dir / "downstream_outputs.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["local_density_threshold"], 0.15)
            self.assertEqual(manifest["components"], 60)


if __name__ == "__main__":
    unittest.main()
