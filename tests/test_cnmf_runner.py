from __future__ import annotations

import subprocess
import sys
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "01_data_processing/00_human_snrna_discovery/02_run_cnmf_discovery.py"
CONFIG = ROOT / "config/cnmf_discovery.yaml"


class CnmfRunnerContractTest(unittest.TestCase):
    def test_dry_run_plans_the_full_100_replicate_chain(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(RUNNER), "--config", str(CONFIG), "--dry-run"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("cnmf prepare", completed.stdout)
        self.assertIn("--n-iter 100", completed.stdout)
        self.assertIn("cnmf factorize", completed.stdout)
        self.assertIn("cnmf combine", completed.stdout)
        self.assertIn("cnmf consensus", completed.stdout)
        self.assertIn("results/cnmf_human_k60", completed.stdout)


if __name__ == "__main__":
    unittest.main()
