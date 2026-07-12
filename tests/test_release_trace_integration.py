from __future__ import annotations

import csv
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ReleaseTraceIntegrationTest(unittest.TestCase):
    def test_final_suppdata_rows_have_documented_finalization_chains(self) -> None:
        with (ROOT / "final_asset_map.tsv").open(newline="", encoding="utf-8") as handle:
            rows = {row["asset_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

        expected = {
            "SuppData1": "02_finalize_numbering.py",
            "SuppData2": "03_finalize_numbering.py",
            "SuppData4": "03_relabel_reembed_450dpi.py",
            "SuppData6": "03_finalize_numbering.py",
        }
        for asset_id, finalizer in expected.items():
            row = rows[asset_id]
            self.assertEqual(row["source_status"], "SOURCE_CHAIN_DOCUMENTED")
            self.assertIn(finalizer, row["source_chain"])
            self.assertTrue((ROOT / row["directory"] / "source" / finalizer).is_file())
            self.assertNotEqual(row["default_output"], "MISSING")

    def test_public_source_contracts_distinguish_manual_access_from_missing_code(self) -> None:
        sources = (ROOT / "DATA_SOURCES.tsv").read_text(encoding="utf-8")
        status = (ROOT / "PUBLIC_SOURCE_STATUS.md").read_text(encoding="utf-8")
        for token in (
            "manual_provider_access",
            "STDF000000000022418",
            "nemo:dat-y5zxh0y",
            "CNP0003837",
            "CNP0002035",
            "1663381185152036865",
        ):
            self.assertIn(token, sources + status)
        self.assertIn("manual provider access", status.lower())
        self.assertNotIn("BLOCKED_PROVIDER_DOWNLOAD", sources + status)


if __name__ == "__main__":
    unittest.main()
