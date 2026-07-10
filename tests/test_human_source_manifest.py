from __future__ import annotations

import csv
import hashlib
import importlib.util
import sys
import tempfile
import types
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
PREPARE = ROOT / "01_data_processing/00_human_snrna_discovery/01_prepare_jorstad_wei.py"
RECORD_WEI = ROOT / "01_data_processing/00_human_snrna_discovery/00_record_wei_stomics_manifest.py"


def load_module():
    spec = importlib.util.spec_from_file_location("prepare_jorstad_wei", PREPARE)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    previous = {name: sys.modules.get(name) for name in ("anndata", "scipy", "scipy.sparse")}
    sys.modules["anndata"] = types.ModuleType("anndata")
    scipy = types.ModuleType("scipy")
    sparse = types.ModuleType("scipy.sparse")
    scipy.sparse = sparse
    sys.modules["scipy"] = scipy
    sys.modules["scipy.sparse"] = sparse
    try:
        spec.loader.exec_module(module)
    finally:
        for name, previous_module in previous.items():
            if previous_module is None:
                del sys.modules[name]
            else:
                sys.modules[name] = previous_module
    return module


def load_wei_recorder():
    spec = importlib.util.spec_from_file_location("record_wei_stomics_manifest", RECORD_WEI)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class HumanSourceManifestTest(unittest.TestCase):
    def test_jorstad_manifest_requires_five_verified_distinct_assets(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = []
            for number in range(5):
                source = root / f"source_{number}.h5ad"
                source.write_bytes(f"fixture-{number}".encode("utf-8"))
                rows.append(
                    {
                        "dataset_id": f"dataset-{number}",
                        "file_name": source.name,
                        "bytes": str(source.stat().st_size),
                        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                    }
                )
            manifest = root / "SOURCE_CHECKSUMS.tsv"
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            self.assertEqual(module.jorstad_paths_from_manifest(manifest), [root / row["file_name"] for row in rows])

    def test_wei_provider_recorder_identifies_the_aggregate_input_and_metadata(self) -> None:
        module = load_wei_recorder()
        spatial_files = [
            {
                "file_id": f"STDF0000000000{number:05d}",
                "file_name": f"section_{number}_cellbin_processed.h5ad",
                "file_path": f"/STSP{number}",
                "file_size": "1MB",
                "data_type": "stomics",
                "file_type": "h5ad",
                "download": False,
            }
            for number in range(22355, 22399)
        ]
        snrna_files = [
            {
                "file_id": f"STDF0000000000{number:05d}",
                "file_name": "snRNA.h5ad" if number == 22418 else f"cell_type_{number}.h5ad",
                "file_path": "/STSP0015456",
                "file_size": "1MB",
                "data_type": "snRNA",
                "file_type": "h5ad",
                "download": False,
            }
            for number in range(22399, 22422)
        ]
        files = spatial_files + snrna_files + [
            {
                "file_id": "STDF000000000022422",
                "file_name": "SnRNA_seurat.RDS",
                "file_path": "/STSP0015456",
                "file_size": "2GB",
                "data_type": "snRNA",
                "file_type": "RDS",
                "download": False,
            },
            {
                "file_id": "STDF000000000022423",
                "file_name": "SnRNA_Meta_with_EdLein_dataset_included.csv",
                "file_path": "/STSP0015456",
                "file_size": "155MB",
                "data_type": "snRNA",
                "file_type": "csv",
                "download": False,
            },
        ]

        rows = module.classify_provider_rows(files)
        self.assertEqual(len(rows), 69)
        aggregate = [row for row in rows if row["role"] == "human_cnmf_input"]
        metadata = [row for row in rows if row["role"] == "human_metadata"]
        self.assertEqual(aggregate, [
            {
                "file_id": "STDF000000000022418",
                "file_name": "snRNA.h5ad",
                "file_path": "/STSP0015456",
                "file_size": "1MB",
                "data_type": "snRNA",
                "provider_download": "false",
                "role": "human_cnmf_input",
            }
        ])
        self.assertEqual(metadata[0]["file_id"], "STDF000000000022423")
        self.assertEqual(metadata[0]["provider_download"], "false")


if __name__ == "__main__":
    unittest.main()
