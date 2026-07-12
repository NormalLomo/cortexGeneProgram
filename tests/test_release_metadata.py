from __future__ import annotations

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ReleaseMetadataTest(unittest.TestCase):
    def test_cff_has_release_identity_and_all_software_creators(self) -> None:
        cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertIn("identifiers:", cff)
        self.assertIn("10.5281/zenodo.21245200", cff)
        self.assertIn("Concept DOI for the collection", cff)
        self.assertIn("version: 2026.07.12", cff)
        self.assertIn("date-released: 2026-07-12", cff)
        self.assertIn("license: MIT", cff)
        self.assertEqual(cff.count("family-names:"), 17)
        self.assertLess(cff.index("family-names: Luo"), cff.index("family-names: Liu"))
        self.assertIn("email: slaket0625@hotmail.com", cff)

    def test_default_conda_environment_uses_the_verified_source_installation_route(self) -> None:
        environment = (ROOT / "environment/environment.yml").read_text(encoding="utf-8")
        requirements = (ROOT / "environment/requirements.txt").read_text(encoding="utf-8")
        self.assertNotIn("r-seuratdisk=", environment)
        self.assertNotIn("r-spacexr=", environment)
        self.assertIn("statsmodels=0.14.2", environment)
        self.assertIn("statsmodels==0.14.2", requirements)
        installer = ROOT / "environment/install_spatial_r_packages.R"
        self.assertTrue(installer.is_file())
        source = installer.read_text(encoding="utf-8")
        self.assertIn("877d4e18ab38c686f5db54f8cd290274ccdbe295", source)
        self.assertIn("9f5dc33c8060f946c6072a138b70e189636e1435", source)

    def test_macosko_bdbag_receipt_is_exact_without_claiming_a_preprocessed_input(self) -> None:
        sources = (ROOT / "DATA_SOURCES.tsv").read_text(encoding="utf-8")
        self.assertIn("nemo:dat-y5zxh0y", sources)
        self.assertIn("DERIVED_INPUT_EDGE_UNRESOLVED", sources)

        fetcher = ROOT / "01_data_processing/00_mouse_snrna_preprocessing/00_fetch_macosko_nemo_bag.sh"
        self.assertTrue(fetcher.is_file())
        source = fetcher.read_text(encoding="utf-8")
        self.assertIn("Raw_data_bag_1_Langlieb_Macosko_WMB_Atlas_2023_Raw_10x.tgz", source)
        self.assertIn("bdbag --resolve-fetch all --validate full", source)
        self.assertIn("does not contain the derived H5AD", source)
        self.assertIn("cell metadata required by 00_fix_subset_mouse.py", source)

    def test_wei_receipt_uses_the_provider_aggregate_without_claiming_downloaded_data(self) -> None:
        sources = (ROOT / "DATA_SOURCES.tsv").read_text(encoding="utf-8")
        self.assertIn("STDF000000000022418", sources)
        self.assertIn("STDF000000000022423", sources)
        self.assertIn("PUBLIC_RECORD_API_VERIFIED_MANUAL_ACCESS_REQUIRED", sources)

        recorder = ROOT / "01_data_processing/00_human_snrna_discovery/00_record_wei_stomics_manifest.py"
        self.assertTrue(recorder.is_file())
        source = recorder.read_text(encoding="utf-8")
        self.assertIn("EXPECTED_PROVIDER_FILES = 69", source)
        self.assertIn("44 spatial and 23 snRNA H5AD files", source)
        self.assertIn("does not claim local materialization", source)


if __name__ == "__main__":
    unittest.main()
