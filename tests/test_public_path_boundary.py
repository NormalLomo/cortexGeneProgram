from __future__ import annotations

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {".py", ".R", ".r", ".sh", ".md", ".txt", ".tsv", ".yaml", ".yml", ".json", ".cff"}
FORBIDDEN = (
    re.compile(r"PROJECT_ROOT"),
    re.compile(r"(?<!CORTEX_PROGRAM_)DATA_ROOT"),
    re.compile(r"/Users/|/home/|/mnt/|192\\.168\\.|fsfy|gpuserver|luomeng@"),
)


class PublicPathBoundaryTest(unittest.TestCase):
    def test_release_sources_have_no_private_or_legacy_root_tokens(self) -> None:
        for path in ROOT.rglob("*"):
            if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
                continue
            if path.name == "SMOKE_TESTS.sh" or path.parts[-2:-1] == ("tests",):
                continue
            content = path.read_text(encoding="utf-8")
            for token in FORBIDDEN:
                self.assertIsNone(token.search(content), path.relative_to(ROOT))

    def test_legacy_root_aliases_are_created_explicitly(self) -> None:
        bootstrap = ROOT / "environment/bootstrap_legacy_path_aliases.sh"
        self.assertTrue(bootstrap.is_file())
        source = bootstrap.read_text(encoding="utf-8")
        self.assertIn("CORTEX_PROGRAM_ROOT", source)
        self.assertIn("CORTEX_PROGRAM_DATA_ROOT", source)


if __name__ == "__main__":
    unittest.main()
