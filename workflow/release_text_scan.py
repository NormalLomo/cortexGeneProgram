#!/usr/bin/env python3
"""Reject prohibited historical terminology from publishable release text."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import re


TEXT_SUFFIXES = {".py", ".R", ".r", ".sh", ".md", ".txt", ".tsv", ".yaml", ".yml", ".json", ".cff"}
CODE_SUFFIXES = {".py", ".R", ".r", ".sh"}
FORBIDDEN_TERMS = ("reco" + "very", "reco" + "vered", "恢" + "復", "恢" + "复")
HISTORY_PATTERNS = (
    re.compile(r"\bcurrent\s+submit\b", re.IGNORECASE),
    re.compile(r"\bcandidate\s+(?:fig(?:ure)?\.?\s*\d+|release)\b", re.IGNORECASE),
    re.compile(r"\b(?:hardening|verdicts(?:\.json)?)\b", re.IGNORECASE),
    re.compile(r"\b(?:submission|submit)[_-](?:final|v\d+)\b", re.IGNORECASE),
    re.compile(r"\b" + "figure" + "_" + "release" + r"\b", re.IGNORECASE),
    re.compile(r"(?<![-_])\b(?:worker|" + "p" + "0" + r")\b", re.IGNORECASE),
    re.compile(
        r"\b(?:fix|patch|update|renumber(?:ed)?|rebuild|font[ _-]?unify|attempt\d+)\b.{0,64}\b20\d{2}-\d{2}-\d{2}\b",
        re.IGNORECASE,
    ),
    re.compile(
        r"\b20\d{2}-\d{2}-\d{2}\b.{0,64}\b(?:fix|patch|update|renumber(?:ed)?|rebuild|font[ _-]?unify|attempt\d+)\b",
        re.IGNORECASE,
    ),
    re.compile(r"\b(?:[a-z]+-)?(?:maj\d+|n\d+|w-[a-z0-9_-]+)\s+(?:fix|patch|update)\b", re.IGNORECASE),
)


def source_chain_paths(root: Path) -> set[Path]:
    """Return concrete code files that define the mapped public asset chains."""
    asset_map = root / "final_asset_map.tsv"
    if not asset_map.is_file():
        return {path for path in root.rglob("*") if path.is_file() and path.suffix in CODE_SUFFIXES}

    paths: set[Path] = set()
    with asset_map.open(newline="", encoding="utf-8") as handle:
        for asset in csv.DictReader(handle, delimiter="\t"):
            asset_directory = root / asset["directory"]
            nodes = []
            for branch in asset["source_chain"].split(";"):
                nodes.extend(token.strip() for token in branch.split("->"))
            nodes.append(asset["entrypoint"])
            for node in nodes:
                if Path(node).suffix not in CODE_SUFFIXES:
                    continue
                candidates = (asset_directory / node, asset_directory / "source" / node, root / node)
                path = next((candidate for candidate in candidates if candidate.is_file()), None)
                if path is not None:
                    paths.add(path)
    return paths


def publication_code_paths(root: Path) -> set[Path]:
    """Return shipped executable sources, including files outside mapped chains."""
    scanner = root / "workflow" / "release_text_scan.py"
    return {
        path
        for path in root.rglob("*")
        if path.is_file()
        and path.suffix in CODE_SUFFIXES
        and "tests" not in path.parts
        and path.name != "SMOKE_TESTS.sh"
        and path != scanner
    }


def scan_release_text(root: Path) -> list[Path]:
    matches: set[Path] = set()
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
            continue
        content = path.read_text(encoding="utf-8")
        if any(term.casefold() in content.casefold() for term in FORBIDDEN_TERMS):
            matches.add(path)

    for path in publication_code_paths(root):
        content = path.read_text(encoding="utf-8")
        if any(pattern.search(content) for pattern in HISTORY_PATTERNS):
            matches.add(path)
    return sorted(matches)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    matches = scan_release_text(args.root)
    if matches:
        raise SystemExit("\n".join(str(path.relative_to(args.root)) for path in matches))
    print("Release text scan passed.")


if __name__ == "__main__":
    main()
