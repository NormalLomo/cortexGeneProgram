#!/usr/bin/env python3
"""Validate the immutable Han section manifest before section-level preprocessing."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t")
    required = {"section_id", "relative_file", "sha256"}
    if not required.issubset(manifest.columns):
        raise ValueError(f"Section manifest requires columns: {', '.join(sorted(required))}")
    if manifest.empty:
        raise ValueError("Han section manifest must declare at least one section")
    if manifest["section_id"].duplicated().any() or manifest["relative_file"].duplicated().any():
        raise ValueError("Han section IDs and relative files must be unique")
    for row in manifest.itertuples(index=False):
        path = args.input_root / row.relative_file
        if not path.is_file():
            raise FileNotFoundError(path)
        if sha256(path) != row.sha256:
            raise ValueError(f"Checksum mismatch for {path}")
    print(f"Validated {len(manifest)} Han sections.")


if __name__ == "__main__":
    main()
