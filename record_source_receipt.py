#!/usr/bin/env python3
"""Record source-file size and SHA256 after retrieval from an official provider."""
from __future__ import annotations

import argparse
import csv
import hashlib
from datetime import date
from pathlib import Path


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--file", type=Path, required=True)
    parser.add_argument("--official-url", required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if not args.file.is_file():
        raise FileNotFoundError(args.file)
    write_header = not args.receipt.exists()
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    with args.receipt.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source_id", "file_name", "bytes", "sha256", "retrieved_on", "official_url"],
            delimiter="\t",
        )
        if write_header:
            writer.writeheader()
        writer.writerow(
            {
                "source_id": args.source_id,
                "file_name": args.file.name,
                "bytes": args.file.stat().st_size,
                "sha256": checksum(args.file),
                "retrieved_on": date.today().isoformat(),
                "official_url": args.official_url,
            }
        )


if __name__ == "__main__":
    main()
