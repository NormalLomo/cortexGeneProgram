#!/usr/bin/env python3
"""Validate exact Chen snRNA files obtained through the provider-controlled route."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


EXPECTED = {
    "snRNA.metadata.2monkeys.rds": (31_466_747, "04e33de308496603119101f88935d493d1c78682078721499b0118b5c0f4dc1c"),
    "snRNA.sparseMatrix_Monkey1.counts.rds": (4_129_681_166, "23c00e17a1b20b1731fb522e4f836819ff6f88315ee52b582a2124d50a7a9db9"),
    "snRNA.sparseMatrix_Monkey2.counts.rds": (3_544_696_257, "57aed275706d1443002758edc6ed323c995e0d38d6c00464cee31aa5e1247f58"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    args = parser.parse_args()
    for name, (expected_size, expected_hash) in EXPECTED.items():
        path = args.input_root / name
        if not path.is_file():
            raise SystemExit(f"Missing provider file: {path}")
        if path.stat().st_size != expected_size:
            raise SystemExit(f"Unexpected byte count for {name}: {path.stat().st_size}")
        if sha256(path) != expected_hash:
            raise SystemExit(f"SHA256 mismatch for {name}")
    print("Validated 3 Chen provider files.")


if __name__ == "__main__":
    main()
