#!/usr/bin/env python3
"""Inspect the authoritative, asset-specific source chains in this release."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ASSET_MAP = ROOT / "final_asset_map.tsv"


def assets() -> dict[str, dict[str, str]]:
    with ASSET_MAP.open(newline="", encoding="utf-8") as handle:
        return {row["asset_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="List final asset IDs and source status.")
    parser.add_argument("--asset", help="Show one asset-specific source chain.")
    parser.add_argument("--validate", action="store_true", help="Validate mapped directories and terminal source files.")
    args = parser.parse_args()
    selected = assets()
    if args.validate:
        failures = []
        for asset in selected.values():
            if not (ROOT / asset["directory"]).is_dir():
                failures.append(f"missing directory: {asset['asset_id']}")
            if not (ROOT / asset["entrypoint"]).is_file():
                failures.append(f"missing terminal source: {asset['asset_id']}")
            if " -> " not in asset["source_chain"]:
                failures.append(f"incomplete chain: {asset['asset_id']}")
        if failures:
            raise SystemExit("\n".join(failures))
        print(f"Validated {len(selected)} mapped asset-specific source chains.")
    elif args.list:
        for asset in selected.values():
            print(f"{asset['asset_id']}\t{asset['source_status']}\t{asset['description']}")
    elif args.asset:
        if args.asset not in selected:
            raise SystemExit(f"Unknown asset: {args.asset}")
        asset = selected[args.asset]
        for key in ("asset_id", "description", "source_status", "entrypoint", "source_chain", "output_evidence", "default_output"):
            print(f"{key}: {asset[key]}")
    else:
        parser.error("Specify --list, --asset, or --validate.")


if __name__ == "__main__":
    main()
