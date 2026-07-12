#!/usr/bin/env python3
"""Inspect the authoritative, asset-specific source chains in this release."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ASSET_MAP = ROOT / "final_asset_map.tsv"
CODE_SUFFIXES = {".py", ".R", ".r", ".sh"}


def assets(asset_map: Path = ASSET_MAP) -> dict[str, dict[str, str]]:
    with asset_map.open(newline="", encoding="utf-8") as handle:
        return {row["asset_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def source_chain_nodes(source_chain: str) -> list[str]:
    nodes = []
    for branch in source_chain.split(";"):
        for token in branch.split("->"):
            node = token.strip()
            if Path(node).suffix in CODE_SUFFIXES:
                nodes.append(node)
    return nodes


def _resolve_source_chain_node(root: Path, asset: dict[str, str], node: str) -> Path | None:
    asset_directory = root / asset["directory"]
    candidates = (
        asset_directory / node,
        asset_directory / "source" / node,
        root / node,
    )
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def validate_asset_map(asset_map: Path = ASSET_MAP, root: Path = ROOT) -> list[str]:
    failures = []
    for asset in assets(asset_map).values():
        asset_id = asset["asset_id"]
        if not (root / asset["directory"]).is_dir():
            failures.append(f"missing directory: {asset_id}")
        if not (root / asset["entrypoint"]).is_file():
            failures.append(f"missing terminal source: {asset_id}")
        if " -> " not in asset["source_chain"]:
            failures.append(f"incomplete chain: {asset_id}")
            continue
        nodes = source_chain_nodes(asset["source_chain"])
        if not nodes:
            failures.append(f"{asset_id}: source chain has no code nodes")
            continue
        for node in nodes:
            if _resolve_source_chain_node(root, asset, node) is None:
                failures.append(f"{asset_id}: missing source-chain node: {node}")
    return failures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="List final asset IDs and source status.")
    parser.add_argument("--asset", help="Show one asset-specific source chain.")
    parser.add_argument("--validate", action="store_true", help="Validate mapped directories and every code node in each source chain.")
    args = parser.parse_args()
    selected = assets()
    if args.validate:
        failures = validate_asset_map()
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
