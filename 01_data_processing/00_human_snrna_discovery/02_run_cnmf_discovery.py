#!/usr/bin/env python3
"""Run the public 100-replicate human cNMF discovery workflow."""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXPORTER = Path(__file__).with_name("03_export_cnmf_outputs.py")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(os.environ.get("CORTEX_PROGRAM_DATA_ROOT", REPOSITORY_ROOT / "data")),
        help="Root containing the configured input data (default: CORTEX_PROGRAM_DATA_ROOT or ./data).",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path(os.environ.get("CORTEX_PROGRAM_RESULTS_ROOT", REPOSITORY_ROOT / "results")),
        help="Directory for generated results (default: CORTEX_PROGRAM_RESULTS_ROOT or ./results).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the commands and paths without reading data or starting cNMF.",
    )
    return parser.parse_args()


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {
        "counts_h5ad",
        "output_dir",
        "components",
        "selected_components",
        "n_iter",
        "seed",
        "num_highvar_genes",
        "local_density_threshold",
    }
    missing = sorted(required.difference(config))
    if missing:
        raise ValueError(f"Missing cNMF configuration keys: {', '.join(missing)}")
    if config["n_iter"] != 100:
        raise ValueError("Public cNMF discovery uses exactly 100 replicates for every K.")
    if config["selected_components"] not in config["components"]:
        raise ValueError("selected_components must be included in components.")
    return config


def resolve(path_text: str, root: Path) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else root / path


def build_plan(config: dict[str, Any], data_root: Path, results_root: Path) -> list[list[str]]:
    counts = resolve(config["counts_h5ad"], data_root)
    output_dir = resolve(config["output_dir"], results_root)
    name = "human_cnmf"
    components = [str(component) for component in config["components"]]
    selected = str(config["selected_components"])
    threshold = str(config["local_density_threshold"])
    return [
        [
            "cnmf",
            "prepare",
            "--output-dir",
            str(output_dir),
            "--name",
            name,
            "-c",
            str(counts),
            "-k",
            *components,
            "--n-iter",
            "100",
            "--seed",
            str(config["seed"]),
            "--numgenes",
            str(config["num_highvar_genes"]),
        ],
        [
            "cnmf",
            "factorize",
            "--output-dir",
            str(output_dir),
            "--name",
            name,
            "--worker-index",
            "0",
            "--total-workers",
            "1",
        ],
        ["cnmf", "combine", "--output-dir", str(output_dir), "--name", name],
        [
            "cnmf",
            "consensus",
            "--output-dir",
            str(output_dir),
            "--name",
            name,
            "--components",
            selected,
            "--local-density-threshold",
            threshold,
        ],
        [
            sys.executable,
            str(EXPORTER),
            "--cnmf-run-dir",
            str(output_dir / name),
            "--output-dir",
            str(output_dir),
            "--name",
            name,
            "--components",
            selected,
            "--local-density-threshold",
            threshold,
        ],
    ]


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    commands = build_plan(config, args.data_root.resolve(), args.results_root.resolve())
    if args.dry_run:
        print("Validated public cNMF configuration: 100 replicates for every configured K.")
        print("Planned commands:")
        for command in commands:
            print(shlex.join(command))
        return

    counts = resolve(config["counts_h5ad"], args.data_root)
    if not counts.is_file():
        raise FileNotFoundError(f"Configured raw-UMI H5AD does not exist: {counts}")
    for command in commands:
        print("+", shlex.join(command), flush=True)
        subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
