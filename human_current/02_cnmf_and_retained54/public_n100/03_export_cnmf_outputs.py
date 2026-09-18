#!/usr/bin/env python3
"""Export canonical K-specific cNMF outputs for downstream public steps."""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cnmf-run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--components", type=int, required=True)
    parser.add_argument("--local-density-threshold", type=float, required=True)
    return parser.parse_args()


def one_file(directory: Path, patterns: list[str]) -> Path:
    matches = sorted({path for pattern in patterns for path in directory.glob(pattern)})
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one cNMF consensus file matching {patterns!r} in {directory}; found {len(matches)}."
        )
    return matches[0]


def main() -> None:
    args = parse_args()
    threshold = f"{args.local_density_threshold:g}"
    suffixes = [
        f"k_{args.components}.dt_{threshold}",
        f"k_{args.components}.dt_{threshold.replace('.', '_')}",
    ]
    def source_patterns(stem: str, consensus: bool = False) -> list[str]:
        end = ".consensus.txt" if consensus else ".txt"
        return [f"{args.name}.{stem}.{suffix}{end}" for suffix in suffixes]

    sources = {
        "human_k60_gene_spectra_score.tsv": one_file(args.cnmf_run_dir, source_patterns("gene_spectra_score")),
        "human_k60_gene_spectra_tpm.tsv": one_file(args.cnmf_run_dir, source_patterns("gene_spectra_tpm")),
        "human_k60_cell_scores.tsv": one_file(args.cnmf_run_dir, source_patterns("usages", consensus=True)),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for destination_name, source in sources.items():
        destination = args.output_dir / destination_name
        shutil.copyfile(source, destination)


if __name__ == "__main__":
    main()
