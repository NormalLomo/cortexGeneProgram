#!/usr/bin/env python3
"""Apply a deterministic torus shift to a two-dimensional derived spatial array."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from stable_seed import stable_seed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Two-dimensional NPY array.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--key", required=True)
    args = parser.parse_args()
    values = np.load(args.input)
    if values.ndim != 2:
        raise ValueError("Input must be two-dimensional")
    generator = np.random.default_rng(stable_seed("torus_shift", args.key))
    shifted = np.roll(values, tuple(generator.integers(0, size, endpoint=False) for size in values.shape), axis=(0, 1))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.save(args.output, shifted)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
