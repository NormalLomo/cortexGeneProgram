"""Resolve the portable analysis root used by figure entrypoints."""
from __future__ import annotations

import argparse
import os
from pathlib import Path


ENVIRONMENT_VARIABLE = "CORTEX_PROGRAM_CANONICAL_ROOT"


def add_canonical_root_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--canonical-root",
        help=f"Analysis root. Overrides ${ENVIRONMENT_VARIABLE}.",
    )


def resolve_canonical_root(value: str | None) -> Path:
    raw_root = value or os.environ.get(ENVIRONMENT_VARIABLE)
    if not raw_root:
        raise ValueError(
            f"--canonical-root or {ENVIRONMENT_VARIABLE} is required for canonical figure inputs"
        )
    root = Path(raw_root).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"canonical root is not a directory: {root}")
    return root
