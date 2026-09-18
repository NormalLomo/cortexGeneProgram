"""Gene-overlap input check used by the retained fixed-H scorer."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.sparse as sp

SEED = 20260822


def assert_raw_counts(matrix) -> None:
    values = matrix.data if sp.issparse(matrix) else np.asarray(matrix).ravel()
    if values.size and (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.equal(values, np.floor(values)).all()
    ):
        raise ValueError("input must contain finite nonnegative raw integer counts")


def assert_unique_axis(values, label: str) -> None:
    text = [str(value) for value in values]
    if len(text) != len(set(text)):
        raise ValueError(f"duplicate {label} identifiers")


def validate_overlap(n_genes: int, identity_mass) -> None:
    """Accept the nonempty reference-query gene intersection, as starCAT does."""
    if int(n_genes) <= 0:
        raise ValueError("reference and query have no overlapping genes")


def refuse_existing_output(path) -> None:
    if Path(path).exists():
        raise FileExistsError(f"output already exists: {path}")


def assert_simplex(values, axis: int = 1, atol: float = 1e-6) -> None:
    array = np.asarray(values, dtype=np.float64)
    if not np.isfinite(array).all() or np.min(array) < -atol:
        raise ValueError("simplex values must be finite and nonnegative")
    if not np.allclose(array.sum(axis=axis), 1.0, atol=atol):
        raise ValueError("simplex rows must sum to one")
