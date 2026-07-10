#!/usr/bin/env python3
"""Validate the streamed Macosko Isocortex subset and write the cNMF-ready H5AD."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--output-h5ad", type=Path, required=True)
    parser.add_argument("--expected-nuclei", type=int, default=761378)
    args = parser.parse_args()
    source = ad.read_h5ad(args.input_h5ad)
    if source.n_obs != args.expected_nuclei:
        raise ValueError(f"Isocortex input has {source.n_obs} nuclei; expected {args.expected_nuclei}")
    matrix = source.X.tocsr() if sparse.issparse(source.X) else sparse.csr_matrix(source.X)
    if np.any(matrix.data < 0) or not np.allclose(matrix.data, np.rint(matrix.data)):
        raise ValueError("Mouse input must contain non-negative integer raw UMI counts")
    result = ad.AnnData(X=matrix.astype(np.float32), obs=source.obs.copy(), var=source.var.copy())
    result.uns["input_contract"] = "Macosko Isocortex raw integer UMI subset"
    result.uns["source_sha256"] = sha256(args.input_h5ad)
    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    result.write(args.output_h5ad)
    print(f"Wrote {args.output_h5ad}: {result.n_obs} nuclei x {result.n_vars} genes")


if __name__ == "__main__":
    main()
