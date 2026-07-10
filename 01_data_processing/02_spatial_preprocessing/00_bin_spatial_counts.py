#!/usr/bin/env python3
"""Aggregate raw spatial counts to bin50 and apply a supplied anatomical-domain filter."""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def coordinates(adata: ad.AnnData) -> np.ndarray:
    if "spatial" in adata.obsm:
        values = np.asarray(adata.obsm["spatial"])
    elif {"x", "y"}.issubset(adata.obs.columns):
        values = adata.obs[["x", "y"]].to_numpy()
    else:
        raise KeyError("Input requires adata.obsm['spatial'] or obs columns x and y")
    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError("Spatial coordinates must have at least two columns")
    return values[:, :2].astype(float)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-h5ad", type=Path, required=True)
    parser.add_argument("--domain-tsv", type=Path, required=True, help="TSV with bin_id and keep columns.")
    parser.add_argument("--output-h5ad", type=Path, required=True)
    parser.add_argument("--bin-size", type=float, default=50.0)
    args = parser.parse_args()
    source = ad.read_h5ad(args.raw_h5ad)
    matrix = source.X.tocsr() if sparse.issparse(source.X) else sparse.csr_matrix(source.X)
    if np.any(matrix.data < 0):
        raise ValueError("Spatial counts must be non-negative")
    xy = coordinates(source)
    grid = np.floor(xy / args.bin_size).astype(int)
    bin_id = pd.Index([f"{x}_{y}" for x, y in grid], name="bin_id")
    codes, unique_bins = pd.factorize(bin_id, sort=True)
    aggregator = sparse.csr_matrix(
        (np.ones(source.n_obs, dtype=np.float32), (codes, np.arange(source.n_obs))),
        shape=(len(unique_bins), source.n_obs),
    )
    binned = aggregator @ matrix
    domain = pd.read_csv(args.domain_tsv, sep="\t")
    if not {"bin_id", "keep"}.issubset(domain.columns):
        raise ValueError("Domain table requires bin_id and keep columns")
    domain = domain.drop_duplicates("bin_id").set_index("bin_id")
    obs = pd.DataFrame(index=pd.Index(unique_bins.astype(str), name="bin_id"))
    obs["x"] = [int(value.split("_")[0]) * args.bin_size for value in obs.index]
    obs["y"] = [int(value.split("_")[1]) * args.bin_size for value in obs.index]
    obs = obs.join(domain, how="left")
    keep = obs["keep"].fillna(False).astype(bool).to_numpy()
    result = ad.AnnData(X=binned[keep].tocsr(), obs=obs.iloc[np.flatnonzero(keep)].copy(), var=source.var.copy())
    result.uns["input_contract"] = "raw spatial counts aggregated to bin50 after external anatomical-domain filtering"
    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    result.write(args.output_h5ad)
    print(f"Wrote {args.output_h5ad}: {result.n_obs} retained bin50 observations")


if __name__ == "__main__":
    main()
