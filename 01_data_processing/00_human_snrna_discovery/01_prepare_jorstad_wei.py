#!/usr/bin/env python3
"""Merge Jorstad and Wei raw-UMI inputs into the public human cNMF input."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


EXPECTED_JORSTAD = 814_034
EXPECTED_WEI = 222_005
EXPECTED_JORSTAD_ASSETS = 5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def jorstad_paths_from_manifest(manifest: Path, expected_assets: int = EXPECTED_JORSTAD_ASSETS) -> list[Path]:
    table = pd.read_csv(manifest, sep="\t")
    required = {"dataset_id", "file_name", "bytes", "sha256"}
    if not required.issubset(table.columns):
        raise ValueError(f"Jorstad manifest requires columns: {', '.join(sorted(required))}")
    if len(table) != expected_assets or table["dataset_id"].nunique() != expected_assets:
        raise ValueError(f"Jorstad manifest must contain exactly {expected_assets} distinct dataset assets")
    paths: list[Path] = []
    for row in table.sort_values("dataset_id").itertuples(index=False):
        path = manifest.parent / row.file_name
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.stat().st_size != int(row.bytes):
            raise ValueError(f"Byte-count mismatch for {path}")
        if sha256(path) != row.sha256:
            raise ValueError(f"Checksum mismatch for {path}")
        paths.append(path)
    return paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jorstad-manifest",
        type=Path,
        required=True,
        help="SOURCE_CHECKSUMS.tsv created by 00_download_jorstad_cellxgene.R for all five Supercluster assets.",
    )
    parser.add_argument("--wei-h5ad", type=Path, required=True)
    parser.add_argument("--blacklist", type=Path, required=True)
    parser.add_argument("--output-h5ad", type=Path, required=True)
    parser.add_argument("--raw-layer", default="counts", help="Raw integer UMI layer, or X when explicitly intended.")
    parser.add_argument("--gene-id-column", default="gene_id")
    parser.add_argument("--min-counts-per-nucleus", type=int, default=500)
    parser.add_argument("--expected-jorstad-nuclei", type=int, default=EXPECTED_JORSTAD)
    parser.add_argument("--expected-wei-nuclei", type=int, default=EXPECTED_WEI)
    return parser.parse_args()


def raw_matrix(adata: ad.AnnData, layer: str, source_name: str) -> sparse.csr_matrix:
    matrix = adata.X if layer == "X" else adata.layers.get(layer)
    if matrix is None:
        raise KeyError(f"{source_name} does not provide the requested raw layer: {layer}")
    matrix = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    if np.any(matrix.data < 0) or not np.allclose(matrix.data, np.rint(matrix.data)):
        raise ValueError(f"{source_name} raw layer is not non-negative integer UMI counts")
    return matrix


def prepare_source(
    path: Path,
    source_name: str,
    layer: str,
    gene_id_column: str,
    min_counts: int,
    expected_nuclei: int | None,
) -> tuple[sparse.csr_matrix, pd.DataFrame, pd.Index]:
    source = ad.read_h5ad(path, backed="r")
    if gene_id_column not in source.var:
        raise KeyError(f"{source_name} is missing the required gene ID column: {gene_id_column}")
    genes = pd.Index(source.var[gene_id_column].astype(str), name="gene_id")
    if genes.has_duplicates:
        raise ValueError(f"{source_name} gene IDs are not unique")
    matrix = raw_matrix(source, layer, source_name)
    keep = np.asarray(matrix.sum(axis=1)).ravel() >= min_counts
    matrix = matrix[keep]
    obs = source.obs.iloc[np.flatnonzero(keep)].copy()
    obs["source_cohort"] = source_name
    if expected_nuclei is not None and matrix.shape[0] != expected_nuclei:
        raise ValueError(
            f"{source_name} nuclei after raw-UMI filtering: {matrix.shape[0]}; expected {expected_nuclei}"
        )
    return matrix, obs, genes


def combine_sources(
    paths: list[Path],
    source_name: str,
    layer: str,
    gene_id_column: str,
    min_counts: int,
    expected_nuclei: int,
) -> tuple[sparse.csr_matrix, pd.DataFrame, pd.Index]:
    prepared = [
        prepare_source(path, f"{source_name}:{path.stem}", layer, gene_id_column, min_counts, None)
        for path in paths
    ]
    common_genes = prepared[0][2]
    for _, _, genes in prepared[1:]:
        common_genes = common_genes.intersection(genes, sort=False)
    if common_genes.empty:
        raise ValueError(f"{source_name} assets have no shared gene IDs")
    matrix = sparse.vstack(
        [reindex_columns(matrix, genes, common_genes) for matrix, _, genes in prepared], format="csr"
    )
    obs = pd.concat([obs for _, obs, _ in prepared], axis=0)
    obs["source_cohort"] = source_name
    if matrix.shape[0] != expected_nuclei:
        raise ValueError(
            f"{source_name} nuclei after merge and raw-UMI filtering: {matrix.shape[0]}; expected {expected_nuclei}"
        )
    return matrix, obs, common_genes


def reindex_columns(matrix: sparse.csr_matrix, source_genes: pd.Index, common_genes: pd.Index) -> sparse.csr_matrix:
    positions = source_genes.get_indexer(common_genes)
    if np.any(positions < 0):
        raise ValueError("Gene-ID intersection is inconsistent")
    return matrix[:, positions]


def main() -> None:
    args = parse_args()
    blacklist = {line.strip() for line in args.blacklist.read_text(encoding="utf-8").splitlines() if line.strip()}
    jorstad_paths = jorstad_paths_from_manifest(args.jorstad_manifest)
    j_matrix, j_obs, j_genes = combine_sources(
        jorstad_paths,
        "Jorstad",
        args.raw_layer,
        args.gene_id_column,
        args.min_counts_per_nucleus,
        args.expected_jorstad_nuclei,
    )
    w_matrix, w_obs, w_genes = prepare_source(
        args.wei_h5ad,
        "Wei",
        args.raw_layer,
        args.gene_id_column,
        args.min_counts_per_nucleus,
        args.expected_wei_nuclei,
    )
    common_genes = j_genes.intersection(w_genes, sort=False)
    common_genes = common_genes[~common_genes.isin(blacklist)]
    if common_genes.empty:
        raise ValueError("No shared genes remain after the configured blacklist")
    matrix = sparse.vstack(
        [reindex_columns(j_matrix, j_genes, common_genes), reindex_columns(w_matrix, w_genes, common_genes)],
        format="csr",
    ).astype(np.float32)
    obs = pd.concat([j_obs, w_obs], axis=0)
    result = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=common_genes))
    result.uns["input_contract"] = "raw integer UMI counts; no batch correction, log normalization, or SCT correction"
    result.uns["source_counts"] = {"Jorstad": int(j_matrix.shape[0]), "Wei": int(w_matrix.shape[0])}
    result.uns["source_sha256"] = {
        "Jorstad": {path.name: sha256(path) for path in jorstad_paths},
        "Wei": sha256(args.wei_h5ad),
    }
    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    result.write(args.output_h5ad)
    print(f"Wrote {args.output_h5ad}: {result.n_obs} nuclei x {result.n_vars} shared genes")


if __name__ == "__main__":
    main()
