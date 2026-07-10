#!/usr/bin/env python
"""
Fig2 panel i (REPLACEMENT): spatial bin50 UMAP in EXPRESSION space + Harmony.

REPLACES the old program-score UMAP (stage_spatial_umap.py). Instead of
embedding bins by their 60 cNMF program-z, we embed bins by their precomputed
SCT Pearson residuals (per-chip vst, depth-controlled, can be negative) ->
HVG (seurat_current on RAW counts) -> PCA -> Harmony integrate by chip, then
scanpy neighbors + UMAP. This shows that ANATOMY (GM/WM/AR) — not chip batch
— organises the integrated expression manifold.

Pipeline (seed=0; params hardcoded + commented):
  1. Per chip h5ad: stratified subsample by majorDomain, ~150k total across 44
     chips (seed=0), tag chip id. Use layers['sct_residuals'] = Pearson
     residuals (per-chip SCT, depth-controlled) as the working matrix X.
  2. Concat all chips -> HVG ~2000 (flavor='seurat_current' on RAW counts; fallback
     = top-N by residual variance) -> PCA(50) directly on residuals (NO
     sc.pp.scale; residuals are already ~mean-0 per gene).
  3. Harmony integrate by batch=chip via harmonypy.run_harmony DIRECTLY
     (scanpy's sc.external.pp.harmony_integrate triggers a heavy scSLAT import
     path that errored before with harmonypy 2.0.0; we bypass it).
  4. sc.pp.neighbors(use_rep='X_pca_harmony', n_neighbors=30, random_state=0)
     -> sc.tl.umap(min_dist=0.3, random_state=0).
  5. Report kNN chip-entropy mixing BEFORE (on raw PCA) vs AFTER (on Harmony)
     to confirm integration worked.
  6. Save embedding csv (bin, UMAP1, UMAP2, majorDomain, region, chip).

Reproducible: seeds + params hardcoded. Read-only on the h5ad source files.

Run:  python scripts/fig2/stage_spatial_umap_expr.py
"""
import os
# cap BLAS / thread pools BEFORE numpy/sklearn import: the default OpenBLAS
# build on this box overflows its precompiled NUM_THREADS during the
# NearestNeighbors path and can hard-crash. Keep threading bounded + deterministic.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "8")
import glob
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp

import anndata as ad
import scanpy as sc
import harmonypy

# ------------------------------------------------------------------ params
SEED = 0
N_TARGET = 150_000          # target subsample size across all chips
N_HVG = 2000                # HVG via seurat flavor on log1p
N_PCS = 50                  # PCA components
N_NEIGHBORS = 50
UMAP_MIN_DIST = 0.1
UMAP_SPREAD = 0.8
METRIC = "euclidean"
KNN_MIX = 30                # k for chip-entropy mixing metric (before/after)

PROJ = "CORTEX_PROGRAM_ROOT"
BIN50_DIR = "CORTEX_PROGRAM_DATA_ROOT/neuropeptide_cortex/data/human/spatial/bin50"
SCORE_PARQUET = f"{PROJ}/results/crossregion_v1/spatial_bin50_program_score_SCT.parquet"

OUT_CSV = f"{PROJ}/figures/fig2/_intermediate/spatial_umap_expr_harmony.csv"

# anatomy ordering for the subsample-composition report / legend
MAJORDOMAIN_ORDER = ["L1", "L2", "L3", "L4", "L5", "L6", "WM", "ARACHNOID"]

np.random.seed(SEED)


def log(msg):
    print(f"[fig2-umap-expr] {msg}", flush=True)


def knn_chip_entropy(rep, chip_codes, k=KNN_MIX, seed=SEED):
    """Mean kNN chip-entropy, normalised by log(n_chips) -> in [0,1].
    1.0 = neighbours are perfectly chip-mixed (batch removed); near 0 = each
    point's neighbours all share its chip (batch-dominated). Uses scanpy's
    pynndescent path via sc.pp.neighbors on a throwaway AnnData to avoid the
    crashing sklearn OpenBLAS NearestNeighbors path on this box.
    """
    n_chips = int(chip_codes.max()) + 1
    tmp = ad.AnnData(np.zeros((rep.shape[0], 1), dtype=np.float32))
    tmp.obsm["_rep"] = np.ascontiguousarray(rep, dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(tmp, n_neighbors=k, use_rep="_rep",
                        metric=METRIC, random_state=seed)
    G = tmp.obsp["distances"].tocsr()
    indptr, indices = G.indptr, G.indices
    n = G.shape[0]
    logn = np.log(n_chips)
    ent = np.empty(n, dtype=np.float64)
    for i in range(n):
        nbr = indices[indptr[i]:indptr[i + 1]]
        if nbr.size == 0:
            ent[i] = np.nan
            continue
        counts = np.bincount(chip_codes[nbr], minlength=n_chips).astype(np.float64)
        p = counts[counts > 0] / counts.sum()
        ent[i] = -(p * np.log(p)).sum() / logn if logn > 0 else np.nan
    return float(np.nanmean(ent))


def main():
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

    files = sorted(glob.glob(f"{BIN50_DIR}/*_bin50.h5ad"))
    log(f"found {len(files)} chip h5ad files")
    cap = max(1, int(round(N_TARGET / len(files))))   # not used directly; per-chip strat below

    rng = np.random.default_rng(SEED)

    # -------------------------------------------------- per-chip stratified subsample
    # Per chip: stratify by majorDomain. Allocate this chip's quota (N_TARGET /
    # n_chips) proportionally across its majorDomains (so every domain present
    # on a chip is represented), seed=0. Keep RAW counts (layers['counts']).
    per_chip_quota = N_TARGET / len(files)
    sub_list = []
    for f in files:
        chip = os.path.basename(f).replace("_bin50.h5ad", "")
        # OOM-FIX: open backed='r' to avoid materialising full-chip dense
        # sct_residuals (~14 GB/chip). Decide subsample indices from obs FIRST,
        # then slice -> to_memory(), which only realises the kept ~3.4k rows.
        a = ad.read_h5ad(f, backed='r')
        if "sct_residuals" not in a.layers:
            raise RuntimeError(f"{chip}: no layers['sct_residuals']")
        if "counts" not in a.layers:
            raise RuntimeError(f"{chip}: no layers['counts']")
        md = a.obs["majorDomain"].astype(str).values
        # proportional-by-domain allocation of this chip's quota (obs-only)
        domains, dcounts = np.unique(md, return_counts=True)
        idx_take = []
        for dom, dc in zip(domains, dcounts):
            dom_idx = np.where(md == dom)[0]
            frac = dc / len(md)
            take = int(round(per_chip_quota * frac))
            take = min(max(take, 1), len(dom_idx))
            idx_take.append(rng.choice(dom_idx, size=take, replace=False))
        idx_take = np.sort(np.concatenate(idx_take))  # sorted for backed slicing
        n_full = a.n_obs
        # materialise ONLY the kept rows (backed slice -> in-memory AnnData)
        s = a[idx_take].to_memory()
        # retarget X to sct_residuals WITHOUT .copy() (avoids dense duplicate)
        s.X = s.layers["sct_residuals"]
        s.obs["chip"] = chip
        # keep only what we need downstream
        keep = ["chip", "x", "y", "domain", "region", "majorDomain", "bin_total_umi"]
        s.obs = s.obs[[c for c in keep if c in s.obs.columns]].copy()
        s.var = s.var[[]].copy()
        # keep both: X = sct_residuals (working), layers['counts'] = raw (HVG via seurat_current)
        if "sct_residuals" in s.layers:
            del s.layers["sct_residuals"]   # already in X; drop to save memory
        if "sct_log1p" in s.layers:
            del s.layers["sct_log1p"]       # not needed downstream
        sub_list.append(s)
        log(f"  {chip}: subsampled {s.n_obs:,} / {n_full:,} bins")
        try:
            a.file.close()
        except Exception:
            pass
        del a

    log("concatenating subsampled chips ...")
    adata = ad.concat(sub_list, join="inner", index_unique=None, label=None)
    adata.obs_names_make_unique()
    del sub_list
    log(f"concatenated: {adata.n_obs:,} bins x {adata.n_vars:,} genes")

    # store the raw-count bin id (= original obs_names from h5ad) for the join
    adata.obs["bin"] = adata.obs_names.astype(str)

    log("subsample composition (n per majorDomain):")
    comp_md = adata.obs["majorDomain"].astype(str).value_counts().reindex(
        MAJORDOMAIN_ORDER).dropna().astype(int)
    for k, v in comp_md.items():
        log(f"    {k:<10s} {v:>7,d}")
    log(f"n chips = {adata.obs['chip'].nunique()}")

    # -------------------------------------------------- HVG (raw counts) + PCA on residuals
    # X holds SCT Pearson residuals (per-chip, depth-controlled, mean-0 ish,
    # can be negative). residuals are NOT log-counts, so HVG must use a method
    # that works on RAW counts: seurat_current (variance-stabilised). Fallback if
    # seurat_current errors (e.g. zero-row chips): top-N genes by residual variance.
    log(f"HVG selection (flavor='seurat_current' on RAW counts, n_top={N_HVG}) ...")
    try:
        sc.pp.highly_variable_genes(
            adata, flavor="seurat_current", n_top_genes=N_HVG, layer="counts")
        hvg_flavor = "seurat_current"
    except Exception as e:
        log(f"  'seurat_current' HVG failed ({e!r}); falling back to top-N residual variance")
        Xall = adata.X
        if sp.issparse(Xall):
            # residuals were stored as dense earlier; this is just safety
            v = np.asarray(Xall.power(2).mean(axis=0)).ravel() - \
                np.asarray(Xall.mean(axis=0)).ravel() ** 2
        else:
            v = Xall.var(axis=0)
        top = np.argsort(-v)[:N_HVG]
        hv = np.zeros(adata.n_vars, dtype=bool); hv[top] = True
        adata.var["highly_variable"] = hv
        hvg_flavor = "residual_var_topN"
    n_hvg = int(adata.var["highly_variable"].sum())
    log(f"  selected {n_hvg} HVGs (flavor={hvg_flavor})")

    # subset to HVGs (X = sct_residuals restricted to HVGs); NO sc.pp.scale —
    # residuals are already ~mean-0 per gene by construction.
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    Xh = adata_hvg.X
    if sp.issparse(Xh):
        Xh = Xh.toarray()
    Xh = np.asarray(Xh, dtype=np.float32)
    n_bad = int((~np.isfinite(Xh)).sum())
    if n_bad:
        log(f"  sanitised {n_bad} non-finite residual entries -> 0")
        Xh = np.nan_to_num(Xh, nan=0.0, posinf=0.0, neginf=0.0)
    # clip extreme residual tails (parity with sc.pp.scale max_value=10) to
    # keep PCA stable; residual range can hit >200 on rare counts.
    np.clip(Xh, -10.0, 10.0, out=Xh)
    adata_hvg.X = Xh
    log(f"  HVG matrix for PCA: {adata_hvg.shape} (sct_residuals, clipped +-10, no scale)")

    log(f"PCA (n_comps={N_PCS}, seed={SEED}) ...")
    sc.pp.pca(adata_hvg, n_comps=N_PCS, random_state=SEED)
    pca = adata_hvg.obsm["X_pca"]
    adata.obsm["X_pca"] = pca

    # -------------------------------------------------- Harmony integrate by chip
    chip_codes = pd.Categorical(adata.obs["chip"]).codes.astype(np.int64)
    log("kNN chip-entropy mixing BEFORE Harmony (on raw PCA) ...")
    mix_before = knn_chip_entropy(pca, chip_codes)
    log(f"  mix_before (PCA) = {mix_before:.4f}  (1.0 = fully chip-mixed)")

    log("Harmony integrate (harmonypy.run_harmony directly, batch=chip) ...")
    meta = adata.obs[["chip"]].copy()
    ho = harmonypy.run_harmony(pca, meta, vars_use=["chip"], random_state=SEED)
    # harmonypy may return Z_corr as (n_pcs, n_cells); transpose to (n_cells, n_pcs)
    Z = np.asarray(ho.Z_corr).astype(np.float32)
    if Z.shape[0] == adata.n_obs:
        harm = Z
    elif Z.shape[1] == adata.n_obs:
        harm = Z.T  # standard harmonypy convention: transpose
    else:
        raise RuntimeError(f"unexpected Z_corr shape {Z.shape} for n_obs={adata.n_obs}")
    adata.obsm["X_pca_harmony"] = harm
    log(f"  harmony embedding shape = {harm.shape}")

    log("kNN chip-entropy mixing AFTER Harmony ...")
    mix_after = knn_chip_entropy(harm, chip_codes)
    log(f"  mix_after (Harmony) = {mix_after:.4f}  (1.0 = fully chip-mixed)")
    log(f"  MIXING DELTA = {mix_after - mix_before:+.4f}  "
        f"(positive = Harmony improved chip mixing)")

    # -------------------------------------------------- neighbors + UMAP on harmony
    log(f"neighbors(use_rep=X_pca_harmony, n_neighbors={N_NEIGHBORS}, seed={SEED}) ...")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_pca_harmony",
                    metric=METRIC, random_state=SEED)
    log(f"umap(min_dist={UMAP_MIN_DIST}, spread={UMAP_SPREAD}, seed={SEED}) ...")
    sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, random_state=SEED)
    emb = adata.obsm["X_umap"]
    adata.obs["UMAP1"] = emb[:, 0]
    adata.obs["UMAP2"] = emb[:, 1]

    # -------------------------------------------------- save embedding csv
    out = adata.obs[["bin", "UMAP1", "UMAP2", "majorDomain", "region", "chip"]].copy()
    out.to_csv(OUT_CSV, index=False)
    log(f"wrote embedding -> {OUT_CSV}  ({len(out):,} rows)")

    # machine-readable summary line for the report
    log(f"SUMMARY_JSON {{'n_sub': {adata.n_obs}, 'n_chips': "
        f"{adata.obs['chip'].nunique()}, 'n_hvg': {n_hvg}, "
        f"'mix_before': {mix_before:.4f}, 'mix_after': {mix_after:.4f}, "
        f"'mix_delta': {mix_after - mix_before:.4f}}}")
    log("DONE")


if __name__ == "__main__":
    main()
