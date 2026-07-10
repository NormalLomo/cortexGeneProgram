#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figure-1 UMAP recomputation using SLAT's Dual PCA (DPCA) to integrate the two
snRNA cohorts (batch = 'us' vs 'edlein') on the cNMF K=60 per-cell USAGE matrix.

This REPLACES the prior Harmony-based UMAP.

WHY DPCA here: SLAT's dual_pca (scSLAT.model.batch.dual_pca) is its feature-embedding
/ batch-correction step. It takes two feature matrices X (cells_X x F) and Y (cells_Y x F)
in the SAME feature space, forms the cross-covariance C = X @ Y.T, and takes its SVD
(randomized_svd, sklearn backend). The left/right singular vectors give per-dataset
embeddings (embd_X, embd_Y) in a shared latent space. This is the CCA / dual-PCA
construction (ref: Xin-Ming Tu blog, cited in SLAT source). It does NOT require spatial
coordinates -- the spatial graph is only used by the *downstream* GNN alignment, not by
the DPCA feature step. On our 60-dim cNMF-usage input this is effectively a CCA-style
cross-cohort alignment in the shared 60-program feature space.

SLAT's DPCA convention (from scSLAT/model/loaddata.py, feature=='dpca' branch):
  - select HVG / normalize_total / log1p   <-- gene-expression specific, SKIPPED here
    because our input is already the cNMF K=60 per-cell USAGE matrix (a derived,
    non-count feature space), not raw gene counts. The 60 programs ARE the shared
    feature space, analogous to the post-HVG feature block SLAT feeds to dual_pca.
  - sc.pp.scale(adata_1); sc.pp.scale(adata_2)   <-- z-score each cohort independently. KEPT.
  - dual_pca(adata_1.X, adata_2.X, dim=50, singular=True, backend='sklearn')   <-- KEPT.
  - SLAT also requires len(adata_1) >= len(adata_2) (check_order). We honor this by
    putting the LARGER cohort (edlein, 117648) as X and us (32352) as Y, then
    reassembling both embeddings back to the original 150k barcode order.

Seed = 0 everywhere. Run on the SAME 150k cells used by the Harmony UMAP (for comparability).

scSLAT 0.3.0. dual_pca imported in ISOLATION (bypassing scSLAT/__init__.py which pulls
heavy/incompatible deps faiss/torch_geometric) -- batch.py only needs numpy/torch/sklearn
and scSLAT.utils.get_free_gpu (pynvml). use_gpu=False: the sklearn backend computes the
SVD on CPU regardless; the only GPU op would be the X@Y.T matmul, and this torch build
(cu117, sm<=86) is not built for the H100 (sm_90), so we keep the matmul on CPU for
determinism/safety.
"""

import importlib.util
import sys
import types

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 0
np.random.seed(SEED)

PROJ = "CORTEX_PROGRAM_ROOT"
UMAP_CSV = f"{PROJ}/figures/fig1/_intermediate/umap_embedding.csv"
CNMF_SCORES = f"{PROJ}/results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_cell_scores.tsv"
OBS_CSV = f"{PROJ}/inputs/snRNA_1M_obs.csv"
DPCA_DIM = 50  # SLAT default; must be <= 60 (input feature dim). cor_var rank <= 60.

# ---------------------------------------------------------------------------
# Isolated import of SLAT's dual_pca (avoid heavy package __init__)
# ---------------------------------------------------------------------------
SP = "PYTHON_SITE_PACKAGES"
_pkg = types.ModuleType("scSLAT"); _pkg.__path__ = [SP + "/scSLAT"]; sys.modules["scSLAT"] = _pkg
_mod = types.ModuleType("scSLAT.model"); _mod.__path__ = [SP + "/scSLAT/model"]; sys.modules["scSLAT.model"] = _mod
_spec_u = importlib.util.spec_from_file_location("scSLAT.utils", SP + "/scSLAT/utils.py")
_um = importlib.util.module_from_spec(_spec_u); sys.modules["scSLAT.utils"] = _um; _spec_u.loader.exec_module(_um)
_spec_b = importlib.util.spec_from_file_location("scSLAT.model.batch", SP + "/scSLAT/model/batch.py")
_bm = importlib.util.module_from_spec(_spec_b); sys.modules["scSLAT.model.batch"] = _bm; _spec_b.loader.exec_module(_bm)
dual_pca = _bm.dual_pca

import scSLAT  # noqa  (works now that submodules pre-stubbed? no -> just print version from metadata)
try:
    from importlib.metadata import version as _ver
    print("scSLAT version:", _ver("scSLAT"))
except Exception as e:
    print("scSLAT version lookup failed:", e)

# ---------------------------------------------------------------------------
# STEP 1: inputs -- reuse SAME 150k cells
# ---------------------------------------------------------------------------
print("\n=== STEP 1: load inputs ===")
umap_old = pd.read_csv(UMAP_CSV, index_col=0)
barcodes = umap_old.index  # preserve order
print("Reference (Harmony) UMAP CSV:", umap_old.shape, "cols:", list(umap_old.columns))
assert umap_old.shape[0] == 150000, "expected 150k cells"

# cNMF K=60 per-cell usages, subset to the 150k barcodes (preserve order)
scores = pd.read_csv(CNMF_SCORES, sep="\t", index_col=0)
scores.columns = [str(c) for c in scores.columns]
print("cNMF scores full:", scores.shape)
assert barcodes.isin(scores.index).all(), "not all 150k barcodes in cNMF scores"
usage = scores.loc[barcodes]  # 150000 x 60, in barcode order
assert (usage.index == barcodes).all(), "barcode order mismatch after subset"
assert usage.shape[1] == 60, f"expected 60 programs, got {usage.shape[1]}"
print("Subset usage matrix:", usage.shape, "(100% barcode match, order preserved)")

# batch per cell from obs.csv (authoritative per spec)
obs = pd.read_csv(OBS_CSV, index_col=0)
assert barcodes.isin(obs.index).all(), "not all barcodes in obs.csv"
batch = obs.loc[barcodes, "batch"].astype(str)
n_us = int((batch == "us").sum()); n_edlein = int((batch == "edlein").sum())
print(f"batch from obs.csv -> us={n_us}  edlein={n_edlein}  total={n_us + n_edlein}")
# cross-check vs the batch column already carried in the Harmony CSV
agree = (batch.values == umap_old["batch"].astype(str).values).mean()
print(f"batch agreement vs existing CSV batch col: {agree:.6f}")

# ---------------------------------------------------------------------------
# STEP 2: SLAT Dual PCA integration
# ---------------------------------------------------------------------------
print("\n=== STEP 2: SLAT dual_pca ===")
# Larger cohort first (SLAT check_order: adata1 >= adata2). edlein(117648) > us(32352).
mask_edlein = (batch == "edlein").values
mask_us = (batch == "us").values
X_edlein = usage.values[mask_edlein].astype("float32")  # (117648, 60)
Y_us = usage.values[mask_us].astype("float32")          # (32352, 60)
print("X (edlein):", X_edlein.shape, " Y (us):", Y_us.shape)

# SLAT convention: z-score each cohort independently (sc.pp.scale) before dual_pca.
ax = sc.AnnData(X_edlein.copy()); ay = sc.AnnData(Y_us.copy())
sc.pp.scale(ax)  # per-program z-score within edlein
sc.pp.scale(ay)  # per-program z-score within us
Xs = np.asarray(ax.X, dtype="float32")
Ys = np.asarray(ay.X, dtype="float32")
print("scaled X mean/std:", float(Xs.mean()), float(Xs.std()),
      " scaled Y mean/std:", float(Ys.mean()), float(Ys.std()))

# dual PCA -> shared embedding. singular=True (SLAT default in load_anndatas),
# backend='sklearn' (randomized_svd, random_state=0 inside dual_pca), use_gpu=False.
Z_x, Z_y = dual_pca(Xs, Ys, dim=DPCA_DIM, singular=True, backend="sklearn", use_gpu=False)
Z_x = np.asarray(Z_x); Z_y = np.asarray(Z_y)
print("DPCA embeddings -> Z_edlein:", Z_x.shape, " Z_us:", Z_y.shape,
      f" (dim={DPCA_DIM}, singular=True, backend=sklearn, use_gpu=False)")

# reassemble into one (150k x dim) matrix in original barcode order
dpca = np.zeros((usage.shape[0], DPCA_DIM), dtype="float32")
dpca[mask_edlein] = Z_x
dpca[mask_us] = Z_y
print("Combined DPCA embedding:", dpca.shape)

# ---------------------------------------------------------------------------
# STEP 3: UMAP
# ---------------------------------------------------------------------------
print("\n=== STEP 3: UMAP ===")
ad = sc.AnnData(np.zeros((usage.shape[0], 1), dtype="float32"))
ad.obs_names = list(barcodes)
ad.obsm["X_dpca"] = dpca
sc.pp.neighbors(ad, use_rep="X_dpca", n_neighbors=30, metric="euclidean", random_state=SEED)
sc.tl.umap(ad, min_dist=0.3, random_state=SEED)
emb = ad.obsm["X_umap"]
print("UMAP coords:", emb.shape,
      "  params: n_neighbors=30, metric=euclidean, min_dist=0.3, random_state=0")

# ---------------------------------------------------------------------------
# STEP 4: write outputs (handled partly by the orchestrator: backup before this)
# ---------------------------------------------------------------------------
print("\n=== STEP 4: write new umap_embedding.csv ===")
out = pd.DataFrame(index=barcodes)
out.index.name = umap_old.index.name  # keep same (unnamed) index header
out["UMAP1"] = emb[:, 0]
out["UMAP2"] = emb[:, 1]
out["subclass"] = umap_old["subclass"].values
out["region"] = umap_old["region"].values
out["class"] = umap_old["class"].values
out["dominant_program"] = umap_old["dominant_program"].values
out["batch"] = batch.values
out.to_csv(UMAP_CSV)
print("Wrote", UMAP_CSV, out.shape, "cols:", list(out.columns))

# stash the dpca embedding + emb for the verify/figure steps
np.save("/tmp/dpca_embedding.npy", dpca)
out.to_csv("/tmp/umap_embedding_dpca.csv")
print("DONE stage1_umap_dualpca")
