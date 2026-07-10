#!/usr/bin/env python
"""
Fig2: Spatial bin50 UMAP in cNMF program-score space (SCT).

Each point = one bin50 spatial bin, embedded by its 60 SCT program z-scores.
Stratified subsample by majorDomain x region (seed=0), StandardScaler the 60
program-z cols, then scanpy neighbors+UMAP. Color by anatomy (majorDomain),
batch (region), and overlay program-score feature plots.

Reproducible: seeds + params hardcoded. Read-only on source parquet.

Run:  python scripts/fig2/stage_spatial_umap.py
"""
import os
# cap BLAS / thread pools BEFORE numpy/sklearn import: on this box the default
# OpenBLAS build overflows its precompiled NUM_THREADS during NearestNeighbors,
# which can hard-crash the process. Keep threading bounded + deterministic.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "8")
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

import scanpy as sc
from sklearn.preprocessing import StandardScaler

# ------------------------------------------------------------------ params
SEED = 0
N_TARGET = 150_000          # target subsample size
N_NEIGHBORS = 30
UMAP_MIN_DIST = 0.3
METRIC = "euclidean"
N_PROG = 60
PROG_COLS = [f"program_{i}" for i in range(1, N_PROG + 1)]

PROJ = "CORTEX_PROGRAM_ROOT"
SCORE_PARQUET = f"{PROJ}/results/crossregion_v1/spatial_bin50_program_score_SCT.parquet"
NAMES_TSV = f"{PROJ}/results/crossregion_v1/program_names.tsv"

OUT_CSV = f"{PROJ}/figures/fig2/_intermediate/spatial_umap_progscore.csv"
OUT_PNG = f"{PROJ}/figures/fig2/_diag/spatial_umap_progscore.png"

# feature programs to overlay: (program_number, short_label_role)
FEATURE_PROGS = [37, 7, 29, 15, 40, 8]

# anatomy ordering for legend / palette
MAJORDOMAIN_ORDER = ["L1", "L2", "L3", "L4", "L5", "L6", "WM", "ARACHNOID"]

np.random.seed(SEED)


def log(msg):
    print(f"[fig2-umap] {msg}", flush=True)


def main():
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)

    # -------------------------------------------------- load
    log("reading parquet (meta + 60 program cols) ...")
    cols = ["bin", "majorDomain", "domain", "region", "bin_total_umi", "x", "y"] + PROG_COLS
    df = pq.read_table(SCORE_PARQUET, columns=cols).to_pandas()
    log(f"loaded {len(df):,} bins x {df.shape[1]} cols")

    # program name map
    names = pd.read_csv(NAMES_TSV, sep="\t")
    name_map = dict(zip(names["program"].astype(int), names["name_short"].astype(str)))

    # -------------------------------------------------- stratified subsample
    # strata = majorDomain x region; cap per stratum so total ~ N_TARGET while
    # guaranteeing every majorDomain and every region is represented.
    df["_stratum"] = df["majorDomain"].astype(str) + "||" + df["region"].astype(str)
    strata = df["_stratum"].unique()
    n_strata = len(strata)
    cap = max(1, int(round(N_TARGET / n_strata)))
    log(f"stratified subsample: {n_strata} strata (majorDomain x region), cap={cap}/stratum")

    rng = np.random.default_rng(SEED)
    parts = []
    for s, g in df.groupby("_stratum", sort=True):
        take = min(cap, len(g))
        idx = rng.choice(g.index.values, size=take, replace=False)
        parts.append(idx)
    sel_idx = np.concatenate(parts)
    sub = df.loc[sel_idx].copy().reset_index(drop=True)
    sub = sub.drop(columns=["_stratum"])
    log(f"subsample size = {len(sub):,} bins")

    log("subsample composition (n per majorDomain):")
    comp_md = sub["majorDomain"].value_counts().reindex(MAJORDOMAIN_ORDER).dropna().astype(int)
    for k, v in comp_md.items():
        log(f"    {k:<10s} {v:>7,d}")
    log("subsample composition (n per region):")
    comp_rg = sub["region"].value_counts().sort_index()
    for k, v in comp_rg.items():
        log(f"    {k:<8s} {v:>7,d}")

    # -------------------------------------------------- UMAP
    X = sub[PROG_COLS].to_numpy(dtype=np.float32)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    Xs = StandardScaler().fit_transform(X).astype(np.float32)
    log(f"StandardScaler applied to {N_PROG} program-z cols; shape={Xs.shape}")

    ad = sc.AnnData(Xs)
    ad.obs["majorDomain"] = pd.Categorical(sub["majorDomain"].values,
                                            categories=MAJORDOMAIN_ORDER)
    ad.obs["region"] = sub["region"].astype(str).values
    log(f"neighbors: n_neighbors={N_NEIGHBORS}, metric={METRIC}, random_state={SEED}")
    sc.pp.neighbors(ad, n_neighbors=N_NEIGHBORS, metric=METRIC,
                    use_rep="X", random_state=SEED)
    log(f"umap: min_dist={UMAP_MIN_DIST}, random_state={SEED}")
    sc.tl.umap(ad, min_dist=UMAP_MIN_DIST, random_state=SEED)

    emb = ad.obsm["X_umap"]
    sub["UMAP1"] = emb[:, 0]
    sub["UMAP2"] = emb[:, 1]

    # -------------------------------------------------- batch check: kNN region mixing
    # For each point's nearest neighbors (in the 60-d scaled program space),
    # fraction whose region differs from the point. Compare observed mean to the
    # expected fraction under random labeling = 1 - sum(p_r^2) (Simpson-style).
    # We REUSE the kNN graph scanpy already built in sc.pp.neighbors (ad.obsp
    # ['distances']) instead of refitting sklearn NearestNeighbors (whose
    # OpenBLAS path crashes on this box).
    log("kNN region-mixing metric (from scanpy neighbor graph) ...")
    reg_codes = pd.Categorical(sub["region"]).codes
    dist = ad.obsp["distances"].tocsr()  # k-1 nonzero neighbors per row (self excluded)
    indptr, indices = dist.indptr, dist.indices
    n = dist.shape[0]
    frac_diff_region = np.empty(n, dtype=np.float64)
    for i in range(n):
        nbr = indices[indptr[i]:indptr[i + 1]]
        if nbr.size == 0:
            frac_diff_region[i] = np.nan
            continue
        frac_diff_region[i] = float(np.mean(reg_codes[nbr] != reg_codes[i]))
    observed_mix = float(np.nanmean(frac_diff_region))
    p = sub["region"].value_counts(normalize=True).to_numpy()
    expected_mix = float(1.0 - np.sum(p ** 2))
    mix_ratio = observed_mix / expected_mix if expected_mix > 0 else float("nan")
    log(f"region mixing: observed={observed_mix:.3f}  expected(random)={expected_mix:.3f}  "
        f"ratio={mix_ratio:.3f}")
    if mix_ratio >= 0.80:
        verdict = ("ANATOMY-DRIVEN: neighbors are well mixed across regions "
                   "(observed ~ random), so the UMAP is NOT dominated by chip/region batch.")
    elif mix_ratio >= 0.55:
        verdict = ("MOSTLY ANATOMY with some regional structure: partial region "
                   "clustering present but not dominant; integration optional.")
    else:
        verdict = ("BATCH-DRIVEN (FLAG): neighbors strongly share region -> UMAP is "
                   "dominated by chip/region batch. Recommend dual_pca/harmony "
                   "integration by chip before interpreting anatomy.")
    log(f"VERDICT: {verdict}")

    # -------------------------------------------------- save embedding csv
    out_emb = sub[["bin", "UMAP1", "UMAP2", "majorDomain", "region", "domain"]].copy()
    out_emb.to_csv(OUT_CSV, index=False)
    log(f"wrote embedding -> {OUT_CSV}  ({len(out_emb):,} rows)")

    # -------------------------------------------------- figure
    log("rendering figure ...")
    plt.rcParams.update({"font.size": 7, "axes.titlesize": 8,
                         "figure.dpi": 150, "savefig.dpi": 150})
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.6))
    x = sub["UMAP1"].to_numpy()
    y = sub["UMAP2"].to_numpy()
    pt = 1.2

    # (a) majorDomain
    ax = axes[0, 0]
    md_cats = [c for c in MAJORDOMAIN_ORDER if c in set(sub["majorDomain"].dropna())]
    cmap_md = plt.get_cmap("tab10")
    md_color = {c: cmap_md(i) for i, c in enumerate(md_cats)}
    for c in md_cats:
        m = (sub["majorDomain"] == c).to_numpy()
        ax.scatter(x[m], y[m], s=pt, c=[md_color[c]], linewidths=0, label=c, rasterized=True)
    ax.set_title("(a) majorDomain (anatomy)")
    ax.legend(markerscale=6, fontsize=5, loc="best", frameon=False, ncol=2)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")

    # (b) region (batch check)
    ax = axes[0, 1]
    rg_cats = sorted(sub["region"].unique())
    cmap_rg = plt.get_cmap("tab20")
    rg_color = {c: cmap_rg(i % 20) for i, c in enumerate(rg_cats)}
    for c in rg_cats:
        m = (sub["region"] == c).to_numpy()
        ax.scatter(x[m], y[m], s=pt, c=[rg_color[c]], linewidths=0, label=c, rasterized=True)
    ax.set_title(f"(b) region (batch check)\nkNN mix ratio={mix_ratio:.2f}")
    ax.legend(markerscale=6, fontsize=4.5, loc="best", frameon=False, ncol=2)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")

    # (c-h) feature overlays
    panel_letters = ["c", "d", "e", "f", "g", "h"]
    flat_axes = [axes[0, 2], axes[0, 3], axes[1, 0], axes[1, 1], axes[1, 2], axes[1, 3]]
    for ax, letter, pnum in zip(flat_axes, panel_letters, FEATURE_PROGS):
        z = sub[f"program_{pnum}"].to_numpy(dtype=np.float32)
        z = np.nan_to_num(z, nan=0.0)
        lo, hi = np.percentile(z, [2, 98])
        vmax = max(abs(lo), abs(hi), 1e-6)
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        order = np.argsort(np.abs(z))  # plot strong-signal points on top
        scat = ax.scatter(x[order], y[order], s=pt, c=z[order], cmap="RdBu_r",
                          norm=norm, linewidths=0, rasterized=True)
        nm = name_map.get(pnum, "")
        ax.set_title(f"({letter}) P{pnum} {nm}")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        cb = fig.colorbar(scat, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=5)
        cb.set_label("program z", fontsize=5)

    fig.suptitle("Fig2: spatial bin UMAP (cNMF program-score space, SCT)",
                 fontsize=11, y=0.995)
    fig.text(0.005, 0.005,
             f"n={len(sub):,} bins | stratified majorDomain x region (cap {cap}) | "
             f"seed={SEED} | neighbors n={N_NEIGHBORS} {METRIC} | umap min_dist={UMAP_MIN_DIST} | "
             f"region kNN mix ratio={mix_ratio:.2f} (obs {observed_mix:.2f} / exp {expected_mix:.2f})",
             fontsize=5, ha="left")
    fig.tight_layout(rect=[0, 0.015, 1, 0.985])
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    log(f"wrote figure -> {OUT_PNG}")

    # machine-readable summary line for the report
    log(f"SUMMARY_JSON {{'n_sub': {len(sub)}, 'cap': {cap}, 'n_strata': {n_strata}, "
        f"'observed_mix': {observed_mix:.4f}, 'expected_mix': {expected_mix:.4f}, "
        f"'mix_ratio': {mix_ratio:.4f}}}")
    log("DONE")


if __name__ == "__main__":
    main()
