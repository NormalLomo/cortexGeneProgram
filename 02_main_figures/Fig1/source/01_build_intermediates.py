#!/usr/bin/env python
# =====================================================================
# Fig.1 stage-0 (REPRODUCIBLE REBUILD of the cached _intermediate CSVs)
#
# WHY: the 9 Fig.1 "_intermediate" CSVs below are read by the LIVE Fig.1
# panel scripts (b / c / f / i / j) but their ONLY producer was lost in a
# project migration (it lived in scripts/fig1/_deprecated and was never
# re-authored -- see _void_20260527/.../REBUILD_NEEDED.txt). This script
# RESTORES an in-tree, self-contained, seeded producer recovered from the
# cached CSV schemas + the surviving stage1.log 7-step recipe.
#
# PRODUCES (into an output dir; default = a _verify dir, NOT the working one):
#   subclass_class.csv               subclass -> class map
#   program_x_subclass_mean.csv      program(1..60) x subclass : mean cNMF usage
#   program_specificity.csv          program, entropy, gini, dom_subclass, dom_class
#   program_program_corr.csv         60x60 Pearson corr of per-cell usages
#   top100_genes_per_program.csv     program, gene, loading, rank (top100 by loading)
#   top5_genes_per_program.csv       program, gene, loading, rank (top5 by loading)
#   marker_on_umap_subsample.csv     150k cells x 25 markers : log1p(CP10k) expr
#   marker_x_subclass_mean.csv       25 markers x subclass : mean log1p(CP10k) expr
#   markers_found.json               {"found": [markers present in atlas]}
# (NOT umap_embedding.csv -- that has its own active producer,
#  scripts/fig1/stage1_umap_dualpca.py.)
#
# INPUTS (all in-project except the intentional external snRNA atlas):
#   IN-PROJECT  results/crossregion_v1/cell_program_region_subclass.parquet
#               (1,036,039 cells x [60 program usages + region + subclass];
#                its 60 usage cols == k60_cell_scores.tsv exactly -- verified)
#   IN-PROJECT  results/cnmf_snrna_joint_full1M_v1/
#                 snrna_joint_full1M_v1_k60_factor_loadings.tsv (60 x 18742 genes)
#                 -- this is the EXACT source the original used: the cached
#                    top-gene loadings match factor_loadings (NOT gene_spectra_tpm,
#                    which does not exist in-project for snRNA k60).
#   IN-PROJECT  inputs/snRNA_1M_obs.csv  (subclass/class/region per cell)
#   IN-PROJECT  inputs/geneInfo_snRNA.csv (gene_name <-> Ensembl gene_id map;
#                    the atlas var_names are Ensembl-93 IDs, not symbols)
#   IN-PROJECT  figures/fig1/_intermediate/umap_embedding.csv  (read ONLY to
#                    lock the SAME 150,000 barcodes for the marker subsample,
#                    so marker_on_umap_subsample stays row-aligned to the UMAP)
#   EXTERNAL    CORTEX_PROGRAM_DATA_ROOT/neuropeptide_cortex/data/human/snrna/
#                    snRNA_1M.h5ad  (raw int counts; intentional shared atlas)
#
# NORMALIZATION (recovered + bit-exact verified against the cached CSV):
#   marker expression = log1p( raw_count / cell_total_counts * 1e4 )
#   i.e. scanpy normalize_total(target_sum=1e4) then log1p. Per-cell totals
#   use the cell's FULL transcriptome counts (all 32,649 atlas genes).
#
# SPECIFICITY metric (recovered, verified ~1e-7 vs cache):
#   on the program x subclass MEAN matrix m[p, :]:
#     entropy = -sum(P*ln P)/ln(n_subclass)   with P = m[p]/sum(m[p])
#     gini    = standard Gini of the subclass-mean vector
#     dominant_subclass = argmax_subclass m[p];  dominant_class = its class
#
# SEED: 0 everywhere. All params hardcoded below. Deterministic given inputs.
#
# NOTE ON program_program_corr.csv: the original was computed on the
# stage-1 150k SUBSAMPLE (stage1.log step [5]). That exact subsample seed
# was lost when umap_embedding.csv was later rebuilt (dualpca, then harmony).
# We therefore recompute corr on the CURRENT umap_embedding.csv 150k
# barcodes -- the canonical, reproducible definition going forward. This is
# expected NOT to bit-match the stale cache (older subsample); see report.
# =====================================================================
import os, sys, json, time, argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp

SEED = 0
np.random.seed(SEED)

PROJ = "CORTEX_PROGRAM_ROOT"
INT = os.path.join(PROJ, "figures/fig1/_intermediate")

PARQUET = os.path.join(PROJ, "results/crossregion_v1/cell_program_region_subclass.parquet")
LOADINGS = os.path.join(PROJ, "results/cnmf_snrna_joint_full1M_v1/"
                              "snrna_joint_full1M_v1_k60_factor_loadings.tsv")
OBS = os.path.join(PROJ, "inputs/snRNA_1M_obs.csv")
GENEINFO = os.path.join(PROJ, "inputs/geneInfo_snRNA.csv")
UMAP_EMB = os.path.join(INT, "umap_embedding.csv")        # read-only: lock 150k barcodes
ATLAS = "CORTEX_PROGRAM_DATA_ROOT/neuropeptide_cortex/data/human/snrna/snRNA_1M.h5ad"  # external

K = 60
PROGCOLS = [str(i) for i in range(1, K + 1)]
N_SUBSAMPLE = 150000
TARGET_SUM = 1e4
# canonical marker panel (order preserved == cached markers_found.json order)
MARKERS = ["SLC17A7", "SATB2", "CUX2", "RORB", "FEZF2", "BCL11B", "TLE4",
           "GAD1", "GAD2", "LAMP5", "PVALB", "SST", "VIP", "AQP4", "GFAP",
           "MBP", "PLP1", "MOBP", "PDGFRA", "OLIG1", "CX3CR1", "P2RY12",
           "CLDN5", "FLT1", "PDGFRB"]


def entropy_norm(v):
    v = np.asarray(v, float)
    p = v / v.sum()
    p = p[p > 0]
    return float(-(p * np.log(p)).sum() / np.log(len(v)))


def gini(v):
    x = np.sort(np.asarray(v, float))
    n = len(x)
    cum = np.cumsum(x)
    return float((n + 1 - 2 * (cum.sum() / cum[-1])) / n)


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    print(f"[stage0] OUTDIR = {outdir}")

    # ---- load parquet (programs + region + subclass), 1M cells ----------
    pq = pd.read_parquet(PARQUET)
    assert list(pq.columns[:K]) == PROGCOLS, "parquet program cols mismatch"
    print(f"[1] parquet: {pq.shape}")

    # ---- subclass -> class map (from obs) -------------------------------
    obs = pd.read_csv(OBS, index_col=0)
    scc = (obs.drop_duplicates("subclass")
              .set_index("subclass")["class"].to_dict())
    subclasses = sorted(pq["subclass"].dropna().unique().tolist())
    sc_df = pd.DataFrame({"subclass": subclasses,
                          "class": [scc[s] for s in subclasses]})
    sc_df.to_csv(os.path.join(outdir, "subclass_class.csv"), index=False)
    print(f"[2] subclass_class.csv: {sc_df.shape}")

    # ---- program x subclass mean usage ----------------------------------
    msub = pq.groupby("subclass")[PROGCOLS].mean()          # subclass x program
    m = msub.T                                              # program x subclass
    m.index = m.index.astype(int)
    m = m[subclasses]                                       # column order = sorted
    m_out = m.reset_index().rename(columns={"index": "program"})
    m_out.columns = ["program"] + subclasses
    m_out.to_csv(os.path.join(outdir, "program_x_subclass_mean.csv"), index=False)
    print(f"[3] program_x_subclass_mean.csv: {m_out.shape}")

    # ---- program specificity (entropy + gini over the mean matrix) ------
    rows = []
    for p in range(1, K + 1):
        v = m.loc[p].values
        dom = m.loc[p].idxmax()
        rows.append({"program": p,
                     "entropy": entropy_norm(v),
                     "gini": gini(v),
                     "dominant_subclass": dom,
                     "dominant_class": scc.get(dom)})
    spec = pd.DataFrame(rows)
    spec.to_csv(os.path.join(outdir, "program_specificity.csv"), index=False)
    print(f"[3b] program_specificity.csv: {spec.shape}")

    # ---- program-program correlation (per-cell usage, 150k subsample) ---
    emb = pd.read_csv(UMAP_EMB, index_col=0)
    bc150 = emb.index.tolist()
    assert len(bc150) == N_SUBSAMPLE, f"expected {N_SUBSAMPLE}, got {len(bc150)}"
    miss = [b for b in bc150 if b not in pq.index]
    assert len(miss) == 0, f"{len(miss)} 150k barcodes missing from parquet"
    usage150 = pq.loc[bc150, PROGCOLS]
    corr = usage150.corr()                                  # 60 x 60 Pearson
    corr.index = PROGCOLS
    corr.columns = PROGCOLS
    corr.to_csv(os.path.join(outdir, "program_program_corr.csv"))
    print(f"[4] program_program_corr.csv: {corr.shape} (on current 150k)")

    # ---- top genes per program (from factor_loadings) -------------------
    fl = pd.read_csv(LOADINGS, sep="\t", index_col=0)       # program x gene
    fl.index = fl.index.astype(int)
    top100_rows, top5_rows = [], []
    for p in range(1, K + 1):
        s = fl.loc[p].sort_values(ascending=False)
        for rank, (gene, loading) in enumerate(s.head(100).items(), start=1):
            rec = {"program": p, "gene": gene, "loading": float(loading), "rank": rank}
            top100_rows.append(rec)
            if rank <= 5:
                top5_rows.append(rec)
    pd.DataFrame(top100_rows).to_csv(
        os.path.join(outdir, "top100_genes_per_program.csv"), index=False)
    pd.DataFrame(top5_rows).to_csv(
        os.path.join(outdir, "top5_genes_per_program.csv"), index=False)
    print(f"[7] top100_genes_per_program.csv: {len(top100_rows)} rows; "
          f"top5: {len(top5_rows)} rows")

    # ---- markers from the snRNA atlas (Ensembl var_names) ---------------
    import anndata as ad
    gi = pd.read_csv(GENEINFO)
    sym2ens = dict(zip(gi["gene_name"], gi["gene_id"]))
    a = ad.read_h5ad(ATLAS, backed="r")
    var_pos = {g: i for i, g in enumerate(a.var_names)}
    found = []
    found_cols = []      # (symbol, ensembl, col_index)
    for sym in MARKERS:
        ens = sym2ens.get(sym)
        if ens is not None and ens in var_pos:
            found.append(sym)
            found_cols.append((sym, ens, var_pos[ens]))
    with open(os.path.join(outdir, "markers_found.json"), "w") as fh:
        json.dump({"found": found}, fh)
    print(f"[6] markers_found.json: {len(found)} / {len(MARKERS)} markers")

    # per-cell total counts (full transcriptome) for ALL 1M cells, then
    # normalize_total(1e4)+log1p on the marker columns only.
    bc_atlas = list(a.obs_names)
    atlas_pos = {b: i for i, b in enumerate(bc_atlas)}
    col_idx = [c[2] for c in found_cols]

    # --- marker_x_subclass_mean : mean log1p(CP10k) over ALL 1M cells ----
    # stream the X matrix in chunks to get totals + marker counts.
    n = a.n_obs
    CHUNK = 100000
    marker_norm_full = np.zeros((n, len(found_cols)), dtype=np.float64)
    for start in range(0, n, CHUNK):
        stop = min(start + CHUNK, n)
        Xc = a.X[start:stop]
        Xc = Xc.tocsr() if sp.issparse(Xc) else sp.csr_matrix(Xc)
        tot = np.asarray(Xc.sum(axis=1)).ravel().astype(np.float64)
        tot[tot == 0] = 1.0
        sub = Xc[:, col_idx].toarray().astype(np.float64)
        marker_norm_full[start:stop] = np.log1p(sub / tot[:, None] * TARGET_SUM)
        print(f"      atlas rows {stop}/{n}")
    mk_df_full = pd.DataFrame(marker_norm_full, index=bc_atlas, columns=found)

    # subclass mean (markers x subclass), aligned to obs subclass
    sub_of = obs["subclass"]
    mk_df_full = mk_df_full.copy()
    mk_df_full["__subclass__"] = sub_of.reindex(mk_df_full.index).values
    mxs = mk_df_full.groupby("__subclass__")[found].mean().T   # marker x subclass
    mxs = mxs[subclasses]
    mxs.index.name = "gene"
    mxs.to_csv(os.path.join(outdir, "marker_x_subclass_mean.csv"))
    print(f"[6b] marker_x_subclass_mean.csv: {mxs.shape}")

    # --- marker_on_umap_subsample : SAME 150k barcodes, row-aligned ------
    sub150 = mk_df_full.loc[bc150, found]
    sub150.index.name = ""
    sub150.to_csv(os.path.join(outdir, "marker_on_umap_subsample.csv"))
    print(f"[5] marker_on_umap_subsample.csv: {sub150.shape} (150k locked)")

    print("[stage0] DONE")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(PROJ, "figures/fig1/_intermediate_verify"),
                    help="output dir (default = _verify dir; does NOT clobber the working one)")
    args = ap.parse_args()
    main(args.outdir)
