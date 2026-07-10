#!/usr/bin/env python
"""
T6: Spatial co-localization (DESCRIPTIVE) + consistency check.

Human cortex cross-region gene-program paper, Fig.6.
- Spatial co-localization matrix: program (60) x cell-type (22) Pearson corr
  between program spatial z-score and RCTD celltype weight, over pass_mask bins.
- snRNA identity reference: program x subclass mean activity -> dominant subclass.
- Consistency: spatial-top-celltype vs snRNA-top-subclass per program + corr of
  the two 22-vectors (consistency_r). Overall match rate + median consistency_r.
- Optional cheap distance-ring profile on one chip for top-3 matched pairs.

DESCRIPTIVE ONLY: no permutation / null statistics.
Memory-safe: chunked single-pass accumulation of sums/sumsq/crossproducts.
"""
import os
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

BASE = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/"
F_PROG = BASE + "spatial_bin50_program_score_SCT.parquet"
F_RCTD = BASE + "spatial_bin50_rctd_weights.parquet"
F_META = BASE + "spatial_bin50_meta.parquet"
F_SNRNA = BASE + "cell_program_region_subclass.parquet"

PROGRAMS = [str(i) for i in range(1, 61)]          # snRNA program cols
PROG_COLS = ["program_%d" % i for i in range(1, 61)]  # spatial program cols
CELLTYPES = ['AST', 'CHANDELIER', 'ENDO', 'ET', 'L2-L3 IT LINC00507',
             'L3-L4 IT RORB', 'L4-L5 IT RORB', 'L6 CAR3', 'L6 CT', 'L6 IT',
             'L6B', 'LAMP5', 'MICRO', 'NDNF', 'NP', 'OLIGO', 'OPC', 'PAX6',
             'PVALB', 'SST', 'VIP', 'VLMC']
NP_, NC = 60, 22


def log(msg):
    print("[T6] " + msg, flush=True)


# --------------------------------------------------------------------------
# 1. Spatial co-localization matrix via chunked single-pass accumulation.
#    Files are row-aligned by `bin` (verified), so read by row-group position.
# --------------------------------------------------------------------------
def compute_coloc():
    pf_prog = pq.ParquetFile(F_PROG)
    n_rg = pf_prog.num_row_groups

    # RCTD weights live in a DIFFERENT row-group layout (6 vs prog's 44) but the
    # global `bin` order is identical (verified). So load RCTD once and slice it
    # by cumulative row offset to stay aligned with prog's per-row-group reads.
    rctd_full = pq.read_table(
        F_RCTD, columns=CELLTYPES + ["rctd_pass_mask"]).to_pandas()
    rctd_Y = rctd_full[CELLTYPES].values.astype(np.float64)        # (N,22)
    rctd_mask = rctd_full["rctd_pass_mask"].values.astype(bool)    # (N,)

    sum_x = np.zeros(NP_)          # Σ program
    sumsq_x = np.zeros(NP_)        # Σ program^2
    sum_y = np.zeros(NC)           # Σ celltype
    sumsq_y = np.zeros(NC)         # Σ celltype^2
    cross = np.zeros((NP_, NC))    # Σ program_p * celltype_c
    n_total = 0
    dropped_nan = 0

    off = 0
    for rg in range(n_rg):
        prog = pf_prog.read_row_group(rg, columns=PROG_COLS).to_pandas()
        nrows = len(prog)
        sl = slice(off, off + nrows)
        off += nrows

        mask = rctd_mask[sl]
        X = prog.values[mask].astype(np.float64)          # (m,60)
        Y = rctd_Y[sl][mask]                              # (m,22)

        # drop any rows with residual NaN/inf in X or Y to keep corr defined
        good = np.isfinite(X).all(axis=1) & np.isfinite(Y).all(axis=1)
        n_drop = int((~good).sum())
        dropped_nan += n_drop
        X = X[good]
        Y = Y[good]

        if X.shape[0] == 0:
            continue
        sum_x += X.sum(axis=0)
        sumsq_x += (X * X).sum(axis=0)
        sum_y += Y.sum(axis=0)
        sumsq_y += (Y * Y).sum(axis=0)
        cross += X.T @ Y
        n_total += X.shape[0]
        log("row-group %d/%d  cum_n=%d (dropped_nan=%d)"
            % (rg + 1, n_rg, n_total, dropped_nan))

    n = float(n_total)
    # Pearson r = (n*Σxy - Σx*Σy) / sqrt((n*Σx2-(Σx)^2)(n*Σy2-(Σy)^2))
    cov = n * cross - np.outer(sum_x, sum_y)             # (60,22)
    var_x = n * sumsq_x - sum_x ** 2                     # (60,)
    var_y = n * sumsq_y - sum_y ** 2                     # (22,)
    denom = np.sqrt(np.outer(var_x, var_y))
    with np.errstate(invalid="ignore", divide="ignore"):
        R = cov / denom
    R[~np.isfinite(R)] = 0.0   # zero-variance column -> define corr=0
    coloc = pd.DataFrame(R, index=["program_%d" % i for i in range(1, 61)],
                         columns=CELLTYPES)
    log("coloc matrix done: shape=%s  n_bins=%d dropped_nan=%d"
        % (str(coloc.shape), n_total, dropped_nan))
    return coloc, n_total


# --------------------------------------------------------------------------
# 2. snRNA program x subclass mean activity (60 x 22) -> dominant subclass.
# --------------------------------------------------------------------------
def compute_snrna_ref():
    df = pq.read_table(F_SNRNA, columns=PROGRAMS + ["subclass"]).to_pandas()
    # mean program activity per subclass
    grp = df.groupby("subclass")[PROGRAMS].mean()        # (subclasses x 60)
    # align to the 22 CELLTYPES order (subclass labels == celltype labels)
    missing = [c for c in CELLTYPES if c not in grp.index]
    if missing:
        log("WARN snRNA missing subclasses: %s" % missing)
    grp = grp.reindex(CELLTYPES)
    # transpose -> program x subclass (60 x 22)
    ref = grp.T
    ref.index = ["program_%d" % i for i in range(1, 61)]
    ref.columns = CELLTYPES
    log("snRNA ref done: shape=%s  cells=%d" % (str(ref.shape), len(df)))
    return ref


# --------------------------------------------------------------------------
# 3. Consistency check.
# --------------------------------------------------------------------------
def consistency(coloc, ref):
    rows = []
    for i, pg in enumerate(["program_%d" % i for i in range(1, 61)]):
        sp = coloc.loc[pg]                 # 22-vector spatial coloc
        sr = ref.loc[pg]                   # 22-vector snRNA mean activity
        sp_top = sp.idxmax()
        sr_valid = sr.dropna()
        sr_top = sr_valid.idxmax() if len(sr_valid) else None
        # consistency r between the two 22-vectors (use shared finite entries)
        a = sp.values.astype(float)
        b = sr.values.astype(float)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() >= 3 and np.std(a[ok]) > 0 and np.std(b[ok]) > 0:
            cr = float(np.corrcoef(a[ok], b[ok])[0, 1])
        else:
            cr = np.nan
        rows.append({
            "program": pg,
            "spatial_top_celltype": sp_top,
            "snrna_top_subclass": sr_top,
            "match": bool(sp_top == sr_top),
            "consistency_r": cr,
        })
    tab = pd.DataFrame(rows)
    return tab


# --------------------------------------------------------------------------
# 4. Optional cheap distance-ring profile for top-3 matched pairs, one chip.
# --------------------------------------------------------------------------
def ring_profile(coloc, tab, n_pairs=3, n_rings=6, hi_quantile=0.95,
                 max_hi=400):
    # pick top-3 matched program x celltype pairs by spatial coloc r
    matched = tab[tab["match"]].copy()
    if matched.empty:
        log("ring: no matched programs -> skip")
        return None
    matched["r"] = [coloc.loc[r["program"], r["spatial_top_celltype"]]
                    for _, r in matched.iterrows()]
    matched = matched.sort_values("r", ascending=False).head(n_pairs)

    # choose ONE representative large chip (most pass bins)
    rctd_chip = pq.read_table(
        F_RCTD, columns=["chip", "rctd_pass_mask"]).to_pandas()
    pass_counts = (rctd_chip[rctd_chip["rctd_pass_mask"]]
                   .groupby("chip").size().sort_values(ascending=False))
    big_chip = pass_counts.index[0]
    log("ring: representative chip=%s (%d pass bins)"
        % (big_chip, int(pass_counts.iloc[0])))

    # load that chip's bins only: x,y + needed programs + celltypes
    meta = pq.read_table(F_META, columns=["bin", "chip", "x", "y"]).to_pandas()
    chip_mask = (meta["chip"] == big_chip).values
    pass_mask = rctd_chip["rctd_pass_mask"].values.astype(bool)
    sel = chip_mask & pass_mask
    idx = np.where(sel)[0]
    log("ring: chip pass bins=%d" % len(idx))

    progs_needed = list(dict.fromkeys(matched["program"].tolist()))
    progs_needed_spatial = [p.replace("program_", "program_") for p in progs_needed]
    cts_needed = list(dict.fromkeys(matched["spatial_top_celltype"].tolist()))

    prog_all = pq.read_table(F_PROG, columns=progs_needed_spatial).to_pandas()
    rctd_all = pq.read_table(F_RCTD, columns=cts_needed).to_pandas()

    xy = meta.loc[idx, ["x", "y"]].values.astype(float)
    out_rows = []
    for _, mr in matched.iterrows():
        pg = mr["program"]; ct = mr["spatial_top_celltype"]
        pscore = prog_all[pg].values[idx]
        cw = rctd_all[ct].values[idx]
        finite = np.isfinite(pscore) & np.isfinite(cw)
        xs = xy[finite]; ps = pscore[finite]; cv = cw[finite]
        if len(ps) < 50:
            continue
        # high-program seed bins
        thr = np.quantile(ps, hi_quantile)
        seed = np.where(ps >= thr)[0]
        if len(seed) > max_hi:
            seed = np.random.RandomState(0).choice(seed, max_hi, replace=False)
        seed_xy = xs[seed]
        # distance from each bin to nearest seed (chunked over bins)
        from scipy.spatial import cKDTree
        tree = cKDTree(seed_xy)
        d, _ = tree.query(xs, k=1)
        # ring edges by distance quantiles
        edges = np.quantile(d, np.linspace(0, 1, n_rings + 1))
        edges = np.unique(edges)
        ring_id = np.clip(np.digitize(d, edges[1:-1]), 0, len(edges) - 2)
        for ridx in range(len(edges) - 1):
            m = ring_id == ridx
            if m.sum() == 0:
                continue
            out_rows.append({
                "program": pg, "celltype": ct, "ring": ridx,
                "dist_lo": float(edges[ridx]), "dist_hi": float(edges[ridx + 1]),
                "n_bins": int(m.sum()),
                "mean_celltype_weight": float(np.mean(cv[m])),
                "mean_program_z": float(np.mean(ps[m])),
            })
    if not out_rows:
        log("ring: no profile rows produced -> skip")
        return None
    rp = pd.DataFrame(out_rows)
    rp.insert(0, "chip", big_chip)
    log("ring: profile rows=%d" % len(rp))
    return rp


def main():
    os.chdir(BASE)
    log("=== compute spatial coloc matrix ===")
    coloc, n_bins = compute_coloc()
    assert coloc.shape == (60, 22), "coloc shape wrong"
    assert not coloc.isna().any().any(), "coloc has NaN"

    log("=== snRNA identity reference ===")
    ref = compute_snrna_ref()

    log("=== consistency check ===")
    tab = consistency(coloc, ref)
    assert len(tab) == 60, "consistency table != 60 rows"

    match_rate = tab["match"].mean()
    med_cr = tab["consistency_r"].median()
    n_match = int(tab["match"].sum())
    log("MATCH RATE = %d/60 = %.3f" % (n_match, match_rate))
    log("MEDIAN consistency_r = %.4f" % med_cr)

    # 3 example matched programs (highest spatial coloc r among matched)
    mt = tab[tab["match"]].copy()
    mt["r"] = [coloc.loc[r["program"], r["spatial_top_celltype"]]
               for _, r in mt.iterrows()]
    mt = mt.sort_values("r", ascending=False)
    log("--- top matched examples ---")
    for _, r in mt.head(5).iterrows():
        log("  %s -> %s  r=%.3f (snrna_top=%s)"
            % (r["program"], r["spatial_top_celltype"], r["r"],
               r["snrna_top_subclass"]))

    # spot-check: report OLIGO column - which program co-localizes most w/ OLIGO
    oligo_top = coloc["OLIGO"].idxmax()
    log("spot-check: program most co-localized with OLIGO weight = %s (r=%.3f)"
        % (oligo_top, coloc.loc[oligo_top, "OLIGO"]))

    # ---- write core outputs ----
    coloc.to_csv(BASE + "program_celltype_coloc.tsv", sep="\t")
    tab.to_csv(BASE + "coloc_vs_composition.tsv", sep="\t", index=False)
    log("wrote program_celltype_coloc.tsv  shape=%s" % str(coloc.shape))
    log("wrote coloc_vs_composition.tsv  shape=%s" % str(tab.shape))

    # ---- optional ring profile ----
    log("=== ring profile (optional) ===")
    try:
        rp = ring_profile(coloc, tab)
    except Exception as e:
        log("ring profile FAILED: %r -> skipping" % e)
        rp = None
    if rp is not None:
        rp.to_csv(BASE + "coloc_ring_profile.tsv", sep="\t", index=False)
        log("wrote coloc_ring_profile.tsv  shape=%s" % str(rp.shape))
    else:
        with open(BASE + "coloc_ring_profile.tsv", "w") as fh:
            fh.write("# SKIPPED: ring profile not computed (see run log)\n")
        log("ring profile SKIPPED (placeholder written)")

    # final summary block for parsing
    log("=== SUMMARY ===")
    log("n_pass_bins_used=%d" % n_bins)
    log("match_rate=%.4f (%d/60)" % (match_rate, n_match))
    log("median_consistency_r=%.4f" % med_cr)


if __name__ == "__main__":
    main()
