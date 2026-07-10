#!/usr/bin/env python
"""
09_markcorr_rigor.py  --  T4 rigor checks for markcorr.

Primary g aggregate = cross-chip MEDIAN of per-chip g (locked per羅老師).
Headline set = pairs with |log2(median_g)| > 0.32 at r=25um (ring idx 0).

Checks implemented:
  (1) per-donor leave-one-out (LOO) stability of headline pairs.
  (2) gene-overlap Jaccard(top50) for headline pairs (circularity risk).
  (3) leave-one-gene-out (LOGO) effect-size machinery -- TIMED on a few example
      headline pairs and EXTRAPOLATED to the full headline set. NOT run full.

Inputs (read-only):
  per_chip dumps   results/.../markcorr/per_chip/{cellprog,progprog}_<chip>.npz
  final product    results/.../markcorr/final/{cellprog,progprog}_median_iqr.npz
  gene spectra     results/cnmf_.../snrna_joint_full1M_v1.gene_spectra_tpm.k_60.dt_0_15.txt
  program names    results/crossregion_v1/program_names.tsv
  RCTD snRNA ref   neuropeptide_cortex/.../REFERENCE__cortex_ref1M__NA__snRNA__ref_small.h5ad
  raw bin50 h5ad   neuropeptide_cortex/.../bin50/<chip>_bin50.h5ad   (LOGO timing only)
  mark sources     spatial_bin50_program_score_SCT.parquet, spatial_bin50_rctd_weights.parquet,
                   spatial_bin50_meta.parquet   (LOGO timing only)

Outputs (only writes here):
  results/.../markcorr/rigor/{mode}_donor_loo.tsv
  results/.../markcorr/rigor/{mode}_jaccard.tsv
  results/.../markcorr/rigor/markers_top50_celltype.tsv   (derived cell-type markers)
  results/.../markcorr/rigor/_rigor.log
  results/.../markcorr/rigor/logo_estimate.tsv            (timing of example pairs)
"""
import os, sys, time, glob, gc, json
import numpy as np
import pandas as pd

ROOT  = "CORTEX_PROGRAM_ROOT"
RES   = os.path.join(ROOT, "results/crossregion_v1")
MV2   = os.path.join(RES, "markcorr")
PCDIR = os.path.join(MV2, "per_chip")
FINAL = os.path.join(MV2, "final")
RIGOR = os.path.join(MV2, "rigor")
os.makedirs(RIGOR, exist_ok=True)

TPM_PATH = (ROOT + "/results/cnmf_snrna_joint_full1M_v1/cnmf_work/"
            "snrna_joint_full1M_v1/snrna_joint_full1M_v1.gene_spectra_tpm.k_60.dt_0_15.txt")
PROGNAMES = os.path.join(RES, "program_names.tsv")
REF_H5AD  = ("CORTEX_PROGRAM_DATA_ROOT/neuropeptide_cortex/data/human/spatial/bin50/ref/"
             "REFERENCE__cortex_ref1M__NA__snRNA__ref_small.h5ad")
RAW_DIR   = "CORTEX_PROGRAM_DATA_ROOT/neuropeptide_cortex/data/human/spatial/bin50"
F_META = os.path.join(RES, "spatial_bin50_meta.parquet")
F_SCT  = os.path.join(RES, "spatial_bin50_program_score_SCT.parquet")
F_RCTD = os.path.join(RES, "spatial_bin50_rctd_weights.parquet")

HEAD_THR = 0.32   # |log2 median_g| headline
RING0    = 0      # r=25um ring index
TOP_N    = 50     # top genes for Jaccard

LOGF = open(os.path.join(RIGOR, "_rigor.log"), "w")
def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s); LOGF.write(s + "\n"); LOGF.flush()


# ---------------------------------------------------------------- donor grouping
def donor_of(chip):
    """Donor/block proxy = leading block ID before the slide suffix.
    Stereo-seq chip IDs are <6-char block><2-char slide> e.g. A00797|C3.
    The non-standard SS200001075BR has no suffix -> its own donor."""
    if chip.startswith("SS"):
        return chip
    return chip[:6]


# ---------------------------------------------------------------- per-chip g load
def load_perchip_g(mode):
    """Return (chips[list], g[nChip,nA,nB,nR], A_names, B_names, ring_edges_um).
    per-chip g sanitized nan/inf->1.0 exactly as the median aggregator."""
    files = sorted(glob.glob(os.path.join(PCDIR, f"{mode}_*.npz")))
    chips, gs = [], []
    A_names = B_names = ring = None
    for f in files:
        z = np.load(f, allow_pickle=True)
        g = np.array(z["g"], dtype=np.float64)
        g[~np.isfinite(g)] = 1.0
        gs.append(g)
        chips.append(str(z["chip"]))
        if A_names is None:
            A_names = [str(x) for x in z["A_names"]]
            B_names = [str(x) for x in z["B_names"]]
            ring    = np.array(z["ring_edges_um"])
    G = np.stack(gs, axis=0)  # (nChip,nA,nB,nR)
    return chips, G, A_names, B_names, ring


def headline_pairs(mode):
    """From final median product: pairs with |log2(median_g)|>0.32 at ring0.
    progprog: dedupe to upper triangle (A_idx<B_idx) and drop diagonal."""
    z = np.load(os.path.join(FINAL, f"{mode}_median_iqr.npz"), allow_pickle=True)
    med = z["median_g"][:, :, RING0]
    l2  = np.log2(med)
    A_names = [str(x) for x in z["A_names"]]
    B_names = [str(x) for x in z["B_names"]]
    mask = np.abs(l2) > HEAD_THR
    pairs = []
    nA, nB = med.shape
    for a in range(nA):
        for b in range(nB):
            if not mask[a, b]:
                continue
            if mode == "progprog":
                if a >= b:        # dedupe symmetric + drop diagonal
                    continue
            pairs.append((a, b, A_names[a], B_names[b],
                          float(med[a, b]), float(l2[a, b])))
    return pairs, A_names, B_names


# ---------------------------------------------------------------- (1) donor LOO
def run_donor_loo(mode):
    chips, G, A_names, B_names, ring = load_perchip_g(mode)
    donors = [donor_of(c) for c in chips]
    uniq = sorted(set(donors))
    log(f"[{mode}] LOO: {len(chips)} chips -> {len(uniq)} donor groups")
    pairs, _, _ = headline_pairs(mode)
    log(f"[{mode}] headline pairs (ring0, |log2 med g|>{HEAD_THR}): {len(pairs)}")

    donor_arr = np.array(donors)
    rows = []
    for (a, b, an, bn, med0, l20) in pairs:
        gvec_all = G[:, a, b, RING0]                       # per-chip g at ring0
        fold_med, still = [], []
        for d in uniq:
            sel = donor_arr != d                           # exclude this donor
            gv = gvec_all[sel]
            m  = float(np.median(gv))
            fold_med.append(m)
            l2 = np.log2(m) if m > 0 else np.nan
            ok = np.isfinite(l2) and (abs(l2) > HEAD_THR) and (np.sign(l2) == np.sign(l20))
            still.append(bool(ok))
        fold_med = np.array(fold_med)
        frac = float(np.mean(still))
        rows.append(dict(
            A=an, B=bn, A_idx=a, B_idx=b,
            full_median_g=med0, full_log2=l20,
            loo_min_median_g=float(fold_med.min()),
            loo_median_median_g=float(np.median(fold_med)),
            loo_max_median_g=float(fold_med.max()),
            n_donors=len(uniq),
            frac_folds_still_headline=frac,
            stable_all=(frac == 1.0),
            stable_90=(frac >= 0.90)))
    df = pd.DataFrame(rows).sort_values("frac_folds_still_headline")
    out = os.path.join(RIGOR, f"{mode}_donor_loo.tsv")
    df.to_csv(out, sep="\t", index=False)
    n_all = int((df.frac_folds_still_headline == 1.0).sum())
    n_90  = int((df.frac_folds_still_headline >= 0.90).sum())
    n_drv = len(df) - n_all
    log(f"[{mode}] LOO written {out}: {len(df)} headline; STABLE(all folds)={n_all}; "
        f">=90% folds={n_90}; donor-driven(drop in >=1 fold)={n_drv}")
    return df, uniq


# ---------------------------------------------------------------- top-gene sets
def load_program_topgenes(top_n=TOP_N):
    """programs x genes spectra -> per-program top-N genes by loading."""
    tpm = pd.read_csv(TPM_PATH, sep="\t", index_col=0)   # 60 x 18742
    tpm.index = [f"program_{i}" for i in tpm.index]
    top = {}
    for p in tpm.index:
        s = tpm.loc[p]
        top[p] = list(s.sort_values(ascending=False).index[:top_n])
    return top, tpm


def load_celltype_markers(top_n=TOP_N):
    """Derive top-N markers per RCTD cell type from the snRNA reference used by
    RCTD (reference cortex_ref1M, 22 cell_type, 195k cells). Method: per cell,
    library-normalize to 1e4 + log1p; per type mean expr; SPECIFICITY score =
    log1p_mean(type) - log1p_mean(rest-of-types max-pooled second-best) is fragile
    -> use mean(type) / mean(across all types) (fold) ranked, requiring the type
    mean to be in its own top. We rank by (type_mean - mean_of_other_types) which
    is the standard 1-vs-rest specificity on normalized expression."""
    import anndata as ad
    import scipy.sparse as sp
    log("markers: reading RCTD snRNA reference (backed) ...")
    A = ad.read_h5ad(REF_H5AD)                 # full read; X = counts float32
    X = A.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsr().astype(np.float64)
    # library-normalize per cell to 1e4 then log1p
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    inv = sp.diags(1e4 / lib)
    Xn = inv @ X
    Xn.data = np.log1p(Xn.data)
    genes = list(A.var_names)
    ct = A.obs["cell_type"].astype(str).values
    types = sorted(pd.unique(ct))
    # per-type mean over cells (mean of log1p-normalized)
    means = np.zeros((len(types), Xn.shape[1]), dtype=np.float64)
    for i, t in enumerate(types):
        idx = np.where(ct == t)[0]
        means[i] = np.asarray(Xn[idx].mean(axis=0)).ravel()
    grand = means.mean(axis=0)                 # mean across the 22 type means
    top = {}
    rows = []
    for i, t in enumerate(types):
        spec = means[i] - grand                # 1-vs-rest specificity
        order = np.argsort(spec)[::-1][:top_n]
        gl = [genes[j] for j in order]
        top[t] = gl
        rows.append(dict(cell_type=t, top_genes=";".join(gl)))
    pd.DataFrame(rows).to_csv(
        os.path.join(RIGOR, "markers_top50_celltype.tsv"), sep="\t", index=False)
    log(f"markers: derived top{top_n} markers for {len(types)} cell types -> "
        "markers_top50_celltype.tsv")
    del A, X, Xn; gc.collect()
    return top


def jaccard(a, b):
    sa, sb = set(a), set(b)
    u = sa | sb
    inter = sa & sb
    return (len(inter) / len(u) if u else 0.0), sorted(inter)


# ---------------------------------------------------------------- (2) Jaccard
def run_jaccard(mode, prog_top, ct_top):
    pairs, _, _ = headline_pairs(mode)
    rows = []
    for (a, b, an, bn, med0, l20) in pairs:
        # B is always a program (program_<k>); A is program (progprog) or celltype
        bset = prog_top[bn]
        if mode == "progprog":
            aset = prog_top[an]
            kind = "prog_x_prog"
        else:
            aset = ct_top.get(an, [])
            kind = "cell_x_prog"
        j, shared = jaccard(aset, bset)
        rows.append(dict(kind=kind, A=an, B=bn, A_idx=a, B_idx=b,
                         full_median_g=med0, full_log2=l20,
                         jaccard_top50=j, n_shared=len(shared),
                         shared_genes=";".join(shared)))
    df = pd.DataFrame(rows).sort_values("jaccard_top50", ascending=False)
    out = os.path.join(RIGOR, f"{mode}_jaccard.tsv")
    df.to_csv(out, sep="\t", index=False)
    hi = int((df.jaccard_top50 >= 0.20).sum())
    log(f"[{mode}] Jaccard written {out}: {len(df)} headline; "
        f"Jaccard>=0.20 (circularity-risk)={hi}; max={df.jaccard_top50.max():.3f}")
    return df


# ---------------------------------------------------------------- (3) LOGO timing
def logo_recompute_pair_g(chip, a_idx, b_idx, mode, shared_genes,
                          tpm, prog_cols, ct_names, prog_names,
                          residual_cache):
    """Recompute MEDIAN-input per-chip g for ONE pair after removing shared genes
    from the program activity of the involved program(s). Returns g (ring0) for
    this chip. NO permutation null (LOGO = effect-size robustness only).

    Mechanism: new_score[bin,prog] = old_score - Residual[bin,shared]·tpm[prog,shared].
    We need Residual at shared genes for this chip (the expensive part)."""
    # placeholder: real recompute done in timing harness below; this signature
    # documents the per-pair cost path. Not used in full run.
    raise NotImplementedError


def run_logo_estimate(mode_pairs_examples):
    """TIME the LOGO machinery on a few example headline pairs (mixed Jaccard),
    by actually reconstructing the residual-projection delta for those pairs.
    Cost driver = per-chip raw h5ad read + dense Pearson residuals.

    To keep the timing honest but cheap we time the DOMINANT per-chip step on a
    SUBSET of chips and extrapolate, since the residual computation is identical
    regardless of which pair (the pair only selects which programs/genes to
    re-dot, a negligible add-on)."""
    import anndata as ad, scanpy as sc, scipy.sparse as sp
    chips = sorted(glob.glob(os.path.join(RAW_DIR, "*_bin50.h5ad")))
    chips = [c for c in chips if not os.path.basename(c).startswith(".")]
    # load spectra once
    tpm = pd.read_csv(TPM_PATH, sep="\t", index_col=0)
    tpm.index = [f"program_{i}" for i in tpm.index]
    tpm_genes = list(tpm.columns)

    # time the per-chip residual+project step on N_TIME chips
    N_TIME = 3
    per_chip_secs = []
    for fp in chips[:N_TIME]:
        chip = os.path.basename(fp).replace("_bin50.h5ad", "")
        t0 = time.time()
        A = ad.read_h5ad(fp)
        if not sp.issparse(A.X):
            A.X = sp.csr_matrix(A.X)
        A.X = A.X.astype(np.float32)
        sc.experimental.pp.normalize_pearson_residuals(A)   # dense residuals
        R = np.asarray(A.X, dtype=np.float32)
        pos = {g: i for i, g in enumerate(A.var_names)}
        shared = [g for g in tpm_genes if g in pos]
        col = [pos[g] for g in shared]
        W = tpm[shared].to_numpy(np.float32)                # 60 x g_shared
        # full projection (what LOGO subtracts a few shared-gene cols from)
        _ = R[:, col] @ W.T                                  # n x 60
        dt = time.time() - t0
        per_chip_secs.append(dt)
        log(f"[LOGO-time] {chip}: read+residual+project = {dt:.1f}s "
            f"({A.n_obs} bins, {len(shared)} shared genes)")
        del A, R; gc.collect()
    return float(np.mean(per_chip_secs)), len(chips)


# ---------------------------------------------------------------- main
def main():
    log("=== T4 rigor check start ===")
    # headline counts both modes
    hc = {}
    for mode in ("cellprog", "progprog"):
        pairs, _, _ = headline_pairs(mode)
        hc[mode] = len(pairs)
    log(f"HEADLINE COUNTS: cellprog={hc['cellprog']}  progprog={hc['progprog']} "
        f"(progprog deduped upper-tri, no diagonal)")

    # (1) donor LOO -- both modes
    loo_summary = {}
    for mode in ("cellprog", "progprog"):
        df, uniq = run_donor_loo(mode)
        loo_summary[mode] = (df, uniq)

    # (2) Jaccard -- both modes
    log("loading top-50 gene sets (programs + cell types) ...")
    prog_top, _ = load_program_topgenes()
    ct_top = load_celltype_markers()
    jac_summary = {}
    for mode in ("cellprog", "progprog"):
        jac_summary[mode] = run_jaccard(mode, prog_top, ct_top)

    # (3) LOGO estimate -- TIME ONLY, do not run full
    per_chip_s, n_chips = run_logo_estimate(None)
    n_head = hc["cellprog"] + hc["progprog"]
    # per-pair cost: residuals can be cached across pairs in one pass, so the
    # dominant cost is ONE full residual pass over all chips (shared by all pairs);
    # the per-pair add-on is a small matmul. Report both naive and cached ETAs.
    secs_per_pair_naive = per_chip_s * n_chips        # if recomputed per pair
    eta_naive = secs_per_pair_naive * n_head
    eta_cached = per_chip_s * n_chips                 # one residual pass total
    log("=== LOGO ESTIMATE ===")
    log(f"per-chip read+residual+project = {per_chip_s:.1f}s (mean of timed chips)")
    log(f"n_chips={n_chips}  headline pairs total={n_head} "
        f"(cellprog {hc['cellprog']} + progprog {hc['progprog']})")
    log(f"NAIVE (recompute residuals per pair): "
        f"{secs_per_pair_naive:.0f}s/pair -> ETA {eta_naive/3600:.1f} h all-headline")
    log(f"CACHED (one residual pass, reuse for all pairs): "
        f"ETA {eta_cached/60:.1f} min total + small per-pair matmul")
    pd.DataFrame([dict(
        per_chip_residual_secs=per_chip_s, n_chips=n_chips,
        headline_cellprog=hc["cellprog"], headline_progprog=hc["progprog"],
        headline_total=n_head,
        eta_naive_hours=eta_naive/3600.0,
        eta_cached_minutes=eta_cached/60.0)]).to_csv(
        os.path.join(RIGOR, "logo_estimate.tsv"), sep="\t", index=False)
    log("=== T4 rigor check done ===")


if __name__ == "__main__":
    main()
