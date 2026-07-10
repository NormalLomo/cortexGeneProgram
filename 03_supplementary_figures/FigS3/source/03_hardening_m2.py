#!/usr/bin/env python
"""TEST M2 -- spatial-autocorrelation-preserving null for Fig7 P40 x DMN.

Original (fig7): Spearman r between P40 region-z (14 ROIs) and Neurosynth
'default mode' z (14 ROIs); null = region-LABEL permutation (10k) -> perm_p=1e-4,
but BH-FDR across the 60xnTerm grid = 0.072 (already only nominal). Region-label
permutation does NOT preserve spatial autocorrelation -> anti-conservative.

Hardening: build a spatial null that preserves the spatial autocorrelation of the
14-ROI map, then re-test P40 x DMN.

Spatial substrate: we derive REAL MNI centroids for each of the 14 cortical ROIs by
loading the same Harvard-Oxford atlas + REGION_MAP used in fig7, taking the
volume centroid (voxel COM -> world mm) of each ROI's label set. Euclidean distance
between centroids -> 14x14 distance matrix D.

Two SA-preserving nulls (both honest about n=14 weakness):
 (1) VARIOGRAM-MATCHED surrogate maps (BrainSMASH algorithm, re-implemented from
     scratch since the brainsmash package is not installed): for the brain map x
     (here P40 region-z), generate surrogates that (a) randomly permute then (b)
     smooth via a distance-kernel and rescale to match x's empirical variogram in
     distance bins. Correlate each surrogate with the DMN map -> null distribution
     of Spearman r; empirical two-sided p.
 (2) MORAN SPECTRAL RANDOMIZATION (MSR): eigen-decompose the row-standardized spatial
     weights W = 1/D (zero diag); generate surrogates of x as random rotations in the
     Moran eigenvector basis preserving the spatial-autocorrelation spectrum; correlate
     with DMN map -> null p. (Singleton, Spatial Statistics 2017.)

We surrogate the P40 map and correlate against the FIXED DMN map (standard direction).
Report observed r, both spatial-null p's, and BH context. n=14 is small; both nulls
have limited resolution -> we state this explicitly. If brainsmash-on-surface were
available it would be preferable; volumetric-centroid variogram is the feasible
approximation here.

Outputs -> results/crossregion_v1/hardening/M2_spatial_null.tsv (+ console verdict).
"""
import warnings; warnings.filterwarnings("ignore")
import os, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

ROOT = "CORTEX_PROGRAM_ROOT"
COG  = f"{ROOT}/results/crossregion_v1/program_cognition"
OUT  = f"{ROOT}/results/crossregion_v1/hardening"
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(0)

REGION_MAP = {
    "FPPFC": [1], "DLPFC": [4], "VLPFC": [5,6], "M1": [7], "ACC": [28,29],
    "S1": [17], "PoCG": [17], "S1E": [17], "STG": [9,10], "ITG": [14,15,16],
    "SMG": [19,20], "SPL": [18], "AG": [21], "V1": [24,47,48],
}

# ---- region maps ----
rt = pd.read_csv(f"{COG}/region_term_neurosynth.tsv", sep="\t", index_col=0)
z  = pd.read_csv(f"{ROOT}/results/crossregion_v1/program_region_zscore.tsv", sep="\t", index_col=0)
z.columns = [str(c) for c in z.columns]
regions = [r for r in REGION_MAP if r in rt.index and r in z.index]
dmn = rt.loc[regions, "default mode"].astype(float).values
p40 = z.loc[regions, "40"].astype(float).values
r_obs = spearmanr(p40, dmn).correlation
print(f"=== TEST M2: spatial null for P40 x DMN (n={len(regions)} ROIs) ===")
print(f"  observed Spearman r = {r_obs:.4f}")

# ---- MNI centroids from Harvard-Oxford ----
print("  building HO atlas centroids ...", flush=True)
from nilearn import datasets, image
ho = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm")
ho_img = ho.maps if hasattr(ho, "maps") else image.load_img(ho["maps"])
data = np.asarray(ho_img.dataobj); aff = ho_img.affine
cent = {}
for reg in regions:
    mask = np.isin(data, REGION_MAP[reg])
    ijk = np.array(np.where(mask)).mean(axis=1)            # voxel COM
    xyz = aff @ np.r_[ijk, 1.0]
    cent[reg] = xyz[:3]
C = np.array([cent[r] for r in regions])
D = np.sqrt(((C[:, None, :] - C[None, :, :]) ** 2).sum(-1))  # 14x14 mm
print("  median inter-ROI dist (mm):", round(float(np.median(D[D > 0])), 1))

# =========================================================
# (1) Variogram-matched surrogates (BrainSMASH-style), surrogate the P40 map
# =========================================================
def variogram(x, D, bins):
    xv = (x - x.mean())
    sv = 0.5 * (xv[:, None] - xv[None, :]) ** 2
    iu = np.triu_indices_from(D, k=1)
    dd = D[iu]; vv = sv[iu]
    out = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (dd >= lo) & (dd < hi)
        out.append(vv[m].mean() if m.any() else np.nan)
    return np.array(out)

def smvar_surrogate(x, D, n=5000):
    """Generate variogram-matched surrogates: permute, distance-weight-smooth over a
    range of kernel bandwidths, pick bandwidth+rescale (regress) minimizing variogram
    mismatch to target. Returns array (n x len(x))."""
    nloc = len(x)
    dmax = D.max()
    bins = np.linspace(0, dmax + 1e-6, 7)
    v_target = variogram(rankdata(x), D, bins)         # match on ranks (Spearman target)
    # candidate kernel bandwidths (fraction of max dist)
    bws = dmax * np.array([0.15, 0.25, 0.4, 0.6, 0.85, 1.2])
    surr = np.empty((n, nloc))
    xr = rankdata(x)
    for i in range(n):
        perm = RNG.permutation(xr)
        best = None; best_err = np.inf
        for bw in bws:
            Kk = np.exp(-(D ** 2) / (2 * bw ** 2)); np.fill_diagonal(Kk, 0)
            sm = (Kk @ perm) / (Kk.sum(1) + 1e-9)
            # regress sm onto a copy to match scale of target ranks
            sm = (sm - sm.mean()) / (sm.std() + 1e-9)
            sm = sm * xr.std() + xr.mean()
            err = np.nansum((variogram(sm, D, bins) - v_target) ** 2)
            if err < best_err:
                best_err = err; best = sm
        surr[i] = best
    return surr

print("  (1) variogram-matched surrogates ...", flush=True)
NSURR = 5000
surr = smvar_surrogate(p40, D, n=NSURR)
null_r1 = np.array([spearmanr(s, dmn).correlation for s in surr])
p_var = (np.sum(np.abs(null_r1) >= abs(r_obs)) + 1) / (NSURR + 1)
print(f"      variogram null: p(|r_null|>=|r_obs|) = {p_var:.4g}  "
      f"(null r: mean={null_r1.mean():.3f} sd={null_r1.std():.3f} "
      f"95pct|r|={np.percentile(np.abs(null_r1),95):.3f})")

# =========================================================
# (2) Moran Spectral Randomization, surrogate the P40 map
# =========================================================
def msr_surrogate(x, W, n=5000):
    # symmetric double-centered weights -> MEM eigenvectors
    Ws = 0.5 * (W + W.T)
    n_ = Ws.shape[0]
    one = np.ones((n_, 1))
    Cmat = np.eye(n_) - one @ one.T / n_
    B = Cmat @ Ws @ Cmat
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(-vals)
    vals = vals[order]; vecs = vecs[:, order]
    keep = np.abs(vals) > 1e-8
    V = vecs[:, keep]                                   # MEM basis (n x m)
    xc = x - x.mean()
    coef = V.T @ xc                                     # projection
    amp = np.abs(coef)
    surr = np.empty((n, n_))
    for i in range(n):
        signs = RNG.choice([-1, 1], size=len(coef))
        # random rotation preserving per-eigenvector amplitude (procrustes-free MSR)
        newc = amp * signs
        surr[i] = x.mean() + V @ newc
    return surr

W = 1.0 / D; np.fill_diagonal(W, 0.0)
print("  (2) Moran spectral randomization ...", flush=True)
surr2 = msr_surrogate(p40, W, n=NSURR)
null_r2 = np.array([spearmanr(s, dmn).correlation for s in surr2])
p_msr = (np.sum(np.abs(null_r2) >= abs(r_obs)) + 1) / (NSURR + 1)
print(f"      MSR null: p(|r_null|>=|r_obs|) = {p_msr:.4g}  "
      f"(null r: mean={null_r2.mean():.3f} sd={null_r2.std():.3f} "
      f"95pct|r|={np.percentile(np.abs(null_r2),95):.3f})")

# original label-perm p + grid FDR for reference
permp = pd.read_csv(f"{COG}/program_term_permp.tsv", sep="\t", index_col=0)
fdr   = pd.read_csv(f"{COG}/program_term_fdr.tsv", sep="\t", index_col=0)
orig_permp = float(permp.loc[40, "default mode"]) if 40 in permp.index else float(permp.loc["40","default mode"])
orig_fdr   = float(fdr.loc[40, "default mode"]) if 40 in fdr.index else float(fdr.loc["40","default mode"])

nom = (p_var < 0.05) and (p_msr < 0.05)
verdict = ("NOMINAL-ONLY (spatial-null p<0.05 but original grid FDR=%.3f>0.05 => not FDR-significant)" % orig_fdr) \
          if nom else "COLLAPSES under spatial null (spatial-null p>=0.05)"
print(f"\n  original label-perm p={orig_permp:.2e}  original grid FDR={orig_fdr:.3f}")
print(f"  VERDICT P40 x DMN: {verdict}")
print("  NOTE: n=14 ROIs => low-resolution spatial null; volumetric-centroid variogram is a")
print("        feasible surrogate (no surface mesh for geodesic BrainSMASH). Interpret as a")
print("        sanity bound, not a definitive spin-test.")

pd.DataFrame([dict(
    pair="P40_x_DMN", n_roi=len(regions), spearman_r=r_obs,
    orig_labelperm_p=orig_permp, orig_grid_FDR=orig_fdr,
    variogram_null_p=p_var, MSR_null_p=p_msr,
    variogram_null_abs_r_95=float(np.percentile(np.abs(null_r1),95)),
    MSR_null_abs_r_95=float(np.percentile(np.abs(null_r2),95)),
    verdict=verdict)]).to_csv(f"{OUT}/M2_spatial_null.tsv", sep="\t", index=False)
print("[write]", f"{OUT}/M2_spatial_null.tsv")
print("\nDONE M2")
