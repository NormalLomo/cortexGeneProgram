#!/usr/bin/env python
"""
Aggregate per-chip lattice-ring program x program spatial mark cross-correlation into
per-REGION tensors, test region-dependence, and make an Extended-Data figure.

Inputs : results/crossregion_v1/spatial_crosscorr/_perchip/<chip>.npz  (C:(10,60,60) cross-cov per ring, npairs, mu, sd)
         results/crossregion_v1/spatial_crosscorr/_chipmap.json
         results/crossregion_v1/program_names.tsv
Outputs: results/crossregion_v1/spatial_crosscorr/region_<REGION>.npz   (per-region C tensor, npairs)
         results/crossregion_v1/spatial_crosscorr/region_tensor_all.npz (stacked: regions x 10 x 60 x 60)
         results/crossregion_v1/spatial_crosscorr/top_region_variable_pairs.tsv
         results/crossregion_v1/spatial_crosscorr/region_similarity.tsv
         figures/extended/ed_spatial_crosscorr_byregion.{pdf,png}
         /tmp/ed_spatial_crosscorr.png
"""
import os, json, glob
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# FONT UNIFY (W-figfont-unify 2026-06-26): Nimbus Sans cross-engine
import matplotlib as _mpl_font
_mpl_font.rcParams["font.family"] = "sans-serif"
_mpl_font.rcParams["font.sans-serif"] = ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"]
_mpl_font.rcParams["pdf.fonttype"] = 42
_mpl_font.rcParams["ps.fonttype"] = 42
_mpl_font.rcParams["svg.fonttype"] = "none"
_mpl_font.rcParams["mathtext.fontset"] = "dejavusans"  # sans math
_mpl_font.rcParams["mathtext.default"] = "regular"  # math uses font.family
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
import matplotlib as mpl
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif", "font.sans-serif": ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6, "axes.titlesize": 7, "axes.labelsize": 6,
    "xtick.labelsize": 5, "ytick.labelsize": 5, "legend.fontsize": 5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
})

ROOT = "CORTEX_PROGRAM_ROOT"
CCDIR = f"{ROOT}/results/crossregion_v1/spatial_crosscorr"
PERCHIP = f"{CCDIR}/_perchip"
FIGDIR = f"{ROOT}/figures/extended"
os.makedirs(FIGDIR, exist_ok=True)
NPROG = 60; NRING = 10; LAT = 50
R_UM = np.arange(1, NRING+1) * LAT  # 50..500

# ---- program names + renumber map (W-tier2 patch 2026-06-26) ----
# Load renumber map: build cNMF(int, 1-based) -> new_P(int) or EXCLUDED
_rmap = pd.read_csv(f"{ROOT}/results/crossregion_v1/program_renumber_map.tsv", sep="\t")
cnmf2new = {}     # cNMF int (1-based) -> new_P int  (kept programs only)
EXCLUDED_CNMF0 = set()  # 0-based cNMF idx for excluded programs
for _, r in _rmap.iterrows():
    old_int = int(r["old_P"])
    new_val = str(r["new_P"])
    if new_val.upper() == "EXCLUDED":
        EXCLUDED_CNMF0.add(old_int - 1)
    else:
        cnmf2new[old_int] = int(new_val)

# program_names: keyed by new_P int -> name_short
pn = pd.read_csv(f"{ROOT}/results/crossregion_v1/program_names.tsv", sep="\t"); pn = pn[pn["new_P"].astype(str).str.startswith("P")]
short = {int(str(r.new_P).lstrip("P")): str(r.name_short) for _, r in pn.iterrows()}
conf  = {int(str(r.new_P).lstrip("P")): str(r.confidence) for _, r in pn.iterrows()}

def plab(i):  # i is 0-based cNMF index (0..59)
    """Translate 0-based cNMF idx -> 'P{new_P} {name_short}'; EXCLUDED -> 'P?(excl)'."""
    cnmf_1based = i + 1
    if i in EXCLUDED_CNMF0 or cnmf_1based not in cnmf2new:
        return f"cNMF{cnmf_1based} (excl)"
    new_p = cnmf2new[cnmf_1based]
    return f"P{new_p} {short.get(new_p, '?')}"

def cnmf0_to_newP(i):
    """0-based cNMF idx -> new_P int (or None if EXCLUDED)."""
    if i in EXCLUDED_CNMF0:
        return None
    return cnmf2new.get(i + 1)

# ---- chip map ----
cmap = json.load(open(f"{CCDIR}/_chipmap.json"))
order = cmap["order"]  # [chip, region, rg]
chip2reg = {c: r for c, r, _ in order}
regions = sorted(set(r for _, r, _ in order))

# ---- aggregate per region (weighted by npairs per ring) ----
reg_chips = {}
for c, r, _ in order:
    reg_chips.setdefault(r, []).append(c)

region_C = {}          # region -> (NRING,60,60) weighted-mean cross-cov
region_npairs = {}     # region -> (NRING,) total pairs
region_nchip = {}
for r in regions:
    acc = np.zeros((NRING, NPROG, NPROG))
    wsum = np.zeros(NRING)
    nch = 0
    for c in reg_chips[r]:
        f = f"{PERCHIP}/{c}.npz"
        if not os.path.exists(f):
            print(f"WARNING missing {f}"); continue
        d = np.load(f, allow_pickle=True)
        C = d["C"].astype(np.float64)        # (NRING,60,60)
        npr = d["npairs"].astype(np.float64) # (NRING,)
        for k in range(NRING):
            if npr[k] > 0:
                acc[k] += C[k] * npr[k]
                wsum[k] += npr[k]
        nch += 1
    for k in range(NRING):
        if wsum[k] > 0:
            acc[k] /= wsum[k]
    region_C[r] = acc
    region_npairs[r] = wsum
    region_nchip[r] = nch
    np.savez_compressed(f"{CCDIR}/region_{r}.npz",
                        C=acc.astype(np.float32), npairs=wsum.astype(np.float64),
                        nchip=nch, region=r)
    print(f"region {r:6s}: nchip={nch}  pairs(r=50)={wsum[0]:.3e}")

# stacked tensor
Rg = regions
T = np.stack([region_C[r] for r in Rg], axis=0)   # (nreg, NRING, 60, 60)
np.savez_compressed(f"{CCDIR}/region_tensor_all.npz",
                    T=T.astype(np.float32), regions=np.array(Rg),
                    r_um=R_UM, nchip=np.array([region_nchip[r] for r in Rg]))
nreg = len(Rg)
print("region tensor:", T.shape)

# single-chip (low replication) regions
singleton = [r for r in Rg if region_nchip[r] == 1]
print("singleton (1-chip) regions:", singleton)

# ================= ANALYSIS =================
iu, ju = np.triu_indices(NPROG, k=0)          # upper triangle incl diagonal -> unique pairs
npair_idx = len(iu)

# (a) per-pair region-variability of g at r=50um (nearest ring, k=0)
M50 = T[:, 0, :, :]                            # (nreg,60,60) cross-cov at 50um
# z-normalize per region by program sd? Use the raw cross-cov (already mark cross-cov of z).
pair_vals_50 = M50[:, iu, ju]                  # (nreg, npair)
pair_var_50 = pair_vals_50.var(axis=0)         # variance across regions per pair
pair_mean_50 = pair_vals_50.mean(axis=0)
pair_range_50 = pair_vals_50.max(axis=0) - pair_vals_50.min(axis=0)
# also at r=150 to be robust
M150 = T[:, 2, :, :]
pair_var_150 = M150[:, iu, ju].var(axis=0)

dfp = pd.DataFrame({
    "progA": iu+1, "progB": ju+1,
    "nameA": [short.get(cnmf2new.get(int(a)+1, -1), "?") for a in iu],
    "nameB": [short.get(cnmf2new.get(int(b)+1, -1), "?") for b in ju],
    "is_self": iu==ju,
    "mean_g50": pair_mean_50,
    "var_g50": pair_var_50,
    "range_g50": pair_range_50,
    "sd_g50": np.sqrt(pair_var_50),
    "var_g150": pair_var_150,
})
# argmin/argmax region per pair
amin = pair_vals_50.argmin(axis=0); amax = pair_vals_50.argmax(axis=0)
dfp["min_region"] = [Rg[i] for i in amin]
dfp["max_region"] = [Rg[i] for i in amax]
dfp["min_g50"] = pair_vals_50.min(axis=0)
dfp["max_g50"] = pair_vals_50.max(axis=0)
# W-tier2 patch 2026-06-26: drop EXCLUDED programs (cNMF 9/18/19/35/52/57) before ranking
_excl_1based = {i + 1 for i in EXCLUDED_CNMF0}  # progA/B in dfp are 1-based cNMF
dfp = dfp[~dfp["progA"].isin(_excl_1based) & ~dfp["progB"].isin(_excl_1based)].reset_index(drop=True)

dfp_sorted = dfp.sort_values("var_g50", ascending=False).reset_index(drop=True)
dfp_sorted.to_csv(f"{CCDIR}/top_region_variable_pairs.tsv", sep="\t", index=False)
print("\n=== TOP 15 region-variable program pairs (var of g(r=50um) across 14 regions) ===")
print(dfp_sorted.head(15)[["progA","nameA","progB","nameB","var_g50","range_g50","min_region","min_g50","max_region","max_g50"]].to_string(index=False))

# off-diagonal only top (cross-program, the interesting co-localization differences)
dfp_off = dfp_sorted[~dfp_sorted["is_self"]].reset_index(drop=True)
print("\n=== TOP 12 CROSS-program (off-diag) region-variable pairs ===")
print(dfp_off.head(12)[["progA","nameA","progB","nameB","var_g50","min_region","min_g50","max_region","max_g50"]].to_string(index=False))

# (b) region x region similarity by g(r=50um) fingerprint (off-diagonal flattened)
io, jo = np.triu_indices(NPROG, k=1)           # strictly off-diagonal
FP = M50[:, io, jo]                            # (nreg, n_offpair) fingerprint
# correlation across regions
Rmat = np.corrcoef(FP)                         # (nreg,nreg)
# distance for clustering
D = 1 - Rmat
np.fill_diagonal(D, 0.0)
D = (D + D.T) / 2
Z = linkage(squareform(D, checks=False), method="average")
dn = dendrogram(Z, labels=Rg, no_plot=True)
leaf_order = dn["ivl"]
leaf_idx = [Rg.index(l) for l in leaf_order]
clusters = fcluster(Z, t=2, criterion="maxclust")
sim_df = pd.DataFrame(Rmat, index=Rg, columns=Rg)
sim_df.to_csv(f"{CCDIR}/region_similarity.tsv", sep="\t")
print("\n=== Region clustering (g(r=50um) fingerprint, 2-cluster cut) ===")
for cl in sorted(set(clusters)):
    print(f"  cluster {cl}: {[Rg[i] for i in range(nreg) if clusters[i]==cl]}")
offdiag = Rmat[np.triu_indices(nreg,1)]
print(f"region-region g(r) fingerprint correlation: mean={offdiag.mean():.3f} min={offdiag.min():.3f} max={offdiag.max():.3f}")
# most/least similar region pairs
ii,jj = np.triu_indices(nreg,1)
o = np.argsort(offdiag)
print(f"  MOST similar regions: {Rg[ii[o[-1]]]}-{Rg[jj[o[-1]]]} r={offdiag[o[-1]]:.3f}")
print(f"  LEAST similar regions: {Rg[ii[o[0]]]}-{Rg[jj[o[0]]]} r={offdiag[o[0]]:.3f}")

# Quantify region-dependence: fraction of pairs whose across-region SD exceeds a threshold
# relative to a within-region noise floor (use range across chips within multi-chip regions as floor)
multi = [r for r in Rg if region_nchip[r] >= 3]
# build per-chip M50 to estimate within-region scatter for floor
chip_M50 = {}
for c, r, _ in order:
    f=f"{PERCHIP}/{c}.npz"
    if os.path.exists(f):
        d=np.load(f,allow_pickle=True); chip_M50[c]=d["C"].astype(np.float64)[0]
within_sd=[]
for r in multi:
    arr=np.stack([chip_M50[c][io,jo] for c in reg_chips[r] if c in chip_M50],axis=0)
    within_sd.append(arr.std(axis=0))
within_floor=np.median(np.concatenate(within_sd)) if within_sd else np.nan
across_sd_off = FP.std(axis=0)
n_var = int((across_sd_off > 2*within_floor).sum())
print(f"\nwithin-region SD floor (median, off-diag pairs)={within_floor:.4f}")
print(f"across-region SD of off-diag pairs: median={np.median(across_sd_off):.4f} max={across_sd_off.max():.4f}")
print(f"# off-diag pairs with across-region SD > 2x within-region floor: {n_var} / {len(io)} "
      f"({100*n_var/len(io):.1f}%)")

# ================= FIGURE =================
fig = plt.figure(figsize=(7.2, 8.6))
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.05, 1.0, 1.05],
              width_ratios=[1.0, 1.0], hspace=0.42, wspace=0.30,
              left=0.085, right=0.965, top=0.945, bottom=0.065)

# ---- Panel a: region x region similarity heatmap (clustered) + dendrogram ----
axA = fig.add_subplot(gs[0, 0])
Rord = Rmat[np.ix_(leaf_idx, leaf_idx)]
im = axA.imshow(Rord, cmap="RdYlBu_r", vmin=np.percentile(offdiag,2), vmax=1.0, aspect="equal")
axA.set_xticks(range(nreg)); axA.set_yticks(range(nreg))
axA.set_xticklabels(leaf_order, rotation=90, fontsize=5)
axA.set_yticklabels(leaf_order, fontsize=5)
axA.set_title("a  Region similarity of g(r=50 µm) fingerprint", loc="left", fontweight="bold")
cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.02)
cb.set_label("Pearson r (cross-program g)", fontsize=5); cb.ax.tick_params(labelsize=4)

# ---- Panel b: dendrogram of regions ----
axB = fig.add_subplot(gs[0, 1])
dendrogram(Z, labels=Rg, ax=axB, color_threshold=0.7*Z[:,2].max(),
           leaf_font_size=5, above_threshold_color="0.6")
axB.set_title("b  Region clustering by spatial program architecture", loc="left", fontweight="bold")
axB.set_ylabel("1 - r")
axB.spines[["top","right"]].set_visible(False)
for lbl in axB.get_xticklabels(): lbl.set_rotation(90)

# ---- Panel c: exemplar cross-program g(r) curves, one line per region (top divergent pairs) ----
# pick top 4 off-diagonal divergent pairs
ex_pairs = dfp_off.head(4)[["progA","progB"]].values.astype(int)
cmap_reg = plt.get_cmap("tab20")
reg_colors = {r: cmap_reg(i % 20) for i, r in enumerate(Rg)}
# 2x2 inside gs[1,:]
gsc = gs[1, :].subgridspec(2, 2, hspace=0.55, wspace=0.28)
for pi, (a, b) in enumerate(ex_pairs):
    ax = fig.add_subplot(gsc[pi // 2, pi % 2])
    ai, bi = a-1, b-1
    for ri, r in enumerate(Rg):
        curve = T[ri, :, ai, bi]
        lw = 1.4 if region_nchip[r] >= 3 else 0.7
        ls = "-" if region_nchip[r] >= 3 else ":"
        ax.plot(R_UM, curve, color=reg_colors[r], lw=lw, ls=ls,
                label=r if region_nchip[r] >= 3 else None)
    ax.set_title(f"{plab(ai)}\nx {plab(bi)}", fontsize=5.2, loc="left")
    ax.set_xlabel("r (µm)", fontsize=5); ax.set_ylabel("g(r)", fontsize=5)
    ax.tick_params(labelsize=4.5)
    ax.axhline(0, color="0.7", lw=0.4, zorder=0)
    ax.spines[["top","right"]].set_visible(False)
    if pi == 0:
        ax.legend(ncol=2, fontsize=4, frameon=False, loc="upper right",
                  handlelength=1.0, columnspacing=0.6, labelspacing=0.2)
# panel-c title (place above the 2x2 block)
fig.text(0.085, 0.635, "c  Exemplar cross-program g(r): one line per region "
         "(solid = ≥3 chips, dotted = ≤2 chips)",
         fontsize=7, fontweight="bold", ha="left")

# ---- Panel d: top region-variable pair matrix (SD of g(r=50) across regions, 60x60) ----
axD = fig.add_subplot(gs[2, 0])
SDmat = np.zeros((NPROG, NPROG))
SDmat[iu, ju] = np.sqrt(pair_var_50)
SDmat[ju, iu] = np.sqrt(pair_var_50)
imd = axD.imshow(SDmat, cmap="magma", aspect="equal",
                 vmax=np.percentile(np.sqrt(pair_var_50), 99))
axD.set_title("d  Across-region SD of g(r=50 µm), 60×60 cNMF pairs (incl. 6 excluded)", loc="left", fontweight="bold")
axD.set_xlabel("program"); axD.set_ylabel("program")
axD.set_xticks([0,14,29,44,59]); axD.set_xticklabels([1,15,30,45,60], fontsize=4.5)
axD.set_yticks([0,14,29,44,59]); axD.set_yticklabels([1,15,30,45,60], fontsize=4.5)
cbd = fig.colorbar(imd, ax=axD, fraction=0.046, pad=0.02)
cbd.set_label("SD across 14 regions", fontsize=5); cbd.ax.tick_params(labelsize=4)

# ---- Panel e: barh of top region-variable cross-program pairs ----
axE = fig.add_subplot(gs[2, 1])
topn = dfp_off.head(12).iloc[::-1]
ylab = [f"P{cnmf2new[int(a)]}×P{cnmf2new[int(b)]}\n{na[:16]} × {nb[:16]}"
        for a,b,na,nb in zip(topn.progA, topn.progB, topn.nameA, topn.nameB)]
yy = np.arange(len(topn))
axE.barh(yy, topn["range_g50"], color="#3b6ea5", height=0.7)
axE.set_yticks(yy); axE.set_yticklabels(ylab, fontsize=4.2)
axE.set_xlabel("max−min g(r=50 µm) across regions", fontsize=5)
axE.set_title("e  Top region-divergent cross-program pairs", loc="left", fontweight="bold")
axE.tick_params(labelsize=4.5)
axE.spines[["top","right"]].set_visible(False)
for i,(rg_min,rg_max) in enumerate(zip(topn.min_region, topn.max_region)):
    axE.text(topn["range_g50"].iloc[i]*1.01, i, f"{rg_min}↓{rg_max}↑",
             va="center", fontsize=3.6, color="0.3")
axE.set_xlim(0, topn["range_g50"].max()*1.28)

fig.savefig(f"{FIGDIR}/ed_spatial_crosscorr_byregion.pdf", dpi=300)
fig.savefig(f"{FIGDIR}/ed_spatial_crosscorr_byregion.png", dpi=220)
fig.savefig("/tmp/ed_spatial_crosscorr.png", dpi=220)
print(f"\nFIGURE saved: {FIGDIR}/ed_spatial_crosscorr_byregion.pdf/.png  size 7.2x8.6 in")
print("DONE_FIGURE")
