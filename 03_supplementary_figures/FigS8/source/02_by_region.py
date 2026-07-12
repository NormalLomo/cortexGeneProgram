#!/usr/bin/env python
"""
Aggregate per-chip lattice-ring program x program spatial mark cross-correlation into
per-REGION tensors, test region-dependence, and make an Extended-Data figure.

Inputs : raw 60-component tensor in results/crossregion_v1/spatial_crosscorr/_perchip/<chip>.npz
         (C:(10,60,60) cross-cov per ring, npairs, mu, sd)
         results/crossregion_v1/spatial_crosscorr/_chipmap.json
         results/crossregion_v1/program_names.tsv
Outputs: retained 54 by 54 outputs after the retained map subsets both program axes.
         results/crossregion_v1/spatial_crosscorr/region_<REGION>.npz   (per-region C tensor, npairs)
         results/crossregion_v1/spatial_crosscorr/region_tensor_all.npz (stacked raw tensor: regions x 10 x 60 x 60)
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
# Use a shared sans-serif font stack across renderers.
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

ROOT = os.environ.get("CORTEX_NMF_ROOT", "CORTEX_PROGRAM_ROOT")
CCDIR = f"{ROOT}/results/crossregion_v1/spatial_crosscorr"
PERCHIP = f"{CCDIR}/_perchip"
OUTDIR = os.environ.get("FIGS8_CACHE_DIR", CCDIR)
FIGDIR = os.environ.get("FIGS8_OUTPUT_DIR", f"{ROOT}/figures/extended")
REGION_TENSOR = os.environ.get("REGION_TENSOR", f"{CCDIR}/region_tensor_all.npz")
RETAIN_MAP = os.environ.get("RETAIN_MAP", f"{ROOT}/tables/TableS3_program_annotation.tsv")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)
NRING = 10; LAT = 50
R_UM = np.arange(1, NRING+1) * LAT  # 50..500

# ---- retained program map ----
# Map source cNMF component IDs to contiguous retained program IDs.
_rmap = pd.read_csv(RETAIN_MAP, sep="\t")
assert len(_rmap) == 54
_rmap["old_int"] = _rmap["cnmf_component"].astype(int)
_rmap["new_int"] = _rmap["new_P"].astype(str).str.removeprefix("P").astype(int)
assert _rmap["new_int"].tolist() == list(range(1, 55))
assert set(range(1, 61)) - set(_rmap["old_int"]) == {9, 18, 19, 35, 52, 57}
KEEP0 = (_rmap["old_int"].to_numpy() - 1).astype(int)
NPROG = len(KEEP0)

# program_names: keyed by new_P int -> name_short
short = dict(zip(_rmap["new_int"], _rmap["functional_name"].astype(str)))
conf  = dict(zip(_rmap["new_int"], _rmap["confidence"].astype(str)))

def plab(i):
    """Translate retained 0-based index to contiguous display identifier."""
    new_p = i + 1
    return f"P{new_p} {short.get(new_p, '?')}"

def cnmf0_to_newP(i):
    return i + 1

# ---- load compact canonical region tensor and subset both program axes ----
tensor = np.load(REGION_TENSOR, allow_pickle=True)
required_keys = {"T", "regions", "r_um", "nchip"}
assert required_keys.issubset(tensor.files), tensor.files
T60 = tensor["T"].astype(np.float64)
assert T60.ndim == 4 and T60.shape[1:] == (10, 60, 60), T60.shape
T = T60[:, :, KEEP0, :][:, :, :, KEEP0]
Rg = tensor["regions"].astype(str).tolist()
R_UM = tensor["r_um"].astype(int)
nchip_values = tensor["nchip"].astype(int)
assert T.shape == (len(Rg), 10, 54, 54), T.shape
region_nchip = dict(zip(Rg, nchip_values))
np.savez_compressed(f"{OUTDIR}/region_tensor_retained54.npz",
                    T=T.astype(np.float32), regions=np.array(Rg),
                    r_um=R_UM, nchip=nchip_values,
                    old_components=_rmap["old_int"].to_numpy(),
                    new_programs=_rmap["new_int"].to_numpy())
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
    "nameA": [short.get(int(a)+1, "?") for a in iu],
    "nameB": [short.get(int(b)+1, "?") for b in ju],
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
dfp_sorted = dfp.sort_values("var_g50", ascending=False).reset_index(drop=True)
dfp_sorted.to_csv(f"{OUTDIR}/top_region_variable_pairs.tsv", sep="\t", index=False)
print("\n=== TOP 15 region-variable program pairs (var of g(r=50um) across 14 regions) ===")
print(dfp_sorted.head(15)[["progA","nameA","progB","nameB","var_g50","range_g50","min_region","min_g50","max_region","max_g50"]].to_string(index=False))

# off-diagonal only top (cross-program, the interesting co-localization differences)
dfp_off = dfp_sorted[~dfp_sorted["is_self"]].reset_index(drop=True)
print("\n=== TOP 12 CROSS-program (off-diag) region-variable pairs ===")
print(dfp_off.head(12)[["progA","nameA","progB","nameB","var_g50","min_region","min_g50","max_region","max_g50"]].to_string(index=False))

# (b) region x region similarity by g(r=50um) fingerprint (off-diagonal flattened)
io, jo = np.triu_indices(NPROG, k=1)           # strictly off-diagonal
FP = M50[:, io, jo]                            # (nreg, n_offpair) fingerprint
assert FP.shape[1] == 1431, FP.shape
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
sim_df.to_csv(f"{OUTDIR}/region_similarity.tsv", sep="\t")
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

# The compact tensor has region aggregates but not per-chip matrices. The retained
# panels are recomputed exactly from those aggregates; the optional within-chip noise
# diagnostic is therefore explicitly unavailable in this tensor-only route.
within_floor = np.nan
across_sd_off = FP.std(axis=0)
print(f"\nwithin-region SD floor (median, off-diag pairs)={within_floor:.4f}")
print(f"across-region SD of off-diag pairs: median={np.median(across_sd_off):.4f} max={across_sd_off.max():.4f}")
print("within-chip threshold count: not computed from region-tensor-only input")

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
    ax.set_title(f"P{a} × P{b}", fontsize=5.2, loc="left")
    ax.set_xlabel("r (µm)", fontsize=5); ax.set_ylabel("g(r)", fontsize=5)
    ax.tick_params(labelsize=4.5)
    ax.axhline(0, color="0.7", lw=0.4, zorder=0)
    ax.spines[["top","right"]].set_visible(False)
    if pi == 0:
        ax.legend(ncol=2, fontsize=4, frameon=False, loc="upper right",
                  handlelength=1.0, columnspacing=0.6, labelspacing=0.2)
# panel-c title (place above the 2x2 block)
fig.text(0.085, 0.655, "c  Exemplar cross-program g(r): one line per region "
         "(solid = ≥3 chips, dotted = ≤2 chips)",
         fontsize=7, fontweight="bold", ha="left")

# ---- Panel d: retained-program SD matrix ----
axD = fig.add_subplot(gs[2, 0])
SDmat = np.zeros((NPROG, NPROG))
SDmat[iu, ju] = np.sqrt(pair_var_50)
SDmat[ju, iu] = np.sqrt(pair_var_50)
assert SDmat.shape == (54, 54)
imd = axD.imshow(SDmat, cmap="magma", aspect="equal",
                 vmax=np.percentile(np.sqrt(pair_var_50), 99))
axD.set_title("d  Across-region SD of g(r=50 µm), retained 54×54", loc="left", fontweight="bold")
axD.set_xlabel("program"); axD.set_ylabel("program")
ticks = [0, 13, 26, 39, 53]
axD.set_xticks(ticks); axD.set_xticklabels([1,14,27,40,54], fontsize=4.5)
axD.set_yticks(ticks); axD.set_yticklabels([1,14,27,40,54], fontsize=4.5)
cbd = fig.colorbar(imd, ax=axD, fraction=0.046, pad=0.02)
cbd.set_label("SD across 14 regions", fontsize=5); cbd.ax.tick_params(labelsize=4)

# ---- Panel e: barh of top region-variable cross-program pairs ----
axE = fig.add_subplot(gs[2, 1])
topn = dfp_off.head(12).iloc[::-1]
ylab = [f"P{int(a)}×P{int(b)}\n{na[:16]} × {nb[:16]}"
        for a,b,na,nb in zip(topn.progA, topn.progB, topn.nameA, topn.nameB)]
yy = np.arange(len(topn))
axE.barh(yy, topn["range_g50"], color="#3b6ea5", height=0.7)
axE.set_yticks(yy); axE.set_yticklabels(ylab, fontsize=4.2)
axE.set_xlabel("max−min g(r=50 µm) across regions", fontsize=5)
axE.set_title("e  Top region-divergent pairs", loc="left", fontweight="bold")
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
