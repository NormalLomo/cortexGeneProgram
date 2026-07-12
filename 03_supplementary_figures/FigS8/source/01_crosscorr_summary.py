#!/usr/bin/env python
"""
Renumber patch for FigS4 (ed8_crosscorr.py):
- Read program_renumber_map.tsv to build old_cNMF → new_P mapping
- Patch: plab(), short dict, Panel D ticks/title, Panel E ylab
- Excluded programs get [FLAG-EXCLUDED] marker in Panel E labels
- Layout / figure structure NOT changed (same GridSpec, same figsize)

RUN_LOG: scripts/extended/RUN_LOG_renumber_figS4.md
"""
import os
import json
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
    "font.size": 6, "axes.titlesize": 7, "axes.labelsize": 6,
    "xtick.labelsize": 5, "ytick.labelsize": 5, "legend.fontsize": 5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
})

ROOT = "CORTEX_PROGRAM_ROOT"
CCDIR = f"{ROOT}/results/crossregion_v1/spatial_crosscorr"
FIGDIR = f"{ROOT}/figures/extended"
NPROG = 60
NRING = 10

# ---- renumber map: old cNMF (1-60) → new P number ----
rmap_df = pd.read_csv(f"{ROOT}/results/crossregion_v1/program_renumber_map.tsv", sep="\t")
# old_P is int, new_P is "EXCLUDED" or int-as-string
old2new = {}   # old_int → new_int or "EXCLUDED"
for _, row in rmap_df.iterrows():
    old = int(row["old_P"])
    new = row["new_P"]
    if str(new).upper() == "EXCLUDED":
        old2new[old] = "EXCLUDED"
    else:
        old2new[old] = int(new)

EXCLUDED_OLD = {k for k, v in old2new.items() if v == "EXCLUDED"}  # {9,18,19,35,52,57}

# ---- program names: key by new_P number ----
pn = pd.read_csv(f"{ROOT}/results/crossregion_v1/program_names.tsv", sep="\t")
# new_P column is "P1".."P54" for kept; EXCLUDED rows use string "EXCLUDED"
# cnmf_component is the old cNMF index (matches program_renumber_map old_P)
short_by_old = {}   # old_cNMF_int → name_short  (for ALL rows incl. EXCLUDED)
for _, row in pn.iterrows():
    comp = int(row["cnmf_component"])
    short_by_old[comp] = str(row["name_short"])

# Also build short_by_new: new_P_int → name_short  (for kept programs)
short_by_new = {}
for old, new in old2new.items():
    if new != "EXCLUDED":
        name = short_by_old.get(old, "?")
        short_by_new[new] = name

def plab_old(old_0based):
    """Given 0-based OLD cNMF index, return label with NEW P number."""
    old1 = old_0based + 1   # 1-based old index
    new = old2new.get(old1, "?")
    name = short_by_old.get(old1, "?")
    if new == "EXCLUDED":
        return f"[FLAG-EXCLUDED-cNMF{old1}] {name}"
    return f"P{new} {name}"

# ---- region tensor ----
dt = np.load(f"{CCDIR}/region_tensor_all.npz", allow_pickle=True)
T = dt["T"].astype(np.float64)            # (nreg, NRING, 60, 60)
Rg = [str(x) for x in dt["regions"]]
R_UM = np.asarray(dt["r_um"], dtype=float)
nchip = {r: int(n) for r, n in zip(Rg, dt["nchip"])}
nreg = len(Rg)
well = [r for r in Rg if nchip[r] >= 3]
well_idx = [Rg.index(r) for r in well]
print("regions:", Rg)
print("well-replicated (>=3 chips):", well)

# ---- region similarity (panels a/b) : recompute from tensor exactly as original ----
io, jo = np.triu_indices(NPROG, k=1)
M50 = T[:, 0, :, :]
FP = M50[:, io, jo]
Rmat = np.corrcoef(FP)
offdiag = Rmat[np.triu_indices(nreg, 1)]
D = 1 - Rmat
np.fill_diagonal(D, 0.0)
D = (D + D.T) / 2
Z = linkage(squareform(D, checks=False), method="average")
dn = dendrogram(Z, labels=Rg, no_plot=True)
leaf_order = dn["ivl"]
leaf_idx = [Rg.index(l) for l in leaf_order]

# ---- panel d data : across-region SD of g(r=50) over all pairs ----
iu, ju = np.triu_indices(NPROG, k=0)
pair_vals_50 = M50[:, iu, ju]
pair_var_50 = pair_vals_50.var(axis=0)

# =====================================================================
# MAGNITUDE-REMOVED exemplar selection for panel c
# =====================================================================
ia, ib = np.triu_indices(NPROG, k=1)
Tw = T[np.ix_(well_idx, range(NRING), range(NPROG), range(NPROG))]
curves = Tw[:, :, ia, ib]
centered = curves - curves.mean(axis=1, keepdims=True)
shape_div = centered.var(axis=0).sum(axis=0)
mag = np.abs(curves).mean(axis=(0, 1))
rel_div = shape_div / (mag + 1e-6)

exdf = pd.DataFrame({
    "progA": ia + 1, "progB": ib + 1,   # 1-based OLD indices
    "shape_div": shape_div, "mag": mag, "rel_div": rel_div,
})
exdf_shape = exdf.sort_values("shape_div", ascending=False).reset_index(drop=True)
print("\n=== TOP 12 by MAGNITUDE-REMOVED shape divergence ===")
for _, row in exdf_shape.head(12).iterrows():
    a, b = int(row.progA), int(row.progB)
    print(f"  old cNMF {a}×{b}  → new {old2new.get(a,'?')}×{old2new.get(b,'?')}  "
          f"shape_div={row.shape_div:.4f}  mag={row.mag:.4f}")

chosen = []
used = set()
for _, row in exdf_shape.iterrows():
    a, b = int(row.progA), int(row.progB)
    if a in used or b in used:
        continue
    chosen.append((a, b))
    used.update([a, b])
    if len(chosen) == 3:
        break
print("\nCHOSEN exemplar pairs (old cNMF, panel c):", chosen)
print("  → new labels:", [(old2new.get(a,'?'), old2new.get(b,'?')) for a,b in chosen])

# Check for EXCLUDED in chosen
for a, b in chosen:
    if a in EXCLUDED_OLD or b in EXCLUDED_OLD:
        print(f"  [FLAG-EXCLUDED] Exemplar pair cNMF {a}×{b} contains EXCLUDED program")

# =====================================================================
# FIGURE  (4 rows: a/b | c | d/e | f-ARI)
# =====================================================================
# Load ARI data for panel f
_ari_json = f"{ROOT}/results/crossregion_v1/cluster_confusion/confusion_summary.json"
with open(_ari_json) as _fh:
    _ari_data = json.load(_fh)
_ari_sub = _ari_data["view3_subclass_x_region"]
_ari_dom = _ari_data["view4_domain_x_region"]

fig = plt.figure(figsize=(7.2, 11.0))
gs = GridSpec(4, 2, figure=fig, height_ratios=[1.05, 1.02, 1.05, 0.55],
              width_ratios=[1.0, 1.0], hspace=0.52, wspace=0.30,
              left=0.085, right=0.965, top=0.948, bottom=0.048)

# ---- Panel a ----
axA = fig.add_subplot(gs[0, 0])
Rord = Rmat[np.ix_(leaf_idx, leaf_idx)]
im = axA.imshow(Rord, cmap="RdYlBu_r", vmin=np.percentile(offdiag, 2), vmax=1.0, aspect="equal")
axA.set_xticks(range(nreg)); axA.set_yticks(range(nreg))
axA.set_xticklabels(leaf_order, rotation=90, fontsize=5)
axA.set_yticklabels(leaf_order, fontsize=5)
axA.set_title("a  Region similarity of g(r=50 µm) fingerprint", loc="left", fontweight="bold")
cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.02)
cb.set_label("Pearson r (cross-program g)", fontsize=5); cb.ax.tick_params(labelsize=4)

# ---- Panel b ----
axB = fig.add_subplot(gs[0, 1])
dendrogram(Z, labels=Rg, ax=axB, color_threshold=0.7 * Z[:, 2].max(),
           leaf_font_size=5, above_threshold_color="0.6")
axB.set_title("b  Region clustering by spatial program architecture", loc="left", fontweight="bold")
axB.set_ylabel("1 - r")
axB.spines[["top", "right"]].set_visible(False)
for lbl in axB.get_xticklabels(): lbl.set_rotation(90)

# ---- Panel c : mean-centered exemplar g(r) curves ----
cmap_reg = plt.get_cmap("tab20")
reg_colors = {r: cmap_reg(i % 20) for i, r in enumerate(Rg)}
gsc = gs[1, :].subgridspec(1, 3, wspace=0.40)
for pi, (a, b) in enumerate(chosen):
    ax = fig.add_subplot(gsc[0, pi])
    ai, bi = a - 1, b - 1   # 0-based
    for ri, r in enumerate(Rg):
        raw = T[ri, :, ai, bi]
        curve = raw - raw.mean()
        lw = 1.4 if nchip[r] >= 3 else 0.6
        ls = "-" if nchip[r] >= 3 else ":"
        z = 3 if nchip[r] >= 3 else 1
        ax.plot(R_UM, curve, color=reg_colors[r], lw=lw, ls=ls, zorder=z,
                label=r if nchip[r] >= 3 else None)
    # Use plab_old with 0-based index
    ax.set_title(f"{plab_old(ai)}\n× {plab_old(bi)}", fontsize=5.0, loc="left", pad=3)
    if pi != 1:
        ax.set_xlabel("r (µm)", fontsize=5)
    if pi == 0:
        ax.set_ylabel("g(r) − mean$_r$ g(r)", fontsize=5)
    ax.tick_params(labelsize=4.5)
    ax.axhline(0, color="0.7", lw=0.4, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

handles = [plt.Line2D([0], [0], color=reg_colors[r], lw=1.4) for r in well]
fig.legend(handles, well, ncol=len(well), fontsize=4.4, frameon=False,
           loc="lower center", bbox_to_anchor=(0.52, 0.385),
           handlelength=1.1, columnspacing=0.8, handletextpad=0.4,
           title="≥3-chip regions (solid; ≤2-chip dotted, unlabeled)",
           title_fontsize=4.6)
fig.text(0.085, 0.640,
         "c  Exemplar cross-program g(r), mean-centered per region",
         fontsize=7, fontweight="bold", ha="left")
fig.text(0.085, 0.628,
         "(shape / decay divergence; magnitude-removed pair ranking)",
         fontsize=5.6, fontweight="normal", ha="left", color="0.25")

# ---- Panel d ----
axD = fig.add_subplot(gs[2, 0])
SDmat = np.zeros((NPROG, NPROG))
SDmat[iu, ju] = np.sqrt(pair_var_50)
SDmat[ju, iu] = np.sqrt(pair_var_50)
imd = axD.imshow(SDmat, cmap="magma", aspect="equal",
                 vmax=np.percentile(np.sqrt(pair_var_50), 99))
axD.set_title("d  Across-region SD of g(r=50 µm), all 54×54 program pairs", loc="left", fontweight="bold")
axD.set_xlabel("program"); axD.set_ylabel("program")
# New tick labels: map old positions to new P numbers
# Positions 0,14,29,44,59 in old 0-based = old cNMF 1,15,30,45,60
# New P numbers: old1→new1, old15→new14, old30→new27, old45→new41, old60→new54
tick_old1 = [1, 15, 30, 45, 60]   # old cNMF 1-based
tick_new = [old2new.get(o, "?") for o in tick_old1]
tick_new_labels = [str(v) if v != "EXCLUDED" else "X" for v in tick_new]
axD.set_xticks([0, 14, 29, 44, 59])
axD.set_xticklabels(tick_new_labels, fontsize=4.5)
axD.set_yticks([0, 14, 29, 44, 59])
axD.set_yticklabels(tick_new_labels, fontsize=4.5)

posD0 = axD.get_position()
axD.set_position([0.070, posD0.y0, 0.300, posD0.height])
posD = axD.get_position()
caxD = fig.add_axes([posD.x1 + 0.012, posD.y0 + 0.012, 0.014, posD.height - 0.024])
cbd = fig.colorbar(imd, cax=caxD)
cbd.set_label("SD across 14 regions", fontsize=5); cbd.ax.tick_params(labelsize=4)

# ---- Panel e : top region-divergent cross-program pairs ----
dfp_full = pd.read_csv(f"{CCDIR}/top_region_variable_pairs.tsv", sep="\t")
dfp_off = dfp_full[~dfp_full["is_self"]].reset_index(drop=True)
axE = fig.add_subplot(gs[2, 1])
pos = axE.get_position()
axE.set_position([pos.x0 + 0.085, pos.y0, pos.width - 0.085, pos.height])
# Re-tally: remove any pair where either end is an excluded program (cNMF 9/18/19/35/52/57)
dfp_clean = dfp_off[
    (~dfp_off["progA"].isin(EXCLUDED_OLD)) & (~dfp_off["progB"].isin(EXCLUDED_OLD))
].reset_index(drop=True)
print(f"Panel E re-tally: {len(dfp_off)} off-diag pairs → {len(dfp_clean)} after removing excluded; taking top 10.")
topn = dfp_clean.head(10).iloc[::-1]

# Build y-labels with NEW P numbers (no excluded → no FLAG)
# Use name_short from program_names.tsv instead of the truncated raw GO term
# names in the nameA/nameB columns.
ylab = []
for _, row in topn.iterrows():
    a_old = int(row["progA"])
    b_old = int(row["progB"])
    na = short_by_old.get(a_old, str(row["nameA"]))  # functional name_short
    nb = short_by_old.get(b_old, str(row["nameB"]))  # functional name_short
    a_new = old2new.get(a_old, "?")
    b_new = old2new.get(b_old, "?")
    a_str = f"P{a_new}" if a_new not in ("EXCLUDED", "?") else f"[ERR-{a_old}]"
    b_str = f"P{b_new}" if b_new not in ("EXCLUDED", "?") else f"[ERR-{b_old}]"
    ylab.append(f"{a_str}×{b_str}  {na} × {nb}")
print("Panel E labels:", ylab[::-1])

yy = np.arange(len(topn))
axE.barh(yy, topn["range_g50"], color="#3b6ea5", height=0.66)
axE.set_yticks(yy); axE.set_yticklabels(ylab, fontsize=5.0)
axE.set_xlabel("max−min g(r=50 µm) across regions", fontsize=5)
axE.set_title("e  Top region-divergent cross-program pairs", loc="left", fontweight="bold")
axE.tick_params(axis="x", labelsize=4.5)
axE.tick_params(axis="y", length=0)
axE.spines[["top", "right"]].set_visible(False)
for i, (rg_min, rg_max) in enumerate(zip(topn.min_region, topn.max_region)):
    axE.text(topn["range_g50"].iloc[i] * 1.01, i, f"{rg_min}↓ {rg_max}↑",
             va="center", fontsize=4.0, color="0.3")
axE.set_xlim(0, topn["range_g50"].max() * 1.30)

# ---- Panel f : ARI cluster-confusion summary (subclass vs region) ----
axF = fig.add_subplot(gs[3, :])
_groups = [
    f"Subclass × region\n(n={_ari_sub['n_rows']}, k={_ari_sub['k']})",
    f"Domain × region\n(n={_ari_dom['n_rows']}, k={_ari_dom['k']})",
]
_ari_id  = [_ari_sub["ARI_identity"],  _ari_dom["ARI_identity"]]
_ari_reg = [_ari_sub["ARI_region"],    _ari_dom["ARI_region"]]
_yy = np.arange(len(_groups))
_hbar = 0.45
_bar_id  = axF.barh(_yy, _ari_id,  _hbar, color="#2c7fb8", label="ARI (identity)")
_bar_reg = axF.barh(_yy + _hbar, _ari_reg, _hbar, color="#d95f0e", label="ARI (region)")
# annotations on bars
for _i, (_vi, _vr) in enumerate(zip(_ari_id, _ari_reg)):
    _xi = max(_vi, 0.02)
    axF.text(_xi + 0.01, _yy[_i], f"{_vi:.2f}", va="center", fontsize=5.5, color="#2c7fb8")
    _xr = max(_vr, 0.0)
    _txt_x = abs(_vr) + 0.01 if _vr < 0 else _vr + 0.01
    axF.text(_txt_x + 0.01, _yy[_i] + _hbar, f"{_vr:.2f}", va="center", fontsize=5.5, color="#d95f0e")
axF.set_yticks(_yy + _hbar / 2)
axF.set_yticklabels(_groups, fontsize=5.5)
axF.set_xlabel("Adjusted Rand Index (ARI)", fontsize=5.5)
axF.axvline(0, color="k", lw=0.5, ls="--")
axF.set_xlim(-0.15, 1.05)
axF.legend(fontsize=5, loc="lower right", frameon=False, ncol=2,
           handlelength=1.0, columnspacing=1.0)
axF.set_title("f  Cluster confusion: cell identity vs. brain region", loc="left", fontweight="bold")
axF.tick_params(labelsize=5)
axF.spines[["top", "right"]].set_visible(False)
print(f"Panel f: ARI subclass/identity={_ari_sub['ARI_identity']:.3f}  ARI subclass/region={_ari_sub['ARI_region']:.3f}")

fig.savefig(f"{FIGDIR}/ed_spatial_crosscorr_byregion.pdf", dpi=300)
fig.savefig(f"{FIGDIR}/ed_spatial_crosscorr_byregion.svg")
fig.savefig(f"{FIGDIR}/ed_spatial_crosscorr_byregion.png", dpi=220)
fig.savefig("/tmp/ed8_crosscorr_renum.png", dpi=220)
print(f"\nFIGURE saved: {FIGDIR}/ed_spatial_crosscorr_byregion.pdf/.png  size 7.2x11.0 in (6 panels a-f)")
print("DONE_FIGURE_RENUMBERED")
