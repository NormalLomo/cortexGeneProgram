#!/usr/bin/env python
"""
Extended Data Fig 3 figure : cNMF K-robustness (native matplotlib, vector PDF).

Panels:
 (a) Spearman concordance of per-program cross-region variability (eta2) ranking,
     each alt-K vs canonical K=60 (all matched pairs + confident-match subset).
 (b) Scatter eta2_Kx vs eta2_K60 for gene-spectra-matched programs, one facet per K.
 (c) Top-driver class recurrence: fraction of top-6 cell-driver subclasses that are
     deep-EX / glia, per K (incl K=60).
 (d) Identity of the single most cross-region-variable program at each K
     (is it laminar-IT / RORB-like?) — annotated grid.

Reads results/crossregion_v1/k_robustness/*. Writes figures/extended/ed_fig3_k_robustness.{pdf,png}
"""
import os
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
from matplotlib.patches import Rectangle
import matplotlib.font_manager as fm

BASE = "CORTEX_PROGRAM_ROOT"
OUTD = f"{BASE}/results/crossregion_v1/k_robustness"
FIGD = f"{BASE}/figures/extended"
os.makedirs(FIGD, exist_ok=True)

# ---- global style: embedded TrueType, >=5pt, journal-ish ----
plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "figure.dpi": 200,
})

ALTK = [40, 50, 65, 70]
KCOL = {40: "#4C72B0", 50: "#55A868", 65: "#C44E52", 70: "#8172B3"}
DEEP = "#2166AC"; GLIA = "#B2182B"; OTHER = "#D9D9D9"

conc = pd.read_csv(f"{OUTD}/concordance_summary.tsv", sep="\t")
drv  = pd.read_csv(f"{OUTD}/topdriver_recurrence.tsv", sep="\t")
tv   = pd.read_csv(f"{OUTD}/top_variable_program_identity.tsv", sep="\t")
scat = {K: pd.read_csv(f"{OUTD}/scatter_eta_k{K}_vs_k60.tsv", sep="\t") for K in ALTK}

# ============================ figure layout ============================
# 2 rows: top = (a)+(c)+(d); bottom = (b) 4 scatter facets
fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 12, height_ratios=[1.0, 1.05], hspace=0.55, wspace=2.4,
                      left=0.075, right=0.965, top=0.955, bottom=0.085)  # editfig_r1: right pad for panel d + reclaim suptitle band

# ---- panel (a): Spearman bars ----
axa = fig.add_subplot(gs[0, 0:5])
x = np.arange(len(ALTK)); w = 0.38
sa = [conc.loc[conc.K == K, "spearman_all"].values[0] for K in ALTK]
sc = [conc.loc[conc.K == K, "spearman_conf"].values[0] for K in ALTK]
b1 = axa.bar(x - w/2, sa, w, color=[KCOL[k] for k in ALTK], edgecolor="black", linewidth=0.5, label="all matched")
b2 = axa.bar(x + w/2, sc, w, color=[KCOL[k] for k in ALTK], edgecolor="black", linewidth=0.5,
             alpha=0.55, hatch="///", label="conf. match (r>0.5)")
axa.axhline(0, color="0.4", lw=0.6)
axa.set_xticks(x); axa.set_xticklabels([f"K{k}" for k in ALTK])
axa.set_ylabel("Spearman ρ  (η² rank vs K60)")
axa.set_ylim(0, 1.0)
axa.set_title("a  Variability-ranking concordance", loc="left", fontweight="bold")
for xi, val in zip(x - w/2, sa):
    if not np.isnan(val):
        axa.text(xi, val + 0.02, f"{val:.2f}", ha="center", va="bottom", fontsize=5.6)
axa.legend(loc="lower right", handlelength=1.2)

# ---- panel (c): top-driver class recurrence (stacked frac of top-6) ----
axc = fig.add_subplot(gs[0, 5:9])
ks = [60] + ALTK
xc = np.arange(len(ks))
nd = [drv.loc[drv.K == k, "n_deepEX"].values[0] / 6 for k in ks]
ng = [drv.loc[drv.K == k, "n_glia"].values[0] / 6 for k in ks]
no = [1 - a - b for a, b in zip(nd, ng)]
axc.bar(xc, nd, 0.62, color=DEEP, edgecolor="black", linewidth=0.5, label="deep-EX")
axc.bar(xc, ng, 0.62, bottom=nd, color=GLIA, edgecolor="black", linewidth=0.5, label="glia")
axc.bar(xc, no, 0.62, bottom=[a+b for a,b in zip(nd,ng)], color=OTHER, edgecolor="black", linewidth=0.5, label="other")
axc.set_xticks(xc); axc.set_xticklabels(["K60*"] + [f"K{k}" for k in ALTK])
axc.set_ylabel("fraction of top-6 drivers")
axc.set_ylim(0, 1.0)
axc.set_title("c  Top cell-driver recurrence", loc="left", fontweight="bold")
axc.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, handlelength=1.0, columnspacing=0.9)

# ---- panel (d): top-variable program identity grid ----
axd = fig.add_subplot(gs[0, 9:12])
axd.set_title("d  Most-variable program", loc="left", fontweight="bold", fontsize=7.0)
axd.axis("off")
ks = [60] + ALTK
nrow = len(ks)
axd.set_xlim(0, 1); axd.set_ylim(0, nrow + 0.6)
# header
axd.text(0.02, nrow + 0.15, "K", fontsize=6.5, fontweight="bold")
axd.text(0.24, nrow + 0.15, "η²", fontsize=6.5, fontweight="bold")
axd.text(0.45, nrow + 0.15, "assoc. subclass", fontsize=6.0, fontweight="bold")  # editfig_r1 current: re-spaced
for i, k in enumerate(ks):
    row = tv.loc[tv.K == k].iloc[0]
    y = nrow - 1 - i + 0.5
    is_lam = bool(row["is_laminar_IT"])
    col = "#2E7D32" if is_lam else "#9E9E9E"
    axd.add_patch(Rectangle((0.0, y - 0.34), 1.0, 0.68, facecolor=col, alpha=0.16,
                            edgecolor="none", zorder=0))
    klab = "K60*" if k == 60 else f"K{k}"
    axd.text(0.02, y, klab, fontsize=6.3, va="center")
    axd.text(0.24, y, f"{row['top_var_eta2']:.3f}", fontsize=6.3, va="center")
    sc_txt = str(row["assoc_subclass"])
    mark = " ✓" if is_lam else ""
    axd.text(0.45, y, sc_txt + mark, fontsize=5.2, va="center",
             color="#1B5E20" if is_lam else "0.25")  # editfig_r1 current: re-spaced + shrink
axd.text(0.0, -0.15, "✓ = laminar-IT / RORB-like", fontsize=5.4, color="#1B5E20")

# ---- panel (b): 4 scatter facets eta2_Kx vs eta2_K60 ----
# determine common axis range
allv = np.concatenate([np.r_[scat[K]["eta_K"].values, scat[K]["eta_60"].values] for K in ALTK])
hi = np.nanpercentile(allv, 99.5) * 1.08
for j, K in enumerate(ALTK):
    axb = fig.add_subplot(gs[1, j*3:(j+1)*3])
    s = scat[K]
    cc = s["corr"].clip(0, 1).values
    sca = axb.scatter(s["eta_60"], s["eta_K"], c=cc, cmap="viridis", vmin=0.3, vmax=1.0,
                      s=14, edgecolor="black", linewidth=0.25, zorder=3)
    axb.plot([0, hi], [0, hi], ls="--", color="0.5", lw=0.7, zorder=1)
    rho = conc.loc[conc.K == K, "spearman_all"].values[0]
    pear = conc.loc[conc.K == K, "pearson_all"].values[0]
    axb.set_xlim(0, hi); axb.set_ylim(0, hi)
    axb.set_aspect("equal", "box")
    axb.set_title(f"K{K} vs K60", loc="left", fontsize=7, color=KCOL[K], fontweight="bold")
    axb.text(0.04, 0.93, f"ρ={rho:.2f}\nr={pear:.2f}", transform=axb.transAxes,
             fontsize=5.8, va="top", ha="left")
    axb.set_xlabel("η²  (K60)")
    if j == 0:
        axb.set_ylabel("η²  (alt-K, matched)")
    axb.tick_params(labelsize=5.8)
# panel-b letter + shared colorbar
fig.text(0.075, 0.49, "b  Per-program variability (η²): alt-K vs K60, gene-spectra matched",
         fontsize=8, fontweight="bold", ha="left")
cax = fig.add_axes([0.40, 0.018, 0.22, 0.013])
cb = fig.colorbar(sca, cax=cax, orientation="horizontal")
cb.set_label("gene-spectra match corr", fontsize=5.6)
cb.ax.tick_params(labelsize=5.2)

# The figure-level title is supplied by the public figure legend.

pdf = f"{FIGD}/ed_fig3_k_robustness.pdf"
png = f"{FIGD}/ed_fig3_k_robustness.png"
fig.savefig(pdf)
fig.savefig(png, dpi=300)
fig.savefig(pdf.replace(".pdf",".svg"))  # editfig_r1: svg for SUBMISSION sync
print("saved:", pdf, png)
print("figsize_in:", fig.get_size_inches(), "AR(W/H):", round(fig.get_size_inches()[0]/fig.get_size_inches()[1], 3))
