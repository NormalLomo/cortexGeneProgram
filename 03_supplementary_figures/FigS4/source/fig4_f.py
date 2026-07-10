#!/usr/bin/env python
# Panel F: UMAP of L3-L4 IT RORB cells, colored by region & by p14 activity
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

RES = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
FIG = "CORTEX_PROGRAM_ROOT/figures/fig4"
plt.rcParams.update({
    "pdf.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.5,
})

# program -> functional label lookup (P{n} {name_short}, "*" if brain-weak)
def _load_prog_labels():
    nm = pd.read_csv(os.path.join(RES, "program_names.tsv"), sep="\t")
    lab = {}
    for _, r in nm.iterrows():
        ns = str(r["name_short"]).strip()
        ns = pd.Series([ns]).str.replace(r"\s+P\d+\s*$", "", regex=True).iloc[0]
        star = "*" if str(r["confidence"]) == "brain-weak" else ""
        lab[int(r["program"])] = f"P{int(r['program'])} {ns}{star}"
    return lab
PROG_LABEL = _load_prog_labels()
def prog_label(p): return PROG_LABEL.get(int(p), f"P{int(p)}")
P14_LAB = prog_label(14)

d = pd.read_csv(os.path.join(RES, "panel_f_umap.tsv"), sep="\t")
# shuffle for even overplot
d = d.sample(frac=1.0, random_state=0).reset_index(drop=True)
regions = sorted(d["region"].unique())
# qualitative palette (distinct, ~14)
cmap = plt.get_cmap("tab20")
rcol = {r: cmap(i % 20) for i, r in enumerate(regions)}

fig = plt.figure(figsize=(7.4, 3.6))
gs = GridSpec(1, 2, figure=fig, wspace=0.28, left=0.04, right=0.93, top=0.86, bottom=0.06)

# --- (i) by region ---
ax0 = fig.add_subplot(gs[0, 0])
for r in regions:
    s = d[d["region"] == r]
    ax0.scatter(s["UMAP1"], s["UMAP2"], s=1.0, c=[rcol[r]], linewidths=0, alpha=0.6, rasterized=True)
ax0.set_title("L3-L4 IT RORB — by region", fontsize=9, loc="left", weight="bold")
ax0.set_xticks([]); ax0.set_yticks([])
for sp in ax0.spines.values(): sp.set_visible(False)
ax0.set_xlabel("UMAP1", fontsize=7); ax0.set_ylabel("UMAP2", fontsize=7)
handles = [Line2D([0],[0], marker='o', linestyle='', markersize=4.5,
                  markerfacecolor=rcol[r], markeredgewidth=0, label=r) for r in regions]
ax0.legend(handles=handles, fontsize=6.5, ncol=2, loc="upper left",
           bbox_to_anchor=(0.0, 1.0), frameon=False, handletextpad=0.3,
           columnspacing=0.6, labelspacing=0.35)

# --- (ii) by p14 ---
ax1 = fig.add_subplot(gs[0, 1])
vmin, vmax = np.percentile(d["p14"], [2, 98])
sc = ax1.scatter(d["UMAP1"], d["UMAP2"], s=1.0, c=d["p14"], cmap="magma",
                 vmin=vmin, vmax=vmax, linewidths=0, alpha=0.7, rasterized=True)
ax1.set_title(f"{P14_LAB} activity", fontsize=9, loc="left", weight="bold")
ax1.set_xticks([]); ax1.set_yticks([])
for sp in ax1.spines.values(): sp.set_visible(False)
ax1.set_xlabel("UMAP1", fontsize=7); ax1.set_ylabel("UMAP2", fontsize=7)
cb = fig.colorbar(sc, ax=ax1, fraction=0.045, pad=0.02)
cb.set_label(P14_LAB, fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
cb.outline.set_linewidth(0.3)

fig.suptitle("Region-structured program gradient within one cell type",
             fontsize=9.5, x=0.04, ha="left", weight="bold", y=0.985)
fig.savefig(os.path.join(FIG, "fig4_f.pdf"), dpi=300, bbox_inches="tight")
print("panel f done")
