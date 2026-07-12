#!/usr/bin/env python
"""
Render the retained-program FigS5 composite from cached producer outputs.

The renderer derives deterministic display orderings from the cached means and
writes vector and raster figure outputs.
"""
import os, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# Use a portable sans-serif font stack for vector and raster output.
import matplotlib as _mpl_font
_mpl_font.rcParams["font.family"] = "sans-serif"
_mpl_font.rcParams["font.sans-serif"] = ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"]
_mpl_font.rcParams["pdf.fonttype"] = 42
_mpl_font.rcParams["ps.fonttype"] = 42
_mpl_font.rcParams["svg.fonttype"] = "none"
_mpl_font.rcParams["mathtext.fontset"] = "dejavusans"  # sans math
_mpl_font.rcParams["mathtext.default"] = "regular"  # math uses font.family
from matplotlib.gridspec import GridSpec
from matplotlib.colors import to_rgba
import matplotlib as mpl
from svgutils import compose as sc
import cairosvg

ROOT = os.environ.get("CORTEX_NMF_ROOT", "CORTEX_PROGRAM_ROOT")
RES  = f"{ROOT}/results/crossregion_v1"
OUTD = os.environ.get("FIGS5_CACHE_DIR", f"{RES}/cluster_confusion")
FIGD = os.environ.get("FIGS5_OUTPUT_DIR", f"{ROOT}/figures/extended")

RETAIN_MAP = os.environ.get("RETAIN_MAP", f"{ROOT}/tables/TableS3_program_annotation.tsv")
MAP = pd.read_csv(RETAIN_MAP, sep="\t")
assert len(MAP) == 54
MAP["new_int"] = MAP["new_P"].astype(str).str.removeprefix("P").astype(int)
MAP["old_int"] = MAP["cnmf_component"].astype(int)
assert MAP["new_int"].tolist() == list(range(1, 55))
assert set(range(1, 61)) - set(MAP["old_int"]) == {9, 18, 19, 35, 52, 57}
PROG = MAP["new_int"].astype(str).tolist()

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif", "font.sans-serif": ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6, "axes.linewidth": 0.4, "xtick.major.width": 0.4,
    "ytick.major.width": 0.4, "xtick.major.size": 2, "ytick.major.size": 2,
})

# ---------------- program short names ----------------
nm = MAP.copy()
nm["new_P"] = nm["new_int"].astype(str)
short = dict(zip(nm["new_P"], nm["functional_name"]))
conf  = dict(zip(nm["new_P"], nm["confidence"]))
def plabel(p):
    s = short.get(p, p)
    star = " *" if str(conf.get(p, "")).endswith("weak") else ""
    return f"P{p} {s}{star}"
PROG_LABELS = [plabel(p) for p in PROG]

# ---------------- helpers ----------------
def corr_ward(M):
    D = pdist(M, metric="correlation")
    D = np.nan_to_num(D, nan=1.0)
    Z = linkage(D, method="ward")
    order = dendrogram(Z, no_plot=True)["leaves"]
    return Z, order

def zscore_cols(M):
    mu = M.mean(0); sd = M.std(0); sd[sd == 0] = 1
    return (M - mu) / sd

def draw_left_dendrogram(ax, Z, n_leaves):
    dd = dendrogram(Z, orientation="left", ax=ax, no_labels=True,
                    color_threshold=0, link_color_func=lambda k: "0.35")
    ycoords = [y for seg in dd.get("icoord", []) for y in seg]
    ymax = (max(ycoords) + 5) if ycoords else (10 * n_leaves)
    ax.set_ylim(ymax, 0)
    return dd

def compose_final_svg(fig, panel_svg, final_svg, final_pdf, final_png):
    width_pt = 7.2 * 72
    height_pt = 7.8 * 72
    fig.savefig(panel_svg)
    panel = sc.MplFigure(fig)
    fig_svg = sc.Figure(f"{width_pt}pt", f"{height_pt}pt", panel)
    fig_svg.save(str(final_svg))
    svg_bytes = fig_svg.tostr()
    if isinstance(svg_bytes, str):
        svg_bytes = svg_bytes.encode("utf-8")
    cairosvg.svg2pdf(bytestring=svg_bytes, write_to=str(final_pdf))
    cairosvg.svg2png(bytestring=svg_bytes, write_to=str(final_png),
                     output_width=1440, output_height=1560)

# ================= LOAD CACHED TABLES (no recompute) =================
sub_mean = pd.read_csv(f"{OUTD}/view1_subclass_mean.tsv", sep="\t", index_col=0)
sub_mean.columns = [str(c) for c in sub_mean.columns]; sub_mean = sub_mean[PROG]

dom_mean = pd.read_csv(f"{OUTD}/view2_domain_mean.tsv", sep="\t", index_col=0)
dom_mean.columns = [str(c) for c in dom_mean.columns]; dom_mean = dom_mean[PROG]

sr_mean = pd.read_csv(f"{OUTD}/view3_subclass_region_mean.tsv", sep="\t",
                      index_col=[0, 1])
sr_mean.columns = [str(c) for c in sr_mean.columns]; sr_mean = sr_mean[PROG]

dr_mean = pd.read_csv(f"{OUTD}/view4_domain_region_mean.tsv", sep="\t",
                      index_col=[0, 1])
dr_mean.columns = [str(c) for c in dr_mean.columns]; dr_mean = dr_mean[PROG]

with open(f"{OUTD}/confusion_summary.json") as f:
    summary = json.load(f)
res3 = summary["view3_subclass_x_region"]
res4 = summary["view4_domain_x_region"]
vp3p = summary["varpart_subclass_region_mean_pct"]
vp4p = summary["varpart_domain_region_mean_pct"]

# dendrograms / orders (display only; deterministic, same metric as producer)
Z1, ord1 = corr_ward(zscore_cols(sub_mean.values))
Z2, ord2 = corr_ward(zscore_cols(dom_mean.values))
Z3, ord3 = corr_ward(zscore_cols(sr_mean.values))
Z4, ord4 = corr_ward(zscore_cols(dr_mean.values))

# ================= FIGURE =================
def cat_colors(cats, cmap):
    cats = list(dict.fromkeys(cats))
    cm = plt.get_cmap(cmap, len(cats))
    return {c: cm(i) for i, c in enumerate(cats)}

all_sub = sorted(set(sub_mean.index) | set(sr_mean.index.get_level_values(0)))
all_dom = sorted(set(dom_mean.index) | set(dr_mean.index.get_level_values(0)))
all_reg = sorted(set(sr_mean.index.get_level_values(1)) | set(dr_mean.index.get_level_values(1)))
SUB_COL = cat_colors(all_sub, "tab20")
DOM_COL = cat_colors(all_dom, "viridis")
REG_COL = cat_colors(all_reg, "tab20b")
CMAP = "RdBu_r"

fig = plt.figure(figsize=(7.2, 7.8))
# top (a,b) shorter; mid (c,d) tall; bottom (e) close behind -> no dead band
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.08, 1.58, 0.46],
              hspace=0.34, wspace=0.42,
              left=0.165, right=0.965, top=0.955, bottom=0.085)

def add_cbar(ax, im, label=None):
    cb = fig.colorbar(im, ax=ax, fraction=0.020, pad=0.012, aspect=14)
    cb.ax.tick_params(labelsize=5, length=1.5)
    cb.outline.set_linewidth(0.3)
    if label:
        cb.set_label(label, fontsize=5)
    return cb

def panel_dendro_heat(subgs, M_df, Z, order, title, vmax=2.5, every=2):
    """Left dendrogram + heatmap (aspect='auto' so rows fill cell and align w/ dendro)."""
    Mz = zscore_cols(M_df.values)
    _, cord = corr_ward(Mz.T)
    Mz = Mz[np.ix_(order, cord)]
    rlabs = [str(M_df.index[i]) for i in order]
    inner = subgs.subgridspec(1, 2, width_ratios=[0.20, 1.0], wspace=0.035)
    axd = fig.add_subplot(inner[0, 0])
    draw_left_dendrogram(axd, Z, len(rlabs))
    axd.axis("off")
    axh = fig.add_subplot(inner[0, 1])
    im = axh.imshow(Mz, aspect="auto", cmap=CMAP, vmin=-vmax, vmax=vmax,
                    interpolation="nearest")
    axh.set_yticks(range(len(rlabs)))
    axh.set_yticklabels(rlabs, fontsize=5)
    # thin program labels to every-other to avoid smear; >=5pt
    cidx = list(range(0, len(cord), every))
    axh.set_xticks(cidx)
    axh.set_xticklabels([f"P{cord[j]+1}" for j in cidx], fontsize=5, rotation=90)
    axh.tick_params(length=1, pad=1.2)
    axh.set_title(title, fontsize=6.5, pad=4, loc="left", fontweight="bold")
    add_cbar(axh, im, "col z")
    return axh

panel_dendro_heat(gs[0, 0], sub_mean, Z1, ord1,
                  "a  Subclass (22) → program profile", every=3)
panel_dendro_heat(gs[0, 1], dom_mean, Z2, ord2,
                  "b  Spatial domain (8) → program profile", every=3)

def panel_confusion(subgs, M_df, Z, order, res, ident_lut, title, letter):
    Mz = zscore_cols(M_df.values)
    _, cord = corr_ward(Mz.T)
    Mz = Mz[np.ix_(order, cord)]
    ident = M_df.index.get_level_values(0).values[order]
    region = M_df.index.get_level_values(1).values[order]
    inner = subgs.subgridspec(1, 4, width_ratios=[0.16, 0.035, 0.035, 1.0], wspace=0.03)
    axd = fig.add_subplot(inner[0, 0])
    draw_left_dendrogram(axd, Z, len(ident))
    axd.axis("off")
    axs1 = fig.add_subplot(inner[0, 1])
    axs1.imshow([[to_rgba(ident_lut[ident[i]])] for i in range(len(ident))],
                aspect="auto"); axs1.set_xticks([]); axs1.set_yticks([])
    axs1.set_xlabel("id", fontsize=5, rotation=90, labelpad=1)
    axs2 = fig.add_subplot(inner[0, 2])
    axs2.imshow([[to_rgba(REG_COL[region[i]])] for i in range(len(region))],
                aspect="auto"); axs2.set_xticks([]); axs2.set_yticks([])
    axs2.set_xlabel("reg", fontsize=5, rotation=90, labelpad=1)
    axh = fig.add_subplot(inner[0, 3])
    im = axh.imshow(Mz, aspect="auto", cmap=CMAP, vmin=-2.5, vmax=2.5,
                    interpolation="nearest")
    ticks = list(range(0, len(PROG), 5))
    axh.set_yticks([]); axh.set_xticks(ticks)
    axh.set_xticklabels([f"P{cord[j]+1}" for j in ticks], fontsize=5,
                        rotation=90)
    axh.tick_params(length=1, pad=1.2)
    ann = (f"ARI(identity)={res['ARI_identity']:.2f}   ARI(region)={res['ARI_region']:.2f}\n"
           f"NMI(id)={res['NMI_identity']:.2f}   NMI(reg)={res['NMI_region']:.2f}")
    axh.set_title(f"{letter}  {title}\n{ann}", fontsize=6, pad=4, loc="left",
                  fontweight="bold")
    add_cbar(axh, im)

panel_confusion(gs[1, 0], sr_mean, Z3, ord3, res3, SUB_COL,
                f"Subclass × region ({res3['n_rows']} groups)", "c")
panel_confusion(gs[1, 1], dr_mean, Z4, ord4, res4, DOM_COL,
                f"Domain × region ({res4['n_rows']} groups)", "d")

# ---- panel e: variance partition bars ----
axe = fig.add_subplot(gs[2, :])
groups = ["Subclass × region", "Domain × region"]
ident_pct = [vp3p["identity"], vp4p["identity"]]
reg_pct   = [vp3p["region"],   vp4p["region"]]
res_pct   = [vp3p["residual"], vp4p["residual"]]
y = np.arange(len(groups)); hbar = 0.55
axe.barh(y, ident_pct, hbar, color="#2c7fb8", label="Identity (subclass/domain)")
axe.barh(y, reg_pct, hbar, left=ident_pct, color="#d95f0e", label="Region")
axe.barh(y, res_pct, hbar, left=np.array(ident_pct)+np.array(reg_pct), color="0.8",
         label="Residual")
for i in range(len(groups)):
    axe.text(ident_pct[i]/2, y[i], f"{ident_pct[i]:.0f}%", ha="center", va="center",
             fontsize=6, color="white", fontweight="bold")
    rx = ident_pct[i] + reg_pct[i]/2
    # region slice can be thin -> place its % just right of the slice if too narrow
    if reg_pct[i] >= 6:
        axe.text(rx, y[i], f"{reg_pct[i]:.0f}%", ha="center", va="center",
                 fontsize=6, color="white", fontweight="bold")
    else:
        axe.text(ident_pct[i]+reg_pct[i]+0.8, y[i], f"{reg_pct[i]:.0f}%", ha="left",
                 va="center", fontsize=6, color="#d95f0e", fontweight="bold")
axe.set_yticks(y); axe.set_yticklabels(groups, fontsize=6)
axe.set_xlabel("mean variance explained across 54 retained programs (%)", fontsize=6)
axe.set_xlim(0, 100)
axe.set_title("e  Variance partition: identity vs region", fontsize=6.5, loc="left",
              fontweight="bold", pad=4)
axe.legend(fontsize=5, loc="lower right", frameon=False, ncol=3,
           handlelength=1.1, columnspacing=1.0)
axe.tick_params(labelsize=6, length=2)
axe.spines[["top", "right"]].set_visible(False)
axe.margins(y=0.18)

Path(FIGD).mkdir(parents=True, exist_ok=True)
panel_svg = Path(FIGD) / "ed_cluster_confusion.panel.svg"
final_svg = Path(FIGD) / "ed_cluster_confusion.svg"
final_pdf = Path(FIGD) / "ed_cluster_confusion.pdf"
final_png = Path(FIGD) / "ed_cluster_confusion.png"
compose_final_svg(fig, panel_svg, final_svg, final_pdf, final_png)
print("WROTE", f"{final_svg}|{final_pdf}|{final_png} via svgutils")
print("KEY NUMS  ARI_sub_ident=%.3f ARI_sub_region=%.3f | ARI_dom_ident=%.3f ARI_dom_region=%.3f | var_sub id=%.1f reg=%.1f"
      % (res3["ARI_identity"], res3["ARI_region"], res4["ARI_identity"],
         res4["ARI_region"], vp3p["identity"], vp3p["region"]))
print("DONE")
