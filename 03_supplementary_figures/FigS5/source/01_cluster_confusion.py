#!/usr/bin/env python
"""
Extended Data: clustering and identity-versus-region analysis.

raw 60-component inputs are subset with the retained map before calculation.
The retained 54-program outputs comprise four views: subclass, spatial domain,
subclass by region, and domain by region. Each view evaluates whether program
profiles group more strongly by cellular identity than by region, with a
per-program identity/region/residual variance partition.
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import to_rgba
import matplotlib as mpl

ROOT = os.environ.get("CORTEX_NMF_ROOT", "CORTEX_PROGRAM_ROOT")
RES  = f"{ROOT}/results/crossregion_v1"
OUTD = os.environ.get("FIGS5_CACHE_DIR", f"{RES}/cluster_confusion")
FIGD = os.environ.get("FIGS5_OUTPUT_DIR", f"{ROOT}/figures/extended")
os.makedirs(OUTD, exist_ok=True); os.makedirs(FIGD, exist_ok=True)

MIN_CELLS = 20
RETAIN_MAP = os.environ.get("RETAIN_MAP", f"{ROOT}/tables/TableS3_program_annotation.tsv")
MAP = pd.read_csv(RETAIN_MAP, sep="\t")
assert len(MAP) == 54
MAP["new_int"] = MAP["new_P"].astype(str).str.removeprefix("P").astype(int)
MAP["old_int"] = MAP["cnmf_component"].astype(int)
assert MAP["new_int"].tolist() == list(range(1, 55))
assert set(range(1, 61)) - set(MAP["old_int"]) == {9, 18, 19, 35, 52, 57}
OLD_PROG = MAP["old_int"].astype(str).tolist()
PROG = MAP["new_int"].astype(str).tolist()
OLD_TO_NEW = dict(zip(OLD_PROG, PROG))
PROG_SP = [f"program_{i}" for i in MAP["old_int"]]

# ---- font / style: >=5pt, vector-friendly ----
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 6, "axes.linewidth": 0.4, "xtick.major.width": 0.4,
    "ytick.major.width": 0.4, "xtick.major.size": 2, "ytick.major.size": 2,
})

# ---------------- program short names ----------------
nm = MAP.copy()
nm["program"] = nm["new_int"].astype(str)
short = dict(zip(nm["program"], nm["functional_name"]))
conf  = dict(zip(nm["program"], nm["confidence"]))
def plabel(p):
    s = short.get(p, p)
    star = " *" if str(conf.get(p, "")).endswith("weak") else ""
    return f"P{p} {s}{star}"
PROG_LABELS = [plabel(p) for p in PROG]

# ---------------- helpers ----------------
def corr_ward(M):
    """rows=objects; correlation distance + ward linkage. Returns Z, leaf order."""
    D = pdist(M, metric="correlation")
    D = np.nan_to_num(D, nan=1.0)
    Z = linkage(D, method="ward")
    order = dendrogram(Z, no_plot=True)["leaves"]
    return Z, order

def zscore_cols(M):
    """z-score each program (column) across rows for display."""
    mu = M.mean(0); sd = M.std(0); sd[sd == 0] = 1
    return (M - mu) / sd

def var_partition(df, group_col, region_col, val_cols):
    """Per-program 2-way (group + region) variance partition via eta^2 (type-I-ish,
    sequential SS group then region). Returns mean % over programs."""
    rows = []
    for p in val_cols:
        y = df[p].values.astype(float)
        gt = y.mean()
        ss_tot = ((y - gt) ** 2).sum()
        if ss_tot <= 0:
            continue
        # group SS
        gmean = df.groupby(group_col)[p].transform("mean").values
        ss_g = ((gmean - gt) ** 2).sum()
        # region SS *after* group (on residuals): fit region on (y - gmean) cell means
        resid = y - gmean
        rmean = pd.Series(resid).groupby(df[region_col].values).transform("mean").values
        ss_r = (rmean ** 2).sum()
        ss_res = ss_tot - ss_g - ss_r
        rows.append({"program": p, "grp": ss_g/ss_tot, "region": ss_r/ss_tot,
                     "resid": max(ss_res, 0)/ss_tot})
    vp = pd.DataFrame(rows)
    return vp

# ================= LOAD CELL-LEVEL (views 1 & 3) =================
print("loading cell-level...", flush=True)
cell = pd.read_parquet(f"{RES}/cell_program_region_subclass.parquet")
cell.columns = [str(c) for c in cell.columns]
# view 1: subclass means
missing_cell = set(OLD_PROG) - set(cell.columns)
assert not missing_cell, f"cell parquet missing retained old components: {sorted(missing_cell)}"
sub_mean = cell.groupby("subclass")[OLD_PROG].mean().rename(columns=OLD_TO_NEW)
assert sub_mean.shape[1] == 54
sub_mean.to_csv(f"{OUTD}/view1_subclass_mean.tsv", sep="\t")
# view 3: subclass x region means, >=20 cells
gsize = cell.groupby(["subclass", "region"]).size()
keep = gsize[gsize >= MIN_CELLS].index
sr_mean = cell.groupby(["subclass", "region"])[OLD_PROG].mean().loc[keep].rename(columns=OLD_TO_NEW)
assert sr_mean.shape[1] == 54
sr_mean.to_csv(f"{OUTD}/view3_subclass_region_mean.tsv", sep="\t")
print(f"  subclass={sub_mean.shape[0]} | subclass*region rows={sr_mean.shape[0]}", flush=True)

# variance partition (subclass x region) -- on raw cells for proper SS
vp3 = var_partition(cell, "subclass", "region", OLD_PROG)
vp3["program"] = vp3["program"].map(OLD_TO_NEW)
assert len(vp3) == 54
vp3.to_csv(f"{OUTD}/varpart_subclass_region.tsv", sep="\t", index=False)

del cell  # free RAM

# ================= LOAD SPATIAL (views 2 & 4) =================
print("loading spatial (SCT)...", flush=True)
sp = pd.read_parquet(f"{RES}/spatial_bin50_program_score_SCT.parquet",
                     columns=["majorDomain", "region"] + PROG_SP)
ren = {f"program_{old}": str(new) for old, new in zip(MAP["old_int"], MAP["new_int"])}
sp = sp.rename(columns=ren)
sp = sp.dropna(subset=["majorDomain", "region"])
# view 2: domain means
dom_mean = sp.groupby("majorDomain")[PROG].mean()
assert dom_mean.shape[1] == 54
dom_mean.to_csv(f"{OUTD}/view2_domain_mean.tsv", sep="\t")
# view 4: domain x region means, >=20 bins
dsize = sp.groupby(["majorDomain", "region"]).size()
dkeep = dsize[dsize >= MIN_CELLS].index
dr_mean = sp.groupby(["majorDomain", "region"])[PROG].mean().loc[dkeep]
assert dr_mean.shape[1] == 54
dr_mean.to_csv(f"{OUTD}/view4_domain_region_mean.tsv", sep="\t")
print(f"  domain={dom_mean.shape[0]} | domain*region rows={dr_mean.shape[0]}", flush=True)

vp4 = var_partition(sp, "majorDomain", "region", PROG)
assert len(vp4) == 54
vp4.to_csv(f"{OUTD}/varpart_domain_region.tsv", sep="\t", index=False)
del sp

# ================= CONFUSION TEST (views 3 & 4) =================
def confusion_test(M_df, identity_name):
    """M_df: rows = (identity, region) MultiIndex, cols = programs.
    Cluster rows (corr/ward), cut at k=n_identity, compute ARI/NMI vs identity & region."""
    ident = M_df.index.get_level_values(0).values
    region = M_df.index.get_level_values(1).values
    k = len(np.unique(ident))
    M = M_df.values
    Z, order = corr_ward(M)
    cl = fcluster(Z, t=k, criterion="maxclust")
    res = {
        "identity_name": identity_name, "n_rows": len(ident), "k": int(k),
        "ARI_identity": float(adjusted_rand_score(ident, cl)),
        "ARI_region":   float(adjusted_rand_score(region, cl)),
        "NMI_identity": float(normalized_mutual_info_score(ident, cl)),
        "NMI_region":   float(normalized_mutual_info_score(region, cl)),
    }
    cont_id  = pd.crosstab(pd.Series(cl, name="cluster"), pd.Series(ident, name=identity_name))
    cont_reg = pd.crosstab(pd.Series(cl, name="cluster"), pd.Series(region, name="region"))
    return res, order, cl, cont_id, cont_reg, Z

res3, ord3, cl3, cont3_id, cont3_reg, Z3 = confusion_test(sr_mean, "subclass")
res4, ord4, cl4, cont4_id, cont4_reg, Z4 = confusion_test(dr_mean, "majorDomain")

cont3_id.to_csv(f"{OUTD}/view3_contingency_cluster_x_subclass.tsv", sep="\t")
cont3_reg.to_csv(f"{OUTD}/view3_contingency_cluster_x_region.tsv", sep="\t")
cont4_id.to_csv(f"{OUTD}/view4_contingency_cluster_x_domain.tsv", sep="\t")
cont4_reg.to_csv(f"{OUTD}/view4_contingency_cluster_x_region.tsv", sep="\t")

summary = {
    "view3_subclass_x_region": res3,
    "view4_domain_x_region": res4,
    "varpart_subclass_region_mean_pct": {
        "identity": float(vp3["grp"].mean()*100), "region": float(vp3["region"].mean()*100),
        "residual": float(vp3["resid"].mean()*100)},
    "varpart_domain_region_mean_pct": {
        "identity": float(vp4["grp"].mean()*100), "region": float(vp4["region"].mean()*100),
        "residual": float(vp4["resid"].mean()*100)},
    "min_cells": MIN_CELLS,
}
with open(f"{OUTD}/confusion_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print("=== SUMMARY ===", flush=True)
print(json.dumps(summary, indent=2), flush=True)

# ================= FIGURE =================
# color maps for identity / region strips
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

fig = plt.figure(figsize=(7.2, 9.6))
# 5 logical panels: a,b (top: subclass & domain dendro+heat), c,d (mid: confusion), e (bars)
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.15, 1.45, 0.7],
              hspace=0.55, wspace=0.45,
              left=0.16, right=0.985, top=0.96, bottom=0.06)

def panel_dendro_heat(subgs, M_df, Z, order, row_label_fmt, title, col_strip=None,
                      strip_lut=None, vmax=2.5):
    """Heatmap with left dendrogram; rows ordered by `order`; columns=programs (clustered)."""
    Mz = zscore_cols(M_df.values)
    # cluster columns too for readability
    _, cord = corr_ward(Mz.T)
    Mz = Mz[np.ix_(order, cord)]
    rlabs = [row_label_fmt(M_df.index[i]) for i in order]
    clabs = [PROG_LABELS[j] for j in cord]
    inner = subgs.subgridspec(1, 3 if col_strip is None else 4,
                              width_ratios=[0.18, (0.0 if col_strip is None else 0.05), 1.0, 0.0][:(3 if col_strip is None else 4)],
                              wspace=0.04)
    # dendrogram
    axd = fig.add_subplot(inner[0, 0])
    dendrogram(Z, orientation="left", ax=axd, no_labels=True,
               color_threshold=0, link_color_func=lambda k: "0.35")
    axd.invert_yaxis(); axd.axis("off")
    hcol = 2 if col_strip is None else 2
    axh = fig.add_subplot(inner[0, hcol])
    im = axh.imshow(Mz, aspect="equal", cmap=CMAP, vmin=-vmax, vmax=vmax,
                    interpolation="nearest")
    axh.set_yticks(range(len(rlabs))); axh.set_yticks(range(len(rlabs)))
    axh.set_yticklabels(rlabs, fontsize=4.6)
    axh.set_xticks(range(len(clabs)))
    axh.set_xticklabels(clabs, fontsize=3.4, rotation=90)
    axh.tick_params(length=1, pad=1)
    axh.set_title(title, fontsize=7, pad=3, loc="left", fontweight="bold")
    cb = fig.colorbar(im, ax=axh, fraction=0.018, pad=0.01)
    cb.ax.tick_params(labelsize=4.5, length=1.5); cb.set_label("col z", fontsize=5)
    return axh

# ---- panel a: subclass ----
Z1, ord1 = corr_ward(zscore_cols(sub_mean.values))
panel_dendro_heat(gs[0, 0], sub_mean, Z1, ord1, lambda idx: str(idx),
                  "a  Subclass (22) -> program profile")
# ---- panel b: domain ----
Z2, ord2 = corr_ward(zscore_cols(dom_mean.values))
panel_dendro_heat(gs[0, 1], dom_mean, Z2, ord2, lambda idx: str(idx),
                  "b  Spatial domain (8) -> program profile")

# ---- panel c: subclass x region confusion ----
def panel_confusion(subgs, M_df, Z, order, res, ident_lut, title, letter):
    Mz = zscore_cols(M_df.values)
    _, cord = corr_ward(Mz.T)
    Mz = Mz[np.ix_(order, cord)]
    ident = M_df.index.get_level_values(0).values[order]
    region = M_df.index.get_level_values(1).values[order]
    # layout: dendro | ident strip | region strip | heatmap
    inner = subgs.subgridspec(1, 4, width_ratios=[0.16, 0.035, 0.035, 1.0], wspace=0.03)
    axd = fig.add_subplot(inner[0, 0])
    dendrogram(Z, orientation="left", ax=axd, no_labels=True, color_threshold=0,
               link_color_func=lambda k: "0.35")
    axd.invert_yaxis(); axd.axis("off")
    # strips
    axs1 = fig.add_subplot(inner[0, 1])
    axs1.imshow([[to_rgba(ident_lut[ident[i]]) for _ in range(1)] for i in range(len(ident))],
                aspect="auto"); axs1.set_xticks([]); axs1.set_yticks([])
    axs1.set_xlabel("id", fontsize=4, rotation=90, labelpad=1)
    axs2 = fig.add_subplot(inner[0, 2])
    axs2.imshow([[to_rgba(REG_COL[region[i]]) for _ in range(1)] for i in range(len(region))],
                aspect="auto"); axs2.set_xticks([]); axs2.set_yticks([])
    axs2.set_xlabel("reg", fontsize=4, rotation=90, labelpad=1)
    axh = fig.add_subplot(inner[0, 3])
    im = axh.imshow(Mz, aspect="auto", cmap=CMAP, vmin=-2.5, vmax=2.5, interpolation="nearest")
    ticks = list(range(0, len(PROG), 5))
    axh.set_yticks([]); axh.set_xticks(ticks)
    axh.set_xticklabels([f"P{cord[j]+1}" for j in ticks], fontsize=3.6, rotation=90)
    axh.tick_params(length=1, pad=1)
    ann = (f"ARI(identity)={res['ARI_identity']:.2f}  ARI(region)={res['ARI_region']:.2f}\n"
           f"NMI(id)={res['NMI_identity']:.2f}  NMI(reg)={res['NMI_region']:.2f}")
    axh.set_title(f"{letter}  {title}\n{ann}", fontsize=6, pad=3, loc="left",
                  fontweight="bold")
    cb = fig.colorbar(im, ax=axh, fraction=0.014, pad=0.01)
    cb.ax.tick_params(labelsize=4.5, length=1.5)

panel_confusion(gs[1, 0], sr_mean, Z3, ord3, res3, SUB_COL,
                f"Subclass x region ({res3['n_rows']} groups)", "c")
panel_confusion(gs[1, 1], dr_mean, Z4, ord4, res4, DOM_COL,
                f"Domain x region ({res4['n_rows']} groups)", "d")

# ---- panel e: variance partition bars ----
axe = fig.add_subplot(gs[2, :])
groups = ["Subclass x region", "Domain x region"]
ident_pct = [vp3["grp"].mean()*100, vp4["grp"].mean()*100]
reg_pct   = [vp3["region"].mean()*100, vp4["region"].mean()*100]
res_pct   = [vp3["resid"].mean()*100, vp4["resid"].mean()*100]
y = np.arange(len(groups)); h = 0.55
axe.barh(y, ident_pct, h, color="#2c7fb8", label="Identity (subclass/domain)")
axe.barh(y, reg_pct, h, left=ident_pct, color="#d95f0e", label="Region")
axe.barh(y, res_pct, h, left=np.array(ident_pct)+np.array(reg_pct), color="0.8",
         label="Residual")
for i in range(len(groups)):
    axe.text(ident_pct[i]/2, y[i], f"{ident_pct[i]:.0f}%", ha="center", va="center",
             fontsize=5.5, color="white", fontweight="bold")
    axe.text(ident_pct[i]+reg_pct[i]/2, y[i], f"{reg_pct[i]:.0f}%", ha="center",
             va="center", fontsize=5.5, color="white", fontweight="bold")
axe.set_yticks(y); axe.set_yticklabels(groups, fontsize=6)
axe.set_xlabel("mean variance explained across 54 retained programs (%)", fontsize=6)
axe.set_xlim(0, 100); axe.set_title("e  Variance partition: identity vs region",
                                    fontsize=7, loc="left", fontweight="bold", pad=3)
axe.legend(fontsize=5, loc="lower right", frameon=False, ncol=3)
axe.spines[["top", "right"]].set_visible(False)

fig.savefig(f"{FIGD}/ed_cluster_confusion.pdf", dpi=600)
fig.savefig(f"{FIGD}/ed_cluster_confusion.png", dpi=300)
fig.savefig("/tmp/ed_cluster_confusion.png", dpi=300)
print(f"WROTE {FIGD}/ed_cluster_confusion.pdf|png", flush=True)
print("DONE", flush=True)
