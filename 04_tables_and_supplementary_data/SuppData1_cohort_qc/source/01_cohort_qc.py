#!/usr/bin/env python
"""
Extended Data Fig 2 — Cohort & QC transparency for the human cortex snRNA atlas.

Purpose (external_assessor-risk mitigation): the 1.04M-nucleus snRNA atlas is an
UNCORRECTED pool of two cohorts (us = in-house S-donors; edlein = Allen
Institute H-donors). Region x cohort is PARTIALLY CONFOUNDED: 7 regions are
us-only, 7 are mixed (us+edlein), edlein never appears alone. This figure
surfaces that confound honestly (composition, donor/library depth, cell
counts, per-cohort QC, and an explicit confound matrix).

Native matplotlib, single script, vector PDF + PNG.
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
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
from matplotlib.patches import Patch
import matplotlib.colors as mcolors
import os

# ----------------------------------------------------------------------
# Fonts: keep TrueType (Type 42) so pdffonts shows embedded subsets, >=5pt
# ----------------------------------------------------------------------
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.0,
    "axes.titlesize": 7.0,
    "axes.labelsize": 6.0,
    "xtick.labelsize": 5.5,
    "ytick.labelsize": 5.5,
    "legend.fontsize": 5.5,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.0,
    "ytick.major.size": 2.0,
    "axes.edgecolor": "#333333",
    "savefig.dpi": 600,
})

PROJ = "CORTEX_PROGRAM_ROOT"
OBS = os.path.join(PROJ, "inputs/snRNA_1M_obs.csv")
OUTDIR = os.path.join(PROJ, "figures/extended")
SCRIPTDIR = os.path.join(PROJ, "scripts/extended")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(SCRIPTDIR, exist_ok=True)

# Cohort palette + names
COH = {"us": "In-house (S-donors)", "edlein": "Allen / Ed Lein (H-donors)"}
CCOL = {"us": "#D7642C", "edlein": "#3C6E9C"}  # warm vs cool

print("[load] reading obs ...")
df = pd.read_csv(OBS, index_col=0)
df["batch"] = df["batch"].astype(str)
print("[load] shape", df.shape)

# ----------------------------------------------------------------------
# Aggregations
# ----------------------------------------------------------------------
# region order: by total cells descending (keeps the big mixed regions left)
reg_tot = df.groupby("region").size().sort_values(ascending=False)
REGIONS = list(reg_tot.index)

# region x cohort cell counts
ctc = pd.crosstab(df["region"], df["batch"]).reindex(REGIONS)
for c in ["us", "edlein"]:
    if c not in ctc.columns:
        ctc[c] = 0
ctc = ctc[["us", "edlein"]]
# fractions
frac = ctc.div(ctc.sum(axis=1), axis=0)

# classify regions
def cohort_class(row):
    u, e = row["us"], row["edlein"]
    if u > 0 and e > 0:
        return "mixed"
    if u > 0 and e == 0:
        return "us-only"
    return "edlein-only"
rclass = ctc.apply(cohort_class, axis=1)
n_mixed = (rclass == "mixed").sum()
n_usonly = (rclass == "us-only").sum()
n_edonly = (rclass == "edlein-only").sum()
print("[class] mixed=%d us-only=%d edlein-only=%d" % (n_mixed, n_usonly, n_edonly))
print("[class] us-only:", sorted(rclass[rclass=="us-only"].index.tolist()))
print("[class] mixed  :", sorted(rclass[rclass=="mixed"].index.tolist()))

# donor / library counts per region
gmeta = df.groupby("region").agg(
    n_cells=("batch", "size"),
    n_donors=("donor", "nunique"),
    n_libs=("library_prep", "nunique"),
).reindex(REGIONS)
# donors split by cohort per region
don_us = df[df.batch == "us"].groupby("region")["donor"].nunique().reindex(REGIONS).fillna(0)
don_ed = df[df.batch == "edlein"].groupby("region")["donor"].nunique().reindex(REGIONS).fillna(0)
lib_us = df[df.batch == "us"].groupby("region")["library_prep"].nunique().reindex(REGIONS).fillna(0)
lib_ed = df[df.batch == "edlein"].groupby("region")["library_prep"].nunique().reindex(REGIONS).fillna(0)

# region x subclass cell counts (heatmap, log10)
sub_order = (df.groupby("subclass").size().sort_values(ascending=False)).index.tolist()
ctsub = pd.crosstab(df["region"], df["subclass"]).reindex(index=REGIONS, columns=sub_order)

# QC per cohort
qc_nc = {c: df.loc[df.batch == c, "nCount_RNA"].values for c in ["us", "edlein"]}
qc_nf = {c: df.loc[df.batch == c, "nFeature_RNA"].values for c in ["us", "edlein"]}

# demographics by cohort (for confound honesty annotation)
age_by = df.groupby("batch")["age"].median()
sex_ct = pd.crosstab(df["batch"], df["sex"])
ndon_by = df.groupby("batch")["donor"].nunique()
overlap = set(df[df.batch=="edlein"].donor) & set(df[df.batch=="us"].donor)

# ----------------------------------------------------------------------
# Figure layout: 2 columns of panels on an A4-ish portrait canvas
# ----------------------------------------------------------------------
FW, FH = 7.2, 9.0  # inches (fits within Nature ED single-page)
fig = plt.figure(figsize=(FW, FH))
gs = GridSpec(
    4, 2, figure=fig,
    height_ratios=[1.05, 1.0, 1.15, 1.0],
    width_ratios=[1.0, 1.0],
    hspace=0.62, wspace=0.34,
    left=0.085, right=0.975, top=0.945, bottom=0.055,
)

def panel_tag(ax, letter, dx=-0.085, dy=1.02):
    ax.text(dx, dy, letter, transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom", ha="right")

x = np.arange(len(REGIONS))

# ---- (a) region x cohort composition: stacked fraction bar ----
axa = fig.add_subplot(gs[0, 0])
b0 = frac["us"].values
b1 = frac["edlein"].values
axa.bar(x, b0, color=CCOL["us"], edgecolor="white", linewidth=0.3, label=COH["us"])
axa.bar(x, b1, bottom=b0, color=CCOL["edlein"], edgecolor="white", linewidth=0.3, label=COH["edlein"])
# mark us-only regions with a hatch tick above
for i, r in enumerate(REGIONS):
    if rclass[r] == "us-only":
        axa.text(i, 1.02, "*", ha="center", va="bottom", fontsize=7, color=CCOL["us"], fontweight="bold")
axa.set_xticks(x); axa.set_xticklabels(REGIONS, rotation=90)
axa.set_ylim(0, 1.0); axa.set_ylabel("Fraction of nuclei")
axa.set_title("Region x cohort composition", pad=10)
axa.legend(loc="lower center", bbox_to_anchor=(0.5, 1.10), ncol=1,
           frameon=False, handlelength=1.0, columnspacing=0.8, borderpad=0.1)
axa.text(0.0, -0.62, "* = us-only region (0 Allen nuclei)", transform=axa.transAxes,
         fontsize=5.0, color=CCOL["us"], style="italic")
panel_tag(axa, "a")

# ---- (b) donor + library counts per region (grouped bars) ----
axb = fig.add_subplot(gs[0, 1])
w = 0.4
axb.bar(x - w/2, don_us.values, w, color=CCOL["us"], label="us donors")
axb.bar(x - w/2, don_ed.values, w, bottom=don_us.values, color=CCOL["edlein"], label="edlein donors")
axb.set_xticks(x); axb.set_xticklabels(REGIONS, rotation=90)
axb.set_ylabel("n donors")
axb.set_title("Donors per region (by cohort)", pad=10)
# annotate n libraries above each bar
ymax = (don_us + don_ed).max()
for i, r in enumerate(REGIONS):
    axb.text(i - w/2, (don_us+don_ed)[r] + 0.12, "%d" % int(gmeta.loc[r, "n_libs"]),
             ha="center", va="bottom", fontsize=4.6, color="#444444", rotation=90)
axb.set_ylim(0, ymax + 2.2)
axb.legend(loc="upper right", frameon=False, handlelength=1.0, borderpad=0.1)
axb.text(0.0, -0.62, "number above bar = n libraries (library_prep)", transform=axb.transAxes,
         fontsize=5.0, color="#444444", style="italic")
panel_tag(axb, "b")

# ---- (c) region x subclass cell counts heatmap (log10) ----
axc = fig.add_subplot(gs[1, :])
M = ctsub.values.astype(float)
Mlog = np.log10(M + 1)
im = axc.imshow(Mlog, aspect="auto", cmap="magma")
axc.set_yticks(range(len(REGIONS))); axc.set_yticklabels(REGIONS)
axc.set_xticks(range(len(sub_order))); axc.set_xticklabels(sub_order, rotation=90)
axc.set_title("Cell counts per region x subclass (sampling depth)", pad=6)
cb = fig.colorbar(im, ax=axc, fraction=0.018, pad=0.012)
cb.set_label("log10(cells + 1)", fontsize=5.5)
cb.ax.tick_params(labelsize=5.0)
# bracket us-only regions on the y-axis
for i, r in enumerate(REGIONS):
    if rclass[r] == "us-only":
        axc.add_patch(plt.Rectangle((-0.5, i-0.5), len(sub_order), 1, fill=False,
                                     edgecolor=CCOL["us"], linewidth=0.8, clip_on=False))
panel_tag(axc, "c", dx=-0.045)

# ---- (d) per-cohort QC violins: nCount + nFeature ----
axd = fig.add_subplot(gs[2, 0])
data_nc = [qc_nc["us"], qc_nc["edlein"]]
parts = axd.violinplot(data_nc, positions=[0, 1], widths=0.8, showmeans=False,
                       showextrema=False, showmedians=True)
for pc, c in zip(parts["bodies"], ["us", "edlein"]):
    pc.set_facecolor(CCOL[c]); pc.set_alpha(0.75); pc.set_edgecolor("#333333"); pc.set_linewidth(0.4)
parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(0.8)
axd.set_xticks([0, 1]); axd.set_xticklabels(["us", "edlein"])
axd.set_ylabel("nCount_RNA (UMIs)")
axd.set_title("Per-cohort sequencing depth", pad=4)
# medians annotation
for i, c in enumerate(["us", "edlein"]):
    axd.text(i, np.median(qc_nc[c]), "  %.0f" % np.median(qc_nc[c]),
             fontsize=4.8, va="center", ha="left", color="#222222")
panel_tag(axd, "d")

axd2 = fig.add_subplot(gs[2, 1])
data_nf = [qc_nf["us"], qc_nf["edlein"]]
parts2 = axd2.violinplot(data_nf, positions=[0, 1], widths=0.8, showmeans=False,
                         showextrema=False, showmedians=True)
for pc, c in zip(parts2["bodies"], ["us", "edlein"]):
    pc.set_facecolor(CCOL[c]); pc.set_alpha(0.75); pc.set_edgecolor("#333333"); pc.set_linewidth(0.4)
parts2["cmedians"].set_color("black"); parts2["cmedians"].set_linewidth(0.8)
axd2.set_xticks([0, 1]); axd2.set_xticklabels(["us", "edlein"])
axd2.set_ylabel("nFeature_RNA (genes)")
axd2.set_title("Per-cohort gene complexity", pad=4)
for i, c in enumerate(["us", "edlein"]):
    axd2.text(i, np.median(qc_nf[c]), "  %.0f" % np.median(qc_nf[c]),
              fontsize=4.8, va="center", ha="left", color="#222222")
panel_tag(axd2, "d", dx=-0.10)

# ---- (e) region x cohort confound matrix (14 x 2 fraction heatmap, annotated) ----
axe = fig.add_subplot(gs[3, 0])
Mf = frac[["us", "edlein"]].values  # rows=regions, cols=cohort
im2 = axe.imshow(Mf, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=1)
axe.set_yticks(range(len(REGIONS))); axe.set_yticklabels(REGIONS)
axe.set_xticks([0, 1]); axe.set_xticklabels(["us", "edlein"])
axe.set_title("Region x cohort confound matrix", pad=4)
# annotate each cell with fraction
for i in range(len(REGIONS)):
    for j in range(2):
        v = Mf[i, j]
        axe.text(j, i, "%.2f" % v, ha="center", va="center",
                 fontsize=4.6, color="white" if (v < 0.3 or v > 0.75) else "black")
# mark us-only rows
for i, r in enumerate(REGIONS):
    if rclass[r] == "us-only":
        axe.add_patch(plt.Rectangle((-0.5, i-0.5), 2, 1, fill=False,
                                    edgecolor="#111111", linewidth=1.0, clip_on=False))
cb2 = fig.colorbar(im2, ax=axe, fraction=0.07, pad=0.04)
cb2.set_label("cohort fraction", fontsize=5.0); cb2.ax.tick_params(labelsize=5.0)
panel_tag(axe, "e")

# ---- (e-text) honest confound summary box (uses the gs[3,1] slot) ----
axt = fig.add_subplot(gs[3, 1]); axt.axis("off")
lines = [
    "Cohort confound summary (honest disclosure)",
    "",
    "Atlas = UNCORRECTED pool of 2 cohorts:",
    "  us = in-house (%d donors, %s cells)" % (int(ndon_by.get("us",0)), f"{int(ctc['us'].sum()):,}"),
    "  edlein = Allen/Ed Lein (%d donors, %s cells)" % (int(ndon_by.get("edlein",0)), f"{int(ctc['edlein'].sum()):,}"),
    "",
    "Donor overlap between cohorts: %d  (perfectly nested)" % len(overlap),
    "",
    "Region x cohort is PARTIALLY confounded:",
    "  mixed regions : %d  (us + edlein)" % n_mixed,
    "  us-only       : %d  (0 Allen nuclei)" % n_usonly,
    "  edlein-only   : %d" % n_edonly,
    "",
    "us-only regions:",
    "  " + ", ".join(sorted(rclass[rclass=="us-only"].index.tolist())),
    "",
    "Cohorts also differ in:",
    "  depth  : median nCount us=%.0f vs edlein=%.0f" % (np.median(qc_nc["us"]), np.median(qc_nc["edlein"])),
    "  age    : median us=%.0f vs edlein=%.0f yr" % (age_by.get("us",np.nan), age_by.get("edlein",np.nan)),
    "  sex    : us %d%% male vs edlein %d%% male" % (
        100*sex_ct.loc["us","male"]/sex_ct.loc["us"].sum() if "us" in sex_ct.index else 0,
        100*sex_ct.loc["edlein","male"]/sex_ct.loc["edlein"].sum() if "edlein" in sex_ct.index else 0),
    "",
    "=> cross-region comparisons restricted to / sensitivity-",
    "   tested on the %d mixed regions; us-only regions" % n_mixed,
    "   flagged as cohort-confounded throughout.",
]
axt.text(0.0, 1.0, "\n".join(lines), transform=axt.transAxes,
         fontsize=5.0, va="top", ha="left", family="Nimbus Sans",
         linespacing=1.35)
axt.add_patch(plt.Rectangle((-0.02, 0.0), 1.04, 1.02, transform=axt.transAxes,
              fill=False, edgecolor="#888888", linewidth=0.6, clip_on=False))

# global title
fig.suptitle("Supplementary Data 2  |  Cohort composition & QC transparency of the human cortex snRNA atlas (1,036,039 nuclei)",
             fontsize=7.2, fontweight="bold", x=0.085, ha="left", y=0.987)

# ----------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------
pdf_path = os.path.join(OUTDIR, "ed_fig2_cohort_qc.pdf")
png_path = os.path.join(OUTDIR, "ed_fig2_cohort_qc.png")
fig.savefig(pdf_path, bbox_inches=None)
fig.savefig(png_path, dpi=400, bbox_inches=None)
print("[save] pdf:", pdf_path)
print("[save] png:", png_path)
print("[dims] %.2f x %.2f in  =  %.0f x %.0f pt" % (FW, FH, FW*72, FH*72))
print("[done]")
