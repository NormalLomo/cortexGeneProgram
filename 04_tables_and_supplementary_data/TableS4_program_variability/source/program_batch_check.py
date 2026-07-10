#!/usr/bin/env python
"""
Extended Data — cNMF program BATCH (cohort) robustness check.

QUESTION (external_assessor-risk): are the 60 cNMF programs strongly driven by the
COHORT confound (us = in-house S-donors vs edlein = Allen H-donors), rather
than by real biology / cross-region structure?

KEY CONFOUND (from ed_fig2_cohort_qc): region x cohort is PARTIALLY confounded
(7 us-only regions, 7 mixed, edlein never alone) AND the two cohorts have
different subclass composition + different sequencing depth. So the *fair*
test of a cohort effect is the cohort effect WITHIN cell identity (subclass).

METHOD (per program, on all 1,036,039 cells):
  Fit additive two-way linear models on standardized usage:
    M_sub  : usage ~ C(subclass)
    M_sc   : usage ~ C(subclass) + C(batch)      (batch = cohort)
    M_sr   : usage ~ C(subclass) + C(region)
  Type-II partial eta^2 for a term = SS(term) / (SS(term) + SS_resid_full),
  where SS(term) = RSS(reduced) - RSS(full).
    cohort_partial_eta2 = (RSS(M_sub) - RSS(M_sc)) / (RSS(M_sub)-RSS(M_sc) + RSS(M_sc))
    region_partial_eta2 = (RSS(M_sub) - RSS(M_sr)) / (RSS(M_sub)-RSS(M_sr) + RSS(M_sr))
  subclass_eta2 (one-way) = 1 - RSS(M_sub)/TSS  (== classic eta^2 of subclass alone)
  Also report the *marginal* one-way cohort eta^2 (usage ~ batch) for contrast
  (this is the NAIVE number that ignores composition).
  Depth check: Pearson corr of each program usage vs log10(nCount_RNA), and the
  cohort mean-usage gap, to flag depth/QC axes.

RSS for an additive categorical model is computed exactly via a cheap
normal-equations solve on the one-hot design (drop reference level), all 60
programs solved at once (X'X is tiny; X'Y is (p x 60)).

Outputs:
  results/crossregion_v1/program_cohort_eta2.tsv
  figures/extended/ed_program_batch_check.{pdf,png}
"""
import os
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

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.0, "axes.titlesize": 7.0, "axes.labelsize": 6.0,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0,
    "axes.edgecolor": "#333333", "savefig.dpi": 600,
})

PROJ = "CORTEX_PROGRAM_ROOT"
SCORES = os.path.join(PROJ, "results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_cell_scores.tsv")
OBS = os.path.join(PROJ, "inputs/snRNA_1M_obs.csv")
NAMES = os.path.join(PROJ, "results/crossregion_v1/program_names.tsv")
VARI = os.path.join(PROJ, "results/crossregion_v1/program_variability.tsv")
RENUMBER = os.path.join(PROJ, "results/crossregion_v1/program_renumber_map.tsv")
OUTDIR = os.path.join(PROJ, "figures/extended")
RESDIR = os.path.join(PROJ, "results/crossregion_v1")
os.makedirs(OUTDIR, exist_ok=True)

CCOL = {"us": "#D7642C", "edlein": "#3C6E9C"}
SIGCOL = {"brain-sig": "#2E6E4E", "brain-weak": "#B0883B"}  # green vs amber

# ----------------------------------------------------------------------
# RENUMBER MAP: old_P (int) -> new label for display
#   kept programs: "P{new_P}"  (new_P is the new integer 1-54)
#   excluded programs (old 9/18/19/35/52/57): "cNMF component {old_P}"
# ----------------------------------------------------------------------
_rmap_df = pd.read_csv(RENUMBER, sep="\t")
# Build dict: old_P_int -> display_label
_OLD_TO_LABEL = {}
for _, row in _rmap_df.iterrows():
    old_p = int(row["old_P"])
    new_p = str(row["new_P"])
    if new_p == "EXCLUDED":
        _OLD_TO_LABEL[old_p] = "cNMF component %d" % old_p
    else:
        _OLD_TO_LABEL[old_p] = "P%s" % new_p

EXCLUDED_OLD = set(int(row["old_P"]) for _, row in _rmap_df.iterrows()
                   if str(row["new_P"]) == "EXCLUDED")

# Load confidence (brain-weak) for asterisk display
# program_names.tsv new_P format: "P4", "P12", etc. (NOT bare int)
_names_df = pd.read_csv(NAMES, sep="\t")
# Key: new_P integer -> confidence string (strip "P" prefix)
_CONF_MAP = {}
for _, row in _names_df.iterrows():
    np_str = str(row["new_P"])
    if np_str not in ("EXCLUDED", "nan") and np_str.startswith("P"):
        try:
            _CONF_MAP[int(np_str[1:])] = str(row["confidence"])
        except (ValueError, TypeError):
            pass

def _star_for_old(old_p):
    """Return '*' if new_P for this old_p is brain-weak, else ''."""
    label = _OLD_TO_LABEL.get(old_p, "")
    if label.startswith("P"):
        try:
            new_int = int(label[1:])
            return "*" if _CONF_MAP.get(new_int, "") == "brain-weak" else ""
        except ValueError:
            pass
    return ""

def prog_label(old_p):
    """Return display label for a program given its old cNMF component integer (+ '*' if brain-weak)."""
    base = _OLD_TO_LABEL.get(old_p, "P%d" % old_p)
    return base + _star_for_old(old_p)

# ----------------------------------------------------------------------
# REPLOT_ONLY: skip the heavy 1M-cell load + ANCOVA, rebuild `tab` from the
# already-saved results table. Use for layout-only re-renders (fix_figS3).
#   REPLOT_ONLY=1 python program_batch_check.py
# ----------------------------------------------------------------------
REPLOT_ONLY = os.environ.get("REPLOT_ONLY", "0") == "1"

def _build_tab_from_data():
    """Heavy path: load 1M-cell scores + obs, fit ANCOVA, assemble + save table."""
    print("[load] scores ...", flush=True)
    S = pd.read_csv(SCORES, sep="\t", index_col=0)
    S.columns = [int(c) for c in S.columns]
    S = S.sort_index(axis=1)
    PROGS = list(S.columns)
    print("  scores:", S.shape, "progs:", PROGS[:3], "...", PROGS[-1], flush=True)

    print("[load] obs ...", flush=True)
    obs = pd.read_csv(OBS, index_col=0,
                      usecols=["Unnamed: 0", "batch", "region", "subclass", "class",
                               "nCount_RNA", "nFeature_RNA"])
    # align
    assert set(obs.index) == set(S.index), "barcode set mismatch"
    obs = obs.reindex(S.index)
    assert (obs.index == S.index).all()
    print("  100%% barcode match. n=%d" % len(obs), flush=True)
    print("  cohorts:", obs["batch"].value_counts().to_dict(), flush=True)

    for c in ["batch", "region", "subclass", "class"]:
        obs[c] = obs[c].astype(str)
    log_nc = np.log10(obs["nCount_RNA"].values.astype(float) + 1.0)

    # ----------------------------------------------------------------------
    # Helpers: exact RSS of additive categorical model via normal equations.
    # Design = intercept + dummies (drop first level) for each factor.
    # Solve beta = (X'X)^-1 X'Y for ALL programs at once; RSS = sum((Y-Xb)^2).
    # ----------------------------------------------------------------------
    def design_matrix(factors):
        """factors: list of pandas Series (categorical). Returns dense float X with
        intercept + one-hot(drop-first) for each factor."""
        n = len(factors[0])
        cols = [np.ones(n)]
        for f in factors:
            cats = pd.Categorical(f)
            codes = cats.codes
            K = len(cats.categories)
            # drop first level (code 0) as reference
            for k in range(1, K):
                cols.append((codes == k).astype(float))
        return np.column_stack(cols)

    def rss_all(X, Y):
        """RSS per column of Y for OLS fit of Y on X. Y: (n x P)."""
        XtX = X.T @ X
        XtY = X.T @ Y
        beta = np.linalg.solve(XtX, XtY)        # (q x P)
        Yhat = X @ beta
        resid = Y - Yhat
        return (resid * resid).sum(axis=0)      # (P,)

    # Standardize each program usage (z-score) so eta^2 are comparable & TSS=N-ish.
    Y = S.values.astype(np.float64)
    Ymean = Y.mean(axis=0, keepdims=True)
    Ystd = Y.std(axis=0, ddof=0, keepdims=True)
    Ystd[Ystd == 0] = 1.0
    Yz = (Y - Ymean) / Ystd
    N = Yz.shape[0]
    TSS = (Yz * Yz).sum(axis=0)   # == N for each (since z-scored)

    sub = obs["subclass"]
    reg = obs["region"]
    bat = obs["batch"]

    print("[fit] design matrices ...", flush=True)
    X_sub = design_matrix([sub])
    X_sc = design_matrix([sub, bat])
    X_sr = design_matrix([sub, reg])
    X_bat = design_matrix([bat])
    X_reg = design_matrix([reg])
    print("  dims: X_sub", X_sub.shape, "X_sc", X_sc.shape, "X_sr", X_sr.shape, flush=True)

    print("[fit] RSS solves ...", flush=True)
    rss_sub = rss_all(X_sub, Yz)
    rss_sc = rss_all(X_sc, Yz)
    rss_sr = rss_all(X_sr, Yz)
    rss_bat = rss_all(X_bat, Yz)
    rss_reg = rss_all(X_reg, Yz)

    # Partial eta^2 (type-II): SS_term / (SS_term + RSS_full)
    ss_cohort = rss_sub - rss_sc                       # added by cohort over subclass
    ss_region = rss_sub - rss_sr                       # added by region over subclass
    cohort_partial_eta2 = ss_cohort / (ss_cohort + rss_sc)
    region_partial_eta2 = ss_region / (ss_region + rss_sr)
    # one-way (marginal) eta^2
    subclass_eta2 = 1.0 - rss_sub / TSS
    cohort_marg_eta2 = 1.0 - rss_bat / TSS
    region_marg_eta2 = 1.0 - rss_reg / TSS

    # ----------------------------------------------------------------------
    # Depth / cohort-gap diagnostics
    # ----------------------------------------------------------------------
    print("[diag] depth corr + cohort gap ...", flush=True)
    # Pearson corr usage(z) vs log10 nCount
    lnc = (log_nc - log_nc.mean()) / log_nc.std()
    depth_corr = (Yz * lnc[:, None]).sum(axis=0) / N    # corr since both ~z
    # cohort mean gap on z-usage (us - edlein), in SD units
    is_us = (bat.values == "us")
    gap = Yz[is_us].mean(axis=0) - Yz[~is_us].mean(axis=0)

    # ----------------------------------------------------------------------
    # Assemble table
    # ----------------------------------------------------------------------
    names_raw = pd.read_csv(NAMES, sep="\t"); names = names_raw[names_raw["new_P"].astype(str).str.startswith("P")].set_index("new_P")
    vari_raw = pd.read_csv(VARI, sep="\t"); vari = vari_raw[vari_raw["new_P"].astype(str).str.startswith("P")].set_index("new_P") if "new_P" in vari_raw.columns else vari_raw.set_index("program")
    F3_VARIABLE = set(vari.index[vari["class"] == "variable"].tolist())
    print("  F3 variable programs (n=%d):" % len(F3_VARIABLE), sorted(F3_VARIABLE), flush=True)

    tab = pd.DataFrame({
        "program": PROGS,
        "name_short": [names.loc[p, "name_short"] if p in names.index else str(p) for p in PROGS],
        "confidence": [names.loc[p, "confidence"] if p in names.index else "NA" for p in PROGS],
        "cohort_partial_eta2": cohort_partial_eta2,
        "region_partial_eta2": region_partial_eta2,
        "subclass_eta2": subclass_eta2,
        "cohort_marginal_eta2": cohort_marg_eta2,
        "region_marginal_eta2": region_marg_eta2,
        "depth_corr_log10nCount": depth_corr,
        "cohort_gap_us_minus_edlein_SD": gap,
        "f3_region_variable": [p in F3_VARIABLE for p in PROGS],
    }).set_index("program")
    tab = tab.sort_values("cohort_partial_eta2", ascending=False)
    OUTTAB = os.path.join(RESDIR, "program_cohort_eta2.tsv")
    tab.to_csv(OUTTAB, sep="\t")
    print("[write]", OUTTAB, flush=True)

    return tab

# ----------------------------------------------------------------------
# Dispatch: replot-only (read saved table) vs full compute
# ----------------------------------------------------------------------
if REPLOT_ONLY:
    OUTTAB = os.path.join(RESDIR, "program_cohort_eta2.tsv")
    print("[replot-only] loading saved table:", OUTTAB, flush=True)
    tab = pd.read_csv(OUTTAB, sep="\t").set_index("program")
else:
    tab = _build_tab_from_data()


# ----------------------------------------------------------------------
# Console summary
# ----------------------------------------------------------------------
ce = tab["cohort_partial_eta2"]
print("\n===== COHORT partial-eta^2 distribution (n=60) =====", flush=True)
print("  median = %.4f" % ce.median(), flush=True)
print("  mean   = %.4f" % ce.mean(), flush=True)
print("  max    = %.4f  (program %s, %s)" % (ce.max(), ce.idxmax(), tab.loc[ce.idxmax(),"name_short"]), flush=True)
print("  n > 0.05 = %d" % (ce > 0.05).sum(), flush=True)
print("  n > 0.10 = %d" % (ce > 0.10).sum(), flush=True)
print("\n  WORST 8 offenders (cohort partial-eta^2):", flush=True)
for p in ce.head(8).index:
    r = tab.loc[p]
    print("   %s  cohort=%.3f  region=%.3f  subclass=%.3f  depthr=%+.2f  gap=%+.2f  %s%s"
          % (prog_label(p), r["cohort_partial_eta2"], r["region_partial_eta2"], r["subclass_eta2"],
             r["depth_corr_log10nCount"], r["cohort_gap_us_minus_edlein_SD"],
             r["name_short"], "  [F3-var]" if r["f3_region_variable"] else ""), flush=True)

print("\n  14 F3 region-VARIABLE programs (cohort vs region partial-eta^2):", flush=True)
sub14 = tab[tab["f3_region_variable"]].sort_values("region_partial_eta2", ascending=False)
for p in sub14.index:
    r = sub14.loc[p]
    ratio = r["region_partial_eta2"] / max(r["cohort_partial_eta2"], 1e-9)
    print("   %s  region=%.3f  cohort=%.3f  (region/cohort=%5.1fx)  depthr=%+.2f  %s"
          % (prog_label(p), r["region_partial_eta2"], r["cohort_partial_eta2"], ratio,
             r["depth_corr_log10nCount"], r["name_short"]), flush=True)

print("\n  Depth flags |corr(log10 nCount)| > 0.30:", flush=True)
dflag = tab[tab["depth_corr_log10nCount"].abs() > 0.30].sort_values("depth_corr_log10nCount")
if len(dflag) == 0:
    print("   (none)", flush=True)
for p in dflag.index:
    r = dflag.loc[p]
    print("   %s  depthr=%+.2f  cohort=%.3f  %s" % (prog_label(p), r["depth_corr_log10nCount"],
          r["cohort_partial_eta2"], r["name_short"]), flush=True)

# ----------------------------------------------------------------------
# FIGURE: 3 panels
# ----------------------------------------------------------------------
print("\n[plot] building figure ...", flush=True)
fig = plt.figure(figsize=(7.2, 7.6))
# fix_figS3: bottom raised 0.075 -> 0.175 so panel-c 40-deg rotated program-name
# x-tick labels have vertical room and their leading words are not clipped by the
# figure bottom edge (was: "ort to cytosol" / "rocyte/myelin").
gs = GridSpec(3, 1, figure=fig, height_ratios=[1.35, 1.0, 1.05],
              hspace=0.62, left=0.10, right=0.965, top=0.945, bottom=0.175)

def sigc(p):
    return SIGCOL.get(tab.loc[p, "confidence"], "#888888")

# ---- (a) cohort vs region partial-eta^2 scatter (60 dots) ----
axa = fig.add_subplot(gs[0, 0])
xr = tab["region_partial_eta2"].values
yc = tab["cohort_partial_eta2"].values
cols = [sigc(p) for p in tab.index]
isf3 = tab["f3_region_variable"].values
# non-F3 dots
axa.scatter(xr[~isf3], yc[~isf3], s=18, c=[cols[i] for i in range(len(cols)) if not isf3[i]],
            alpha=0.78, linewidths=0.3, edgecolors="white", zorder=3)
# F3-variable dots: ring them
axa.scatter(xr[isf3], yc[isf3], s=34, facecolors=[cols[i] for i in range(len(cols)) if isf3[i]],
            edgecolors="#111111", linewidths=0.9, alpha=0.95, zorder=4)
mx = max(xr.max(), yc.max()) * 1.08
axa.plot([0, mx], [0, mx], ls="--", lw=0.6, c="#999999", zorder=1)
# y=x reference + guide lines at 0.05/0.10 cohort
for thr, lab in [(0.05, "0.05"), (0.10, "0.10")]:
    axa.axhline(thr, ls=":", lw=0.5, c="#C24C4C", zorder=1)
    axa.text(mx*0.995, thr, " cohort=%s" % lab, fontsize=4.6, color="#C24C4C",
             va="bottom", ha="right")
# label offenders (top cohort) + all F3 dots
# fix_figS3 current (2026-06-27): use adjustText to prevent P9/P3 (and similar
# near-collisions) from rendering as overlapping "P9 P3" duplicate-label visual artifacts.
to_label = set(ce.head(6).index) | set(tab.index[isf3])
try:
    from adjustText import adjust_text
    texts = []
    for p in to_label:
        texts.append(axa.text(tab.loc[p, "region_partial_eta2"],
                              tab.loc[p, "cohort_partial_eta2"],
                              prog_label(p),
                              fontsize=4.7, color="#222222", zorder=5))
    adjust_text(texts, ax=axa, only_move={"text": "xy"},
                expand_text=(1.05, 1.20), expand_points=(1.05, 1.20),
                force_text=(0.40, 0.55), force_points=(0.25, 0.30),
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.25, alpha=0.7),
                max_move=8.0)
except Exception as _e:
    print("  [warn] adjustText failed (%s); falling back to fixed offsets" % _e, flush=True)
    for p in to_label:
        axa.annotate(prog_label(p),
                     (tab.loc[p, "region_partial_eta2"], tab.loc[p, "cohort_partial_eta2"]),
                     fontsize=4.7, xytext=(2.2, 2.2), textcoords="offset points",
                     color="#222222", zorder=5)
axa.set_xlabel("Region partial-$\\eta^2$  (usage ~ subclass + region)")
axa.set_ylabel("Cohort partial-$\\eta^2$\n(usage ~ subclass + cohort)")
axa.set_xlim(-mx*0.02, mx); axa.set_ylim(-mx*0.02, mx)
axa.set_title("(a)  Per-program cohort vs region effect, controlling for subclass identity",
              loc="left", pad=4)
leg = [Patch(fc=SIGCOL["brain-sig"], label="brain-sig"),
       Patch(fc=SIGCOL["brain-weak"], label="brain-weak"),
       plt.Line2D([0],[0], marker="o", ls="", mfc="#ccc", mec="#111", mew=0.9, ms=5,
                  label="F3 region-variable (n=14)")]
axa.legend(handles=leg, loc="upper left", frameon=False, fontsize=5.0, handletextpad=0.4,
           borderaxespad=0.2)

# ---- (b) ranked bar of cohort partial-eta^2 for all 60 ----
axb = fig.add_subplot(gs[1, 0])
order = tab.index.tolist()  # already sorted desc by cohort eta2
vals = tab["cohort_partial_eta2"].values
barcols = [sigc(p) for p in order]
xpos = np.arange(len(order))
bars = axb.bar(xpos, vals, color=barcols, width=0.78, linewidth=0.2, edgecolor="white")
# hatch / outline the F3-variable ones
for i, p in enumerate(order):
    if tab.loc[p, "f3_region_variable"]:
        bars[i].set_edgecolor("#111111"); bars[i].set_linewidth(0.8)
axb.axhline(0.05, ls=":", lw=0.6, c="#C24C4C")
axb.axhline(0.10, ls="--", lw=0.6, c="#8B2E2E")
axb.text(len(order)-0.5, 0.05, " 0.05", fontsize=4.8, color="#C24C4C", va="bottom", ha="right")
axb.text(len(order)-0.5, 0.10, " 0.10", fontsize=4.8, color="#8B2E2E", va="bottom", ha="right")
axb.set_xticks(xpos)
axb.set_xticklabels([prog_label(p) for p in order], rotation=90, fontsize=4.2)
axb.set_xlim(-0.7, len(order)-0.3)
axb.set_ylabel("Cohort partial-$\\eta^2$")
axb.set_title("(b)  All 54 programs ranked by cohort partial-$\eta^2$  "
              "(outlined = F3 region-variable; color = brain-sig/weak)", loc="left", pad=4)
axb.margins(y=0.12)

# ---- (c) the 14 F3 region-variable programs: region vs cohort grouped bars ----
axc = fig.add_subplot(gs[2, 0])
s14 = tab[tab["f3_region_variable"]].sort_values("region_partial_eta2", ascending=False)
p14 = s14.index.tolist()
xp = np.arange(len(p14))
w = 0.40
axc.bar(xp - w/2, s14["region_partial_eta2"].values, width=w, color="#3C6E9C",
        label="region partial-$\\eta^2$", linewidth=0.2, edgecolor="white")
axc.bar(xp + w/2, s14["cohort_partial_eta2"].values, width=w, color="#D7642C",
        label="cohort partial-$\\eta^2$", linewidth=0.2, edgecolor="white")
axc.axhline(0.05, ls=":", lw=0.5, c="#555555")
axc.set_xticks(xp)
# Panel (c): use label + name_short for x-axis
# For excluded programs: "cNMF component N name"
# For kept programs: "P{new} name"
lbls = ["%s %s" % (prog_label(p), s14.loc[p, "name_short"]) for p in p14]
axc.set_xticklabels(lbls, rotation=40, ha="right", fontsize=4.6)
axc.set_ylabel("partial-$\\eta^2$")
# how many of the 14 are region-dominated vs cohort-dominated?
n_regdom = int((s14["region_partial_eta2"] > s14["cohort_partial_eta2"]).sum())
n_cohdom = len(s14) - n_regdom
axc.set_title("(c)  F3 region-variable programs: %d/14 region-dominated, "
              "%d/14 cohort-dominated (orange>blue = batch-sensitive, flag)"
              % (n_regdom, n_cohdom), loc="left", pad=4)
# mark the cohort-dominated (batch-sensitive) ones with a red dot above the pair
ymax_c = max(s14["region_partial_eta2"].max(), s14["cohort_partial_eta2"].max())
for i, p in enumerate(p14):
    if s14.loc[p, "cohort_partial_eta2"] > s14.loc[p, "region_partial_eta2"]:
        axc.plot(i, ymax_c * 1.04, marker="v", ms=3.2, color="#8B2E2E", zorder=5)
axc.legend(loc="upper right", frameon=False, fontsize=5.2, handletextpad=0.4)
axc.margins(y=0.12)

# fix_figS3 / editfig_r1: verdict footnote removed (relegated to the figure_release
# Figure legend). It previously spanned the full width at the figure bottom and
# was the cause of the right-edge clip; keeping it out also frees the bottom band
# for the panel-c rotated labels.

OUTPDF = os.path.join(OUTDIR, "ed_program_batch_check.pdf")
OUTPNG = os.path.join(OUTDIR, "ed_program_batch_check.png")
OUTSVG = os.path.join(OUTDIR, "ed_program_batch_check.svg")
fig.savefig(OUTPDF)
fig.savefig(OUTPNG, dpi=300)
fig.savefig(OUTSVG)  # fix_figS3: svg for SUBMISSION sync
print("[write]", OUTPDF, flush=True)
print("[write]", OUTPNG, flush=True)
print("[write]", OUTSVG, flush=True)
print("[done]", flush=True)
