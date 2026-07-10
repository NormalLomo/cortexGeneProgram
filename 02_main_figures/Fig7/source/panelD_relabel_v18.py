#!/usr/bin/env python
# v18 (single-source fix, 2026-06-15): regenerate ONLY panel d (== panels_v6/panel_f.svg, the
# mouse areal-variability lollipop) using the AUTHORITY naming table directly:
#   results/crossregion_v1/program_names.tsv  ->  name_short (+ "*" for brain-weak programs)
# Previously this read an intermediate curated tsv (program_names_curated_forfig_v18.tsv) whose
# P56 display differed from authority ("Blood vessel morphog." vs "Blood vessel morphogenesis").
# Now there is exactly ONE name source. Pure DISPLAY-STRING swap: same data, same ranks, same
# colours, same layout, same figure code. We import make_panels_v6 (defines functions + constants,
# no __main__ side effects) and override its NM map, then call panel_f() -> overwrites panel_f.svg.
#
# fig54甲 (2026-06-20): in ADDITION to the name swap, panel d is now ranked over the 54
# biologically-interpreted programs (RANK OF 54), not the full 60. The six cohort-technical
# programs (P9, P18, P19, P35, P52, P57) are dropped from the mouse areal-variability ranking
# denominator and the ranks are recomputed by descending area_median_sd (same metric/method as
# the of-60 ranks, validated rank-desc match). The top-14 colouring / dashed cut-off threshold
# stays at 14 (now read as "top 14 of 54"), the x-axis label reads "of 54", and the x-limit is
# tightened from 64 to 58 (54 * 64/60 ~= 58) to keep the same right-margin proportion. Everything
# else in the lollipop (figsize, lw, scatter size, #rank text, top-14 dashed line + caption box,
# y-labels = AUTHORITY brain names, subplots_adjust, save_svg) is BYTE-IDENTICAL to the original
# panel_f(). Only panel d (panel_f.svg) is regenerated; no other panels_v6 svg is touched.
import os, sys
import numpy as np
import pandas as pd

FIGDIR = "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1/spatial_xspecies/figures/Fig_spatial_univ"
AUTHORITY = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/program_names.tsv"
sys.path.insert(0, FIGDIR)
import make_panels_v6 as M  # noqa: E402  (top-level only defines; safe to import)

# AUTHORITY naming table is the single source of truth for display names + brain-weak star marks.
# v6_WorkerA fix (2026-06-24): program_names.tsv columns = new_P + cnmf_component (post-renumber).
# Build map keyed by "P{cnmf_component}" (the legacy display tag used in D_mouse_areal_universality.tsv
# and panel_f.svg) -> display label "{new_P} {name_short}{star}".
auth = pd.read_csv(AUTHORITY, sep="\t")
# keep only the 54 biologically-interpreted programs (drop EXCLUDED)
auth = auth[auth["new_P"].astype(str).str.upper() != "EXCLUDED"].copy()
def _disp(row):
    star = " *" if str(row["confidence"]).strip() == "brain-weak" else ""
    newP = str(row["new_P"])
    nm = str(row["name_short"])
    return newP + " " + nm + star
# keyed by LEGACY "P{cnmf_component}" tag (what D_mouse table + panel_f.svg use); display = new_P + name
# e.g. legacy "P14" (cnmf_component=14) -> "P13 Pos. reg. cation channel"
M.NM = {f"P{int(row['cnmf_component'])}": _disp(row) for _, row in auth.iterrows()}

# ----------------------------------------------------------------------------------------------
# fig54甲: RANK OF 54 override for panel d.
# Drop the six cohort-technical programs from the mouse areal-variability ranking denominator and
# recompute the rank over the remaining 54 by descending area_median_sd (the exact metric/method
# the of-60 mouse_var_rank column was built from). Then redefine M.panel_f to plot rank-of-54.
DROP6 = ["program_9", "program_18", "program_19", "program_35", "program_52", "program_57"]

def panel_f_of54():
    D = pd.read_csv(f"{M.AGG}/D_mouse_areal_universality.tsv", sep="\t")
    # rank of 54: drop the six cohort-technical programs, re-rank by descending area_median_sd
    sub = D[~D["program"].isin(DROP6)].copy()
    assert len(sub) == 54, f"expected 54 after drop, got {len(sub)}"
    sub["rank54"] = sub["area_median_sd"].rank(ascending=False, method="average")
    # the seven human robustly region-variable programs, ordered by their rank-of-54
    rob = sub[sub["human_robust_variable"] == True].copy().sort_values("rank54")
    assert len(rob) == 7, f"expected 7 human-robust programs, got {len(rob)}"
    # v6_WorkerA: M.NM[P_legacy] already returns "{new_P} {name_short}{star}", so use it directly
    # (avoid double-prefixing the legacy P-tag onto the already-prefixed new_P display string).
    rob["lab"] = rob["P"].map(lambda P: M.NM.get(P, P))
    fig, ax = M.plt.subplots(figsize=(3.0, 1.98))
    y = np.arange(len(rob))[::-1]; ranks = rob["rank54"].values
    top14 = ranks <= 14; cols = [M.C_MOUSE if t else M.C_GREY for t in top14]
    ax.hlines(y, 0, ranks, color=cols, lw=1.3, alpha=0.8, zorder=1)
    ax.scatter(ranks, y, s=34, c=cols, edgecolor="white", linewidth=0.6, zorder=3)
    for yi, r in zip(y, ranks):
        ax.text(r + 1.4, yi, f"#{int(r)}", va="center", fontsize=M.FS_SM, color="#333")
    ax.axvline(14, color="#c0492f", lw=0.7, ls="--", zorder=0)
    ax.text(14, len(rob)-0.35, "top-14\nmost area-variable", fontsize=M.FS_CAV, color="#c0492f",
            ha="center", va="bottom")
    ax.set_yticks(y); ax.set_yticklabels(rob["lab"], fontsize=M.FS_CAV)
    ax.set_xlabel("mouse areal-variability rank (1 = most variable, of 54)", fontsize=M.FS_SM)
    ax.set_xlim(0, 58); ax.set_ylim(-0.7, len(rob)+0.2)
    ax.set_title("Human region-variable programs in mouse areas", fontsize=M.FS_TITLE,
                 fontweight="bold", pad=12)
    n_top = int(top14.sum())
    ax.text(0.985, 0.78, f"{n_top}/7 land in mouse top-14\nall 7 area-significant (ANOVA F≫1)",
            transform=ax.transAxes, fontsize=M.FS_CAV, ha="right", va="center",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#999", lw=0.4))
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    fig.subplots_adjust(left=0.50, right=0.97, top=0.85, bottom=0.18)
    M.save_svg(fig, "panel_f")
    print(f"panel d rank-of-54: P8 rank={rob.loc[rob.P=='P8','rank54'].iloc[0]:.0f}, "
          f"n_top14={n_top}/7 ({list(rob.loc[top14,'P'])})")

# regenerate ONLY panel_f.svg (figure panel d) with AUTHORITY brain labels + rank of 54.
panel_f_of54()
print("DONE: regenerated panels_v6/panel_f.svg (panel d) with AUTHORITY brain labels + RANK OF 54")
