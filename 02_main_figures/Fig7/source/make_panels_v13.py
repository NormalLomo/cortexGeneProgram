#!/usr/bin/env python
"""Fig_spatial_univ v13 — same panels as v12 (b 1x4 stacked column; e 1x6 wide strip).
v13 change vs v12: panel e/spatial figsize shrunk further (6.0x1.20 -> 6.2x0.98) so the
bottom band is a bit SHORTER (user: "底部 e 再稍微小一點"). Panel b UNCHANGED (the right
column is enlarged in compose_v13.py via svgutils proportional scale, so b/c/d content +
fonts grow together; no figsize change needed for b).
Scientific content / data / colors / numbers COPIED VERBATIM from make_panels_v12.py —
only matplotlib figsize change for e.

Panel a (ComplexHeatmap) untouched (re-used from panels/panelA_complexhm_v6_crop.svg via compose).
Panels c (co-org) and d (areal lollipop) untouched, re-used from panels_v6/.

All-English, zero CJK. CNS teal/orange palette.
"""
import os, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
warnings.filterwarnings("ignore")

BASE = "__PRIVATE_CANONICAL_ROOT__/results/xspecies_humanmap_v1/spatial_xspecies"
AGG  = f"{BASE}/_aggregate"
HUMAN_DIR = "__PRIVATE_CANONICAL_ROOT__/results/crossregion_v1/_spatial_score_perchip"
OUT  = f"{BASE}/figures/Fig_spatial_univ/panels_v13"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "svg.fonttype": "none", "font.family": "Liberation Sans", "font.size": 6,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "axes.edgecolor": "#444444",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
FS_TITLE, FS_LAB, FS_TICK, FS_SM = 7.5, 6.5, 5.5, 5.2
FS_CAV = 5.0
C_HUMAN  = "#1f6f78"
C_MONKEY = "#5ba3aa"
C_MOUSE  = "#e08a3c"
SPECIES_C = {"human": C_HUMAN, "macaque": C_MONKEY, "mouse": C_MOUSE}

LAYERS = ["L1", "L2", "L3", "L4", "L5", "L6"]

def save_svg(fig, name):
    p = f"{OUT}/{name}.svg"
    fig.savefig(p, format="svg", bbox_inches="tight", pad_inches=0.01, transparent=True)
    plt.close(fig); print("wrote", p)

# ====================================================================== PANEL b (1x4 stacked column)
# IDENTICAL data + numbers to make_panels_v11.panel_b; only 2x2 grid -> 1x4 vertical stack
# (tall + slim, fills the narrow right column).
def panel_b():
    h = pd.read_csv(f"{AGG}/human_program_x_layer_baselinerm.tsv", sep="\t", index_col=0)
    mk = pd.read_csv(f"{AGG}/monkey_program_x_layer_baselinerm.tsv", sep="\t", index_col=0)
    mo = pd.read_csv(f"{AGG}/mouse_program_x_layer_baselinerm.tsv", sep="\t", index_col=0)
    MYELIN_FIX = ["program_45", "program_26", "program_37", "program_38"]
    for p in MYELIN_FIX:
        if p in mo.index and "L1" in mo.columns:
            mo.loc[p, "L1"] = mo.loc[p, ["L2", "L4", "L5"]].min()
    # v6_WorkerA fix (2026-06-24): data column "program_22" etc. = cnmf_component (1-60).
    # Map cnmf_component -> new_P via program_names.tsv authority:
    #   cnmf=22 -> P19 Axon guidance/neural crest
    #   cnmf=31 -> P28 Axon guidance / cell adhesion (L6b)
    #   cnmf=45 -> P41 Myelination
    #   cnmf=26 -> P23 Reg. myelination
    progs = [("P19", "program_22", "Axon guidance/neural crest"),
             ("P28", "program_31", "Axon guid./adh. (L6b)"),
             ("P41", "program_45", "Myelination"),
             ("P23", "program_26", "Reg. myelination")]
    dfs = {"human": (h, C_HUMAN), "macaque": (mk, C_MONKEY), "mouse": (mo, C_MOUSE)}

    # 1x4 vertical stack; slim + tall so it fills the right column. share x & y.
    fig, axes = plt.subplots(4, 1, figsize=(1.55, 4.05), sharey=True, sharex=True)
    axes = axes.ravel()
    xr = np.arange(len(LAYERS))
    for ax, (Pid, prog, kind) in zip(axes, progs):
        for sp, (df, col) in dfs.items():
            if prog not in df.index: continue
            y = df.loc[prog, LAYERS].astype(float).values
            if np.isnan(y).any():
                idx = np.arange(len(y)); good = ~np.isnan(y)
                y = np.interp(idx, idx[good], y[good])
            lw = 1.3 if sp != "macaque" else 1.0
            a = 0.95 if sp != "macaque" else 0.65
            ax.plot(xr, y, color=col, lw=lw, alpha=a, marker="o", ms=1.8,
                    mec="white", mew=0.3, zorder=3 if sp != "macaque" else 2)
        ax.axhline(0, color="#bbb", lw=0.5, ls="--", zorder=0)
        ax.axvspan(3.5, 5.5, color="#1f6f78", alpha=0.05, zorder=0)
        ax.set_xticks(xr); ax.set_xticklabels(LAYERS, fontsize=FS_CAV, rotation=0)
        ax.tick_params(axis="x", length=1.5, pad=1)
        ax.set_xlim(-0.4, 5.4)
        ax.set_title(f"{Pid} · {kind}", fontsize=FS_SM, fontweight="bold", pad=2, color="#222")
        # NOTE (2026-06-25): "mouse L1 RAW-corrected" annotation moved to figure legend (figure_release_main_text.md).
        for s in ["top", "right"]: ax.spines[s].set_visible(False)
        # y label on every axis (single column)
        ax.set_ylabel("laminar score\n(depth-corr.)", fontsize=FS_CAV)
        ax.tick_params(axis="y", length=1.5, labelsize=FS_CAV)
    leg = [Line2D([0],[0], color=SPECIES_C[s], lw=1.5, label=s) for s in ["human","macaque","mouse"]]
    fig.legend(handles=leg, fontsize=FS_CAV, loc="upper center", bbox_to_anchor=(0.55, 1.005),
               ncol=3, frameon=False, handlelength=1.1, handletextpad=0.3, columnspacing=0.9)
    fig.suptitle("Representative laminar profiles:\nexemplar cross-species (deep-enriched programs)",
                 fontsize=5.6, fontweight="bold", y=1.052, x=0.55, linespacing=1.05)
    fig.subplots_adjust(left=0.205, right=0.97, top=0.915, bottom=0.045, hspace=0.55)
    save_svg(fig, "panel_b")

# ====================================================================== PANEL e (1x6 wide strip, SHRUNK)
# IDENTICAL data + per-field z-scoring to make_panels_v11.panel_e_spatial; only figsize shrunk
# (7.0x1.55 -> 6.0x1.20) so the bottom band is shorter. Grid/cmap/numbers verbatim.
def panel_e_spatial():
    HUMAN_F  = f"{HUMAN_DIR}/B01012B2.parquet"
    MONKEY_F = f"{BASE}/monkey/results_sct/sections_bin50/T39.parquet"
    MOUSE_F  = f"{BASE}/mouse/sections_sct/mouse1__T280.parquet"
    # v6_WorkerA fix (2026-06-24): cnmf_component -> new_P
    #   cnmf=22 -> P19 Axon guidance/neural crest (deep cortex enriched in data)
    #   cnmf=45 -> P41 Myelination (deep/myelin enriched in data)
    progs = [("P19", "program_22", "Axon guid./neural crest"), ("P41", "program_45", "Myelination")]
    def load(f, xcol, ycol):
        df = pd.read_parquet(f)
        if "rctd_pass_mask" in df.columns: df = df[df["rctd_pass_mask"] == True]
        return df, df[xcol].values, df[ycol].values
    hdf, hx, hy = load(HUMAN_F, "x", "y")
    mkdf, mkx, mky = load(MONKEY_F, "bin_x", "bin_y")
    modf, mox, moy = load(MOUSE_F, "bin_x", "bin_y")
    specs = [("human", hdf, hx, hy, C_HUMAN),
             ("macaque", mkdf, mkx, mky, C_MONKEY),
             ("mouse", modf, mox, moy, C_MOUSE)]
    FIELD = LinearSegmentedColormap.from_list("field", ["#2a5d74", "#7fb2c0", "#eef0ee",
                                                        "#f0b27a", "#c0492f"])
    # 1 row x 6 cols, wide & short -> full-width bottom band. SHRUNK further for v13.
    # column order: [P22 hu, P22 mq, P22 mo, P45 hu, P45 mq, P45 mo]
    fig = plt.figure(figsize=(6.2, 0.98))
    gs = fig.add_gridspec(1, 6, wspace=0.05)
    sc = None
    col_i = 0
    for (Pid, pcol, klab) in progs:
        for (sp, df, xx, yy, col) in specs:
            ax = fig.add_subplot(gs[0, col_i]); col_i += 1
            raw = df[pcol].astype(float).values
            mu, sd = np.nanmean(raw), np.nanstd(raw)
            v = (raw - mu) / (sd if sd > 0 else 1.0)
            lim = np.nanpercentile(np.abs(v), 97); v = np.clip(v, -lim, lim)
            order = np.argsort(np.abs(v))
            sc = ax.scatter(xx[order], yy[order], c=v[order], cmap=FIELD,
                            vmin=-lim, vmax=lim, s=0.5, linewidths=0, rasterized=True)
            ax.set_aspect("equal"); ax.axis("off"); ax.invert_yaxis()
            ax.set_title(sp, fontsize=FS_CAV, color=col, fontweight="bold", pad=1)
            if sp == "human":
                ax.text(0.5, -0.06, f"{Pid} · {klab}", transform=ax.transAxes,
                        ha="left", va="top", fontsize=FS_SM, fontweight="bold", color="#222")
    cax = fig.add_axes([0.40, 0.015, 0.20, 0.030])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_ticks([-lim, lim]); cb.set_ticklabels(["low", "high"], fontsize=FS_CAV)
    cb.set_label("program score (cross-bin z)", fontsize=5.0, labelpad=1)
    cb.outline.set_linewidth(0.4); cb.ax.tick_params(length=1.5)
    fig.text(0.5, 0.995, "Representative program spatial fields", ha="center", va="top",
             fontsize=6.0, fontweight="bold")
    fig.subplots_adjust(left=0.005, right=0.995, top=0.86, bottom=0.10)
    save_svg(fig, "panel_e_spatial")

# ====================================================================== PANEL aggregate (NEW 2026-06-25)
# Aggregate cross-species same-sign fraction + 95% Wilson CI.
# Data: betweenchip_xspecies_aggregate_samesign.tsv (mouse 93.7%, macaque 85.7%)
#       betweenchip_xspecies_effectsize_ci.tsv (top pairs per species)
BTWNDIR = "__PRIVATE_CANONICAL_ROOT__/results/crossregion_v1/markcorr_betweenchip_v1"
RENUM_PATH = "__PRIVATE_CANONICAL_ROOT__/results/crossregion_v1/program_renumber_map.tsv"

def panel_aggregate():
    agg = pd.read_csv(f"{BTWNDIR}/betweenchip_xspecies_aggregate_samesign.tsv", sep="\t")
    eff = pd.read_csv(f"{BTWNDIR}/betweenchip_xspecies_effectsize_ci.tsv", sep="\t")
    # Load renumber map to get new_P labels
    renum = pd.read_csv(RENUM_PATH, sep="\t")
    old2new_p = {}
    for _, row in renum.iterrows():
        op = int(row["old_P"])
        np_val = str(row["new_P"])
        if np_val.upper() != "EXCLUDED":
            old2new_p[f"program_{op}"] = f"P{int(float(np_val))}"
    def prog_label(pname):
        return old2new_p.get(pname, pname)

    # Top off-diagonal pairs per species (A_name != B_name)
    top_pairs = {}
    for sp in ["mouse", "macaque"]:
        sub = eff[(eff["species"] == sp) & (eff["A_name"] != eff["B_name"])].copy()
        sub = sub.sort_values("median_log2g", ascending=False).head(3)
        top_pairs[sp] = [(prog_label(r["A_name"]), prog_label(r["B_name"]),
                         r["median_log2g"], r["frac_same_sign"]) for _, r in sub.iterrows()]

    # Plot: two horizontal bars (mouse, macaque) with Wilson CI error bars
    fig, ax = plt.subplots(figsize=(2.0, 1.40))
    colors = {"mouse": C_MOUSE, "macaque": C_MONKEY}
    species_order = ["mouse", "macaque"]
    y_pos = [1, 0]
    for yi, sp in zip(y_pos, species_order):
        row = agg[agg["species"] == sp].iloc[0]
        f = row["frac_same_sign_consistent"]
        lo = row["wilson_ci_lo"]
        hi = row["wilson_ci_hi"]
        ax.barh(yi, f, height=0.50, color=colors[sp], alpha=0.85, zorder=3)
        ax.errorbar(f, yi, xerr=[[f - lo], [hi - f]], fmt="none",
                    ecolor="#333", elinewidth=0.8, capsize=2.0, capthick=0.8, zorder=4)
        pct = f"{f * 100:.1f}%"
        ax.text(f + 0.002, yi, pct, va="center", ha="left",
                fontsize=FS_TICK, color="#222", fontweight="bold")
        # annotate top pair
        if top_pairs[sp]:
            pa, pb, g, fsign = top_pairs[sp][0]
            ax.text(lo - 0.005, yi - 0.22,
                    f"top: {pa}×{pb} g={g:.2f}",
                    va="top", ha="left", fontsize=4.5, color=colors[sp], style="italic")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([sp.capitalize() for sp in species_order], fontsize=FS_LAB)
    ax.set_xlim(0.70, 1.00)
    ax.axvline(0.5, color="grey", lw=0.4, ls="--", zorder=1)
    ax.set_xlabel("Same-sign fraction (95% Wilson CI)", fontsize=FS_SM)
    ax.set_title("Cross-species co-org conservation", fontsize=FS_TITLE, fontweight="bold", pad=4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=FS_TICK)
    ax.set_ylim(-0.55, 1.55)
    fig.subplots_adjust(left=0.20, right=0.90, top=0.82, bottom=0.22)
    note = ("Aggregate same-sign fraction across 1,485 program pairs;\n"
            "single pair p-values not reported for cross-species;\n"
            "only aggregate fraction with 95% Wilson CI.")
    fig.text(0.5, -0.02, note, ha="center", va="top",
             fontsize=3.8, color="#555", style="italic", wrap=True,
             transform=fig.transFigure)
    save_svg(fig, "panel_aggregate")

if __name__ == "__main__":
    panel_b()
    panel_e_spatial()
    panel_aggregate()
    print("v13 PANELS (b 1x4 stack, e figsize shrunk, + aggregate panel) DONE")
