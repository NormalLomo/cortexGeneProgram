#!/usr/bin/env python
"""ED Fig 1 (renumbered 2026-06-20): cross-chip spatial reproducibility of 14 representative programs (SCT smoothed).

Layout: 14 programs (ROWS) x 14 chips (COLS) of per-bin x,y scatters
colored by per-chip KNN-median smoothed program-z.

Per-program shared color limits (robust 2-98 pct across all chips).
Each cell ~25x25 mm. Title with method.

Renumber: display label P{old} → P{new} per results/crossregion_v1/program_renumber_map.tsv
Parquet keys stay as program_{cNMF_component} = program_{old_P}.
Excluded cNMF components: 9/18/19/35/52/57 — none appear in this figure's 14 programs.
old_P → new_P: 5→5, 8→8, 34→31, 7→7, 20→17, 29→26, 37→33, 43→39, 15→14, 51→47, 53→48, 40→36, 54→49, 56→51
"""
import os, sys, time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import matplotlib
matplotlib.use("Agg")
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
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from sklearn.neighbors import NearestNeighbors

T0 = time.time()
ROOT = "CORTEX_PROGRAM_ROOT"
PQ = f"{ROOT}/results/crossregion_v1/spatial_bin50_program_score_SCT.parquet"
OUT_DIR = f"{ROOT}/figures/extended"
os.makedirs(OUT_DIR, exist_ok=True)

# ===== 14 representative programs (row order: EX laminar -> INH -> glia/vascular) =====
# Selected for breadth across cell classes + high cross-chip spatial reproducibility
# (spatial_score from fig2b) + distinct spatial signature (program-program z-corr).
#
# RENUMBER 2026-06-20: display labels updated (P{old} → P{new}) per program_renumber_map.tsv.
# Parquet keys = program_{cNMF_component} = program_{old_P} — unchanged.
# old_P(cNMF) → new_P(display): 5→5, 8→8, 34→31, 7→7, 20→17, 29→26, 37→33, 43→39,
#                                 15→14, 51→47, 53→48, 40→36, 54→49, 56→51
PROGS = [
    # --- excitatory layer programs (superficial -> deep) ---
    ("program_5",  "P5 L2/3-IT (LINC00507)",      "EX"),    # old P5  → new P5  (unchanged)
    ("program_8",  "P8 L4/5-IT RORB",             "EX"),    # old P8  → new P8  (unchanged)
    ("program_34", "P31 L6-IT / activity-IEG",    "EX"),    # old P34 → new P31
    # --- inhibitory subclasses ---
    ("program_7",  "P7 PVALB interneuron",        "INH"),   # old P7  → new P7  (unchanged)
    ("program_20", "P17 SST interneuron",         "INH"),   # old P20 → new P17
    ("program_29", "P26 LAMP5 interneuron*",      "INH"),   # old P29 → new P26
    # --- oligo / OPC ---
    ("program_37", "P33 Oligodendrocyte/myelin",  "OLIGO"), # old P37 → new P33
    ("program_43", "P39 OPC / ECM proteoglycan",  "OPC"),   # old P43 → new P39
    # --- astrocyte ---
    ("program_15", "P14 Astrocyte Glu transport", "AST"),   # old P15 → new P14
    ("program_51", "P47 Reactive astro/vascular", "AST"),   # old P51 → new P47
    ("program_53", "P48 Astrocyte metal homeost.","AST"),   # old P53 → new P48
    # --- microglia ---
    ("program_40", "P36 Microglial activation",   "MG"),    # old P40 → new P36
    ("program_54", "P49 Microglial complement/MHC","MG"),   # old P54 → new P49
    # --- endothelial / vascular ---
    ("program_56", "P51 Blood vessel morphog.",   "ENDO"),  # old P56 → new P51
]
PROG_KEYS = [p[0] for p in PROGS]
ROW_LABELS = [p[1] for p in PROGS]

# ===== chip picks (one per 14 regions) =====
CHIPS = [
    ("D00865A2",     "ACC"),
    ("B02222C3",     "AG"),
    ("SS200001075BR","DLPFC"),
    ("B01012B5",     "FPPFC"),
    ("B02223C3",     "ITG"),
    ("C02137B1",     "M1"),
    ("B02221E3",     "PoCG"),
    ("B02111D1",     "S1"),
    ("B02203D5",     "S1E"),
    ("A01186B3",     "SMG"),
    ("D00865C3",     "SPL"),
    ("B02203E1",     "STG"),
    ("B01012C5",     "V1"),
    ("B01012B1",     "VLPFC"),
]
CHIP_IDS = [c[0] for c in CHIPS]
COL_LABELS = [f"{c[0]}\n{c[1]}" for c in CHIPS]
N_ROW, N_COL = len(PROGS), len(CHIPS)
K = 25

# ===== 1) Load only what we need =====
print(f"[{time.time()-T0:6.1f}s] Loading parquet (bin,x,y + {len(PROG_KEYS)} progs)...", flush=True)
cols = ["bin","x","y"] + PROG_KEYS
df_all = pq.read_table(PQ, columns=cols).to_pandas()
df_all["chip"] = df_all["bin"].astype(str).str.split("_").str[0]
df = df_all[df_all["chip"].isin(CHIP_IDS)].copy()
print(f"   loaded {len(df_all):,} bins, kept {len(df):,} in {df['chip'].nunique()} chips",
      flush=True)
del df_all

# ===== 2) per-chip KNN-median smoothing =====
print(f"[{time.time()-T0:6.1f}s] Per-chip KNN(k={K}) median smoothing...", flush=True)
sm_parts = []
chip_dims = {}  # chip -> (xspan, yspan)
for chip in CHIP_IDS:
    g = df[df["chip"]==chip].reset_index(drop=True)
    if len(g) < K+5:
        print(f"   SKIP chip {chip}: only {len(g)} bins (too sparse)", flush=True)
        continue
    xy = g[["x","y"]].to_numpy(float)
    chip_dims[chip] = (np.ptp(xy[:,0]), np.ptp(xy[:,1]))
    nn = NearestNeighbors(n_neighbors=K, n_jobs=-1).fit(xy)
    _, idx = nn.kneighbors(xy)
    sm = g[["bin","chip","x","y"]].copy()
    for p in PROG_KEYS:
        v = g[p].to_numpy(float)
        sm[p] = np.median(v[idx], axis=1)
    sm_parts.append(sm)
    print(f"   chip {chip}: {len(g):,} bins, xspan={chip_dims[chip][0]:.0f} yspan={chip_dims[chip][1]:.0f}",
          flush=True)
sm = pd.concat(sm_parts, ignore_index=True)
print(f"[{time.time()-T0:6.1f}s] smoothed table: {len(sm):,} rows", flush=True)

# ===== 3) per-program color limits (robust 2-98 pct over all kept chips) =====
prog_lim = {}
for p in PROG_KEYS:
    v = sm[p].to_numpy()
    lo, hi = np.percentile(v, [2, 98])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < 1e-9:
        lo, hi = float(v.min()), float(v.max())
    prog_lim[p] = (lo, hi)
    print(f"   limits {p}: [{lo:.3f}, {hi:.3f}]", flush=True)

# ===== 4) figure =====
# 14 cols x 14 rows, each panel ~25mm = 0.984 in. Add headers/labels space.
CELL_IN = 25.0/25.4
LEFT_LABEL_IN = 1.95   # row labels (longer 'P{n} short-name' strings)
TOP_HEAD_IN   = 0.85   # column headers + title
RIGHT_CBAR_IN = 0.85
BOT_PAD_IN    = 0.35
FIG_W = LEFT_LABEL_IN + N_COL*CELL_IN + RIGHT_CBAR_IN
FIG_H = TOP_HEAD_IN + N_ROW*CELL_IN + BOT_PAD_IN
print(f"[{time.time()-T0:6.1f}s] page = {FIG_W:.2f}in x {FIG_H:.2f}in (AR={FIG_W/FIG_H:.2f})",
      flush=True)

fig = plt.figure(figsize=(FIG_W, FIG_H))

# layout via fractional axes (cleaner than gridspec for fixed-mm cells)
def cell_rect(ci, ri):
    x0 = (LEFT_LABEL_IN + ci*CELL_IN) / FIG_W
    y0 = 1 - (TOP_HEAD_IN + (ri+1)*CELL_IN) / FIG_H
    w  = CELL_IN / FIG_W
    h  = CELL_IN / FIG_H
    return [x0, y0, w, h]

# Title
fig.text(LEFT_LABEL_IN/FIG_W, 1 - 0.10/FIG_H,
         "ED Fig 1: cross-chip spatial reproducibility of representative programs (SCT smoothed)",
         ha="left", va="top", fontsize=8, weight="bold")

# Column headers
for ci, (chip, region) in enumerate(CHIPS):
    rect = cell_rect(ci, 0)
    fig.text(rect[0] + rect[2]/2, 1 - 0.30/FIG_H,
             f"{chip}\n{region}", ha="center", va="top", fontsize=5.5)

# Row labels
for ri, lab in enumerate(ROW_LABELS):
    rect = cell_rect(0, ri)
    fig.text((LEFT_LABEL_IN - 0.10)/FIG_W, rect[1] + rect[3]/2,
             lab, ha="right", va="center", fontsize=5.5, weight="bold")

# Colorbar axes (one per row, on the right)
CBAR_W = 0.18 / FIG_W
for ri, p in enumerate(PROG_KEYS):
    rect = cell_rect(N_COL-1, ri)
    cax_x = rect[0] + rect[2] + 0.10/FIG_W
    cax = fig.add_axes([cax_x, rect[1] + rect[3]*0.10,
                        CBAR_W, rect[3]*0.80])
    lo, hi = prog_lim[p]
    sm_map = ScalarMappable(norm=Normalize(vmin=lo, vmax=hi), cmap="magma")
    sm_map.set_array([])
    cb = fig.colorbar(sm_map, cax=cax)
    cb.ax.tick_params(labelsize=4.5, length=1.2, pad=1)
    cb.outline.set_linewidth(0.3)
    cb.set_ticks([lo, (lo+hi)/2, hi])

# Plot cells
for ci, chip in enumerate(CHIP_IDS):
    g = sm[sm["chip"]==chip]
    if len(g)==0:
        for ri in range(N_ROW):
            ax = fig.add_axes(cell_rect(ci, ri))
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values(): s.set_visible(False)
            ax.text(0.5,0.5,"n/a",ha="center",va="center",fontsize=5,transform=ax.transAxes)
        continue
    x = g["x"].to_numpy()
    y = g["y"].to_numpy()
    # decide point size from local density (cells are 25mm; use small s)
    s_pt = 0.30
    for ri, p in enumerate(PROG_KEYS):
        ax = fig.add_axes(cell_rect(ci, ri))
        v = g[p].to_numpy()
        lo, hi = prog_lim[p]
        # flip y so dorsal up matches typical orientation
        ax.scatter(x, -y, c=v, s=s_pt, cmap="magma", vmin=lo, vmax=hi,
                   marker="s", linewidths=0, rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal")
        for s in ax.spines.values():
            s.set_linewidth(0.3); s.set_color("#444")

# ===== 5) save =====
pdf = f"{OUT_DIR}/ed_fig1_cross_chip_atlas_current_renum.pdf"
png = f"{OUT_DIR}/ed_fig1_cross_chip_atlas_current_renum.png"
print(f"[{time.time()-T0:6.1f}s] writing {pdf}", flush=True)
fig.savefig(pdf, dpi=300, bbox_inches=None)
fig.savefig(png, dpi=150, bbox_inches=None)
plt.close(fig)
print(f"[{time.time()-T0:6.1f}s] DONE. AR={FIG_W/FIG_H:.2f}  size={FIG_W:.2f}x{FIG_H:.2f} in",
      flush=True)
