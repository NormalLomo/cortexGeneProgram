#!/usr/bin/env python3
"""
ED13 | cNMF K-selection diagnostic.

Polished, CNS-grade rendering of the native cNMF k_selection_stats:
  panel a — solution stability (consensus silhouette) vs K
  panel b — reconstruction (prediction) error vs K
K = 60 (the chosen rank) is highlighted in both panels.

Reads the cNMF native stats npz (numpy-2 written) with the cellist env.
Outputs vector PDF + PNG to figures/extended/ed_kselection.{pdf,png}.

Run:  python scripts/extended/ed_kselection.py
"""
import pathlib
import numpy as np
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
from matplotlib import font_manager as fm

ROOT = pathlib.Path("CORTEX_PROGRAM_ROOT")
NPZ = ROOT / ("results/cnmf_snrna_joint_full1M_v1/cnmf_work/"
              "snrna_joint_full1M_v1/snrna_joint_full1M_v1.k_selection_stats.df.npz")
OUT = ROOT / "figures/extended"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- data
d = np.load(NPZ, allow_pickle=True)
cols = list(d["columns"])
data = d["data"]
col = {c: data[:, i] for i, c in enumerate(cols)}
K = col["k"].astype(int)
sil = col["silhouette"].astype(float)            # = consensus stability
err = col["prediction_error"].astype(float)
err_b = err / 1e9                                 # display in 1e9 units
K_CHOSEN = 60
ci = int(np.where(K == K_CHOSEN)[0][0])

# ---------------------------------------------------------------- style
# Liberation Sans is metric-compatible with Arial; clean journal sans.
avail = {f.name for f in fm.fontManager.ttflist}
fam = "Nimbus Sans" if "Nimbus Sans" in avail else (
    "Liberation Sans" if "Liberation Sans" in avail else "DejaVu Sans")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "pdf.fonttype": 42,   # embed TrueType (editable vector text)
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

INK = "#1a1a1a"
C_STAB = "#2c6fbb"      # stability blue
C_ERR = "#c0392b"       # error red
C_HL = "#e8a800"        # highlight gold for K=60
GRID = "#d9d9d9"

# square-ish balanced panels, two side by side
fig, (axA, axB) = plt.subplots(
    1, 2, figsize=(6.6, 3.05), dpi=300,
    gridspec_kw=dict(wspace=0.34, left=0.085, right=0.985, top=0.86, bottom=0.155))


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK)
    ax.tick_params(colors=INK)
    ax.set_xticks(K)
    ax.set_xlim(K.min() - 5, K.max() + 5)
    ax.set_xlabel("Number of programs (K)")
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def highlight(ax, x, y, color, txt, dy, ha="center", xoff=0.0, va="bottom"):
    ax.axvline(K_CHOSEN, color=C_HL, lw=1.1, ls=(0, (4, 2)), zorder=1.5)
    ax.scatter([x], [y], s=70, facecolor="none", edgecolor=C_HL,
               linewidth=1.6, zorder=6)
    ax.annotate(txt, (x, y), xytext=(x + xoff, y + dy), ha=ha, va=va,
                fontsize=6.5, fontweight="bold", color="#8a6500", zorder=7)


# ---- panel a: stability (silhouette) ----
axA.plot(K, sil, "-", color=C_STAB, lw=1.4, zorder=3)
axA.scatter(K, sil, s=24, color=C_STAB, zorder=4, edgecolor="white", linewidth=0.6)
style_axis(axA)
axA.set_ylabel("Consensus stability\n(silhouette)", color=C_STAB)
axA.tick_params(axis="y", colors=C_STAB)
axA.spines["left"].set_color(C_STAB)
yr = sil.max() - sil.min()
axA.set_ylim(sil.min() - 0.10 * yr, sil.max() + 0.26 * yr)
highlight(axA, K_CHOSEN, sil[ci], C_STAB,
          f"K = 60 (chosen)\nstability {sil[ci]:.3f}", 0.085 * yr,
          ha="left", xoff=2.0)
# note the K=30 global max so the choice is defended honestly
imax = int(np.argmax(sil))
axA.annotate(f"global max at K = {K[imax]}\n(lower resolution)",
             (K[imax], sil[imax]), xytext=(K[imax] + 4.5, sil[imax] + 0.004),
             ha="left", va="center", fontsize=5.6, color="#555",
             arrowprops=dict(arrowstyle="-", color="#999", lw=0.6,
                             shrinkA=2, shrinkB=2))
axA.set_title("a   Solution stability vs K", loc="left",
              fontweight="bold", fontsize=8.5, color=INK)

# ---- panel b: prediction error ----
axB.plot(K, err_b, "-", color=C_ERR, lw=1.4, zorder=3)
axB.scatter(K, err_b, s=24, color=C_ERR, zorder=4, edgecolor="white", linewidth=0.6)
style_axis(axB)
axB.set_ylabel(r"Reconstruction error ($\times10^{9}$)", color=C_ERR)
axB.tick_params(axis="y", colors=C_ERR)
axB.spines["left"].set_color(C_ERR)
er = err_b.max() - err_b.min()
axB.set_ylim(err_b.min() - 0.12 * er, err_b.max() + 0.18 * er)
highlight(axB, K_CHOSEN, err_b[ci], C_ERR,
          f"K = 60 (chosen)\nerror {err_b[ci]:.3f}", -0.14 * er,
          ha="left", xoff=2.0, va="top")
axB.annotate("monotonically decreasing\n(no elbow; favours larger K)",
             (K[1], err_b[1]), xytext=(K[2] + 1, err_b[1] + 0.012),
             ha="left", va="bottom", fontsize=5.6, color="#555")
axB.set_title("b   Reconstruction error vs K", loc="left",
              fontweight="bold", fontsize=8.5, color=INK)

fig.savefig(OUT / "ed_kselection.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig(OUT / "ed_kselection.png", bbox_inches="tight", pad_inches=0.02, dpi=300)
print("WROTE", OUT / "ed_kselection.pdf")
print("WROTE", OUT / "ed_kselection.png")
print("font:", fam)
print("K     :", K.tolist())
print("sil   :", [round(x, 3) for x in sil])
print("err/1e9:", [round(x, 3) for x in err_b])
