#!/usr/bin/env python
"""
Extended Data Fig 4 / Supplementary Data 03:
54-program GO:BP annotation reference table (excluded: old P9/P18/P19/P35/P52/P57).

Identical layout to the 60-program original; filters to kept=54 programs
and uses new_P numbering (P1-P54) from program_names.tsv.

Outputs:
  - figures/extended/ed_fig4_program_table_54.{pdf,png}
"""
import os
import re
import numpy as np
import pandas as pd
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
from matplotlib.patches import Rectangle
import matplotlib.font_manager as fm


ROOT = "CORTEX_PROGRAM_ROOT"
CR = os.path.join(ROOT, "results", "crossregion_v1")
FIGDIR = os.path.join(ROOT, "figures", "extended")
os.makedirs(FIGDIR, exist_ok=True)

# ---------------------------------------------------------------- load
names = pd.read_csv(os.path.join(CR, "program_names.tsv"), sep="\t")
# bio-MAJ1 fix 2026-06-25: program_annotation.tsv P01 stale (inh/SST);
# use TableS1_program_annotation.tsv which has correct dominant_class/dominant_subclass
SUPP_DIR = "CORTEX_PROGRAM_ROOT/figure_release/SUBMISSION_final/supplementary"
ann   = pd.read_csv(os.path.join(SUPP_DIR, "TableS1_program_annotation.tsv"), sep="	")
gobp  = pd.read_csv(os.path.join(CR, "program_annotation_gobp.tsv"), sep="\t")

# ---------------------------------------------------------------- filter: keep 54
# new_P == "EXCLUDED" marks the 6 cohort-technical programs
names_kept = names[names["new_P"] != "EXCLUDED"].copy().reset_index(drop=True)
assert len(names_kept) == 54, f"expected 54 kept programs, got {len(names_kept)}"
excluded_old = names[names["new_P"] == "EXCLUDED"]["cnmf_component"].tolist()
# excluded_old should be [9, 18, 19, 35, 52, 57]
print(f"Excluded old-component ids: {sorted(excluded_old)}")

# ---------------------------------------------------------------- join on cnmf_component (=old integer id)
# ann has program column in "P01" format (old numbering); convert to integer
def old_pid(x):
    return int(re.sub(r"[^0-9]", "", str(x)))

# TableS1 already has cnmf_component column; old conversion not needed (bio-MAJ1 fix)
# gobp has integer program column
gobp["cnmf_component"] = gobp["program"].astype(int)

# filter excluded from ann and gobp
ann_kept  = ann[~ann["cnmf_component"].isin(excluded_old)].copy()
gobp_kept = gobp[~gobp["cnmf_component"].isin(excluded_old)].copy()
assert len(ann_kept)  == 54, f"expected 54 ann rows, got {len(ann_kept)}"
assert len(gobp_kept) == 54, f"expected 54 gobp rows, got {len(gobp_kept)}"

# ---------------------------------------------------------------- split GO:BP term text
def split_term(s):
    s = str(s)
    m = re.search(r"\(GO:(\d+)\)", s)
    go = ("GO:" + m.group(1)) if m else ""
    txt = re.sub(r"\s*\(GO:\d+\)\s*$", "", s).strip()
    return txt, go

_split = [split_term(s) for s in gobp_kept["top_BP_term"]]
gobp_kept = gobp_kept.copy()
gobp_kept["gobp_term"] = [t for (t, g) in _split]
gobp_kept["gobp_id"]   = [g for (t, g) in _split]

def top8_genes(s):
    parts = re.split(r"[;,]", str(s))
    parts = [p.strip() for p in parts if p.strip()]
    return ", ".join(parts[:8])

gobp_kept["top8_loading_genes"] = gobp_kept["top10_loading_genes"].map(top8_genes)

# ---------------------------------------------------------------- merge
m = names_kept.merge(
    ann_kept[["cnmf_component", "dominant_subclass", "dominant_class"]].rename(columns={"dominant_class": "class"}), on="cnmf_component", how="left"
).merge(
    gobp_kept[["cnmf_component", "gobp_term", "gobp_id", "top8_loading_genes"]],
    on="cnmf_component", how="left"
)
# sort by new P number (extract integer from new_P like "P1", "P54")
m["new_pid_int"] = m["new_P"].map(lambda x: int(re.sub(r"[^0-9]", "", str(x))))
m = m.sort_values("new_pid_int").reset_index(drop=True)
assert len(m) == 54, f"expected 54 rows after merge, got {len(m)}"
assert m["new_pid_int"].tolist() == list(range(1, 55)), f"program ids not 1..54 contiguous: {m['new_pid_int'].tolist()}"

# ---------------------------------------------------------------- build supp dataframe for rendering
supp = pd.DataFrame({
    "program_id":        m["new_P"],
    "functional_name":   m["name_short"],
    "dominant_class":    m["class"],
    "dominant_subclass": m["dominant_subclass"],
    "top_GOBP_term":     m["gobp_term"],
    "GOBP_id":           m["gobp_id"],
    "NES":               m["brain_term_NES"].round(3),
    "FDR":               m["fdr"],
    "confidence":        m["confidence"],
    "top8_loading_genes":m["top8_loading_genes"],
})
print(f"brain-sig: {(supp['confidence']=='brain-sig').sum()}, brain-weak: {(supp['confidence']=='brain-weak').sum()}")

# ---------------------------------------------------------------- rendered table
CLASS_COL = {
    "exc":       "#3B6FB6",
    "inh":       "#C25E5E",
    "glia":      "#5E9E7A",
    "nonneuron": "#8A6FB0",
    "vascular":  "#C99A3B",
}
CLASS_LABEL = {
    "exc": "Excitatory", "inh": "Inhibitory", "glia": "Glia",
    "nonneuron": "Immune/Other", "vascular": "Vascular",
}

def fdr_fmt(v):
    try:
        v = float(v)
    except Exception:
        return str(v)
    if v < 1e-3:
        return f"{v:.1e}"
    return f"{v:.3f}"

def truncate(s, n):
    s = "" if pd.isna(s) else str(s)
    return s if len(s) <= n else s[: n - 1] + "…"

# Column layout
COLS = [
    ("P",          0.000, 0.034, "left"),
    ("Program",    0.034, 0.150, "left"),
    ("Class",      0.150, 0.215, "left"),
    ("Subclass",   0.215, 0.305, "left"),
    ("Top GO:BP term", 0.305, 0.560, "left"),
    ("NES",        0.560, 0.605, "right"),
    ("FDR",        0.605, 0.670, "right"),
    ("Top-8 loading genes", 0.670, 1.000, "left"),
]
TRUNC = {"Program": 26, "Subclass": 17, "Top GO:BP term": 40, "Top-8 loading genes": 52}

N = len(supp)     # 54
# Split 54 into two blocks: 27 per block (matching 60/2=30 original; use 27 for 54)
HALF = 27
blocks = [supp.iloc[:HALF].reset_index(drop=True),
          supp.iloc[HALF:].reset_index(drop=True)]

FIG_W, FIG_H = 17.0, 11.0   # same canvas as original
fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)

HEADER_FS = 6.4
CELL_FS   = 5.4
TITLE_FS  = 11.0
SUB_FS    = 7.0

fig.text(0.5, 0.975,
         "Supplementary Data 3  |  Cross-region cNMF gene programs: GO:BP functional annotation (P1–P54)",
         ha="center", va="top", fontsize=TITLE_FS, fontweight="bold")
fig.text(0.5, 0.952,
         "Functional name, dominant cell class/subclass, top GO:BP term, brain-term enrichment (NES, FDR) and top-8 loading genes "
         "for each of the 54 biologically-interpreted programs.  "
         "Bold = brain-significant (FDR<0.05);  grey = brain-weak (suggestive).  "
         "Six cohort-technical programs excluded (Methods).",
         ha="center", va="top", fontsize=SUB_FS - 0.5, color="#222222")

# class legend swatches
lx = 0.30
for cls in ["exc", "inh", "glia", "nonneuron", "vascular"]:
    fig.patches.append(Rectangle((lx, 0.928), 0.011, 0.011, transform=fig.transFigure,
                                 facecolor=CLASS_COL[cls], edgecolor="none", zorder=5))
    fig.text(lx + 0.014, 0.9335, CLASS_LABEL[cls], ha="left", va="center",
             fontsize=SUB_FS - 1.0, color="#222222")
    lx += 0.014 + 0.011 + 0.009 * len(CLASS_LABEL[cls]) + 0.010

PANEL_TOP  = 0.915
PANEL_BOT  = 0.018
PANEL_GAP  = 0.030
PANEL_W    = (1.0 - 2 * 0.012 - PANEL_GAP) / 2.0
PANEL_LEFTS= [0.012, 0.012 + PANEL_W + PANEL_GAP]

for bi, (df, ax_left) in enumerate(zip(blocks, PANEL_LEFTS)):
    ax = fig.add_axes([ax_left, PANEL_BOT, PANEL_W, PANEL_TOP - PANEL_BOT])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    nrows   = len(df)
    header_h = 1.0 / (nrows + 1.6)
    row_h    = (1.0 - header_h) / nrows

    def ytop_of_row(i):
        return 1.0 - header_h - i * row_h

    ax.add_patch(Rectangle((0, 1.0 - header_h), 1.0, header_h,
                           facecolor="#2A2A2A", edgecolor="none", zorder=1))
    for (label, x0, x1, align) in COLS:
        if align == "right":
            tx, ha = x1 - 0.004, "right"
        else:
            tx, ha = x0 + 0.004, "left"
        ax.text(tx, 1.0 - header_h / 2.0, label, ha=ha, va="center",
                fontsize=HEADER_FS, color="white", fontweight="bold", zorder=3)

    for i in range(nrows):
        r    = df.iloc[i]
        y1   = ytop_of_row(i)
        y0   = y1 - row_h
        cls  = r["dominant_class"]
        ccol = CLASS_COL.get(cls, "#999999")
        sig  = (r["confidence"] == "brain-sig")

        if i % 2 == 0:
            ax.add_patch(Rectangle((0, y0), 1.0, row_h, facecolor="#F2F2F2",
                                   edgecolor="none", zorder=0.5))
        ax.add_patch(Rectangle((0, y0), 0.006, row_h, facecolor=ccol,
                               edgecolor="none", zorder=2))

        txt_col = "#111111" if sig else "#8A8A8A"
        weight  = "bold" if sig else "normal"
        yc      = (y0 + y1) / 2.0

        cells = {
            "P":           r["program_id"],
            "Program":     truncate(r["functional_name"], TRUNC["Program"]),
            "Class":       CLASS_LABEL.get(cls, str(cls)),
            "Subclass":    truncate(r["dominant_subclass"], TRUNC["Subclass"]),
            "Top GO:BP term": truncate(r["top_GOBP_term"], TRUNC["Top GO:BP term"]),
            "NES":         f"{float(r['NES']):.2f}",
            "FDR":         fdr_fmt(r["FDR"]),
            "Top-8 loading genes": truncate(r["top8_loading_genes"], TRUNC["Top-8 loading genes"]),
        }
        for (label, x0c, x1c, align) in COLS:
            val = cells[label]
            if label == "Class":
                cc = ccol if sig else "#9AA0A6"
                ax.text(x0c + 0.010, yc, val, ha="left", va="center",
                        fontsize=CELL_FS, color=cc, fontweight=weight,
                        fontstyle="italic", zorder=3)
                continue
            if align == "right":
                tx, ha = x1c - 0.004, "right"
            else:
                tx, ha = (x0c + 0.010 if label == "P" else x0c + 0.004), "left"
            fst = "italic" if label == "Top-8 loading genes" else "normal"
            ax.text(tx, yc, val, ha=ha, va="center", fontsize=CELL_FS,
                    color=txt_col, fontweight=weight, fontstyle=fst, zorder=3)

    ax.add_patch(Rectangle((0, 0), 1.0, 1.0, fill=False, edgecolor="#BBBBBB",
                           linewidth=0.4, zorder=4))

pdf_path = os.path.join(FIGDIR, "ed_fig4_program_table_54.pdf")
png_path = os.path.join(FIGDIR, "ed_fig4_program_table_54.png")
fig.savefig(pdf_path)
fig.savefig(png_path, dpi=300)
print(f"wrote {pdf_path}")
print(f"wrote {png_path}")
print(f"page_inches={FIG_W}x{FIG_H}")
print(f"min_font_pt={min(CELL_FS, HEADER_FS, SUB_FS - 1.0)}")
