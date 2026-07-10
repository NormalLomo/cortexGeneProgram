#!/usr/bin/env python
"""
Extended Data Fig 4 / Supplementary Table:
60-program GO:BP annotation reference table for the cross-region cNMF programs.

Outputs:
  - results/crossregion_v1/supp_table_program_annotation.tsv  (machine-readable, all cols)
  - figures/extended/ed_fig4_program_table.{pdf,png}          (rendered figure-table, vector)

Sources (joined on integer program id 1..60):
  - results/crossregion_v1/program_names.tsv           -> name_short, confidence(brain-sig/weak), brain_term_NES, fdr
  - results/crossregion_v1/program_annotation.tsv      -> dominant_subclass, class
  - results/crossregion_v1/program_annotation_gobp.tsv -> top GO:BP term + term_id, top10 loading genes
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.font_manager as fm

# deterministic, embeddable fonts (Type42/TrueType so pdffonts shows embedded subsets)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "DejaVu Sans"

ROOT = "CORTEX_PROGRAM_ROOT"
CR = os.path.join(ROOT, "results", "crossregion_v1")
FIGDIR = os.path.join(ROOT, "figures", "extended")
os.makedirs(FIGDIR, exist_ok=True)

# ---------------------------------------------------------------- load + merge
names = pd.read_csv(os.path.join(CR, "program_names.tsv"), sep="\t")
ann = pd.read_csv(os.path.join(CR, "program_annotation.tsv"), sep="\t")
gobp = pd.read_csv(os.path.join(CR, "program_annotation_gobp.tsv"), sep="\t")


def pid(x):
    """normalise any program id form ('1','P01','P1') to int."""
    return int(re.sub(r"[^0-9]", "", str(x)))


names["pid"] = names["program"].map(pid)
ann["pid"] = ann["program"].map(pid)
gobp["pid"] = gobp["program"].map(pid)

# split GO:BP term text and GO id from "Name (GO:nnnnnnn)"
def split_term(s):
    s = str(s)
    m = re.search(r"\(GO:(\d+)\)", s)
    go = ("GO:" + m.group(1)) if m else ""
    txt = re.sub(r"\s*\(GO:\d+\)\s*$", "", s).strip()
    return txt, go


_split = [split_term(s) for s in gobp["top_BP_term"]]
gobp["gobp_term"] = [t for (t, g) in _split]
gobp["gobp_id"] = [g for (t, g) in _split]


def top8_genes(s):
    parts = re.split(r"[;,]", str(s))
    parts = [p.strip() for p in parts if p.strip()]
    return ", ".join(parts[:8])


gobp["top8_loading_genes"] = gobp["top10_loading_genes"].map(top8_genes)

# brain_term_NES / fdr / confidence come from program_names.tsv because that file's
# NES & fdr are what its brain-sig / brain-weak call is derived from (gobp.tsv carries
# a different NES/fdr for some programs; keep the confidence-consistent pair).
m = names[["pid", "name_short", "confidence", "brain_term_NES", "fdr"]].merge(
    ann[["pid", "dominant_subclass", "class"]], on="pid", how="left"
).merge(
    gobp[["pid", "gobp_term", "gobp_id", "top8_loading_genes"]], on="pid", how="left"
)
m = m.sort_values("pid").reset_index(drop=True)
assert len(m) == 60, f"expected 60 programs, got {len(m)}"
assert m["pid"].tolist() == list(range(1, 61)), "program ids not 1..60 contiguous"

# tidy supplementary table (machine-readable, all cols)
supp = pd.DataFrame({
    "program_id": "P" + m["pid"].astype(str).str.zfill(2),
    "functional_name": m["name_short"],
    "dominant_class": m["class"],
    "dominant_subclass": m["dominant_subclass"],
    "top_GOBP_term": m["gobp_term"],
    "GOBP_id": m["gobp_id"],
    "NES": m["brain_term_NES"].round(3),
    "FDR": m["fdr"],
    "confidence": m["confidence"],
    "top8_loading_genes": m["top8_loading_genes"],
})
supp_path = os.path.join(CR, "supp_table_program_annotation.tsv")
supp.to_csv(supp_path, sep="\t", index=False)
print("wrote", supp_path, supp.shape)
print("brain-sig:", (supp["confidence"] == "brain-sig").sum(),
      "brain-weak:", (supp["confidence"] == "brain-weak").sum())

# ------------------------------------------------------------- rendered table
# class -> colour (left accent stripe + tag); journal-ish muted palette
CLASS_COL = {
    "exc":       "#3B6FB6",  # excitatory neuron - blue
    "inh":       "#C25E5E",  # inhibitory neuron - red
    "glia":      "#5E9E7A",  # glia (astro/oligo/opc) - green
    "nonneuron": "#8A6FB0",  # microglia / immune - purple
    "vascular":  "#C99A3B",  # vascular - amber
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


# column layout (axis fraction within one 30-row sub-table). widths sum ~1.
# id | name | class | subcl | GO:BP term | NES | FDR | top-8 genes
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

N = len(supp)
HALF = 30
blocks = [supp.iloc[:HALF].reset_index(drop=True),
          supp.iloc[HALF:].reset_index(drop=True)]

# figure geometry: 2 side-by-side blocks of 30 rows. A3-landscape-ish canvas.
FIG_W, FIG_H = 17.0, 11.0          # inches
fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)

HEADER_FS = 6.4
CELL_FS = 5.4
TITLE_FS = 11.0
SUB_FS = 7.0

# overall title + legend band at top
fig.text(0.5, 0.975,
         "Extended Data Fig. 4  |  Cross-region cNMF gene programs: GO:BP functional annotation (P1–P60)",
         ha="center", va="top", fontsize=TITLE_FS, fontweight="bold")
fig.text(0.5, 0.952,
         "Functional name, dominant cell class/subclass, top GO:BP term, brain-term enrichment (NES, FDR) and top-8 loading genes for each of the 60 programs.  "
         "Bold = brain-significant (FDR<0.05);  grey = brain-weak (suggestive).",
         ha="center", va="top", fontsize=SUB_FS - 0.5, color="#222222")

# class legend swatches
lx = 0.30
for cls in ["exc", "inh", "glia", "nonneuron", "vascular"]:
    fig.patches.append(Rectangle((lx, 0.928), 0.011, 0.011, transform=fig.transFigure,
                                 facecolor=CLASS_COL[cls], edgecolor="none", zorder=5))
    fig.text(lx + 0.014, 0.9335, CLASS_LABEL[cls], ha="left", va="center",
             fontsize=SUB_FS - 1.0, color="#222222")
    lx += 0.014 + 0.011 + 0.009 * len(CLASS_LABEL[cls]) + 0.010

# two axes panels (left block, right block)
PANEL_TOP = 0.915
PANEL_BOT = 0.018
PANEL_GAP = 0.030
PANEL_W = (1.0 - 2 * 0.012 - PANEL_GAP) / 2.0
PANEL_LEFTS = [0.012, 0.012 + PANEL_W + PANEL_GAP]

for bi, (df, ax_left) in enumerate(zip(blocks, PANEL_LEFTS)):
    ax = fig.add_axes([ax_left, PANEL_BOT, PANEL_W, PANEL_TOP - PANEL_BOT])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    nrows = len(df)
    header_h = 1.0 / (nrows + 1.6)       # a bit taller header
    row_h = (1.0 - header_h) / nrows

    def ytop_of_row(i):
        return 1.0 - header_h - i * row_h

    # header background
    ax.add_patch(Rectangle((0, 1.0 - header_h), 1.0, header_h,
                           facecolor="#2A2A2A", edgecolor="none", zorder=1))
    for (label, x0, x1, align) in COLS:
        if align == "right":
            tx, ha = x1 - 0.004, "right"
        else:
            tx, ha = x0 + 0.004, "left"
        ax.text(tx, 1.0 - header_h / 2.0, label, ha=ha, va="center",
                fontsize=HEADER_FS, color="white", fontweight="bold", zorder=3)

    # rows
    for i in range(nrows):
        r = df.iloc[i]
        y1 = ytop_of_row(i)
        y0 = y1 - row_h
        cls = r["dominant_class"]
        ccol = CLASS_COL.get(cls, "#999999")
        sig = (r["confidence"] == "brain-sig")

        # zebra background
        if i % 2 == 0:
            ax.add_patch(Rectangle((0, y0), 1.0, row_h, facecolor="#F2F2F2",
                                   edgecolor="none", zorder=0.5))
        # left class accent stripe
        ax.add_patch(Rectangle((0, y0), 0.006, row_h, facecolor=ccol,
                               edgecolor="none", zorder=2))

        txt_col = "#111111" if sig else "#8A8A8A"
        weight = "bold" if sig else "normal"
        yc = (y0 + y1) / 2.0

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
        for (label, x0, x1, align) in COLS:
            val = cells[label]
            if label == "Class":
                # colour the class word itself
                cc = ccol if sig else "#9AA0A6"
                ax.text(x0 + 0.010, yc, val, ha="left", va="center",
                        fontsize=CELL_FS, color=cc, fontweight=weight,
                        fontstyle="italic", zorder=3)
                continue
            if align == "right":
                tx, ha = x1 - 0.004, "right"
            else:
                tx, ha = (x0 + 0.010 if label == "P" else x0 + 0.004), "left"
            fst = "italic" if label == "Top-8 loading genes" else "normal"
            ax.text(tx, yc, val, ha=ha, va="center", fontsize=CELL_FS,
                    color=txt_col, fontweight=weight, fontstyle=fst, zorder=3)

    # thin header underline + outer frame
    ax.add_patch(Rectangle((0, 0), 1.0, 1.0, fill=False, edgecolor="#BBBBBB",
                           linewidth=0.4, zorder=4))

pdf_path = os.path.join(FIGDIR, "ed_fig4_program_table.pdf")
png_path = os.path.join(FIGDIR, "ed_fig4_program_table.png")
fig.savefig(pdf_path)
fig.savefig(png_path, dpi=300)
print("wrote", pdf_path)
print("wrote", png_path)
print(f"page_inches={FIG_W}x{FIG_H}  (~{FIG_W*25.4:.0f}x{FIG_H*25.4:.0f} mm)")
print("min_font_pt=", min(CELL_FS, HEADER_FS, SUB_FS - 1.0))
