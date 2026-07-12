#!/usr/bin/env python3
"""Compose the Fig. 7 vector panels from the current source tables."""
import argparse
import os
import sys
from pathlib import Path
import svgutils.transform as sg
from svgutils.compose import Unit
from lxml import etree

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPOSITORY_ROOT))
from workflow.root_contract import add_canonical_root_argument, resolve_canonical_root

parser = argparse.ArgumentParser(description=__doc__)
add_canonical_root_argument(parser)
parser.epilog = "Uses --canonical-root or CORTEX_PROGRAM_CANONICAL_ROOT."
args = parser.parse_args()
canonical_root = resolve_canonical_root(args.canonical_root)
FIG = str(canonical_root / "results/xspecies_humanmap_v1/spatial_xspecies/figures/Fig_spatial_univ")
SUPPORTING_PANELS = f"{FIG}/panels_supporting"
PANELS = f"{FIG}/panels"
A_SVG = f"{PANELS}/panel_a_complex_heatmap_crop.svg"
OUT = FIG

PT_PER_MM = 72.0 / 25.4
GUT    = 2.0 * PT_PER_MM
GUT_AR = 2.2 * PT_PER_MM
GUT_V  = 2.0 * PT_PER_MM
L_MARGIN = 5.0 * PT_PER_MM
R_MARGIN = 3.0 * PT_PER_MM
T_MARGIN = 3.0 * PT_PER_MM
RIGHT_W_TARGET = 40.0 * PT_PER_MM
LET_DX, LET_DY = -1.0, 9.0
LET_FS = 13

def load_path(path):
    fig = sg.fromfile(path)
    w = Unit(fig.width).to("pt").value if isinstance(fig.width, str) else float(fig.width)
    h = Unit(fig.height).to("pt").value if isinstance(fig.height, str) else float(fig.height)
    return fig, w, h

SRC = {
    "a": A_SVG,
    "b": f"{PANELS}/panel_b.svg",
    "c": f"{SUPPORTING_PANELS}/panel_e.svg",
    "d": f"{SUPPORTING_PANELS}/panel_f.svg",
    "f": f"{PANELS}/panel_aggregate.svg",   # aggregate cross-species fraction + Wilson CI
    "e": f"{PANELS}/panel_e_spatial.svg",
}
P = {k: load_path(p) for k, p in SRC.items()}

def placed(name, x, y, s):
    fig, w, h = P[name]
    root = fig.getroot(); root.moveto(x, y, scale_x=s, scale_y=s)
    return root, w*s, h*s

def letter(ch, x, y):
    return sg.TextElement(x, y, ch, size=LET_FS, weight="bold", font="Liberation Sans")

PAGE_W = 180.0 * PT_PER_MM
inner_w = PAGE_W - L_MARGIN - R_MARGIN
x0 = L_MARGIN; y0 = T_MARGIN

elems = []; letters = []

aw0, ah0 = P["a"][1], P["a"][2]
bw0, bh0 = P["b"][1], P["b"][2]
cw0, ch0 = P["c"][1], P["c"][2]
dw0, dh0 = P["d"][1], P["d"][2]
fw0, fh0 = P["f"][1], P["f"][2]
ew0, eh0 = P["e"][1], P["e"][2]

# ---------- TOP-LEFT: a, width = inner - right_col - gutter ----------
right_w = RIGHT_W_TARGET
aw = inner_w - right_w - GUT_AR
s_a = aw / aw0
ah = ah0 * s_a
ra, _, _ = placed("a", x0, y0, s_a)
elems.append(ra); letters.append(("a", x0+LET_DX, y0+LET_DY))

# ---------- RIGHT COLUMN: b / c / d stacked, LEFT-ALIGNED, fills a's height ----------
right_x = x0 + aw + GUT_AR
def hw(w0, h0, w): return h0 * (w / w0)
bh = hw(bw0, bh0, right_w)
ch = hw(cw0, ch0, right_w)
dh = hw(dw0, dh0, right_w)
fh = hw(fw0, fh0, right_w)
nat_stack = bh + ch + dh + fh + 3*GUT
stack_scale = ah / nat_stack
bh *= stack_scale; ch *= stack_scale; dh *= stack_scale; fh *= stack_scale
def place_right(name, y, target_h):
    fig, w0, h0 = P[name]
    s = target_h / h0
    w = w0 * s
    if w > right_w:
        s = right_w / w0; w = right_w; target_h = h0 * s
    x = right_x + (right_w - w) / 2.0   # center within the narrowed slot
    r, _, _ = placed(name, x, y, s)
    return r, x, w, h0 * s

ry = y0
rb, bx, bw_, bh_ = place_right("b", ry, bh)
elems.append(rb); letters.append(("b", right_x+LET_DX, ry+LET_DY))
ry += bh_ + GUT
rc, cx, cw_, ch_ = place_right("c", ry, ch)
elems.append(rc); letters.append(("c", right_x+LET_DX, ry+LET_DY))
ry += ch_ + GUT
rd, dx, dw_, dh_ = place_right("d", ry, dh)
elems.append(rd); letters.append(("d", right_x+LET_DX, ry+LET_DY))
ry += dh_ + GUT
rf, fx, fw_, fh_ = place_right("f", ry, fh)
elems.append(rf); letters.append(("f", right_x+LET_DX, ry+LET_DY))
ry += fh_

# bottom of the top region = the LOWER of (a bottom, right stack bottom)
top_bottom = max(y0 + ah, ry)
ycur = top_bottom + GUT_V               # tight vertical gutter

# ---------- BOTTOM BAND: e full width (spatial fields 1x6) ----------
s_e = inner_w / ew0
ew = inner_w * 1.0; eh = eh0 * s_e
re_, _, _ = placed("e", x0, ycur, s_e)
elems.append(re_); letters.append(("e", x0+LET_DX, ycur+LET_DY))
ycur += eh

content_bottom = ycur

# ---------- finalize: unified macaque caveat footer (tightened) ----------
FOOT_GAP = 1.2 * PT_PER_MM
FOOT_FS = 5.8
B_MARGIN = 2.0 * PT_PER_MM
PAGE_H = content_bottom + FOOT_GAP + FOOT_FS + B_MARGIN

fig = sg.SVGFigure()
fig.set_size((f"{PAGE_W}pt", f"{PAGE_H}pt"))
fig.root.set("viewBox", f"0 0 {PAGE_W} {PAGE_H}")
bg = etree.Element("rect", {"x": "0", "y": "0", "width": str(PAGE_W),
                            "height": str(PAGE_H), "fill": "white"})
fig.root.append(bg)
for e in elems: fig.append(e)
for ch, lx, ly in letters: fig.append(letter(ch, lx, ly))

foot_y = content_bottom + FOOT_GAP + FOOT_FS
foot = sg.TextElement(
    PAGE_W/2.0, foot_y,
    "Macaque: cellbin-sparse data with a global laminar offset; macaque panels are qualitative and technically limited.  "
    "|  Cross-species: single pair p-values not reported; only aggregate same-sign fraction with 95% Wilson CI reported (panel f).",
    size=FOOT_FS, font="Liberation Sans", color="#777")
foot.root.set("text-anchor", "middle")
foot.root.set("font-style", "italic")
fig.append(foot)

svg_out = f"{OUT}/Figure_7.svg"
fig.save(svg_out)
print("wrote", svg_out, f"page {PAGE_W/PT_PER_MM:.1f} x {PAGE_H/PT_PER_MM:.1f} mm")
print(f"a: {aw/PT_PER_MM:.1f}x{ah/PT_PER_MM:.1f}mm  = {aw/PAGE_W*100:.0f}% W x {ah/PAGE_H*100:.0f}% H")
print(f"right col w={right_w/PT_PER_MM:.1f}mm start_x={right_x/PT_PER_MM:.1f}mm  b={bh_/PT_PER_MM:.0f} c={ch_/PT_PER_MM:.0f} d={dh_/PT_PER_MM:.0f}mm")
print(f"b placed width={bw_/PT_PER_MM:.1f}mm  a right edge={(x0+aw)/PT_PER_MM:.1f}mm")
print(f"e band: {ew/PT_PER_MM:.0f}x{eh/PT_PER_MM:.0f}mm  page_h={PAGE_H/PT_PER_MM:.1f}mm")
