#!/usr/bin/env python3
"""Ink-crop the complex heatmap SVG to a tight bounding box.
Render to hi-res PNG via rsvg-convert, find non-white bbox (PIL/numpy),
map px -> svg pt, rewrite viewBox + width/height. Pure geometry crop -> NO visual change.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path
import numpy as np
from PIL import Image
from lxml import etree

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPOSITORY_ROOT))
from workflow.root_contract import add_canonical_root_argument, resolve_canonical_root

parser = argparse.ArgumentParser(description=__doc__)
add_canonical_root_argument(parser)
args = parser.parse_args()
canonical_root = resolve_canonical_root(args.canonical_root)
FIG = canonical_root / "results/xspecies_humanmap_v1/spatial_xspecies/figures/Fig_spatial_univ"
SRC = str(FIG / "panels/panel_a_complex_heatmap.svg")
OUT = str(FIG / "panels/panel_a_complex_heatmap_crop.svg")
RSVG = "rsvg-convert"
DPI = 300
PAD_MM = 0.3  # tiny breathing pad so border line not shaved
PT_PER_MM = 72.0 / 25.4

tree = etree.parse(SRC)
root = tree.getroot()
vb = root.get("viewBox")
if vb:
    _, _, W_pt, H_pt = [float(x) for x in vb.split()]
else:
    def topt(s):
        m = re.match(r"([0-9.]+)", s); return float(m.group(1))
    W_pt = topt(root.get("width")); H_pt = topt(root.get("height"))

png = "/tmp/fig7_panel_a_crop_probe.png"
subprocess.run([RSVG, "-d", str(DPI), "-p", str(DPI), SRC, "-o", png], check=True)
im = Image.open(png).convert("RGBA")
arr = np.asarray(im)
H_px, W_px = arr.shape[:2]
alpha = arr[..., 3]
rgb = arr[..., :3]
nonwhite = (alpha > 8) & ((rgb < 250).any(axis=2))
ys, xs = np.where(nonwhite)
if len(xs) == 0:
    sys.exit("no ink found")
x0p, x1p = xs.min(), xs.max() + 1
y0p, y1p = ys.min(), ys.max() + 1

sx = W_pt / W_px
sy = H_pt / H_px
pad = PAD_MM * PT_PER_MM
x0 = max(0.0, x0p * sx - pad)
y0 = max(0.0, y0p * sy - pad)
x1 = min(W_pt, x1p * sx + pad)
y1 = min(H_pt, y1p * sy + pad)
new_w = x1 - x0
new_h = y1 - y0

root.set("viewBox", f"{x0:.3f} {y0:.3f} {new_w:.3f} {new_h:.3f}")
root.set("width", f"{new_w:.3f}pt")
root.set("height", f"{new_h:.3f}pt")
tree.write(OUT, xml_declaration=True, encoding="UTF-8")
print(f"orig {W_pt:.1f}x{H_pt:.1f}pt ({W_pt/PT_PER_MM:.1f}x{H_pt/PT_PER_MM:.1f}mm)")
print(f"crop {new_w:.1f}x{new_h:.1f}pt ({new_w/PT_PER_MM:.1f}x{new_h/PT_PER_MM:.1f}mm) AR={new_w/new_h:.3f}")
print(f"trimmed L={x0/PT_PER_MM:.1f} T={y0/PT_PER_MM:.1f} "
      f"R={(W_pt-x1)/PT_PER_MM:.1f} B={(H_pt-y1)/PT_PER_MM:.1f} mm")
print("wrote", OUT)
