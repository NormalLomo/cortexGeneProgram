#!/usr/bin/env python3
"""Compose the Fig. 6 vector panels using the documented page geometry."""
import os
import re
from lxml import etree
from svgutils.transform import SVGFigure, TextElement, fromfile

BASE = "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1/figures/Fig10_v2"
PANELS = os.environ.get("FIG10_PANELS_DIR", f"{BASE}/svg_panels")
OUT = os.environ.get("FIG10_OUT_SVG", f"{BASE}/Fig10_v2.svg")

PAGE_W, PAGE_H = 532.9134, 724.3710

BOX = {
    "a": (11.34, 11.34, 245.32, 193.13),
    "b": (262.94, 11.34, 251.07, 193.13),
    "c": (11.25, 212.97, 170.98, 155.56),
    "d": (187.57, 175.33, 337.74, 196.10),
    "e": (1.57, 377.04, 521.97, 74.07),
    "f": (93.85, 430.80, 344.79, 178.55),
    "g": (10.43, 616.54, 511.14, 115.76),
}

LETTERS = {
    "a": (11.34, 9.84),
    "b": (265.16, 9.84),
    "c": (11.34, 211.47),
    "d": (194.99, 226.22),
    "e": (15.59, 390.29),
    "f": (113.39, 429.30),
    "g": (15.59, 630.43),
}

FILL_WIDTH = {"e", "g"}
BOTTOM_ANCHOR = {"g": 712.0}


def parse_len(value: str) -> float:
    match = re.match(r"^\s*([-0-9.eE]+)\s*([a-z%]*)\s*$", str(value))
    number = float(match.group(1))
    unit = (match.group(2) or "pt").lower()
    if unit in ("pt", ""):
        return number
    if unit == "px":
        return number * 72.0 / 96.0
    if unit == "in":
        return number * 72.0
    if unit == "mm":
        return number * 72.0 / 25.4
    if unit == "cm":
        return number * 72.0 / 2.54
    return number


fig = SVGFigure()
fig.set_size((f"{PAGE_W}pt", f"{PAGE_H}pt"))
fig.root.set("viewBox", f"0 0 {PAGE_W} {PAGE_H}")

report = []
for name in ["a", "b", "c", "d", "e", "f", "g"]:
    bx, by, bw, bh = BOX[name]
    svg = fromfile(f"{PANELS}/{name}.svg")
    width = parse_len(svg.get_size()[0])
    height = parse_len(svg.get_size()[1])
    root = svg.getroot()

    if name in FILL_WIDTH:
        scale = bw / width
        drawn_w, drawn_h = width * scale, height * scale
        ox = bx
        oy = BOTTOM_ANCHOR[name] - drawn_h if name in BOTTOM_ANCHOR else by + (bh - drawn_h) / 2.0
    else:
        scale = min(bw / width, bh / height)
        drawn_w, drawn_h = width * scale, height * scale
        ox = bx + (bw - drawn_w) / 2.0
        oy = by + (bh - drawn_h) / 2.0

    root.moveto(ox, oy, scale_x=scale, scale_y=scale)
    fig.append(root)
    report.append(
        f"{name}: raw {width:.1f}x{height:.1f}pt scale={scale:.4f} "
        f"-> {drawn_w:.1f}x{drawn_h:.1f} at ({ox:.1f},{oy:.1f})"
    )

for name, (lx, ly) in LETTERS.items():
    fig.append(TextElement(lx, ly, name, size=11, weight="bold", font="Helvetica", color="black"))

out_dir = os.path.dirname(OUT)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)
fig.save(OUT)
print("WROTE", OUT)
for item in report:
    print(" ", item)
