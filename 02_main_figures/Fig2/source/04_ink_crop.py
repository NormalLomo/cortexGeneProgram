#!/usr/bin/env python3
# =====================================================================
# Fig.2 panel INK-CROP (copied from scripts/fig1/fig1_ink_crop.py; the
# fig1 version hardcodes the fig1 dir + "abcdefghij" so this is a fig2
# sibling pointed at figures/fig3/svg_panels and looping "abcdefgh").
# For each panel SVG (figures/fig3/svg_panels/fig3_{a..h}.svg):
#   1. render SVG -> PDF (rsvg; 1 SVG user-unit == 1 PDF pt, 1:1)
#   2. gs -sDEVICE=bbox -> tight HiRes ink BoundingBox (PDF pts, y up)
#   3. rewrite SVG viewBox + width/height to that tight bbox (y-flip into
#      SVG top-left space). Pure vector: only the header changes.
# Writes cropped SVGs IN PLACE and a JSON sidecar of aspects.
# =====================================================================
import os, re, json, subprocess, tempfile
from lxml import etree

PROJ = "CORTEX_PROGRAM_ROOT"
SVGD = os.path.join(PROJ, "figures/fig3/svg_panels")
RSCRIPT = "/usr/local/bin/Rscript"
PAD = 0.5  # pt safety pad so anti-aliased strokes/text are not clipped

SVG_NS = "http://www.w3.org/2000/svg"


def native_size_pt(root):
    vb = root.get("viewBox")
    if vb:
        parts = [float(x) for x in re.split(r"[\s,]+", vb.strip())]
        return parts[2], parts[3]
    w = float(re.sub(r"[^\d.]", "", root.get("width", "0")))
    h = float(re.sub(r"[^\d.]", "", root.get("height", "0")))
    return w, h


def ink_bbox_pt(svg_path):
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tf:
        pdf = tf.name
    try:
        r = subprocess.run(
            [RSCRIPT, "-e", f'rsvg::rsvg_pdf("{svg_path}", "{pdf}")'],
            capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"rsvg_pdf failed: {r.stderr.strip()}")
        g = subprocess.run(
            ["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=bbox", pdf],
            capture_output=True, text=True)
        out = g.stderr + g.stdout
        m = re.search(
            r"%%HiResBoundingBox:\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+"
            r"([\d.eE+-]+)\s+([\d.eE+-]+)", out)
        if not m:
            raise RuntimeError(f"no HiResBoundingBox in gs output:\n{out}")
        return tuple(float(m.group(i)) for i in range(1, 5))
    finally:
        if os.path.exists(pdf):
            os.remove(pdf)


def crop_one(pid):
    svg_path = os.path.join(SVGD, f"fig3_{pid}.svg")
    tree = etree.parse(svg_path)
    root = tree.getroot()
    cw, ch = native_size_pt(root)
    x0, y0, x1, y1 = ink_bbox_pt(svg_path)

    x0 = max(0.0, x0 - PAD); y0 = max(0.0, y0 - PAD)
    x1 = min(cw, x1 + PAD); y1 = min(ch, y1 + PAD)
    w = x1 - x0
    h = y1 - y0
    svg_x = x0
    svg_y = ch - y1

    root.set("viewBox", f"{svg_x:.3f} {svg_y:.3f} {w:.3f} {h:.3f}")
    root.set("width", f"{w:.3f}pt")
    root.set("height", f"{h:.3f}pt")
    tree.write(svg_path, xml_declaration=True, encoding="UTF-8")

    return {
        "panel": pid,
        "canvas_pt": [round(cw, 2), round(ch, 2)],
        "ink_bbox_pt": [round(x0, 2), round(y0, 2), round(x1, 2), round(y1, 2)],
        "cropped_pt": [round(w, 2), round(h, 2)],
        "aspect_w_h": round(w / h, 4),
        "trim_frac": round(1.0 - (w * h) / (cw * ch), 4),
    }


def main():
    results = {}
    for pid in "abcdefghi":
        r = crop_one(pid)
        results[pid] = r
        print(f"panel {pid}: canvas {r['canvas_pt']} -> cropped {r['cropped_pt']} "
              f"pt  AR {r['aspect_w_h']}  (trimmed {r['trim_frac']*100:.1f}% area)")
    side = os.path.join(SVGD, "_cropped_aspects.json")
    with open(side, "w") as fh:
        json.dump(results, fh, indent=2)
    print("WROTE", side)


if __name__ == "__main__":
    main()
