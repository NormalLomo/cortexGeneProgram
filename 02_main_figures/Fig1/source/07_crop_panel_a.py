#!/usr/bin/env python3
# Ink-crop the fig1 panel-A graphical-abstract SVG to its tight content bbox
# (rsvg->PDF->gs bbox route, same logic as figcrossspecies/ink_crop.py).
# NON-destructive: reads fig1_a_GA.svg, writes fig1_a_GA_crop.svg (keeps original).
import os, re, json, subprocess, tempfile
from lxml import etree

SVGD = "CORTEX_PROGRAM_ROOT/figures/fig1/svg_panels"
SRC = os.path.join(SVGD, "fig1_a_GA.svg")
DST = os.path.join(SVGD, "fig1_a_GA_crop.svg")
RSCRIPT = "/usr/local/bin/Rscript"
PAD = 0.5


def native_size_pt(root):
    vb = root.get("viewBox")
    if vb:
        parts = [float(x) for x in re.split(r"[\s,]+", vb.strip())]
        return parts[0], parts[1], parts[2], parts[3]
    w = float(re.sub(r"[^\d.]", "", root.get("width", "0")))
    h = float(re.sub(r"[^\d.]", "", root.get("height", "0")))
    return 0.0, 0.0, w, h


def ink_bbox_pt(svg_path):
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tf:
        pdf = tf.name
    try:
        r = subprocess.run([RSCRIPT, "-e", f'rsvg::rsvg_pdf("{svg_path}", "{pdf}")'],
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"rsvg_pdf failed: {r.stderr.strip()}")
        g = subprocess.run(["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=bbox", pdf],
                           capture_output=True, text=True)
        out = g.stderr + g.stdout
        m = re.search(r"%%HiResBoundingBox:\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+"
                      r"([\d.eE+-]+)\s+([\d.eE+-]+)", out)
        if not m:
            raise RuntimeError(f"no HiResBoundingBox:\n{out}")
        return tuple(float(m.group(i)) for i in range(1, 5))
    finally:
        if os.path.exists(pdf):
            os.remove(pdf)


def main():
    tree = etree.parse(SRC)
    root = tree.getroot()
    vx, vy, cw, ch = native_size_pt(root)
    # rsvg renders the SVG to a PDF whose page == the SVG viewBox size (cw x ch),
    # so gs bbox is in that same pt space with PDF origin bottom-left.
    x0, y0, x1, y1 = ink_bbox_pt(SRC)
    x0 = max(0.0, x0 - PAD); y0 = max(0.0, y0 - PAD)
    x1 = min(cw, x1 + PAD); y1 = min(ch, y1 + PAD)
    w = x1 - x0; h = y1 - y0
    # SVG viewBox uses top-left origin; PDF bbox uses bottom-left -> min-y = ch - y1.
    # add back the original viewBox offset (vx,vy) which is 0 here but keep general.
    new_minx = vx + x0
    new_miny = vy + (ch - y1)
    root.set("viewBox", f"{new_minx:.3f} {new_miny:.3f} {w:.3f} {h:.3f}")
    root.set("width", f"{w:.3f}pt"); root.set("height", f"{h:.3f}pt")
    tree.write(DST, xml_declaration=True, encoding="UTF-8")
    res = {"src_canvas_pt": [round(cw, 2), round(ch, 2)],
           "ink_bbox_pt": [round(x0, 2), round(y0, 2), round(x1, 2), round(y1, 2)],
           "cropped_pt": [round(w, 2), round(h, 2)],
           "aspect_w_h": round(w / h, 4),
           "trim_frac": round(1.0 - (w * h) / (cw * ch), 4)}
    print(json.dumps(res, indent=2))
    with open(os.path.join(SVGD, "_cropped_aspect_a.json"), "w") as fh:
        json.dump(res, fh, indent=2)


if __name__ == "__main__":
    main()
