#!/usr/bin/env python3
"""ink_crop.py — crop every SVG in a dir to its tight drawn-ink bounding box.

WHY: even a "square" panel rendered zero-margin still carries edge whitespace
(axis tick-label gutters, legend padding, expand()). Composition is only tight
if each SVG's box == its drawn content. This rewrites each SVG's viewBox +
width/height to the ink bbox so downstream svgutils places real content, not air.

HOW (pure-vector, header-only edit):
  SVG --(rsvg::rsvg_pdf, 1 user-unit == 1 pt)--> PDF
  PDF --(Ghostscript -sDEVICE=bbox)--> HiRes ink bbox in PDF pt
  rewrite the SVG <svg> viewBox + width/height to that bbox
    (+pad, y-flipped from PDF bottom-left origin to SVG top-left origin)

Only Ghostscript + R(rsvg) are required (both usually present on the render host).
inkscape --export-area-drawing or pdfcrop also work but are often NOT installed.

USAGE
  python ink_crop.py PANELDIR [--out OUTDIR] [--pad PT] [--rscript RSCRIPT] [--gs GS]
    PANELDIR   dir of *.svg panels (zero-margin renders from stage 1)
    --out      where to write cropped SVGs (default: PANELDIR/cropped)
    --pad      extra padding in pt around the ink bbox (default 0.5)
    --rscript  path to Rscript (default: Rscript on PATH)
    --gs       path to ghostscript (default: gs on PATH)

OUTPUT: OUTDIR/<name>.svg for each input, with corrected viewBox; prints a TSV
(name, orig_w, orig_h, crop_w, crop_h, pct_trimmed) so you can see how much air
each panel carried (12-26% trim was typical).

NOTE: ink-crop removes OUTER whitespace only. INTERIOR dead space (sparse
schematic with a far-below subtitle) survives — fix that at panel design time.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

SVG_DIMS_RE = {
    "viewBox": re.compile(r'''viewBox\s*=\s*["']([\-\d.eE ]+)["']'''),
    "width": re.compile(r'''\bwidth\s*=\s*["']([\d.eE]+)(?:pt|px)?["']'''),
    "height": re.compile(r'''\bheight\s*=\s*["']([\d.eE]+)(?:pt|px)?["']'''),
}


def svg_to_pdf(svg_path, pdf_path, rscript):
    """Render SVG->PDF at 1 user-unit == 1pt via R rsvg::rsvg_pdf."""
    r_code = (
        f'rsvg::rsvg_pdf("{svg_path}", "{pdf_path}")'
    )
    res = subprocess.run(
        [rscript, "-e", r_code],
        capture_output=True, text=True,
    )
    if res.returncode != 0 or not os.path.exists(pdf_path):
        raise RuntimeError(
            f"rsvg_pdf failed for {svg_path}:\n{res.stderr or res.stdout}"
        )


def pdf_bbox(pdf_path, gs):
    """Return (x0, y0, x1, y1) HiRes ink bbox in PDF points via gs -sDEVICE=bbox."""
    res = subprocess.run(
        [gs, "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=bbox", pdf_path],
        capture_output=True, text=True,
    )
    out = res.stderr + res.stdout
    m = re.search(r"%%HiResBoundingBox:\s*([\d.\-eE]+) ([\d.\-eE]+) ([\d.\-eE]+) ([\d.\-eE]+)", out)
    if not m:
        m = re.search(r"%%BoundingBox:\s*([\d.\-eE]+) ([\d.\-eE]+) ([\d.\-eE]+) ([\d.\-eE]+)", out)
    if not m:
        raise RuntimeError(f"gs bbox parse failed for {pdf_path}:\n{out}")
    return tuple(float(g) for g in m.groups())


def read_svg_header(svg_text):
    """Read width/height/viewBox from the <svg> tag. Returns (w, h, vb) where
    vb=(minx,miny,vw,vh). Falls back gracefully if some attrs absent."""
    w = h = None
    vb = None
    m = SVG_DIMS_RE["width"].search(svg_text)
    if m:
        w = float(m.group(1))
    m = SVG_DIMS_RE["height"].search(svg_text)
    if m:
        h = float(m.group(1))
    m = SVG_DIMS_RE["viewBox"].search(svg_text)
    if m:
        parts = [float(x) for x in m.group(1).split()]
        if len(parts) == 4:
            vb = tuple(parts)
    if vb is None and w is not None and h is not None:
        vb = (0.0, 0.0, w, h)
    if w is None and vb is not None:
        w = vb[2]
    if h is None and vb is not None:
        h = vb[3]
    return w, h, vb


def rewrite_svg(svg_text, new_vb, new_w, new_h):
    """Replace viewBox + width + height in the first <svg ...> tag."""
    vb_str = f'viewBox="{new_vb[0]:.4f} {new_vb[1]:.4f} {new_vb[2]:.4f} {new_vb[3]:.4f}"'
    text = svg_text
    if SVG_DIMS_RE["viewBox"].search(text):
        text = SVG_DIMS_RE["viewBox"].sub(vb_str, text, count=1)
    else:
        text = re.sub(r"<svg\b", "<svg " + vb_str, text, count=1)
    text = SVG_DIMS_RE["width"].sub(f'width="{new_w:.4f}pt"', text, count=1)
    text = SVG_DIMS_RE["height"].sub(f'height="{new_h:.4f}pt"', text, count=1)
    return text


def crop_one(svg_path, out_path, pad, rscript, gs):
    with open(svg_path, "r", encoding="utf-8") as f:
        svg_text = f.read()
    w, h, vb = read_svg_header(svg_text)
    if vb is None:
        raise RuntimeError(f"cannot read viewBox/size from {svg_path}")
    minx, miny, vw, vh = vb

    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tf:
        pdf_path = tf.name
    try:
        svg_to_pdf(svg_path, pdf_path, rscript)
        x0, y0, x1, y1 = pdf_bbox(pdf_path, gs)
    finally:
        if os.path.exists(pdf_path):
            os.remove(pdf_path)

    # gs bbox is in PDF pt with bottom-left origin; SVG user-units (rsvg 1:1) match
    # the rendered page so page height == vh (in pt). Convert y to SVG top-left.
    page_h = vh  # rsvg renders 1 user-unit -> 1 pt, page height == viewBox height
    crop_w = (x1 - x0) + 2 * pad
    crop_h = (y1 - y0) + 2 * pad
    # SVG-space top-left of the bbox: x stays; y flips against page height
    new_minx = minx + (x0 - pad)
    new_miny = miny + (page_h - y1 - pad)
    new_vb = (new_minx, new_miny, crop_w, crop_h)

    new_text = rewrite_svg(svg_text, new_vb, crop_w, crop_h)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(new_text)

    pct = 100.0 * (1.0 - (crop_w * crop_h) / (vw * vh)) if vw * vh else 0.0
    return (vw, vh, crop_w, crop_h, pct)


def main():
    ap = argparse.ArgumentParser(description="Crop SVGs to tight ink bbox (gs route).")
    ap.add_argument("paneldir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--pad", type=float, default=0.5)
    ap.add_argument("--rscript", default="Rscript")
    ap.add_argument("--gs", default="gs")
    args = ap.parse_args()

    outdir = args.out or os.path.join(args.paneldir, "cropped")
    os.makedirs(outdir, exist_ok=True)

    svgs = sorted(
        f for f in os.listdir(args.paneldir)
        if f.lower().endswith(".svg")
    )
    if not svgs:
        print(f"no .svg in {args.paneldir}", file=sys.stderr)
        sys.exit(1)

    print("name\torig_w\torig_h\tcrop_w\tcrop_h\tpct_trimmed")
    for name in svgs:
        src = os.path.join(args.paneldir, name)
        dst = os.path.join(outdir, name)
        try:
            vw, vh, cw, ch, pct = crop_one(src, dst, args.pad, args.rscript, args.gs)
            print(f"{name}\t{vw:.1f}\t{vh:.1f}\t{cw:.1f}\t{ch:.1f}\t{pct:.1f}%")
        except Exception as e:
            print(f"{name}\tERROR\t{e}", file=sys.stderr)


if __name__ == "__main__":
    main()
