#!/usr/bin/env python3
"""compose_svgutils.py — place ink-cropped panels per template -> composite SVG -> vector PDF.

WHY: svgutils lets you place mixed-engine SVG panels (ggplot + ComplexHeatmap +
ggraph + matplotlib) into ONE figure with precise scale/position — something
patchwork cannot do across engines. This reads layout_template.json (from
make_template.py), loads each cropped panel SVG, scales it UNIFORMLY to its box
(no stretch), positions it, adds bold panel tags, and writes the composite SVG,
then converts to a VECTOR PDF.

TWO GOTCHAS THIS HANDLES:
  1. Cropped viewBox offset: ink-crop sets viewBox min-x/min-y != 0, but
     svgutils' moveto() ignores the reference viewBox origin. We compensate by
     translating by -(minx,miny)*scale before moveto-ing the box origin, so the
     content lands exactly in the box (this was the "panel D compressed/shifted"
     bug).
  2. Tag overlapping the panel's own corner title: zero-margin panels put their
     title flush in the top-left corner, so a tag placed there collides. We nudge
     each tag OUT into the gutter/margin band (left+up of the box corner).

SVG->PDF: use R rsvg::rsvg_pdf (reliable). cairosvg mis-scales the page by 0.75
(96<->72 dpi confusion) — do NOT use it. Verify vector with `pdffonts` after.

USAGE
  python compose_svgutils.py CROPDIR TEMPLATE_JSON --out-svg OUT.svg --out-pdf OUT.pdf \
      [--tag-size 11] [--tag-dx -3] [--tag-dy -1] [--no-tags] \
      [--title "Figure 1"] [--rscript RSCRIPT]
    CROPDIR        dir of ink-cropped panel SVGs (id = stem, matches template)
    TEMPLATE_JSON  layout_template.json from make_template.py (mm units)
    --out-svg      composite SVG path
    --out-pdf      composite PDF path (vector, via rsvg)
    --tag-size     tag font pt (default 11, bold)
    --tag-dx/dy    tag nudge in mm relative to box top-left (default -3,-1 => out)
    --no-tags      skip a,b,c tags (e.g. tags already baked into panels)
    --title        optional figure title drawn in the top margin band
    --rscript      Rscript path (remote often /usr/local/bin/Rscript)

Requires: python `svgutils`, R `rsvg`.
"""
import argparse
import json
import os
import re
import subprocess
import sys

from svgutils import compose as sc

MM_PER_PT = 25.4 / 72.0   # 1 pt = 0.3528 mm
PT_PER_MM = 72.0 / 25.4   # 1 mm = 2.8346 pt

VIEWBOX_RE = re.compile(r"""viewBox\s*=\s*["']([\-\d.eE,\s]+)["']""")
WIDTH_RE = re.compile(r"""\bwidth\s*=\s*["']([\d.eE]+)(?:pt|px)?["']""")
HEIGHT_RE = re.compile(r"""\bheight\s*=\s*["']([\d.eE]+)(?:pt|px)?["']""")


def read_svg_box(svg_path):
    """Return (minx, miny, w, h) of the cropped panel viewBox, in the SVG's own
    user units (== pt, since panels are rendered 1 user-unit==1pt)."""
    with open(svg_path, "r", encoding="utf-8") as f:
        head = f.read(4000)
    m = VIEWBOX_RE.search(head)
    if m:
        parts = [float(x) for x in re.split(r"[\s,]+", m.group(1).strip()) if x]
        if len(parts) == 4:
            return tuple(parts)
    mw, mh = WIDTH_RE.search(head), HEIGHT_RE.search(head)
    if mw and mh:
        return (0.0, 0.0, float(mw.group(1)), float(mh.group(1)))
    raise RuntimeError(f"cannot read viewBox/size from {svg_path}")


def compose(cropdir, template, out_svg, tag_size, tag_dx, tag_dy, with_tags, title,
            tag_sequential=False, tag_start=0, skip_tag_ids=()):
    skip_tag_ids = set(skip_tag_ids)
    pg = template["page"]
    # work the whole figure in pt (svgutils default unit)
    page_w_pt = pg["w_mm"] * PT_PER_MM
    page_h_pt = pg["h_mm"] * PT_PER_MM

    elements = []
    tags = "abcdefghijklmnopqrstuvwxyz"
    # If every panel id is already a single letter (a,b,c...), the ids ARE the
    # tags — use them, so a nested/tangram traversal order (not alphabetical)
    # still labels each panel with its own letter. Else fall back to positional.
    # tag_sequential forces consecutive a,b,c... by reading order — use this when
    # the panel set is a curated subset whose ids skip letters (e.g. a,b,c,f,h,j
    # after merging/dropping panels) so the printed tags don't jump.
    tag_by_id = (not tag_sequential) and all(
        len(str(p["id"])) == 1 and str(p["id"]).isalpha()
        for p in template["panels"]
    )

    seq_counter = 0  # sequential-tag index that advances ONLY on non-skipped panels
    for i, p in enumerate(template["panels"]):
        svg_path = os.path.join(cropdir, p["id"] + ".svg")
        minx, miny, src_w, src_h = read_svg_box(svg_path)

        box_w_pt = p["w_mm"] * PT_PER_MM
        box_h_pt = p["h_mm"] * PT_PER_MM
        box_x_pt = p["x_mm"] * PT_PER_MM
        box_y_pt = p["y_mm"] * PT_PER_MM

        # uniform scale to fit the box (src is already the right AR after crop,
        # so sx ~= sy; take min to be safe against rounding -> no stretch)
        sx = box_w_pt / src_w
        sy = box_h_pt / src_h
        scale = min(sx, sy)

        panel = sc.SVG(svg_path)
        panel.scale(scale)
        # compensate the cropped viewBox origin (svgutils moveto ignores it),
        # THEN move the (now origin-aligned) content to the box top-left.
        panel.moveto(box_x_pt - minx * scale, box_y_pt - miny * scale)
        elements.append(panel)

        if with_tags and str(p["id"]) not in skip_tag_ids:
            tx = (p["x_mm"] + tag_dx) * PT_PER_MM
            ty = (p["y_mm"] + tag_dy) * PT_PER_MM
            # clamp into the page so tags in the first column/row stay visible
            tx = max(tx, 1.0)
            ty = max(ty, tag_size)
            if tag_by_id:
                label = str(p["id"])
            else:
                # sequential: count only the panels that actually get a tag, so a
                # skipped title band ('t') does not consume a letter.
                idx = seq_counter + tag_start
                label = tags[idx] if idx < len(tags) else f"#{idx}"
                seq_counter += 1
            elements.append(
                sc.Text(label, tx, ty, size=tag_size, weight="bold", font="sans-serif")
            )

    if title:
        elements.append(
            sc.Text(title, 4 * PT_PER_MM, 0.6 * tag_size + 2,
                    size=tag_size, weight="bold", font="sans-serif")
        )

    fig = sc.Figure(f"{page_w_pt}pt", f"{page_h_pt}pt", *elements)
    fig.save(out_svg)
    return out_svg


def svg_to_pdf(svg_path, pdf_path, rscript):
    r_code = f'rsvg::rsvg_pdf("{svg_path}", "{pdf_path}")'
    res = subprocess.run([rscript, "-e", r_code], capture_output=True, text=True)
    if res.returncode != 0 or not os.path.exists(pdf_path):
        raise RuntimeError(f"rsvg_pdf failed:\n{res.stderr or res.stdout}")


def main():
    ap = argparse.ArgumentParser(description="Compose cropped panels per template -> SVG -> vector PDF.")
    ap.add_argument("cropdir")
    ap.add_argument("template_json")
    ap.add_argument("--out-svg", required=True)
    ap.add_argument("--out-pdf", default=None)
    ap.add_argument("--tag-size", type=float, default=11.0)
    ap.add_argument("--tag-dx", type=float, default=-3.0, help="mm, relative to box top-left")
    ap.add_argument("--tag-dy", type=float, default=-1.0, help="mm, relative to box top-left")
    ap.add_argument("--no-tags", action="store_true")
    ap.add_argument("--tag-sequential", action="store_true",
                    help="force consecutive a,b,c... by reading order (use when the "
                         "panel-id set skips letters after merging/dropping panels)")
    ap.add_argument("--tag-start", default="a",
                    help="first tag letter for sequential mode (e.g. 'b' to reserve 'a' "
                         "for a separately-made banner). Implies --tag-sequential.")
    ap.add_argument("--title", default=None)
    ap.add_argument("--skip-tag-ids", default="",
                    help="comma-sep panel ids to NOT draw a tag for (e.g. a title band 't')")
    ap.add_argument("--rscript", default="Rscript")
    args = ap.parse_args()

    with open(args.template_json, "r", encoding="utf-8") as f:
        template = json.load(f)

    tag_start = ord(args.tag_start.lower()) - ord("a")
    out_svg = compose(
        args.cropdir, template, args.out_svg,
        args.tag_size, args.tag_dx, args.tag_dy,
        not args.no_tags, args.title,
        tag_sequential=args.tag_sequential or tag_start != 0,
        tag_start=tag_start,
        skip_tag_ids=[s.strip() for s in args.skip_tag_ids.split(",") if s.strip()],
    )
    print(f"wrote composite SVG: {out_svg}")

    if args.out_pdf:
        svg_to_pdf(out_svg, args.out_pdf, args.rscript)
        print(f"wrote vector PDF: {args.out_pdf}")
        print("verify vector with:  pdffonts " + args.out_pdf + "  (embedded fonts > 0)")


if __name__ == "__main__":
    main()
