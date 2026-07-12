#!/usr/bin/env python3
"""compose_svgutils.py — place ink-cropped panels per template -> composite SVG -> vector PDF.

WHY: svgutils lets you place mixed-engine SVG panels (ggplot + ComplexHeatmap +
ggraph + matplotlib) into ONE figure with precise scale/position — something
patchwork cannot do across engines. This reads layout_template.json (from
make_template.py), loads each cropped panel SVG, scales it UNIFORMLY to its box
(no stretch), positions it, adds bold panel tags, and writes the composite SVG,
then converts to a VECTOR PDF.

SCHEMA: reads BOTH flat and nested layout_template.json shapes:
  FLAT (older): page.{w_mm,h_mm}, panels[].{id, x_mm, y_mm, w_mm, h_mm}
  NESTED (current builders):
    page.{w,h}, panels[].{panel, box:{x,y,w,h}}, plus optional cropped_aspects
The reader dispatches on field presence (see page_size_mm / panel_box_mm).

THREE GOTCHAS THIS HANDLES:
  1. Cropped viewBox offset: ink-crop sets viewBox min-x/min-y != 0, but
     svgutils' moveto() ignores the local viewBox origin. We compensate by
     translating by -(minx,miny)*scale before moveto-ing the box origin, so the
     content lands exactly in the box (this was the "panel D compressed/shifted"
     bug).
  2. Tag overlapping the panel's own corner title: zero-margin panels put their
     title flush in the top-left corner, so a tag placed there collides. We nudge
     each tag OUT into the gutter/margin band (left+up of the box corner).
  3. Page-size pt unit: sc.Figure(width, height) writes width/height verbatim;
     this composer passes the "pt" suffix so cairo/rsvg preserves the intended
     millimetre page dimensions.

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
    --rscript      Rscript path

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

VIEWBOX_RE = re.compile(r'viewBox\s*=\s*["\']([\-\d.eE ]+)["\']')
WIDTH_RE = re.compile(r'\bwidth\s*=\s*["\']([\d.eE]+)(?:pt|px)?["\']')
HEIGHT_RE = re.compile(r'\bheight\s*=\s*["\']([\d.eE]+)(?:pt|px)?["\']')


def read_svg_box(svg_path):
    """Return (minx, miny, w, h) of the cropped panel viewBox, in the SVG's own
    user units (== pt, since panels are rendered 1 user-unit==1pt)."""
    with open(svg_path, "r", encoding="utf-8") as f:
        head = f.read(4000)
    m = VIEWBOX_RE.search(head)
    if m:
        parts = [float(x) for x in m.group(1).split()]
        if len(parts) == 4:
            return tuple(parts)
    mw, mh = WIDTH_RE.search(head), HEIGHT_RE.search(head)
    if mw and mh:
        return (0.0, 0.0, float(mw.group(1)), float(mh.group(1)))
    raise RuntimeError(f"cannot read viewBox/size from {svg_path}")


# ----- schema readers ------------------------------------------------------
# Two on-disk JSON schemas are supported:
#   FLAT (older skill schema):   panel = {id, x_mm, y_mm, w_mm, h_mm, ...},
#                                page  = {w_mm, h_mm, ...}
#   NESTED (current builders):   panel = {panel, box:{x,y,w,h}, ...},
#                                page  = {w, h, ...}   (mm units implicit)
# A NESTED template may also carry `cropped_aspects: {id: aspect}` — we read it
# only as a sanity reference (panel boxes already encode the AR via box.w/box.h).

def page_size_mm(page):
    """Return (w_mm, h_mm) regardless of flat-vs-nested page schema."""
    if "w_mm" in page and "h_mm" in page:
        return float(page["w_mm"]), float(page["h_mm"])
    if "w" in page and "h" in page:
        return float(page["w"]), float(page["h"])
    raise KeyError("page block missing w_mm/h_mm or w/h")


def panel_box_mm(p):
    """Return (id, x_mm, y_mm, w_mm, h_mm) for either schema."""
    pid = p.get("id", p.get("panel"))
    if pid is None:
        raise KeyError("panel record missing 'id' or 'panel'")
    if "box" in p and isinstance(p["box"], dict):
        b = p["box"]
        return str(pid), float(b["x"]), float(b["y"]), float(b["w"]), float(b["h"])
    return (str(pid),
            float(p["x_mm"]), float(p["y_mm"]),
            float(p["w_mm"]), float(p["h_mm"]))


def inject_bg_rect(svg_path, page_w_pt, page_h_pt, color="#ffffff"):
    """Insert a full-page opaque background rect as the FIRST child of <svg>,
    so it sits UNDER every panel. WHY: ink-cropped panels have transparent
    backgrounds; without a page rect the transparent areas render BLACK in
    PDF/PNG (this was the 'black bottom/band' bug). Regex-insert right after the
    opening <svg ...> tag (robust; avoids lxml namespace round-trip of svgutils
    output). Idempotent-ish: only inserts if no bg-rect marker is present."""
    with open(svg_path, "r", encoding="utf-8") as f:
        svg = f.read()
    if 'data-bg-rect="1"' in svg:
        return
    rect = ('<rect data-bg-rect="1" x="0" y="0" width="%.4f" height="%.4f" '
            'fill="%s"/>') % (page_w_pt, page_h_pt, color)
    m = re.search(r"<svg\b[^>]*>", svg)
    if not m:
        raise RuntimeError("no <svg> opening tag in composed output")
    pos = m.end()
    svg = svg[:pos] + rect + svg[pos:]
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg)


def compose(cropdir, template, out_svg, tag_size, tag_dx, tag_dy, with_tags, title,
            tag_sequential=False, tag_start=0, bg_color="#ffffff"):
    pg = template["page"]
    # work the whole figure in pt (svgutils default unit)
    page_w_mm, page_h_mm = page_size_mm(pg)
    page_w_pt = page_w_mm * PT_PER_MM
    page_h_pt = page_h_mm * PT_PER_MM

    # Sanity reference: if cropped_aspects is present, compare to box AR; this
    # surfaces stale-template drift (panel-d-was-redesigned-but-template-wasn't).
    cropped_aspects = template.get("cropped_aspects") or {}

    elements = []
    tags = "abcdefghijklmnopqrstuvwxyz"
    # If every panel id is already a single letter (a,b,c...), the ids ARE the
    # tags — use them, so a nested/tangram traversal order (not alphabetical)
    # still labels each panel with its own letter. Else fall back to positional.
    # tag_sequential forces consecutive a,b,c... by reading order — use this when
    # the panel set is a curated subset whose ids skip letters (e.g. a,b,c,f,h,j
    # after merging/dropping panels) so the printed tags don't jump.
    panel_ids = [panel_box_mm(p)[0] for p in template["panels"]]
    tag_by_id = (not tag_sequential) and all(
        len(pid) == 1 and pid.isalpha() for pid in panel_ids
    )

    for i, p in enumerate(template["panels"]):
        pid, x_mm, y_mm, w_mm, h_mm = panel_box_mm(p)
        svg_path = os.path.join(cropdir, pid + ".svg")
        if not os.path.exists(svg_path):
            # also accept fig1_<id>.svg naming (matches current svg_panels/ layout)
            alt = os.path.join(cropdir, f"fig1_{pid}.svg")
            if os.path.exists(alt):
                svg_path = alt
            else:
                raise FileNotFoundError(f"panel SVG not found: {svg_path} (or {alt})")
        minx, miny, src_w, src_h = read_svg_box(svg_path)

        box_w_pt = w_mm * PT_PER_MM
        box_h_pt = h_mm * PT_PER_MM
        box_x_pt = x_mm * PT_PER_MM
        box_y_pt = y_mm * PT_PER_MM

        # uniform scale to fit the box (src is already the right AR after crop,
        # so sx ~= sy; take min to be safe against rounding -> no stretch).
        # If the box AR drifted from the panel AR (stale template), the scale
        # collapses to the tighter axis and we letterbox into the box — emit a
        # warning so the upstream template can be rebuilt rather than silently
        # producing a thin strip inside a fat box (this was the panel-d bug).
        sx = box_w_pt / src_w
        sy = box_h_pt / src_h
        scale = min(sx, sy)

        src_ar = src_w / src_h
        box_ar = box_w_pt / box_h_pt
        ar_gap = abs(src_ar - box_ar)
        cropped_ref = cropped_aspects.get(pid)
        if ar_gap > 0.05:
            ref_msg = f"  cropped_aspects[{pid}]={cropped_ref}" if cropped_ref else ""
            print(f"WARN: panel {pid} src AR {src_ar:.3f} != box AR {box_ar:.3f} "
                  f"(gap {ar_gap:.3f}); will letterbox.{ref_msg} "
                  f"Rebuild layout_template.json from current SVGs.", file=sys.stderr)

        panel = sc.SVG(svg_path)
        panel.scale(scale)
        # compensate the cropped viewBox origin (svgutils moveto ignores it),
        # THEN move the (now origin-aligned) content to the box top-left.
        panel.moveto(box_x_pt - minx * scale, box_y_pt - miny * scale)
        elements.append(panel)

        if with_tags:
            tx = (x_mm + tag_dx) * PT_PER_MM
            ty = (y_mm + tag_dy) * PT_PER_MM
            # clamp into the page so tags in the first column/row stay visible
            tx = max(tx, 1.0)
            ty = max(ty, tag_size)
            if tag_by_id:
                label = pid
            else:
                idx = i + tag_start
                label = tags[idx] if idx < len(tags) else f"#{idx}"
            elements.append(
                sc.Text(label, tx, ty, size=tag_size, weight="bold", font="sans-serif")
            )

    if title:
        elements.append(
            sc.Text(title, 4 * PT_PER_MM, 0.6 * tag_size + 2,
                    size=tag_size, weight="bold", font="sans-serif")
        )

    # svgutils.compose.Figure writes width/height with the literal suffix you
    # pass — keep the "pt" suffix so cairo/rsvg interprets the MediaBox in pt
    # (524.41pt = 185mm). Without the unit, cairo treats the number as user
    # units == pt, which still works here but is fragile if a caller passes
    # a bare number; the explicit "pt" guarantees the intended physical size.
    fig = sc.Figure(f"{page_w_pt}pt", f"{page_h_pt}pt", *elements)
    fig.save(out_svg)
    # BUG1 fix: lay an opaque page background UNDER all panels so transparent
    # ink-cropped panels don't render as black. Skip with bg_color="none".
    if bg_color and bg_color.lower() != "none":
        inject_bg_rect(out_svg, page_w_pt, page_h_pt, bg_color)
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
    ap.add_argument("--bg-color", default="#ffffff",
                    help="full-page background fill placed UNDER all panels "
                         "(default #ffffff; use 'none' to keep transparent)")
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
        bg_color=args.bg_color,
    )
    print(f"wrote composite SVG: {out_svg}")

    if args.out_pdf:
        svg_to_pdf(out_svg, args.out_pdf, args.rscript)
        print(f"wrote vector PDF: {args.out_pdf}")
        print("verify vector with:  pdffonts " + args.out_pdf + "  (embedded fonts > 0)")


if __name__ == "__main__":
    main()
