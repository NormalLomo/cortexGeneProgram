#!/usr/bin/env python3
"""make_template.py — measure cropped panel aspects -> gutter-0 row-pack template.

WHY: after ink-crop each panel SVG's box == its content, so its true aspect
(w/h) is known. A clean multi-panel figure packs panels so that EVERY box's
aspect == its panel's natural aspect (no stretch, no letterbox) AND panels ABUT
with GUTTER 0 (the owner's strong preference — same as patchwork `design`).
This is a row-based packer: you choose which panels share each row; within a
row, panels abut to a common content width and the row height is solved so no
panel is stretched.

ROW MATH (gutter 0):
  Each panel i has natural aspect a_i = w_i / h_i.
  In a row all panels share height H_row and abut horizontally.
  panel width = a_i * H_row.  Sum of widths == content_width =>
    H_row = content_width / Σ(a_i).
  Each box: w = a_i * H_row, h = H_row  -> box aspect == a_i exactly (no stretch).
Pair WIDE panels (banner/network/ridgeline) with tall/square ones in the same
row so no row goes over-tall; aim page AR ~0.7-0.8 (A4 portrait).

USAGE
  python make_template.py CROPDIR --rows "a,b;c,d,e;f" \
      [--page-width-mm 185] [--margin-mm 4] [--out OUTDIR] [--rscript RSCRIPT] [--gs GS]
    CROPDIR        dir of ink-cropped *.svg (panel id = filename stem)
    --rows         row spec: panels separated by ',', rows separated by ';'.
                   e.g. "a,b;c,d,e;f"  -> row1=[a,b] row2=[c,d,e] row3=[f]
                   IDs must match SVG stems in CROPDIR.
    --page-width-mm  content width (excl. outer margin). Default 185 (~A4 usable).
    --margin-mm    outer margin for tags + figure title only. Default 4.
    --out          output dir (default CROPDIR). Writes layout_template.json
                   + layout_template_wireframe.png.
    --row-weights  optional ';'-separated multipliers to bias a row taller/shorter
                   (e.g. "1;1.3;1" makes the middle row 1.3x its solved height).
                   Use to emphasise an important data row. Default all 1.

It measures each SVG's aspect by reading its (cropped) viewBox directly — no
render needed. Verifies max_fill_gap_from_1 == 0 (box AR == panel AR for all).

OUTPUT layout_template.json:
  { page: {w_mm, h_mm, margin_mm, content_w_mm},
    panels: [ {id, x_mm,y_mm,w_mm,h_mm, x_rel,y_rel,w_rel,h_rel,
               natural_aspect, box_aspect, tag_pos:{x_mm,y_mm}} ... ] }
plus a labeled wireframe PNG (matplotlib) to SHOW THE OWNER before composing.
"""
import argparse
import json
import os
import re
import sys

VIEWBOX_RE = re.compile(r'viewBox\s*=\s*"([\-\d.eE ]+)"')
WIDTH_RE = re.compile(r'\bwidth\s*=\s*"([\d.eE]+)(?:pt|px)?"')
HEIGHT_RE = re.compile(r'\bheight\s*=\s*"([\d.eE]+)(?:pt|px)?"')


def svg_aspect(svg_path):
    """Return (w, h, aspect) from the SVG header (cropped viewBox)."""
    with open(svg_path, "r", encoding="utf-8") as f:
        head = f.read(4000)
    m = VIEWBOX_RE.search(head)
    if m:
        parts = [float(x) for x in m.group(1).split()]
        if len(parts) == 4:
            w, h = parts[2], parts[3]
            return w, h, w / h
    mw, mh = WIDTH_RE.search(head), HEIGHT_RE.search(head)
    if mw and mh:
        w, h = float(mw.group(1)), float(mh.group(1))
        return w, h, w / h
    raise RuntimeError(f"cannot read size/viewBox from {svg_path}")


def parse_rows(rows_spec):
    rows = []
    for row in rows_spec.split(";"):
        ids = [s.strip() for s in row.split(",") if s.strip()]
        if ids:
            rows.append(ids)
    return rows


# ---------------------------------------------------------------------------
# NESTED (tangram) layout: a slicing tree of H/V nodes over leaf panels.
#   H[...]  horizontal split: children ABUT side-by-side, share one HEIGHT.
#   V[...]  vertical split:   children STACK, share one WIDTH.
#   leaf    a panel id.
# This generalises flat rows: a tall panel can be a column SPANNING the height
# of two stacked shorter panels beside it ("拼七巧板", not rigid rows). With
# fixed leaf aspects every leaf keeps its natural aspect (no stretch), gutter 0.
#
# Aspect (w/h) of any node, bottom-up:
#   leaf:   a_i
#   H node: Σ child_aspect            (shared height -> widths add)
#   V node: 1 / Σ(1/child_aspect)     (shared width  -> heights add)
# Sizes top-down: root width = content_width; height = width / root_aspect.
# Example --layout:  "V[ a , H[ V[b,c] , d ] , H[g,j] ]"
# ---------------------------------------------------------------------------

def parse_layout(spec):
    """Parse 'H[a,b,V[c,d]]' into a tree of dicts.
    leaf -> {'t':'leaf','id':str};  node -> {'t':'H'/'V','kids':[...]}."""
    pos = 0
    s = spec

    def skip_ws():
        nonlocal pos
        while pos < len(s) and s[pos].isspace():
            pos += 1

    def parse_node():
        nonlocal pos
        skip_ws()
        # H[...] or V[...] ?
        if pos < len(s) and s[pos] in "HV" and pos + 1 < len(s) and s[pos + 1:pos + 2] == "[":
            typ = s[pos]
            pos += 2  # consume 'H[' / 'V['
            kids = []
            while True:
                kids.append(parse_node())
                skip_ws()
                if pos < len(s) and s[pos] == ",":
                    pos += 1
                    continue
                if pos < len(s) and s[pos] == "]":
                    pos += 1
                    break
                raise ValueError(f"layout parse error near pos {pos}: {s[pos:pos+12]!r}")
            return {"t": typ, "kids": kids}
        # leaf identifier
        start = pos
        while pos < len(s) and (s[pos].isalnum() or s[pos] in "_-."):
            pos += 1
        ident = s[start:pos].strip()
        if not ident:
            raise ValueError(f"layout parse error: empty leaf near pos {pos}")
        return {"t": "leaf", "id": ident}

    tree = parse_node()
    skip_ws()
    if pos != len(s):
        raise ValueError(f"layout trailing chars: {s[pos:]!r}")
    return tree


def node_aspect(node, aspect_map):
    if node["t"] == "leaf":
        if node["id"] not in aspect_map:
            raise KeyError(f"panel '{node['id']}' has no SVG in cropdir")
        return aspect_map[node["id"]]
    kids_a = [node_aspect(k, aspect_map) for k in node["kids"]]
    if node["t"] == "H":
        return sum(kids_a)
    # V
    return 1.0 / sum(1.0 / a for a in kids_a)


def assign_layout(node, x, y, w, h, aspect_map, out, tag_dx=-3.0, tag_dy=-1.0):
    """Recursively place leaves; append leaf boxes to out (mm)."""
    if node["t"] == "leaf":
        a = aspect_map[node["id"]]
        out.append({
            "id": node["id"],
            "x_mm": round(x, 3), "y_mm": round(y, 3),
            "w_mm": round(w, 3), "h_mm": round(h, 3),
            "natural_aspect": round(a, 4),
            "box_aspect": round(w / h, 4),
            "tag_pos": {"x_mm": round(x + tag_dx, 3), "y_mm": round(y + tag_dy, 3)},
        })
        return
    if node["t"] == "H":
        # children share height h; widths = a_child * h
        cx = x
        for k in node["kids"]:
            ka = node_aspect(k, aspect_map)
            kw = ka * h
            assign_layout(k, cx, y, kw, h, aspect_map, out, tag_dx, tag_dy)
            cx += kw
    else:  # V: children share width w; heights = w / a_child
        cy = y
        for k in node["kids"]:
            ka = node_aspect(k, aspect_map)
            kh = w / ka
            assign_layout(k, x, cy, w, kh, aspect_map, out, tag_dx, tag_dy)
            cy += kh


def build_template_nested(cropdir, layout_spec, page_width_mm, margin_mm):
    tree = parse_layout(layout_spec)
    # collect leaf ids, measure aspects
    aspect = {}

    def walk(n):
        if n["t"] == "leaf":
            path = os.path.join(cropdir, n["id"] + ".svg")
            if not os.path.exists(path):
                raise FileNotFoundError(f"panel SVG not found: {path}")
            aspect[n["id"]] = svg_aspect(path)[2]
        else:
            for k in n["kids"]:
                walk(k)

    walk(tree)
    content_w = page_width_mm
    root_a = node_aspect(tree, aspect)
    content_h = content_w / root_a

    panels = []
    assign_layout(tree, margin_mm, margin_mm, content_w, content_h, aspect, panels)

    page_w = page_width_mm + 2 * margin_mm
    page_h = content_h + 2 * margin_mm
    max_gap = max(abs(p["box_aspect"] - p["natural_aspect"]) for p in panels) if panels else 0.0
    for p in panels:
        p["x_rel"] = round(p["x_mm"] / page_w, 5)
        p["y_rel"] = round(p["y_mm"] / page_h, 5)
        p["w_rel"] = round(p["w_mm"] / page_w, 5)
        p["h_rel"] = round(p["h_mm"] / page_h, 5)
    return {
        "page": {
            "w_mm": round(page_w, 3), "h_mm": round(page_h, 3),
            "margin_mm": margin_mm, "content_w_mm": content_w,
            "aspect": round(page_w / page_h, 4),
        },
        "layout_spec": layout_spec,
        "panels": panels,
        "max_fill_gap_from_1": round(max_gap, 6),
    }


def build_template(cropdir, rows, page_width_mm, margin_mm, row_weights):
    # measure aspects
    aspect = {}
    for ids in rows:
        for pid in ids:
            path = os.path.join(cropdir, pid + ".svg")
            if not os.path.exists(path):
                raise FileNotFoundError(f"panel SVG not found: {path}")
            w, h, a = svg_aspect(path)
            aspect[pid] = a

    content_w = page_width_mm
    if row_weights is None:
        row_weights = [1.0] * len(rows)
    if len(row_weights) != len(rows):
        raise ValueError("--row-weights count must match number of rows")

    panels = []
    y = margin_mm
    max_gap = 0.0
    for ridx, ids in enumerate(rows):
        sum_a = sum(aspect[p] for p in ids)
        row_h = (content_w / sum_a) * row_weights[ridx]
        x = margin_mm
        for pid in ids:
            a = aspect[pid]
            w = a * row_h
            h = row_h
            box_a = w / h
            max_gap = max(max_gap, abs(box_a - a))
            panels.append({
                "id": pid,
                "x_mm": round(x, 3),
                "y_mm": round(y, 3),
                "w_mm": round(w, 3),
                "h_mm": round(h, 3),
                "natural_aspect": round(a, 4),
                "box_aspect": round(box_a, 4),
                # tag nudged into the top-left margin band, OUT of panel content
                "tag_pos": {"x_mm": round(x - 3.0, 3), "y_mm": round(y - 1.0, 3)},
            })
            x += w  # gutter 0: abut
        y += row_h  # gutter 0 between rows
    page_w = page_width_mm + 2 * margin_mm
    page_h = y + margin_mm

    # relative coords (fraction of page) for renderers that want 0-1
    for p in panels:
        p["x_rel"] = round(p["x_mm"] / page_w, 5)
        p["y_rel"] = round(p["y_mm"] / page_h, 5)
        p["w_rel"] = round(p["w_mm"] / page_w, 5)
        p["h_rel"] = round(p["h_mm"] / page_h, 5)

    template = {
        "page": {
            "w_mm": round(page_w, 3),
            "h_mm": round(page_h, 3),
            "margin_mm": margin_mm,
            "content_w_mm": content_w,
            "aspect": round(page_w / page_h, 4),
        },
        "panels": panels,
        "max_fill_gap_from_1": round(max_gap, 6),
    }
    return template


def draw_wireframe(template, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    pg = template["page"]
    fig, ax = plt.subplots(figsize=(pg["w_mm"] / 25.4, pg["h_mm"] / 25.4))
    ax.set_xlim(0, pg["w_mm"])
    ax.set_ylim(0, pg["h_mm"])
    ax.invert_yaxis()  # top-left origin like SVG/page
    ax.set_aspect("equal")
    tags = "abcdefghijklmnopqrstuvwxyz"
    for i, p in enumerate(template["panels"]):
        ax.add_patch(Rectangle(
            (p["x_mm"], p["y_mm"]), p["w_mm"], p["h_mm"],
            fill=False, edgecolor="black", linewidth=1.0,
        ))
        ax.text(
            p["x_mm"] + p["w_mm"] / 2, p["y_mm"] + p["h_mm"] / 2,
            f'{p["id"]}\nAR {p["natural_aspect"]:.2f}',
            ha="center", va="center", fontsize=8,
        )
        ax.text(
            p["tag_pos"]["x_mm"], p["tag_pos"]["y_mm"],
            tags[i] if i < len(tags) else f"#{i}",
            ha="left", va="bottom", fontsize=11, fontweight="bold", color="red",
        )
    ax.set_title(
        f'wireframe  page {pg["w_mm"]:.0f}x{pg["h_mm"]:.0f}mm  AR {pg["aspect"]:.2f}'
        f'  fill_gap={template["max_fill_gap_from_1"]:.4f}',
        fontsize=9,
    )
    ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Gutter-0 template from cropped panels (row-pack or nested tangram).")
    ap.add_argument("cropdir")
    ap.add_argument("--rows", default=None, help='flat rows, e.g. "a,b;c,d,e;f"')
    ap.add_argument("--layout", default=None,
                    help='nested tangram, e.g. "V[a,H[V[b,c],d],H[g,j]]" '
                         '(H=abut/share-height, V=stack/share-width). Overrides --rows.')
    ap.add_argument("--page-width-mm", type=float, default=185.0)
    ap.add_argument("--margin-mm", type=float, default=4.0)
    ap.add_argument("--row-weights", default=None, help='row mode only, e.g. "1;1.3;1"')
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.layout:
        template = build_template_nested(
            args.cropdir, args.layout, args.page_width_mm, args.margin_mm
        )
    elif args.rows:
        rows = parse_rows(args.rows)
        if not rows:
            print("empty --rows", file=sys.stderr); sys.exit(1)
        row_weights = None
        if args.row_weights:
            row_weights = [float(x) for x in args.row_weights.split(";")]
        template = build_template(
            args.cropdir, rows, args.page_width_mm, args.margin_mm, row_weights
        )
    else:
        print("need --layout (nested) or --rows (flat)", file=sys.stderr); sys.exit(1)

    outdir = args.out or args.cropdir
    os.makedirs(outdir, exist_ok=True)
    json_path = os.path.join(outdir, "layout_template.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(template, f, indent=2)
    png_path = os.path.join(outdir, "layout_template_wireframe.png")
    try:
        draw_wireframe(template, png_path)
    except Exception as e:
        print(f"wireframe render skipped: {e}", file=sys.stderr)
        png_path = None

    print(f"page: {template['page']['w_mm']:.0f} x {template['page']['h_mm']:.0f} mm  "
          f"AR {template['page']['aspect']:.3f}")
    print(f"max_fill_gap_from_1: {template['max_fill_gap_from_1']}  "
          f"(0 == every box AR matches its panel, no stretch)")
    print(f"wrote: {json_path}")
    if png_path:
        print(f"wrote: {png_path}  <- SHOW THIS TO THE OWNER before composing")


if __name__ == "__main__":
    main()
