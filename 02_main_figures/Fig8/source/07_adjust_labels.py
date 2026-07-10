#!/usr/bin/env python3
"""Add minimal aging-source hierarchy annotations to Fig8 v17 SVG and export v18."""

from __future__ import annotations

import subprocess
import sys
import xml.etree.ElementTree as ET
import os
from pathlib import Path


SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"

ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)

WORKDIR = Path(os.environ.get("FIG8_WORKDIR", Path(__file__).resolve().parents[1]))
OUTDIR = WORKDIR / "outputs"
SRC_SVG = OUTDIR / "Fig8_imagegen_layout_v17.svg"
DST_STEM = OUTDIR / "Fig8_imagegen_layout_v18"
CAIRO_ENV = os.environ | {"DYLD_FALLBACK_LIBRARY_PATH": "/opt/homebrew/lib"}


def svg_el(tag: str, **attrs: str) -> ET.Element:
    return ET.Element(f"{{{SVG_NS}}}{tag}", attrs)


def add_text(
    parent: ET.Element,
    *,
    gid: str,
    x: float,
    y: float,
    text: str | None = None,
    style: str,
    transform: str | None = None,
    tspans: list[tuple[float, float, str]] | None = None,
) -> ET.Element:
    group = svg_el("g", id=gid)
    attrs = {"style": style, "x": f"{x:.6f}", "y": f"{y:.6f}"}
    if transform:
        attrs["transform"] = transform
    text_el = svg_el("text", **attrs)
    if text is not None:
        text_el.text = text
    for tx, ty, body in tspans or []:
        tspan = svg_el("tspan", x=f"{tx:.6f}", y=f"{ty:.6f}")
        tspan.text = body
        text_el.append(tspan)
    group.append(text_el)
    parent.append(group)
    return group


def add_path(parent: ET.Element, *, gid: str, d: str, style: str) -> ET.Element:
    group = svg_el("g", id=gid)
    group.append(svg_el("path", d=d, style=style))
    parent.append(group)
    return group


def main() -> int:
    if not SRC_SVG.exists():
        raise FileNotFoundError(SRC_SVG)

    tree = ET.parse(SRC_SVG)
    root = tree.getroot()
    figure = root.find(f"./{{{SVG_NS}}}g[@id='figure_1']")
    if figure is None:
        raise RuntimeError("figure_1 group not found in source SVG")

    ann = svg_el("g", id="v18_aging_hierarchy_annotations")

    # Panel f: subtle evidence split between AA/DA and curated references.
    add_path(
        ann,
        gid="v18_panel_f_split",
        d="M 130.110236 343.400000 L 130.110236 462.529134",
        style="fill: none; stroke: #c7c7c2; stroke-width: 0.55; stroke-dasharray: 1.4,1.4; stroke-linecap: butt",
    )
    add_text(
        ann,
        gid="v18_panel_f_note",
        x=74.200000,
        y=492.600000,
        style="font-size: 4px; font-family: 'Liberation Sans', 'Nimbus Sans', 'Arial', 'DejaVu Sans', sans-serif; fill: #666666",
        tspans=[
            (74.200000, 492.600000, "AA/DA: main evidence from Cell 2026 aging-pace proteins"),
            (74.200000, 497.200000, "SenMayo, CellAge+, CellAge-, GenAge: supporting curated references"),
        ],
    )

    # Panel h: make the left-side source hierarchy explicit without changing node geometry.
    add_path(
        ann,
        gid="v18_panel_h_split",
        d="M 21.940157 556.900000 L 58.200000 556.900000",
        style="fill: none; stroke: #d2d2cd; stroke-width: 0.50; stroke-linecap: butt",
    )
    add_text(
        ann,
        gid="v18_panel_h_main",
        x=22.800000,
        y=528.200000,
        text="main evidence: Cell 2026 aging-pace proteins",
        style="font-size: 3.7px; font-weight: 700; font-family: 'Liberation Sans', 'Nimbus Sans', 'Arial', 'DejaVu Sans', sans-serif; fill: #666666",
    )
    add_text(
        ann,
        gid="v18_panel_h_support",
        x=22.800000,
        y=558.800000,
        text="supporting curated senescence/aging references",
        style="font-size: 3.7px; font-weight: 700; font-family: 'Liberation Sans', 'Nimbus Sans', 'Arial', 'DejaVu Sans', sans-serif; fill: #666666",
    )

    figure.append(ann)

    dst_svg = DST_STEM.with_suffix(".svg")
    tree.write(dst_svg, encoding="utf-8", xml_declaration=True)

    subprocess.run(
        ["cairosvg", str(dst_svg), "-o", str(DST_STEM.with_suffix(".pdf"))],
        check=True,
        stdout=subprocess.DEVNULL,
        env=CAIRO_ENV,
    )
    subprocess.run(
        ["cairosvg", str(dst_svg), "-o", str(DST_STEM.with_suffix(".png")), "-d", "450"],
        check=True,
        stdout=subprocess.DEVNULL,
        env=CAIRO_ENV,
    )

    print(dst_svg)
    print(DST_STEM.with_suffix(".pdf"))
    print(DST_STEM.with_suffix(".png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
