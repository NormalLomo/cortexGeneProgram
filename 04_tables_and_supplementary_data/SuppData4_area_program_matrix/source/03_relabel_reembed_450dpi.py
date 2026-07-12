#!/usr/bin/env python3
"""Rasterize and re-embed the producer-rendered Supplementary Data 4 PDF."""
from __future__ import annotations

import argparse
from pathlib import Path

import fitz


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-pdf", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=450)
    args = parser.parse_args()

    source = fitz.open(args.input_pdf)
    if len(source) != 1:
        raise RuntimeError("Supplementary Data 4 finalization expects one page")
    page = source[0]
    pixmap = page.get_pixmap(dpi=args.dpi, alpha=False)
    output = fitz.open()
    rendered = output.new_page(width=page.rect.width, height=page.rect.height)
    rendered.insert_image(rendered.rect, pixmap=pixmap)
    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
    output.save(args.output_pdf, garbage=4, deflate=True)


if __name__ == "__main__":
    main()
