#!/usr/bin/env python3
"""Relabel Supplementary Data 4, then rasterize and re-embed its one-page PDF at 450 dpi."""
from __future__ import annotations

import argparse
from pathlib import Path

import fitz


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-pdf", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--source-title", default="Fig. S6")
    parser.add_argument("--final-title", default="Supplementary Data 4")
    args = parser.parse_args()

    source = fitz.open(args.input_pdf)
    if len(source) != 1:
        raise RuntimeError("Supplementary Data 4 finalization expects one page")
    page = source[0]
    matches = page.search_for(args.source_title)
    if not matches:
        raise RuntimeError(f"Expected title not found: {args.source_title!r}")
    for rect in matches:
        page.add_redact_annot(rect, fill=(1, 1, 1))
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE)
    for rect in matches:
        page.insert_text((rect.x0, rect.y1 - 1), args.final_title, fontsize=max(7, rect.height * 0.82), fontname="helv")

    pixmap = page.get_pixmap(dpi=450, alpha=False)
    output = fitz.open()
    rendered = output.new_page(width=page.rect.width, height=page.rect.height)
    rendered.insert_image(rendered.rect, pixmap=pixmap)
    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
    output.save(args.output_pdf, garbage=4, deflate=True)


if __name__ == "__main__":
    main()
