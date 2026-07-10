#!/usr/bin/env python3
"""Apply the documented old-2 to final Supplementary Data 1 title correction."""
from __future__ import annotations

import argparse
from pathlib import Path

import fitz


def replace_title(page: fitz.Page, old: str, new: str) -> None:
    matches = page.search_for(old)
    if not matches:
        raise RuntimeError(f"Expected title not found: {old!r}")
    for rect in matches:
        page.add_redact_annot(rect, fill=(1, 1, 1))
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE)
    for rect in matches:
        page.insert_text((rect.x0, rect.y1 - 1), new, fontsize=max(7, rect.height * 0.82), fontname="helv")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-pdf", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--source-title", default="Supplementary Data 2")
    parser.add_argument("--final-title", default="Supplementary Data 1")
    args = parser.parse_args()
    document = fitz.open(args.input_pdf)
    if len(document) != 1:
        raise RuntimeError("Supplementary Data 1 finalization expects one page")
    replace_title(document[0], args.source_title, args.final_title)
    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
    document.save(args.output_pdf, garbage=4, deflate=True)


if __name__ == "__main__":
    main()
