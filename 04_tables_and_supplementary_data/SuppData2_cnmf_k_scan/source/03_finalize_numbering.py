#!/usr/bin/env python3
"""Apply the documented old-6 to final Supplementary Data 2 corrections."""
from __future__ import annotations

import argparse
from pathlib import Path

import fitz


def find_required_text(page: fitz.Page, text: str) -> list[fitz.Rect]:
    matches = page.search_for(text)
    if not matches:
        raise RuntimeError(f"Expected text not found: {text!r}")
    return matches


def redact_and_replace(page: fitz.Page, replacements: list[tuple[list[fitz.Rect], str]]) -> None:
    for matches, _ in replacements:
        for rect in matches:
            page.add_redact_annot(rect, fill=(1, 1, 1))
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE)
    for matches, new in replacements:
        for rect in matches:
            page.insert_text((rect.x0, rect.y1 - 1), new, fontsize=max(7, rect.height * 0.82), fontname="helv")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-pdf", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--source-title", default="Supplementary Data 6")
    parser.add_argument("--final-title", default="Supplementary Data 2")
    args = parser.parse_args()
    document = fitz.open(args.input_pdf)
    if len(document) != 1:
        raise RuntimeError("Supplementary Data 2 finalization expects one page")
    page = document[0]
    redact_and_replace(
        page,
        [
            (find_required_text(page, args.source_title), args.final_title),
            (find_required_text(page, "K = 60 (chosen)"), "60"),
        ],
    )
    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
    document.save(args.output_pdf, garbage=4, deflate=True)


if __name__ == "__main__":
    main()
