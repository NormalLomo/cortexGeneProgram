#!/usr/bin/env python3
"""Rebuild the 20-page supplementary-figure PDF from source figure PDFs.

The existing English legend pages are preserved byte-for-byte at the page-object
level; only figure pages are refreshed from FigS1.pdf through FigS10.pdf.
"""

from pathlib import Path
import os

from pypdf import PdfReader, PdfWriter


ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = ROOT / "source_figure_pdfs/supplementary_figures_pdf"
TARGET = ROOT / "supplementary_data/Supplementary_Figures.pdf"


def main() -> None:
    current = PdfReader(str(TARGET))
    if len(current.pages) != 20:
        raise RuntimeError(f"expected 20 current pages, found {len(current.pages)}")

    sources = [PdfReader(str(SOURCE_DIR / f"FigS{i}.pdf")) for i in range(1, 11)]
    if any(len(reader.pages) != 1 for reader in sources):
        raise RuntimeError("each FigS1-FigS10 source PDF must contain exactly one page")

    writer = PdfWriter()
    for index, source in enumerate(sources):
        writer.add_page(source.pages[0])
        writer.add_page(current.pages[2 * index + 1])
    writer.add_metadata(
        {
            "/Title": "Supplementary Figures S1-S10 with legends",
            "/Subject": "Self-contained supplementary figure and legend file",
        }
    )

    temp = TARGET.with_name(f".{TARGET.name}.figure-recompute.tmp")
    with temp.open("wb") as handle:
        writer.write(handle)
    rebuilt = PdfReader(str(temp))
    if len(rebuilt.pages) != 20:
        temp.unlink(missing_ok=True)
        raise RuntimeError(f"expected 20 rebuilt pages, found {len(rebuilt.pages)}")
    os.replace(temp, TARGET)
    print(f"rebuilt {TARGET} with {len(rebuilt.pages)} pages")


if __name__ == "__main__":
    main()
