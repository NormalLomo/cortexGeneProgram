#!/usr/bin/env python3
"""Produce the current regional supplementary figures and table package.

This entry point only edits the retained S4--S6 presentation and removes the
retired nucleus-screen labels/fields.  It consumes existing figure PDFs and
stored workbook values; it does not fit models or recompute statistics.
"""
from __future__ import annotations

import io
import re
import subprocess
import sys
from pathlib import Path

import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[3]
FIG_DIR = ROOT / "source_figure_pdfs" / "supplementary_figures_pdf"
PNG_DIR = ROOT / "figures_png" / "supplementary_figures"
TABLE_S4 = ROOT / "tables" / "TableS4_program_variability_validity_human.xlsx"
TABLES = ROOT / "tables" / "Supplementary_Tables_S1-S6_human.xlsx"
SUPP_DIR = ROOT / "supplementary_data"
ATTACH_DIR = SUPP_DIR / "gigascience_supplementary_material"

NEUTRAL = (0.55, 0.55, 0.55)
OLD_COLORS = (
    (0.769, 0.306, 0.322),  # cohort-robust marker
    (0.886, 0.635, 0.173),  # cohort-sensitive marker
)


def _same_rgb(a, b, tol=0.002):
    return a is not None and all(abs(float(a[i]) - b[i]) <= tol for i in range(3))


def _draw_neutralized(page, drawing):
    shape = page.new_shape()
    for item in drawing.get("items", []):
        kind = item[0]
        if kind == "l":
            shape.draw_line(item[1], item[2])
        elif kind == "c":
            shape.draw_bezier(item[1], item[2], item[3], item[4])
        elif kind == "re":
            shape.draw_rect(item[1])
        elif kind == "qu":
            shape.draw_quad(item[1])
    if drawing.get("fill") is not None:
        fill = NEUTRAL
    else:
        fill = None
    if drawing.get("color") is not None:
        color = NEUTRAL
    else:
        color = None
    shape.finish(
        color=color,
        fill=fill,
        width=float(drawing.get("width", 0.5)),
        closePath=bool(drawing.get("closePath", False)),
        fill_opacity=float(drawing.get("fill_opacity", 1.0)),
        stroke_opacity=float(drawing.get("stroke_opacity", 1.0)),
    )
    shape.commit()


def _remove_text(page, predicate):
    pending = []
    for block in page.get_text("dict").get("blocks", []):
        if block.get("type") != 0:
            continue
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                text = span.get("text", "")
                if predicate(text):
                    pending.append(span)
    for span in pending:
        page.add_redact_annot(fitz.Rect(span["bbox"]), fill=(1, 1, 1), cross_out=False)
    if pending:
        page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=0, text=0)


def _mask_text(page, predicate):
    """Cover retired legend words without invoking recursive PDF redaction."""
    for block in page.get_text("dict").get("blocks", []):
        if block.get("type") != 0:
            continue
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                if predicate(span.get("text", "")):
                    page.draw_rect(fitz.Rect(span["bbox"]), color=None,
                                   fill=(1, 1, 1), overlay=True)


def clean_s4(source: Path, output: Path):
    doc = fitz.open(source)
    page = doc[0]
    spans = []
    for block in page.get_text("dict").get("blocks", []):
        if block.get("type") != 0:
            continue
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                if "†" in span.get("text", ""):
                    spans.append(span)
    for span in spans:
        page.add_redact_annot(fitz.Rect(span["bbox"]), fill=(1, 1, 1), cross_out=False)
    if spans:
        page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=0, text=0)
        for span in spans:
            replacement = span["text"].replace("†", "")
            if replacement.strip():
                page.insert_text(
                    fitz.Point(span["origin"]), replacement,
                    fontname="helv", fontsize=float(span["size"]),
                    color=(0.12, 0.12, 0.12), overlay=True,
                )
    _remove_text(page, lambda text: "Historical eight" in text)
    output.write_bytes(doc.tobytes(garbage=4, deflate=True))
    doc.close()


def _table_s4_values(path: Path):
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb["Table S4"]
    headers = [ws.cell(1, c).value for c in range(1, ws.max_column + 1)]
    headers = [str(x) if x is not None else "" for x in headers]
    region_col = next(
        (i for i, h in enumerate(headers)
         if "nucleus screen" in h.lower() and "eta-squared" in h.lower()),
        None,
    )
    if region_col is None:
        region_col = next(i for i, h in enumerate(headers)
                          if "nucleus-level regional eta-squared" in h.lower())
    cohort_col = next(i for i, h in enumerate(headers)
                      if h.lower() == "cohort partial eta-squared, sensitivity")
    conf_col = headers.index("Program confidence")
    program_col = headers.index("Program")
    rows = []
    for r in range(2, ws.max_row + 1):
        program = ws.cell(r, program_col + 1).value
        if not program:
            continue
        rows.append({
            "program": str(program),
            "region_eta": float(ws.cell(r, region_col + 1).value),
            "cohort_eta": float(ws.cell(r, cohort_col + 1).value),
            "confidence": str(ws.cell(r, conf_col + 1).value),
        })
    return rows


def clean_workbook(path: Path, rows_for_plot=None):
    wb = load_workbook(path)
    ws = wb["Table S4"]
    headers = [ws.cell(1, c).value for c in range(1, ws.max_column + 1)]
    remove = {
        "Historical nucleus screen F",
        "Historical nucleus screen P",
        "Historical nucleus screen BH-FDR",
        "Historical nucleus screen decision",
        "Historical nucleus screen eight",
        "Historical OLS decision, not primary",
        "Historical permutation regional decision",
        "Historical cohort sensitivity description",
    }
    for c in range(ws.max_column, 0, -1):
        if ws.cell(1, c).value in remove:
            ws.delete_cols(c, 1)
    for c in range(1, ws.max_column + 1):
        if ws.cell(1, c).value == "Historical nucleus screen eta-squared":
            ws.cell(1, c).value = "Nucleus-level regional eta-squared (descriptive)"
    for name in ["S4 Human target8 evidence"]:
        if name in wb.sheetnames:
            del wb[name]
    if "S4 Human donor robustness" in wb.sheetnames:
        donor = wb["S4 Human donor robustness"]
        header_row = 3
        for c in range(donor.max_column, 0, -1):
            if str(donor.cell(header_row, c).value) == "historical_nucleus_screen_eight":
                donor.delete_cols(c, 1)
    wb.save(path)


def render_s5(rows, output: Path):
    rows = sorted(rows, key=lambda x: int(x["program"][1:]))
    programs = [x["program"] for x in rows]
    region = np.array([x["region_eta"] for x in rows], float)
    cohort = np.array([x["cohort_eta"] for x in rows], float)
    high = np.array([x["confidence"].lower().startswith("higher") for x in rows])
    colors = np.where(high, "#2E6E4E", "#B0883B")
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(7.3, 9.8), gridspec_kw={"height_ratios": (1.0, 1.2)})
    ax_a.scatter(region, cohort, c=colors, s=22, edgecolor="white", linewidth=0.35)
    for x, y, p in zip(region, cohort, programs):
        if y >= np.quantile(cohort, 0.75) or x >= np.quantile(region, 0.75):
            ax_a.text(x, y, p, fontsize=6, ha="left", va="bottom")
    ax_a.set_xlabel("Nucleus-level regional eta² (descriptive)")
    ax_a.set_ylabel("Cohort partial eta² (sensitivity)")
    ax_a.set_title("a  Regional and cohort contributions across all 54 programs", loc="left")
    ax_a.grid(True, color="#e6e8eb", linewidth=0.5)
    ax_a.legend(
        handles=[
            plt.Line2D([], [], marker="o", linestyle="", color="#2E6E4E", label="Higher annotation confidence"),
            plt.Line2D([], [], marker="o", linestyle="", color="#B0883B", label="Lower annotation confidence"),
        ], frameon=False, loc="upper left", fontsize=7,
    )
    order = np.argsort(cohort)[::-1]
    ax_b.bar(np.arange(len(rows)), cohort[order], color=colors[order], width=0.8)
    ax_b.set_xticks(np.arange(len(rows)))
    ax_b.set_xticklabels(np.array(programs)[order], rotation=90, fontsize=5.2)
    ax_b.set_ylabel("Cohort partial eta² (sensitivity)")
    ax_b.set_title("b  Cohort contribution ranking across all 54 programs", loc="left")
    ax_b.grid(axis="y", color="#e6e8eb", linewidth=0.5)
    fig.tight_layout(h_pad=2.0)
    fig.savefig(output, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _neutralize_old_categories(page):
    # The old cohort classes use two exact vector RGB triplets.  Replacing
    # those operators in reachable content streams preserves every path,
    # curve, label and coordinate while removing only the retired encoding.
    doc = page.parent
    patterns = [
        re.compile(rb"(?<![0-9.])(?:0)?\\.769\\d*\\s+(?:0)?\\.306\\d*\\s+(?:0)?\\.322\\d*"),
        re.compile(rb"(?<![0-9.])(?:0)?\\.886\\d*\\s+(?:0)?\\.635\\d*\\s+(?:0)?\\.173\\d*"),
    ]
    # Replace only the two exact retired cohort hues in every raw stream.
    # Enumerating streams is cheap; parsing nested vector drawings is not.
    for xref in range(1, doc.xref_length()):
        try:
            stream = doc.xref_stream(xref)
        except Exception:
            continue
        if not stream:
            continue
        updated = stream
        for pattern in patterns:
            updated = pattern.sub(b"0.549 0.549 0.549", updated)
        if updated != stream:
            doc.update_stream(xref, updated)
    # The retired cohort legend occupies this fixed native slot in panel a.
    # Mask it directly; extracting/redacting text would recursively traverse
    # the million-object nested figure source.
    page.draw_rect(fitz.Rect(914, 205, 1020, 292), color=None,
                   fill=(1, 1, 1), overlay=True)
    page.draw_rect(fitz.Rect(130, 728, 310, 800), color=None,
                   fill=(1, 1, 1), overlay=True)


def clean_s6(source: Path, output: Path):
    source_doc = fitz.open(source)
    source_page = source_doc[0]
    _neutralize_old_categories(source_page)
    width, height = source_page.rect.width, source_page.rect.height
    result = fitz.open()
    page = result.new_page(width=width, height=height)
    crops = [
        (fitz.Rect(0, 0, 1040, 470), fitz.Rect(0, 0, 1040, 470), None),
        (fitz.Rect(0, 480, 340, 820), fitz.Rect(0, 480, 340, 820), None),
        (fitz.Rect(260, 480, 660, 820), fitz.Rect(260, 480, 700, 820), None),
        (fitz.Rect(0, 840, 400, 1195), fitz.Rect(0, 840, 400, 1195), "d"),
        (fitz.Rect(400, 840, 1040, 1195), fitz.Rect(400, 840, 1040, 1195), "e"),
        (fitz.Rect(0, 1215, 340, 1510), fitz.Rect(0, 1215, 340, 1510), "f"),
        (fitz.Rect(340, 1215, 1040, 1510), fitz.Rect(340, 1215, 1040, 1510), "g"),
        (fitz.Rect(340, 1510, 1040, 1810), fitz.Rect(340, 1510, 1040, 1810), "h"),
    ]
    for clip, dest, label in crops:
        page.show_pdf_page(dest, source_doc, 0, clip=clip, keep_proportion=True, overlay=True)
        if label:
            page.draw_rect(fitz.Rect(dest.x0, dest.y0, dest.x0 + 28, dest.y0 + 30), color=None, fill=(1, 1, 1), overlay=True)
            page.insert_text(fitz.Point(dest.x0 + 15, dest.y0 + 20), label, fontname="helv", fontsize=12, color=(0.05, 0.05, 0.05), overlay=True)
    caption = (
        "Fig. S6. Regional program profiles and variability. a, Regional standardized activity with "
        "annotation and eta-squared descriptors. b, Program ranking by regional eta squared. c, Mean "
        "activity and cross-region coefficient of variation, with point size indicating regional eta squared. "
        "d, Regional PCA, with 47% and 19% of variation on PC1 and PC2. e, Regional profiles of P6, P1, "
        "P13, P4, P8 and P9. f, P13, P6 and P1 activity across nuclei in each region. g, Within-program "
        "regional z scores for 38 programs, with frames marking each region's largest positive value. h, "
        "Relative rankings of ten programs across five lobes. These profiles describe regional expression; "
        "the donor-adjusted 21-program analysis is specified in Fig. 2a,c and Supplementary Table S4."
    )
    page.insert_textbox(fitz.Rect(70, 1810, 970, 1995), caption, fontname="helv", fontsize=10, lineheight=1.25, color=(0.05, 0.05, 0.05), overlay=True)
    # Keep the native form resources intact; full garbage collection traverses
    # the million-object source PDF and is unnecessary for this retained crop.
    output.write_bytes(result.tobytes(garbage=0, deflate=True))
    result.close(); source_doc.close()


def render_png(pdf: Path, png: Path):
    if pdf.name == "FigS6.pdf":
        # Poppler handles the nested native vector display list without
        # materializing the full MuPDF display list in memory.
        subprocess.run([
            "pdftocairo", "-png", "-singlefile", "-r", "150",
            str(pdf), str(png.with_suffix("")),
        ], check=True)
        return
    doc = fitz.open(pdf)
    pix = doc[0].get_pixmap(dpi=300, alpha=False)
    pix.save(png)
    doc.close()


def patch_s6_preview_png(png: Path):
    """Apply the same retired-legend masks to the retained S6 preview."""
    from PIL import Image, ImageDraw
    image = Image.open(png).convert("RGB")
    draw = ImageDraw.Draw(image)
    sx, sy = image.width / 1040.0, image.height / 2062.1975

    def rect(x0, y0, x1, y1):
        draw.rectangle((round(x0 * sx), round(y0 * sy),
                        round(x1 * sx), round(y1 * sy)), fill=(255, 255, 255))

    rect(914, 205, 1020, 292)
    rect(130, 728, 310, 800)
    rect(660, 480, 700, 820)
    # Neutralize the two retired cohort hues within the b/c region only;
    # the primary vector PDF carries the exact source colors and geometry.
    px = image.load()
    x0, x1 = round(0 * sx), round(700 * sx)
    y0, y1 = round(480 * sy), round(820 * sy)
    for y in range(y0, min(y1 + 1, image.height)):
        for x in range(x0, min(x1 + 1, image.width)):
            r, g, b = px[x, y]
            if (r >= 170 and g <= 135 and b <= 150 and r - g >= 45) or \
               (r >= 185 and 105 <= g <= 195 and b <= 125 and r - b >= 70):
                px[x, y] = (140, 140, 140)
    image.save(png)


def combine_supplements():
    combined = fitz.open()
    for number in range(1, 21):
        source = fitz.open(FIG_DIR / f"FigS{number}.pdf")
        combined.insert_pdf(source)
        source.close()
    combined.save(SUPP_DIR / "Supplementary_Figures.pdf", garbage=0, deflate=True)
    combined.close()


def finalize_existing():
    """Finalize tables and attachment copies from already-rendered S4--S6."""
    clean_workbook(TABLE_S4)
    clean_workbook(TABLES)
    (ATTACH_DIR / "Additional_file_1_supplementary material_Tables_S1-S6.xlsx").write_bytes(TABLES.read_bytes())
    for n in (4, 5, 6):
        (ATTACH_DIR / f"Additional_file_{n+1}_supplementary material_FigS{n}.pdf").write_bytes((FIG_DIR / f"FigS{n}.pdf").read_bytes())
    patch_s6_preview_png(PNG_DIR / "FigS6.png")
    combine_supplements()


def main():
    rows = _table_s4_values(TABLE_S4)
    clean_s4(FIG_DIR / "FigS4.pdf", FIG_DIR / "FigS4.pdf")
    render_s5(rows, FIG_DIR / "FigS5.pdf")
    clean_s6(FIG_DIR / "FigS6.pdf", FIG_DIR / "FigS6.pdf")
    for n in (4, 5, 6):
        render_png(FIG_DIR / f"FigS{n}.pdf", PNG_DIR / f"FigS{n}.png")
    patch_s6_preview_png(PNG_DIR / "FigS6.png")
    clean_workbook(TABLE_S4)
    clean_workbook(TABLES)
    (ATTACH_DIR / "Additional_file_1_supplementary material_Tables_S1-S6.xlsx").write_bytes(TABLES.read_bytes())
    for n in (4, 5, 6):
        (ATTACH_DIR / f"Additional_file_{n+1}_supplementary material_FigS{n}.pdf").write_bytes((FIG_DIR / f"FigS{n}.pdf").read_bytes())
    combine_supplements()


if __name__ == "__main__":
    if "--finalize-existing" in sys.argv[1:]:
        finalize_existing()
    else:
        main()
