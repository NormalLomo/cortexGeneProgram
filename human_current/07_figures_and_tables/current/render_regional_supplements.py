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
TABLE_S3 = ROOT / "tables" / "TableS3_program_variability_validity_human.xlsx"
TABLES = ROOT / "tables" / "Supplementary_Tables_S1-S5_human.xlsx"
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


def _table_s3_values(path: Path):
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb["Table S3"]
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
    ws = wb["Table S3"]
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
    for name in ["S3 Human target8 evidence"]:
        if name in wb.sheetnames:
            del wb[name]
    if "S3 Human donor robustness" in wb.sheetnames:
        donor = wb["S3 Human donor robustness"]
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
        re.compile(rb"(?<![0-9.])(?:0)?\.769\d*\s+(?:0)?\.306\d*\s+(?:0)?\.322\d*"),
        re.compile(rb"(?<![0-9.])(?:0)?\.886\d*\s+(?:0)?\.635\d*\s+(?:0)?\.173\d*"),
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
        # Remove the retired cohort-robustness legend as native PDF operators
        # (not a white preview mask).  The regional lobe and z-score legends
        # remain part of the retained descriptive panel.
        start = updated.find(
            b"0 0 0 rg\nBT\n4.876 0 0 4.876 492.5195 576.1533 Tm\n(cohort robustness) Tj"
        )
        for marker in (
            b"0.549 0.549 0.549 rg\n492.52 571.398",
            b"0.698 0.094 0.169 rg\n492.52 571.398",
        ):
            candidate = updated.rfind(marker, 0, start if start >= 0 else len(updated))
            if candidate >= 0:
                start = candidate
                break
        end = updated.find(
            b"0 0 0 rg\nBT\n4.876 0 0 4.876 492.5195 531.1846 Tm\n(z\\055score) Tj"
        )
        if start >= 0 and end > start:
            updated = updated[:start] + updated[end:]
        # The same retired cohort key is repeated beneath the ranking panel;
        # remove its native text and marker operators from that compact source
        # region while retaining the ranking labels and axes.
        start_b = updated.find(
            b"3.418 0 0 3.418 60.377 335.1646 Tm\n(c) Tj"
        )
        end_b = updated.find(
            b"4.1894 0 0 4.1894 13.9141 442.3643 Tm",
            start_b if start_b >= 0 else 0,
        )
        if start_b >= 0 and end_b > start_b:
            updated = updated[:start_b] + updated[end_b:]
        # The retired old-d panel occupies one contiguous native clip region
        # in the shared compact source form.  Remove that operator span before
        # any retained panel form is imported; this avoids carrying the old
        # gradient/classification objects behind a crop boundary.
        start_old_d = updated.find(
            b"q\n0 661 547 -661 re\nW\nn\nq\n1 0 0 1 240.8359"
        )
        end_old_d = updated.find(
            b"q\n0 661 547 -661 re\nW\nn\nq\n1 0 0 1 397.7383",
            start_old_d if start_old_d >= 0 else 0,
        )
        if start_old_d >= 0 and end_old_d > start_old_d:
            updated = updated[:start_old_d] + updated[end_old_d:]
        if updated != stream:
            doc.update_stream(xref, updated)


def _build_clean_source_page(doc: fitz.Document) -> fitz.Page:
    """Expose the retained source panels without the retired old-d form.

    The current source PDF stores the original nine-panel composition as one
    compact form (``fzFrm0``--``fzFrm8``).  Build a same-document helper page
    that calls a new composition stream containing a--c and e--i only.  The
    omitted ``fzFrm3`` is therefore absent from the retained page content,
    rather than covered by a later white rectangle.
    """
    layout_xref = None
    layout_stream = None
    for xref in range(1, doc.xref_length()):
        try:
            stream = doc.xref_stream(xref)
        except Exception:
            continue
        if stream and b"/fzFrm8 Do" in stream and b"/fzFrm0 Do" in stream:
            layout_xref = xref
            layout_stream = stream
            break
    if layout_xref is None or layout_stream is None:
        raise RuntimeError("current S6 source composition form was not found")
    rtype, rval = doc.xref_get_key(layout_xref, "Resources")
    if rtype != "xref":
        raise RuntimeError("current S6 source composition has no resource dictionary")
    layout_res = int(re.search(r"\d+", rval).group(0))
    resource_obj = doc.xref_object(layout_res)
    panel_refs = {}
    for n in (0, 1, 2, 4, 5, 6, 7, 8):
        match = re.search(rf"/fzFrm{n}\s+(\d+)\s+0\s+R", resource_obj)
        if not match:
            raise RuntimeError(f"current S6 panel form fzFrm{n} is missing")
        panel_refs[n] = int(match.group(1))

    # Keep the original vector panel forms and their native transforms.  The
    # fourth source panel (old d) is intentionally not called.
    calls = []
    for n in (0, 1, 2, 4, 5, 6, 7, 8):
        calls.append(f"q /fzFrm{n} Do Q")
    composition = ("\n".join(calls) + "\n").encode()

    # Copy the source resource dictionary but remove the retired form key so
    # it is not reachable from the helper page's content graph.
    clean_resource = re.sub(r"\n\s*/fzFrm3\s+\d+\s+0\s+R", "", resource_obj)
    clean_res_xref = doc.get_new_xref()
    doc.update_object(clean_res_xref, clean_resource)
    clean_form_xref = doc.get_new_xref()
    doc.update_object(
        clean_form_xref,
        "<< /Type /XObject /Subtype /Form /BBox [0 0 1040 1810] "
        f"/Matrix [1 0 0 1 0 252.1975] /Resources {clean_res_xref} 0 R >>",
    )
    doc.update_stream(clean_form_xref, composition)

    helper = doc.new_page(width=1040, height=2062.1976)
    helper_res_xref = doc.get_new_xref()
    doc.update_object(helper_res_xref, f"<< /XObject << /base {clean_form_xref} 0 R >> >>")
    helper_contents_xref = doc.get_new_xref()
    doc.update_object(helper_contents_xref, "<< >>")
    doc.update_stream(helper_contents_xref, b"q /base Do Q")
    doc.xref_set_key(helper.xref, "Resources", f"{helper_res_xref} 0 R")
    doc.xref_set_key(helper.xref, "Contents", f"{helper_contents_xref} 0 R")
    return helper


def clean_s6(source: Path, output: Path):
    source_doc = fitz.open(source)
    _neutralize_old_categories(source_doc[0])
    source_page = _build_clean_source_page(source_doc)
    source_page_index = source_doc.page_count - 1
    width, height = 1040, 2062.1976
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
        page.show_pdf_page(dest, source_doc, source_page_index, clip=clip,
                           keep_proportion=True, overlay=True)
        if label:
            page.draw_rect(fitz.Rect(dest.x0, dest.y0, dest.x0 + 28, dest.y0 + 30), color=None, fill=(1, 1, 1), overlay=True)
            page.insert_text(fitz.Point(dest.x0 + 15, dest.y0 + 20), label, fontname="helv", fontsize=12, color=(0.05, 0.05, 0.05), overlay=True)
    page.insert_text(fitz.Point(18, 30), "Figure S6", fontname="hebo", fontsize=12,
                     color=(0.08, 0.12, 0.19), overlay=True)
    caption = (
        "Fig. S6. Regional program profiles and variability. a, Regional standardized activity with "
        "annotation and eta-squared descriptors. b, Program ranking by regional eta squared. c, Mean "
        "activity and cross-region coefficient of variation, with point size indicating regional eta squared. "
        "d, Regional PCA, with 47% and 19% of variation on PC1 and PC2. e, Regional profiles of P6, P1, "
        "P13, P4, P8 and P9. f, P13, P6 and P1 activity across nuclei in each region. g, Within-program "
        "regional z scores for 38 programs, with frames marking each region's largest positive value. h, "
        "Relative rankings of ten programs across five lobes. These profiles describe regional expression; "
        "the donor-adjusted 21-program analysis is specified in Fig. 2a,c and Supplementary Table S3."
    )
    page.insert_textbox(fitz.Rect(70, 1810, 970, 1995), caption, fontname="helv", fontsize=10, lineheight=1.25, color=(0.05, 0.05, 0.05), overlay=True)
    # Keep the native form resources intact; full garbage collection traverses
    # the million-object source PDF and is unnecessary for this retained crop.
    output.write_bytes(result.tobytes(garbage=0, deflate=True))
    result.close(); source_doc.close()


def render_png(pdf: Path, png: Path):
    if pdf.name == "FigS6.pdf":
        # Quartz/PDFKit draws this nested native page without the legacy
        # Poppler/MuPDF display-list stall.  The retained PNG is therefore a
        # direct raster companion of the exact S6 PDF just written above.
        swift_source = r'''
import Foundation
import AppKit
import PDFKit

let pdfURL = URL(fileURLWithPath: CommandLine.arguments[1])
let pngURL = URL(fileURLWithPath: CommandLine.arguments[2])
guard let document = PDFDocument(url: pdfURL), let page = document.page(at: 0) else {
    exit(2)
}
let box = page.bounds(for: .mediaBox)
let scale: CGFloat = 150.0 / 72.0
let size = NSSize(width: box.width * scale, height: box.height * scale)
let image = NSImage(size: size)
image.lockFocus()
NSColor.white.setFill()
NSRect(origin: .zero, size: size).fill()
if let context = NSGraphicsContext.current?.cgContext {
    context.saveGState()
    context.scaleBy(x: scale, y: scale)
    page.draw(with: .mediaBox, to: context)
    context.restoreGState()
}
image.unlockFocus()
guard let tiff = image.tiffRepresentation,
      let bitmap = NSBitmapImageRep(data: tiff),
      let png = bitmap.representation(using: .png, properties: [:]) else {
    exit(3)
}
try png.write(to: pngURL)
'''
        subprocess.run(["swift", "-", str(pdf), str(png)], input=swift_source.encode(), check=True)
        return
    doc = fitz.open(pdf)
    pix = doc[0].get_pixmap(dpi=300, alpha=False)
    pix.save(png)
    doc.close()


def repair_s6_font_resources(pdf: Path):
    """Attach standard Helvetica aliases to reachable S6 form resources.

    The retained nested page contains a few legacy form streams that refer to
    ``/hebo`` and ``/helv`` without carrying those aliases in their local
    resource dictionaries.  Adding the aliases to the existing font operators
    repairs reader rendering without touching paths, images, curves or text.
    """
    doc = fitz.open(pdf)
    bold = doc.get_new_xref()
    regular = doc.get_new_xref()
    doc.update_object(bold, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>")
    doc.update_object(regular, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>")
    page = doc[0]
    resource_refs = set()
    rtype, rval = doc.xref_get_key(page.xref, "Resources")
    if rtype == "xref":
        resource_refs.add(int(re.search(r"\d+", rval).group(0)))
    for xref in range(doc.xref_length()):
        try:
            stream = doc.xref_stream(xref)
        except Exception:
            continue
        if not stream or (b"/hebo" not in stream and b"/helv" not in stream):
            continue
        rtype, rval = doc.xref_get_key(xref, "Resources")
        if rtype == "xref":
            resource_refs.add(int(re.search(r"\d+", rval).group(0)))
        else:
            resource_refs.update(resource_refs)
    for resource in resource_refs:
        doc.xref_set_key(resource, "Font/hebo", f"{bold} 0 R")
        doc.xref_set_key(resource, "Font/helv", f"{regular} 0 R")
    payload = doc.tobytes(garbage=0, deflate=True)
    doc.close()
    pdf.write_bytes(payload)


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
    clean_workbook(TABLE_S3)
    clean_workbook(TABLES)
    (ATTACH_DIR / "Additional_file_1_supplementary material_Tables_S1-S5.xlsx").write_bytes(TABLES.read_bytes())
    repair_s6_font_resources(FIG_DIR / "FigS6.pdf")
    render_png(FIG_DIR / "FigS6.pdf", PNG_DIR / "FigS6.png")
    for n in (4, 5, 6):
        (ATTACH_DIR / f"Additional_file_{n+1}_supplementary material_FigS{n}.pdf").write_bytes((FIG_DIR / f"FigS{n}.pdf").read_bytes())
    # The retained S6 PNG is already the designated companion to the current
    # native S6 PDF.  Re-rasterizing its 72k-object display list through
    # Poppler is unnecessary and can stall on the legacy embedded font tags.
    # Reuse the retained pair while finalizing tables and attachments.
    combine_supplements()


def restore_retained_s6_distribution():
    """Reveal the existing lower P1 plot without rescaling any native object."""
    source = FIG_DIR / "FigS6.pdf"
    doc = fitz.open(source)
    page = doc[0]
    target = None
    for xref, name, parent, _ in page.get_xobjects():
        if name != "fzFrm12":
            continue
        kind, value = doc.xref_get_key(xref, "BBox")
        if kind != "array":
            continue
        box = [float(v) for v in re.findall(r"-?\d+(?:\.\d+)?", value)]
        if (len(box) == 4 and abs(box[0]) < .01
                and abs(box[2] - 340) < .01
                and abs(box[3] - (page.rect.height - 1215)) < .02):
            target = xref
            break
    if target is None:
        raise RuntimeError("The retained S6f native crop could not be located.")
    # P13/P6 occupy the upper row; P1 already exists directly beneath P13.
    # Only the lower clipping boundary hid P1. The right-hand g/h clips stay put.
    doc.xref_set_key(target, "BBox",
                     f"[0 {page.rect.height - 1690:.6f} 340 {page.rect.height - 1215:.6f}]")
    payload = doc.tobytes(garbage=0, clean=False, deflate=True, no_new_id=True)
    doc.close()
    source.write_bytes(payload)
    print("Produced", source, "retained P1 crop restored; fonts unchanged", flush=True)
    render_png(source, PNG_DIR / "FigS6.png")
    print("Produced", PNG_DIR / "FigS6.png", flush=True)
    attachment = ATTACH_DIR / "Additional_file_7_supplementary material_FigS6.pdf"
    attachment.write_bytes(payload)
    print("Produced", attachment, flush=True)


def main():
    rows = _table_s3_values(TABLE_S3)
    clean_s4(FIG_DIR / "FigS4.pdf", FIG_DIR / "FigS4.pdf")
    render_s5(rows, FIG_DIR / "FigS5.pdf")
    clean_s6(FIG_DIR / "FigS6.pdf", FIG_DIR / "FigS6.pdf")
    for n in (4, 5):
        render_png(FIG_DIR / f"FigS{n}.pdf", PNG_DIR / f"FigS{n}.png")
    # Reuse the retained S6 PNG companion; its native PDF display list is too
    # large for a redundant Poppler raster pass during package finalization.
    clean_workbook(TABLE_S3)
    clean_workbook(TABLES)
    (ATTACH_DIR / "Additional_file_1_supplementary material_Tables_S1-S5.xlsx").write_bytes(TABLES.read_bytes())
    for n in (4, 5, 6):
        (ATTACH_DIR / f"Additional_file_{n+1}_supplementary material_FigS{n}.pdf").write_bytes((FIG_DIR / f"FigS{n}.pdf").read_bytes())
    combine_supplements()


if __name__ == "__main__":
    if "--s6-layout-only" in sys.argv[1:]:
        restore_retained_s6_distribution()
    elif "--finalize-existing" in sys.argv[1:]:
        finalize_existing()
    else:
        main()
