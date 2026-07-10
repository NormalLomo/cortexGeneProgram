#!/usr/bin/env bash
# NAMEFIX 20260614 — re-render Fig.8 (figB_large) full compose chain after
# program-name source fix (program_names.tsv name_short authoritative).
# Chain: R emit per-panel SVGs -> ink-crop -> uniq ids -> nested template -> compose vector PDF.
set -euo pipefail

ROOT=CORTEX_PROGRAM_ROOT
SC=$ROOT/scripts/figmarkcorr_B
SVGD=$SC/svg_panels
OUT=$ROOT/figures/markcorr_B
RSCRIPT=/opt/R/4.5.0/bin/Rscript
PY=python

LAYOUT="V[t,H[a,b],H[c,d],H[e,f,g],H[h,i]]"
BAND_TITLE="Low-gene-overlap spatial co-organization of cortical gene programs"
BAND_SUB="*=weak functional annotation (brain-weak).  Co-org = spatial mark-correlation g(r=25 um), 44 chips; bin50 = 25 um (multi-cell, not single cells)."

echo "=== [1/5] Rscript figB_large.R (emit per-panel SVGs) ==="
$RSCRIPT $SC/figB_large.R

echo "=== [2/5] ink-crop ==="
$PY $SC/figB_ink_crop.py

echo "=== [3/5] uniq svg ids ==="
$PY $SC/uniq_svg_ids.py

echo "=== [4/5] make nested template (page-width 180, margin 4) ==="
$PY $SC/make_template_nested.py "$SVGD" \
    --layout "$LAYOUT" --page-width-mm 180 --margin-mm 4 --out "$SVGD"

echo "=== [5/5] compose -> vector PDF + SVG ==="
$PY $SC/compose_svgutils_skill.py "$SVGD" "$SVGD/layout_template.json" \
    --out-svg $OUT/figB_large.svg --out-pdf $OUT/figB_large.pdf \
    --skip-tag-ids t \
    --text-band-id t --band-title "$BAND_TITLE" --band-subtitle "$BAND_SUB" \
    --rscript $RSCRIPT

echo "=== render PNG from composite SVG ==="
$RSCRIPT -e "library(rsvg); rsvg_png('$OUT/figB_large.svg','$OUT/figB_large.png',width=2126)"

echo "=== DONE; fonts: ==="
pdffonts $OUT/figB_large.pdf || true
ls -la $OUT/figB_large.{svg,pdf,png}
