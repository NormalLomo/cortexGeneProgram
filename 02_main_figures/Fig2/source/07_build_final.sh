#!/usr/bin/env bash
# =====================================================================
# Re-render Fig.3 -> fig3_full_v12.{svg,pdf}  (panel-d shrink, OPTION B)
# Approach B (per 羅老師): shrink panel d by changing only its reference
# rendering (line-end label = capped short name, bottom 1-row lobe legend,
# smaller right x-expand, squarer canvas) so its cropped natural_aspect
# drops ~1.24 -> ~1.0-1.15.  The tangram LAYOUT-SPEC is UNCHANGED from v10:
#   V[a,H[b,c,d,e],f,H[g,h,i]]   (d stays in row-2 original slot; a-i order kept)
# With d's aspect lowered, the packer naturally gives d a normal width and
# b/c/e recover normal widths too.
# NO data/science/program-name/robust-7-highlight change; only panel d render.
# Does NOT overwrite v10 or v11.
# =====================================================================
set -euo pipefail
PROJ=__PRIVATE_CANONICAL_ROOT__
F3=$PROJ/scripts/fig3
COMP=$PROJ/figures/fig3/_compose
RSCRIPT=/opt/R/4.5.0/bin/Rscript
LAYOUT='V[a,H[b,c,d,e],f,H[g,h,i]]'        # SAME original layout as v10

echo "[1/4] panels (R)"
"$RSCRIPT" "$F3/fig3_svg_panels.R"

echo "[2/4] ink-crop (py)"
python3 "$F3/fig3_ink_crop.py"

echo "[3/4] gutter-0 template (py)"
python3 "$F3/make_template_nested.py" "$COMP" --layout "$LAYOUT" --out "$COMP"

echo "[4/4] compose -> vector PDF (py)"
python3 "$F3/compose_svgutils_skill.py" "$COMP" "$COMP/layout_template.json" \
    --out-svg "$PROJ/figures/fig3/fig3_full_v12.svg" \
    --out-pdf "$PROJ/figures/fig3/fig3_full_v12.pdf" \
    --rscript "$RSCRIPT"

echo "[png] contact sheet (rsvg, 200dpi)"
"$RSCRIPT" -e "rsvg::rsvg_png('$PROJ/figures/fig3/fig3_full_v12.svg', '$PROJ/figures/fig3/fig3_full_v12.png', width = round(193/25.4*200))"

echo "DONE fig3_full_v12"
