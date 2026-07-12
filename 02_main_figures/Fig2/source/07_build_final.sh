#!/usr/bin/env bash
# Build the Fig. 2 regional program matrix from an external analysis root.
set -euo pipefail

PROJ="${CORTEX_PROGRAM_CANONICAL_ROOT:-}"
while [ "$#" -gt 0 ]; do
  case "$1" in
    --canonical-root)
      PROJ="$2"
      shift 2
      ;;
    --help)
      printf 'Usage: %s [--canonical-root PATH]\n' "$0"
      printf 'Uses --canonical-root or CORTEX_PROGRAM_CANONICAL_ROOT.\n'
      exit 0
      ;;
    *)
      printf 'Unknown argument: %s\n' "$1" >&2
      exit 2
      ;;
  esac
done

if [ -z "$PROJ" ]; then
  printf '%s\n' 'Set --canonical-root or CORTEX_PROGRAM_CANONICAL_ROOT.' >&2
  exit 2
fi

export CORTEX_PROGRAM_CANONICAL_ROOT="$PROJ"
F2="$PROJ/scripts/fig2"
F3="$PROJ/scripts/fig3"
COMP="$PROJ/figures/fig2/_compose"
RSCRIPT="${RSCRIPT:-Rscript}"
LAYOUT='V[a,H[b,c,d,e],f,H[g,h,i]]'

if [ ! -d "$F2" ] || [ ! -d "$F3" ]; then
  printf '%s\n' 'The canonical root must contain scripts/fig2 and scripts/fig3.' >&2
  exit 2
fi

echo "[1/4] panels (R)"
"$RSCRIPT" "$F2/fig2_svg_panels.R"

echo "[2/4] ink-crop (py)"
python3 "$F2/fig2_ink_crop.py"

echo "[3/4] gutter-0 template (py)"
python3 "$F3/make_template_nested.py" "$COMP" --layout "$LAYOUT" --out "$COMP"

echo "[4/4] compose -> vector PDF (py)"
python3 "$F3/compose_svgutils_skill.py" "$COMP" "$COMP/layout_template.json" \
    --out-svg "$PROJ/figures/fig2/Figure_2.svg" \
    --out-pdf "$PROJ/figures/fig2/Figure_2.pdf" \
    --rscript "$RSCRIPT"

echo "[png] contact sheet (rsvg, 200dpi)"
"$RSCRIPT" -e "rsvg::rsvg_png('$PROJ/figures/fig2/Figure_2.svg', '$PROJ/figures/fig2/Figure_2.png', width = round(193/25.4*200))"

echo "DONE Figure_2"
