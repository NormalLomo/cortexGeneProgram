#!/usr/bin/env Rscript
# v18: rasterize/vectorize the composed SVG to PDF (vector) + PNG (preview).
suppressMessages(library(rsvg))
FIG <- "__PRIVATE_CANONICAL_ROOT__/results/xspecies_humanmap_v1/spatial_xspecies/figures/Fig_spatial_univ"
svg <- file.path(FIG, "Fig_spatial_univ_v18.svg")
rsvg::rsvg_pdf(svg, file.path(FIG, "Fig_spatial_univ_v18.pdf"))
rsvg::rsvg_png(svg, file.path(FIG, "Fig_spatial_univ_v18.png"), width = 2600)
cat("wrote v18 pdf + png\n")
