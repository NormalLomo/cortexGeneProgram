#!/usr/bin/env Rscript
# Render the composed SVG to vector PDF and PNG preview.
suppressMessages(library(rsvg))
args <- commandArgs(trailingOnly=TRUE)
root_idx <- match("--canonical-root", args)
canonical_root <- if (!is.na(root_idx) && length(args) > root_idx) args[[root_idx + 1]] else Sys.getenv("CORTEX_PROGRAM_CANONICAL_ROOT")
if (!nzchar(canonical_root)) stop("Set --canonical-root or CORTEX_PROGRAM_CANONICAL_ROOT.")
FIG <- file.path(canonical_root, "results/xspecies_humanmap_v1/spatial_xspecies/figures/Fig_spatial_univ")
svg <- file.path(FIG, "Figure_7.svg")
rsvg::rsvg_pdf(svg, file.path(FIG, "Figure_7.pdf"))
rsvg::rsvg_png(svg, file.path(FIG, "Figure_7.png"), width = 2600)
cat("wrote pdf + png\n")
