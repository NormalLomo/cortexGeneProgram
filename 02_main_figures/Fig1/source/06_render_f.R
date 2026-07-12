#!/usr/bin/env Rscript
# =====================================================================
# Fig.1 -- render EACH of 10 panels (a-j) to a STANDALONE VECTOR SVG.
# Panel mapping (owner spec):
#   a = graphical abstract banner       (wide)
#   b = program x subclass heatmap      (ComplexHeatmap, vector)
#   c = UMAP dominant program + subclass contours/labels (ggplot)
#   d = top-5 gene loadings dotplot     (ggplot)
#   e = program functional-identity lollipop (ggplot)
#   f = subclass program-sharing chord  (circlize -> svglite, vector)
#   g = specificity ridgeline           (ggridges)
#   h = marker x subclass heatmap       (ComplexHeatmap, vector)
#   i = marker UMAP 2x2                  (ggplot + ggrastr)
#   j = program co-activity network     (ggraph)
# Dense point layers (c, i) rasterized IN-PLACE via ggrastr so the SVG
# stays a true vector wrapper (text/axes/legends vector). Heatmaps drawn
# directly to svglite (vector). Chord embeds a hi-dpi PNG (base graphics).
# Each SVG has ZERO plot.margin; x-y panels use coord_fixed(ratio=1).
# Ink-crop (fig1_ink_crop.py) trims to true content downstream.
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(ComplexHeatmap)
  library(circlize); library(ggrepel); library(ggridges)
  library(ggraph); library(igraph); library(tidygraph)
  library(dplyr); library(tidyr); library(grid)
  library(scales); library(RColorBrewer)
  library(svglite); library(ggrastr); library(png)
})
set.seed(42)

PROJ <- "CORTEX_PROGRAM_ROOT"
RES  <- file.path(PROJ, "results/crossregion_v1")
INT  <- file.path(PROJ, "figures/fig1/_intermediate")
OUT  <- file.path(PROJ, "figures/fig1")
SVGD <- file.path(OUT, "svg_panels")
dir.create(SVGD, recursive = TRUE, showWarnings = FALSE)

# ---- per-panel device size (mm). Sized to each panel's natural content
# aspect; ink-crop trims to true drawn extent so exact h is non-critical. -
BOX <- list(
  a = c(w = 173.0,   h = 40.0),    # wide graphical-abstract banner
  b = c(w = 105.861, h = 102.553), # program x subclass heatmap
  c = c(w = 62.0,    h = 62.0),    # UMAP (square, coord_fixed)
  d = c(w = 220.0,   h = 52.0),    # Manhattan top-gene loadings (FLAT banner, AR~4.2)
  e = c(w = 64.0,    h = 92.0),    # functional lollipop (60 rows -> tall)
  f = c(w = 90.0,    h = 90.0),    # subclass program-sharing chord (square)
  g = c(w = 78.0,    h = 46.0),    # ridgeline (2 facets, wide-ish)
  h = c(w = 55.139,  h = 42.414),  # marker x subclass heatmap
  i = c(w = 61.923,  h = 61.923),  # marker UMAP 2x2 (square)
  j = c(w = 99.077,  h = 61.923)   # co-activity network
)
mm2in <- function(x) x / 25.4

# ---- palettes -------------------------------------------------------
CLASS_COL <- c(excitatory = "#3B6FB6", inhibitory = "#C0392B",
               `non-neuronal` = "#27865E")
CLASS_LV  <- c("excitatory", "inhibitory", "non-neuronal")

theme_nat <- function(base = 6.6) {
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      axis.line   = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.35, colour = "black"),
      axis.title  = element_text(size = base),
      axis.text   = element_text(size = base - 0.6, colour = "black"),
      legend.title= element_text(size = base - 0.4),
      legend.text = element_text(size = base - 1.0),
      legend.key.size = unit(3.2, "mm"),
      strip.text  = element_text(size = base - 0.2, face = "bold"),
      strip.background = element_blank(),
      plot.title  = element_text(size = base + 1, face = "bold"),
      plot.subtitle = element_text(size = base - 0.6, colour = "grey30"),
      panel.grid  = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}
theme_set(theme_nat())

# FRESH functional names from program_names.tsv (authoritative; replaces stale
# program_labels.tsv). Convention matches panel b/j: name_short, '*' for brain-weak
# (confidence column == "brain-weak" = the fixed 24-program whitelist; not an fdr cutoff).
pn <- read.delim(file.path(RES, "program_names.tsv"), stringsAsFactors = FALSE,
                 quote = "", comment.char = "", check.names = FALSE)
pn$program <- as.integer(pn$program)
pn$fdr     <- suppressWarnings(as.numeric(pn$fdr))
# star = brain-weak (fixed 24-program whitelist via the authoritative confidence
# column), NOT a dynamic fdr cutoff. brain-weak == 24 programs; brain-sig == 36.
pn$lab_short <- ifelse(pn$confidence == "brain-weak",
                       paste0(pn$name_short, "*"), pn$name_short)
lab_map   <- setNames(pn$name_full,  pn$program)
short_map <- setNames(pn$lab_short,  pn$program)
.lab  <- function(p) ifelse(as.character(p) %in% names(lab_map),
                            lab_map[as.character(p)], paste0("P", p))
.lsh  <- function(p) ifelse(as.character(p) %in% names(short_map),
                            short_map[as.character(p)], paste0("P", p))

scc <- read.csv(file.path(INT, "subclass_class.csv"), stringsAsFactors = FALSE)
class_map <- setNames(scc$class, scc$subclass)
spec <- read.csv(file.path(INT, "program_specificity.csv"))
spec$program <- as.integer(spec$program)

# helper: save a ggplot to a vector svg sized to a box aspect
save_gg_svg <- function(p, panel, scale_in = 1.0) {
  b <- BOX[[panel]]
  w <- mm2in(b["w"]) * scale_in
  h <- mm2in(b["h"]) * scale_in
  f <- file.path(SVGD, sprintf("fig1_%s.svg", panel))
  svglite(f, width = as.numeric(w), height = as.numeric(h), bg = "white")
  print(p); dev.off()
  cat(sprintf("panel %s SVG -> %s (%.2f x %.2f in)\n", panel, f, w, h))
}


## =====================================================================
## PANEL f -- SUBCLASS program-sharing CHORD (circlize -> svglite, vector).
## 22 subclasses arranged around the circle, GROUPED by major class
## (excitatory -> inhibitory -> non-neuronal). Ribbons = how strongly two
## subclasses SHARE programs. Sharing metric S[i,j]: for each of the 56
## (technical-dropped) programs, row-normalise its 22 subclass loadings to
## weights w (sum 1) and accumulate S += outer(w, w); the diagonal is zeroed.
## So S[i,j] = sum_p w_p[i]*w_p[j] = co-loading strength of subclasses i,j
## across programs. Only the strongest links (>= upper-quintile of off-diag
## S) are drawn (thin, transparent) so it is not a hairball; ribbon colour =
## source sector's class. 22 labels are legible (vs the old 56-program ring).
## =====================================================================
{
  m <- read.csv(file.path(INT, "program_x_subclass_mean.csv"), check.names = FALSE)
  progs_all <- as.integer(m$program)
  mat <- as.matrix(m[, -1]); rownames(mat) <- progs_all
  # drop technical programs (consistent with panel b + network) -> 56 programs
  tech <- c(4, 18, 35, 52)
  mat  <- mat[!(progs_all %in% tech), , drop = FALSE]
  n_prog_f <- nrow(mat)
  frac <- mat / rowSums(mat)                       # row-normalised loadings (w)

  subs <- colnames(frac)                           # 22 subclasses
  # --- SUBCLASS x SUBCLASS sharing matrix S = sum_p outer(w_p, w_p) ----
  S <- matrix(0, length(subs), length(subs), dimnames = list(subs, subs))
  for (i in seq_len(nrow(frac))) { w <- frac[i, ]; S <- S + outer(w, w) }
  diag(S) <- 0                                     # no self-loops

  # --- order sectors by major class: excitatory -> inhibitory -> non-neur -
  sub_cls <- class_map[subs]
  ord     <- order(match(sub_cls, CLASS_LV), subs)
  subs_o  <- subs[ord]
  cls_o   <- sub_cls[ord]
  S <- S[subs_o, subs_o]

  # --- THRESHOLD: keep only strong sharing links (>= 80th pctile off-diag) -
  offdiag <- S[upper.tri(S)]
  thr <- as.numeric(quantile(offdiag, 0.80))
  Sthr <- S; Sthr[Sthr < thr] <- 0

  # grid (sector) colours by class; ribbon colour by SOURCE sector's class
  grid.col <- setNames(unname(CLASS_COL[cls_o]), subs_o)

  # build a SOURCE-class colour matrix for the links (alpha for legibility)
  col_mat <- matrix(NA_character_, nrow(Sthr), ncol(Sthr),
                    dimnames = dimnames(Sthr))
  for (i in seq_len(nrow(Sthr))) {
    cc <- grDevices::adjustcolor(unname(CLASS_COL[cls_o[i]]), alpha.f = 0.45)
    col_mat[i, ] <- cc
  }
  col_mat[Sthr == 0] <- "#00000000"

  bf <- BOX[["f"]]
  fsvg <- file.path(SVGD, "fig1_f.svg")
  svglite(fsvg, width = mm2in(bf["w"]), height = mm2in(bf["h"]), bg = "white")
  par(mar = c(0, 0, 0, 0), xpd = NA)
  circlize::circos.clear()
  # gaps: small within a class, larger between the 3 class blocks. Named by
  # sector so circlize never tries to "reduce" (drop) tiny sectors.
  cls_run <- rle(as.character(cls_o))$lengths
  gvec <- unlist(lapply(seq_along(cls_run), function(k)
    c(rep(1.6, cls_run[k] - 1), 7)))
  names(gvec) <- subs_o
  # expand canvas so radial labels + title + legend are NOT clipped
  circlize::circos.par(gap.after = gvec, start.degree = 90,
                       points.overflow.warning = FALSE,
                       cell.padding = c(0, 0, 0, 0),
                       canvas.xlim = c(-1.45, 1.45),
                       canvas.ylim = c(-1.45, 1.45))
  circlize::chordDiagram(
    Sthr, order = subs_o, grid.col = grid.col, col = col_mat,
    transparency = 0.5, directional = 0, symmetric = TRUE,
    reduce = 0,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.18),
    link.lwd = 0.3, link.border = NA, scale = FALSE)

  # sector LABELS just outside the grid track, radial, >=5pt -----------------
  circlize::circos.trackPlotRegion(
    track.index = 1, bg.border = NA, panel.fun = function(x, y) {
      s   <- circlize::get.cell.meta.data("sector.index")
      xlm <- circlize::get.cell.meta.data("xlim")
      circlize::circos.text(mean(xlm), circlize::get.cell.meta.data("ylim")[2] + 0.4,
        labels = s, facing = "clockwise", niceFacing = TRUE,
        adj = c(0, 0.5), cex = 0.46, col = "black")    # cex 0.46 ~= 5.3pt
    })
  circlize::circos.clear()

  # title + class legend (in expanded -1.45..1.45 canvas space) -------------
  text(0, 1.40, "Program sharing", cex = 0.80, font = 2, adj = c(0.5, 1))
  ly <- -1.36; lx <- -0.66
  for (k in seq_along(CLASS_LV)) {
    points(lx + (k - 1) * 0.66, ly, pch = 22, cex = 0.9,
           bg = unname(CLASS_COL[CLASS_LV[k]]), col = "grey30", lwd = 0.4)
    text(lx + (k - 1) * 0.66 + 0.05, ly, CLASS_LV[k], adj = c(0, 0.5), cex = 0.46)
  }
  dev.off()
  n_links <- sum(Sthr[upper.tri(Sthr)] > 0)
  cat(sprintf("panel f SVG -> %s (subclass chord: 22 sectors, %d programs, %d links, thr=%.4f)\n",
              fsvg, n_prog_f, n_links, thr))
}
