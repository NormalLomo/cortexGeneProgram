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
# program_labels.tsv). Convention matches panel b/j: name_short, '*' for the 24
# brain-weak programs (authority confidence == "brain-weak").
pn <- read.delim(file.path(RES, "program_names.tsv"), stringsAsFactors = FALSE,
                 quote = "", comment.char = "", check.names = FALSE)
# program_names.tsv has cols: new_P, cnmf_component, name_full, name_short, confidence, ...
# Use cnmf_component (old cNMF number, integer) as lookup key for panels that need it.
# Filter EXCLUDED rows (new_P == "EXCLUDED") before building maps.
pn <- pn[pn$new_P != "EXCLUDED", ]
pn$program <- as.integer(pn$cnmf_component)
pn$fdr     <- suppressWarnings(as.numeric(pn$fdr))
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
## PANEL a -- GRAPHICAL ABSTRACT (rich banner conveying the paper essence)

## === FOCUSED RE-RENDER: panels c + i only (DPCA (SLAT dual_pca) coords) ===
cat("[render] re-rendering panels c and i from new embedding\n")

## PANEL c -- UMAP dominant program + SUBCLASS CONTOURS + LABELS
## program coloring underneath; per-subclass density contour + centroid
## label on top (ggrepel) so cell-type identity is legible.
## =====================================================================
{
  df <- read.csv(file.path(INT, "umap_embedding.csv"), row.names = 1,
                 check.names = FALSE)
  df$dominant_program <- as.integer(df$dominant_program)
  df$subclass <- as.character(df$subclass)
  progs_u <- sort(unique(df$dominant_program))
  base_cols <- c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"),
                 brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"),
                 brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
  prog_pal <- setNames(rep(base_cols, length.out = length(progs_u)), progs_u)

  # subclass centroids for labels (medoid-ish: median of coords)
  sub_cent <- df %>% group_by(subclass) %>%
    summarise(x = median(UMAP1), y = median(UMAP2), n = n(), .groups = "drop") %>%
    filter(n >= 200)
  # drop tiny/degenerate subclasses from contour set to avoid noise
  keep_sub <- sub_cent$subclass
  dfc <- df[df$subclass %in% keep_sub, ]

  # PROGRAM centroids for "P{n}" labels (median UMAP over cells whose
  # dominant_program == n); renumbered to 54 new_P; 6 excluded old programs dropped.
  # OLD2NEW: old cNMF component -> new P number (EXCLUDED -> NA)
  EXCLUDED_OLD <- c(9L, 18L, 19L, 35L, 52L, 57L)
  OLD2NEW <- c(
    `1`=1,`2`=2,`3`=3,`4`=4,`5`=5,`6`=6,`7`=7,`8`=8,
    `10`=9,`11`=10,`12`=11,`13`=12,`14`=13,`15`=14,`16`=15,`17`=16,
    `20`=17,`21`=18,`22`=19,`23`=20,`24`=21,`25`=22,`26`=23,`27`=24,
    `28`=25,`29`=26,`30`=27,`31`=28,`32`=29,`33`=30,`34`=31,
    `36`=32,`37`=33,`38`=34,`39`=35,`40`=36,`41`=37,`42`=38,`43`=39,
    `44`=40,`45`=41,`46`=42,`47`=43,`48`=44,`49`=45,`50`=46,`51`=47,
    `53`=48,`54`=49,`55`=50,`56`=51,`58`=52,`59`=53,`60`=54
  )
  prog_cent <- df %>%
    filter(!dominant_program %in% EXCLUDED_OLD) %>%
    group_by(dominant_program) %>%
    summarise(x = median(UMAP1), y = median(UMAP2), .groups = "drop") %>%
    mutate(new_p = OLD2NEW[as.character(dominant_program)],
           plab  = paste0("P", new_p))

  x0 <- min(df$UMAP1); y0 <- min(df$UMAP2)
  p_c <- ggplot(df, aes(UMAP1, UMAP2)) +
    # program coloring underneath (rasterized point cloud)
    ggrastr::rasterise(
      geom_point(aes(colour = factor(dominant_program)),
                 size = 0.18, alpha = 0.55, show.legend = FALSE), dpi = 400) +
    scale_colour_manual(values = prog_pal, guide = "none") +
    # per-subclass density contours ON TOP (single mid-level contour line)
    stat_density_2d(data = dfc, aes(UMAP1, UMAP2, group = subclass),
                    colour = "grey15", linewidth = 0.22, alpha = 0.6,
                    breaks = NULL, bins = 4, contour = TRUE,
                    inherit.aes = FALSE) +
    # subclass text labels at centroids, repelled to avoid overlap
    ggrepel::geom_text_repel(
      data = sub_cent, aes(x, y, label = subclass), inherit.aes = FALSE,
      size = 1.7, fontface = "bold", colour = "black",
      bg.color = "white", bg.r = 0.12,
      max.overlaps = 50, segment.size = 0.18, segment.colour = "grey40",
      min.segment.length = 0, box.padding = 0.18, point.padding = 0,
      seed = 42) +
    # PROGRAM "P{n}" labels at program centroids, repelled (54 kept, new_P labels)
    ggrepel::geom_text_repel(
      data = prog_cent, aes(x, y, label = plab), inherit.aes = FALSE,
      size = 1.5, fontface = "plain", colour = "grey10",
      bg.color = "white", bg.r = 0.10,
      max.overlaps = Inf, segment.size = 0.12, segment.colour = "grey55",
      min.segment.length = 0, box.padding = 0.12, point.padding = 0,
      seed = 7) +
    # mini UMAP axis arrows (bottom-left)
    annotate("segment", x = x0, y = y0, xend = x0 + 4, yend = y0,
             arrow = arrow(length = unit(1.2, "mm")), colour = "#555", linewidth = 0.4) +
    annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + 4,
             arrow = arrow(length = unit(1.2, "mm")), colour = "#555", linewidth = 0.4) +
    annotate("text", x = x0 + 2, y = y0 - 1.6, label = "UMAP1", size = 1.6, colour = "#555") +
    annotate("text", x = x0 - 1.6, y = y0 + 2, label = "UMAP2", size = 1.6, colour = "#555",
             angle = 90) +
    labs(title = "Programs on UMAP") +
    coord_fixed(ratio = 1, clip = "off") +
    theme_void(base_family = "Helvetica") +
    theme(plot.title = element_text(size = 6.6, face = "bold", hjust = 0.5),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(0, 0, 0, 0))
  save_gg_svg(p_c, "c")
}

## =====================================================================

## PANEL i -- marker UMAP 2x2 (square, ggplot + ggrastr point layers)
## =====================================================================
{
  df  <- read.csv(file.path(INT, "umap_embedding.csv"), row.names = 1, check.names = FALSE)
  mk  <- read.csv(file.path(INT, "marker_on_umap_subsample.csv"), row.names = 1,
                  check.names = FALSE)
  common <- intersect(rownames(df), rownames(mk))
  df <- df[common, ]; mk <- mk[common, ]
  sel <- intersect(c("SLC17A7", "GAD1", "AQP4", "PLP1"), colnames(mk))
  ml <- do.call(rbind, lapply(sel, function(g) {
    v <- mk[[g]]; cap <- quantile(v, 0.995)
    data.frame(UMAP1 = df$UMAP1, UMAP2 = df$UMAP2,
               expr = pmin(v, cap), marker = g)
  }))
  ml$marker <- factor(ml$marker, levels = sel)
  ml <- ml[order(ml$marker, ml$expr), ]
  p_i <- ggplot(ml, aes(UMAP1, UMAP2, colour = expr)) +
    ggrastr::rasterise(geom_point(size = 0.16, alpha = 0.75), dpi = 400) +
    facet_wrap(~ marker, nrow = 2, ncol = 2) +
    scale_colour_gradientn(colours = rev(viridisLite::magma(256)), name = "expr") +
    coord_fixed(ratio = 1, clip = "off") +
    labs(title = "Canonical markers (UMAP)") +
    theme_void(base_family = "Helvetica") +
    theme(plot.title = element_text(size = 7.2, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 6.4, face = "bold.italic"),
          aspect.ratio = 1,
          panel.spacing = unit(2.2, "mm"),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.key.width = unit(6, "mm"), legend.key.height = unit(2, "mm"),
          legend.title = element_text(size = 5.5), legend.text = element_text(size = 5),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(0, 0, 0, 0))
  save_gg_svg(p_i, "i")
}

## =====================================================================

cat("[render] DONE panels c + i\n")
