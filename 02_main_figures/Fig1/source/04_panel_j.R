#!/usr/bin/env Rscript
# PARTIAL runner: re-render ONLY panel j (program co-activity network) with
# refreshed functional names from program_names.tsv and the 4 technical
# programs (4,18,35,52) dropped, to mirror panel b's 56-program set.
# Shared header copied verbatim from fig1_svg_panels.R (libs/BOX/helpers),
# EXCEPT panel-j node labels now come from program_names.tsv (not the stale
# program_labels.tsv), with `*` appended for fdr >= 0.25 (panel b convention).
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

BOX <- list(j = c(w = 99.077,  h = 61.923))   # co-activity network
mm2in <- function(x) x / 25.4

CLASS_COL <- c(excitatory = "#2166AC", inhibitory = "#B2182B",
               `non-neuronal` = "#1B7837")

# ---- FRESH functional names from program_names.tsv (replaces stale short_map)
pn <- read.delim(file.path(RES, "program_names.tsv"), quote = "", comment.char = "",
                 stringsAsFactors = FALSE, check.names = FALSE)
pn$program <- as.integer(pn$program)
pn$fdr     <- suppressWarnings(as.numeric(pn$fdr))
# panel b convention: weak programs (fdr >= 0.25) get a trailing '*'
pn$lab_short <- ifelse(!is.na(pn$fdr) & pn$fdr >= 0.25,
                       paste0(pn$name_short, "*"), pn$name_short)
short_map <- setNames(pn$lab_short, pn$program)
.lsh <- function(p) ifelse(as.character(p) %in% names(short_map),
                           short_map[as.character(p)], paste0("P", p))

scc <- read.csv(file.path(INT, "subclass_class.csv"), stringsAsFactors = FALSE)
class_map <- setNames(scc$class, scc$subclass)
spec <- read.csv(file.path(INT, "program_specificity.csv"))
spec$program <- as.integer(spec$program)

save_gg_svg <- function(p, panel, scale_in = 1.0) {
  b <- BOX[[panel]]
  w <- mm2in(b["w"]) * scale_in
  h <- mm2in(b["h"]) * scale_in
  f <- file.path(SVGD, sprintf("fig1_%s.svg", panel))
  svglite(f, width = as.numeric(w), height = as.numeric(h), bg = "white")
  print(p); dev.off()
  cat(sprintf("panel %s SVG -> %s (%.2f x %.2f in)\n", panel, f, w, h))
}

theme_nat <- function(base = 6.6) {
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(panel.grid = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
}

## =====================================================================
## PANEL j -- program co-activity network (ggraph; hubs = functional names)
## =====================================================================
DROP <- c(4L, 18L, 35L, 52L)   # technical-artifact programs (match panel b)
{
  cm <- read.csv(file.path(INT, "program_program_corr.csv"), check.names = FALSE,
                 row.names = 1)
  cm <- as.matrix(cm); diag(cm) <- 0
  progs_all <- as.integer(colnames(cm))
  # ---- drop the 4 technical programs (rows+cols) BEFORE thresholding so edges
  # to/from them vanish and the threshold reflects the 56-program set ----
  keep <- !(progs_all %in% DROP)
  cm <- cm[keep, keep, drop = FALSE]
  progs_j <- as.integer(colnames(cm))
  cat(sprintf("panel j: %d program nodes after dropping %s\n",
              length(progs_j), paste(DROP, collapse = ",")))
  # ---- CLEAN edge backbone (no hairball): keep each node's TOP-K strongest
  # co-activity links, gated by a minimum correlation. Union of per-node top
  # edges -> sparse, interpretable graph. Width/alpha scale with weight. ----
  nn   <- nrow(cm)
  KTOP <- 3L          # strongest links kept per node
  MINW <- 0.10        # minimum co-activity correlation to draw an edge
  ekeep <- matrix(FALSE, nn, nn)
  for (i in seq_len(nn)) {
    ord <- order(cm[i, ], decreasing = TRUE)
    ord <- ord[ord != i][seq_len(KTOP)]
    ord <- ord[cm[i, ord] > MINW]
    ekeep[i, ord] <- TRUE
  }
  ekeep <- ekeep | t(ekeep)            # symmetrise (top-of-either-endpoint)
  idx <- which(upper.tri(ekeep) & ekeep, arr.ind = TRUE)
  cat(sprintf("panel j: edge backbone top-%d per node, corr > %.2f -> %d edges\n",
              KTOP, MINW, nrow(idx)))
  # reference edges by NODE NAME (char program id) so tidygraph maps by name,
  # not row position -- required now that program ids are non-contiguous after
  # dropping the 4 technical programs (else add_vertices negative-count error).
  edges <- data.frame(from = as.character(progs_j[idx[, 1]]),
                      to   = as.character(progs_j[idx[, 2]]),
                      weight = cm[idx], stringsAsFactors = FALSE)
  nodes <- data.frame(program = progs_j) %>%
    left_join(spec[, c("program", "dominant_class", "gini")], by = "program")
  nodes$name <- as.character(nodes$program)
  g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE, node_key = "name") %>%
    activate(nodes) %>% mutate(deg = centrality_degree(weights = NULL))
  nd <- g %>% activate(nodes) %>% as_tibble()
  # ---- HUB labelling: only the highest-degree nodes get a functional name,
  # so the panel is not crowded with 56 labels. Top ~12 by degree (ties incl).
  N_HUB <- 12L
  deg_cut <- sort(nd$deg, decreasing = TRUE)[min(N_HUB, length(nd$deg))]
  is_hub  <- nd$deg >= deg_cut & nd$deg > 0
  cat(sprintf("panel j: labelling %d hub nodes (degree >= %d)\n",
              sum(is_hub), deg_cut))
  sl <- .lsh(nd$program)
  g <- g %>% activate(nodes) %>%
    mutate(hub_lab = ifelse(is_hub & nzchar(sl), sl, NA_character_),
           is_hub  = is_hub)
  set.seed(7)
  FF_AR <- as.numeric(BOX[["j"]]["w"] / BOX[["j"]]["h"])  # box width:height
  # Fruchterman-Reingold force-directed layout (good cluster separation).
  lay <- create_layout(g, layout = "fr", niter = 2000)
  # fit the layout to the panel box aspect (keep relative geometry, just scale)
  lay$x <- scales::rescale(lay$x, to = c(0, FF_AR))
  lay$y <- scales::rescale(lay$y, to = c(0, 1.0))
  p_j <- ggraph(lay) +
    geom_edge_link(aes(width = weight, alpha = weight), colour = "grey60") +
    geom_node_point(aes(fill = dominant_class, size = deg + 1),
                    shape = 21, colour = "white", stroke = 0.3) +
    ggrepel::geom_text_repel(
      data = subset(lay, !is.na(hub_lab)),
      aes(x = x, y = y, label = hub_lab),
      size = 1.95, fontface = "bold", colour = "grey10",
      bg.color = "white", bg.r = 0.12,
      box.padding = 0.18, point.padding = 0.10,
      min.segment.length = 0, segment.size = 0.2, segment.colour = "grey55",
      force = 0.4, force_pull = 1.2, max.overlaps = 20, max.iter = 4000,
      seed = 7, na.rm = TRUE) +
    scale_edge_width(range = c(0.12, 0.9), guide = "none") +
    scale_edge_alpha(range = c(0.12, 0.6), guide = "none") +
    scale_fill_manual(values = CLASS_COL, name = "dominant class") +
    scale_size_continuous(range = c(1.2, 5), name = "degree") +
    scale_x_continuous(expand = expansion(mult = 0.06)) +
    scale_y_continuous(expand = expansion(mult = 0.08)) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(title = "Co-activity network") +
    theme_void(base_size = 7, base_family = "Helvetica") +
    theme(plot.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 5.5), legend.title = element_text(size = 6),
          legend.key.size = unit(3.2, "mm"), legend.position = "right",
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(0, 0, 0, 0))
  save_gg_svg(p_j, "j")
}
cat("PANEL j RE-RENDER COMPLETE\n")
