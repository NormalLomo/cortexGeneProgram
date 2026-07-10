#!/usr/bin/env Rscript
# =====================================================================
# Fig.4  Cellular system x region driver of cross-region program diffs
# (human cortex, cNMF K=60).  Independent ZERO-margin vector SVG panels
# (a-h) for svgutils template-first composition (mirrors fig3_svg_panels.R).
#
# Reuses the 8 build_* logics from scripts/fig4/fig4_full.R, but each panel
# is rendered to its OWN standalone svglite SVG with:
#   * plot.margin = 0   (ink hugs content; ink_crop trims the rest)
#   * heatmaps a,d = SQUARE CELLS (width=ncol*unit(s), height=nrow*unit(s))
#   * continuous x-y panels (f UMAP, h network) coord_fixed / void+balanced
#   * panel g legibility fix: facet strip + y labels >=5pt, top-8 programs
#   * panel f legibility fix: ggplot twin (vector), bigger titles/legend,
#       no tiny arrow annotations (the old f had unreadable micro-arrows)
# + ggraph).  data unchanged (cell-level; no bin50 dependency).
# =====================================================================
source("CORTEX_PROGRAM_ROOT/figures/fig4/fig4_theme.R")
suppressMessages({
  library(patchwork); library(grid); library(ggridges); library(ggrepel)
  library(ComplexHeatmap); library(circlize)
  library(ggraph); library(tidygraph); library(igraph); library(svglite)
})
set.seed(42)

OUT  <- "CORTEX_PROGRAM_ROOT/figures/fig4"
SVGD <- file.path(OUT, "svg_panels")
dir.create(SVGD, recursive = TRUE, showWarnings = FALSE)
svgf <- function(id) file.path(SVGD, sprintf("fig4_%s.svg", id))

# zero out plot.margin on any ggplot
nomar <- function(p) p + theme(plot.margin = margin(0, 0, 0, 0))

## ===================================================================
## panel a : subclass x program eta2 heatmap  (SQUARE CELLS -> banner)
## ===================================================================
build_a_grob <- function() {
  eta <- read.delim(file.path(RES, "within_subclass_region_eta2.tsv"))
  # RENUMBER 2026-06-20: filter out excluded programs (old cNMF indices)
  eta <- eta[!eta$program %in% .EXCLUDED_OLD, ]
  eta$program <- as.character(eta$program)
  M <- eta %>% select(subclass, program, eta2) %>%
    pivot_wider(names_from = program, values_from = eta2) %>% as.data.frame()
  rownames(M) <- M$subclass; M$subclass <- NULL
  M <- M[, order(as.integer(colnames(M)))]; M <- as.matrix(M)
  # Rename columns from old cNMF index to new P number
  colnames(M) <- as.character(.old2new_int[colnames(M)])
  drv <- read.delim(file.path(RES, "subclass_driver_rank.tsv")); M <- M[drv$subclass, ]
  cls <- sc_class(rownames(M))
  col_fun <- colorRamp2(quantile(M, c(0, 0.5, 0.8, 0.95, 1)),
                        c("#FCFDBF", "#FEC287", "#F1605D", "#8C2981", "#2D1160"))
  ra <- rowAnnotation(Class = cls, col = list(Class = pal_class),
        annotation_name_gp = gpar(fontsize = 6),
        annotation_legend_param = list(Class = list(labels_gp = gpar(fontsize = 6),
            title_gp = gpar(fontsize = 7))), simple_anno_size = unit(2.5, "mm"))
  # top 12 by eta2, among kept programs (old index used, map to new for mark matrix)
  eta_kept <- eta[!eta$program %in% as.character(.EXCLUDED_OLD), ]
  top <- eta_kept %>% arrange(desc(eta2)) %>% head(12)
  mk <- matrix("", nrow(M), ncol(M), dimnames = dimnames(M))
  for (i in seq_len(nrow(top))) {
    s <- top$subclass[i]
    p_old <- as.character(top$program[i])
    p_new <- as.character(.old2new_int[p_old])
    if (s %in% rownames(mk) && p_new %in% colnames(mk)) mk[s, p_new] <- "*"
  }
  # highlight programs: use new P numbers (colnames are now new P ints as strings)
  hi_progs_new <- sort(as.integer(colnames(M)[apply(mk == "*", 2, any)]))
  # prog_label_selective expects old cNMF int, but colnames are now new P ints
  # Build col labels directly from new P number + PROG_LABEL (via old key mapping)
  new2old <- setNames(.rmap$old_P[.rmap$new_P != "EXCLUDED"],
                       as.character(.old2new_int))
  col_labs <- sapply(colnames(M), function(np) {
    old_p <- new2old[np]
    if (!is.na(old_p) && as.character(old_p) %in% names(PROG_LABEL)) {
      if (as.integer(np) %in% hi_progs_new) PROG_LABEL[as.character(old_p)]
      else paste0("P", np)
    } else paste0("P", np)
  })
  s_mm <- 2.7   # square cell side
  ht <- Heatmap(M, name = "eta2", col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE,
    width  = ncol(M) * unit(s_mm, "mm"),
    height = nrow(M) * unit(s_mm, "mm"),
    show_row_dend = TRUE, show_column_dend = FALSE,  # drop col-dend height -> wider banner AR
    row_dend_width = unit(5, "mm"),
    row_names_gp = gpar(fontsize = 6), column_labels = col_labs, column_names_gp = gpar(fontsize = 5),
    column_title = "Programs (K=60, 6 cohort-technical excluded = 54 shown)", column_title_gp = gpar(fontsize = 7.5),
    row_title = "Subclass", row_title_gp = gpar(fontsize = 7.5), left_annotation = ra,
    heatmap_legend_param = list(title = expression(eta^2), labels_gp = gpar(fontsize = 6.5),
      title_gp = gpar(fontsize = 8), legend_height = unit(16, "mm"), grid_width = unit(3.5, "mm")),
    cell_fun = function(j, i, x, y, w, h, fill) { if (mk[i, j] == "*")
      grid.text("*", x, y, gp = gpar(fontsize = 6, col = "white", fontface = "bold")) },
    rect_gp = gpar(col = NA))
  grid.grabExpr(draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
                     padding = unit(c(0, 0, 0, 0), "mm")))
}

## ===================================================================
## panel b : region-driver ranking ridgeline (tall-narrow)
## ===================================================================
build_b <- function() {
  eta <- read.delim(file.path(RES, "within_subclass_region_eta2.tsv"))
  # RENUMBER 2026-06-20: filter excluded
  eta <- eta[!eta$program %in% .EXCLUDED_OLD, ]
  drv <- read.delim(file.path(RES, "subclass_driver_rank.tsv"))
  ord <- drv %>% arrange(median_eta2) %>% pull(subclass)
  eta$subclass <- factor(eta$subclass, levels = ord); eta$Class <- sc_class(as.character(eta$subclass))
  ggplot(eta, aes(eta2, subclass, fill = Class)) +
    geom_density_ridges(scale = 2.0, alpha = 0.78, linewidth = 0.22, colour = "white",
        quantile_lines = TRUE, quantiles = 2, vline_colour = "grey25", vline_size = 0.25) +
    scale_fill_manual(values = pal_class, name = "Class") +
    scale_x_continuous(expand = expansion(c(0.01, 0.05))) +
    labs(x = expression(paste("Within-subclass region ", eta^2)), y = NULL,
         title = "Region-driver ranking", subtitle = "over 54 programs") +
    theme_fig4(7) + theme(legend.position = c(0.78, 0.20),
      legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
      panel.grid.major.y = element_blank())
}

## ===================================================================
## panel c : variance sources stacked bars (tall-narrow)
## ===================================================================
build_c <- function() {
  pc <- read.delim(file.path(RES, "panel_c_partition.tsv"))
  # RENUMBER 2026-06-20: filter excluded
  pc <- pc[!pc$program %in% .EXCLUDED_OLD, ]
  pc$program <- factor(prog_label(pc$program), levels = prog_label(pc$program[order(pc$eta2_region)]))
  long <- pc %>% select(program, eta2_region, cell_autonomous_frac, compositional_frac) %>%
    pivot_longer(c(cell_autonomous_frac, compositional_frac), names_to = "component", values_to = "frac") %>%
    mutate(component = recode(component, cell_autonomous_frac = "Cell-autonomous",
           compositional_frac = "Compositional"), contrib = frac * eta2_region)
  ggplot(long, aes(contrib, program, fill = component)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.18) +
    scale_fill_manual(values = c(`Cell-autonomous` = "#8E44AD", Compositional = "#F39C12"), name = NULL) +
    scale_x_continuous(expand = expansion(c(0, 0.04))) +
    labs(x = expression(paste("Region ", eta^2, " partitioned")), y = NULL,
         title = "Variance sources", subtitle = "top 15 programs") +
    theme_fig4(7) + theme(legend.position = "top", legend.justification = "left",
      panel.grid.major.y = element_blank())
}

## ===================================================================
## panel d : spotlight L3-L4 IT RORB heatmap (SQUARE CELLS)
## ===================================================================
build_d <- function() {
  m <- read.delim(file.path(RES, "region_subclass_program_mean.tsv"))
  SC <- "L3-L4 IT RORB"; sub <- m %>% filter(subclass == SC); sub$program <- as.integer(sub$program)
  # RENUMBER 2026-06-20: filter excluded
  sub <- sub[!sub$program %in% .EXCLUDED_OLD, ]
  rng <- sub %>% group_by(program) %>% summarise(rng = max(mean) - min(mean), .groups = "drop") %>%
    arrange(desc(rng)) %>% head(20)
  sub2 <- sub %>% filter(program %in% rng$program) %>% group_by(program) %>%
    mutate(z = (mean - mean(mean)) / (sd(mean) + 1e-9)) %>% ungroup()
  sub2$program_f <- factor(prog_label(sub2$program), levels = prog_label(rev(rng$program)))
  p14_lab <- prog_label(14)
  reg_ord <- sub %>% filter(program == 14) %>% arrange(mean) %>% pull(region)
  sub2$region <- factor(sub2$region, levels = reg_ord)
  ggplot(sub2, aes(region, program_f, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    coord_fixed(ratio = 1) +     # SQUARE cells (14 regions x 20 programs)
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, name = "z") +
    annotate("rect", xmin = 0.5, xmax = length(reg_ord) + 0.5,
      ymin = which(levels(sub2$program_f) == p14_lab) - 0.5, ymax = which(levels(sub2$program_f) == p14_lab) + 0.5,
      fill = NA, colour = "#111111", linewidth = 0.7) +
    labs(x = paste0("Region (by ", p14_lab, ")"), y = NULL, title = "Spotlight L3-L4 IT RORB",
         subtitle = paste0(p14_lab, " highlighted")) +
    theme_fig4(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5), panel.grid = element_blank())
}

## ===================================================================
## panel e : top driver pairs lollipop/dotplot (tall-narrow)
## ===================================================================
build_e <- function() {
  eta <- read.delim(file.path(RES, "within_subclass_region_eta2.tsv")); eta$Class <- sc_class(eta$subclass)
  # RENUMBER 2026-06-20: filter excluded
  eta <- eta[!eta$program %in% .EXCLUDED_OLD, ]
  top <- eta %>% filter(fdr < 0.05) %>% arrange(desc(eta2)) %>% head(25)
  top$lab <- factor(paste0(top$subclass, " · ", prog_label(top$program)),
                    levels = rev(paste0(top$subclass, " · ", prog_label(top$program))))
  ggplot(top, aes(eta2, lab, colour = Class)) +
    geom_segment(aes(x = 0, xend = eta2, yend = lab), linewidth = 0.35, colour = "grey75") +
    geom_point(aes(size = n_cells)) +
    scale_colour_manual(values = pal_class, name = "Class") +
    scale_size_continuous(range = c(1, 4.5), name = "n cells", breaks = c(10000, 50000, 100000),
      labels = c("10k", "50k", "100k")) + scale_x_continuous(expand = expansion(c(0, 0.06))) +
    labs(x = expression(paste("region ", eta^2)), y = NULL, title = "Top driver pairs", subtitle = "FDR<0.05") +
    theme_fig4(7) + theme(axis.text.y = element_text(size = 5), legend.position = "right",
      legend.box = "vertical", panel.grid.major.y = element_blank())
}

## ===================================================================
## panel f : UMAP twin (region | p14)  -- ggplot, coord_fixed, vector
##   legibility fix: bigger titles (9pt) + legend text (6pt); NO tiny arrows
## ===================================================================
build_f <- function() {
  d <- read.delim(file.path(RES, "panel_f_umap.tsv"))
  d <- d[sample(nrow(d)), ]
  p_reg <- ggplot(d, aes(UMAP1, UMAP2, colour = region)) +
    ggrastr::rasterise(geom_point(size = 0.18, alpha = 0.6), dpi = 400) +
    coord_fixed(ratio = 1) +
    guides(colour = guide_legend(override.aes = list(size = 2.4), ncol = 2)) +
    labs(title = "L3-L4 IT RORB · region", x = NULL, y = NULL) +
    theme_void(base_size = 7) + theme(legend.position = "right",
      legend.text = element_text(size = 6), legend.key.size = unit(3.2, "mm"),
      legend.title = element_blank(), plot.title = element_text(face = "bold", size = 9, hjust = 0),
      plot.margin = margin(0, 0, 0, 0))
  p14_lab <- prog_label(14)
  p_p14 <- ggplot(d, aes(UMAP1, UMAP2, colour = p14)) +
    ggrastr::rasterise(geom_point(size = 0.18, alpha = 0.7), dpi = 400) +
    coord_fixed(ratio = 1) +
    scale_colour_gradientn(colours = c("#000004", "#51127C", "#B63679", "#FB8861", "#FCFDBF"),
      name = p14_lab, limits = quantile(d$p14, c(.02, .98)), oob = scales::squish) +
    labs(title = paste0(p14_lab, " activity"), x = NULL, y = NULL) +
    theme_void(base_size = 7) + theme(legend.position = "right",
      legend.key.height = unit(7, "mm"), legend.key.width = unit(2.6, "mm"),
      legend.text = element_text(size = 6), legend.title = element_text(size = 7),
      plot.title = element_text(face = "bold", size = 9, hjust = 0),
      plot.margin = margin(0, 0, 0, 0))
  (p_reg | p_p14) + plot_layout(widths = c(1, 1))
}

## ===================================================================
## panel g : bootstrap robustness boxplot grid
##   LEGIBILITY REBUILD (2026-06-10): the 8-facet 4x2 grid shrank to ~0.53x
##   at submission size so all tick/strip/label fell <5pt (unreadable).
##   Fix: (1) trim to 6 REPRESENTATIVE driver subclasses (3 top excitatory
##   drivers + NP + the L3-L4 IT RORB spotlight + AST as the non-neuronal
##   representative), 3x2 grid -> each facet ~1.33x wider; (2) top-6 programs
##   per subclass (was 8) -> fewer, taller rows per facet; (3) base font 7->10,
##   strip 6.5->10.5, axis text 5->10 (compose shrink ~0.50x => strip 5.3pt,
##   axis 5.05pt, all >=5pt at the final 180 mm submission width). Banner aspect
##   kept wide so the H[g,h] row height (and page AR) is unchanged.
## ===================================================================
build_g <- function() {
  b <- read.delim(file.path(RES, "panel_g_bootstrap.tsv")); b$Class <- sc_class(b$subclass)
  # RENUMBER 2026-06-20: filter excluded
  b <- b[!b$program %in% .EXCLUDED_OLD, ]
  drv <- read.delim(file.path(RES, "subclass_driver_rank.tsv"))
  # representative subset: top excitatory drivers + spotlight + non-neuronal rep
  keep_sc <- c("L6 CT", "L6 IT", "ET", "NP", "L3-L4 IT RORB", "AST")
  b <- b %>% filter(subclass %in% keep_sc)
  sc_levels <- intersect(drv$subclass, keep_sc); b$subclass <- factor(b$subclass, levels = sc_levels)
  # keep top-6 programs per subclass (by full_eta2) to de-clutter the dense grid
  keep <- b %>% distinct(subclass, program, full_eta2) %>% group_by(subclass) %>%
    arrange(desc(full_eta2)) %>% slice_head(n = 6) %>% ungroup() %>%
    mutate(k = paste(subclass, program)) %>% pull(k)
  b <- b %>% mutate(k = paste(subclass, program)) %>% filter(k %in% keep)
  prog_order <- b %>% distinct(subclass, program, full_eta2) %>% arrange(subclass, desc(full_eta2)) %>%
    mutate(pl = prog_label(program))
  b$pl <- prog_label(b$program); b$key <- paste(b$subclass, b$pl, sep = "::")
  key_levels <- prog_order %>% mutate(key = paste(subclass, pl, sep = "::")) %>% pull(key)
  b$key <- factor(b$key, levels = rev(key_levels))
  ggplot(b, aes(eta2, key, fill = Class)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.28, width = 0.62, outlier.alpha = 0.3) +
    geom_point(aes(x = full_eta2), shape = 23, size = 1.3, fill = "white", colour = "black", stroke = 0.4) +
    facet_wrap(~subclass, scales = "free", ncol = 3,
               labeller = labeller(subclass = function(x) gsub("L3-L4 IT RORB", "L3-L4 IT\nRORB", x))) +
    scale_fill_manual(values = pal_class, name = "Class") +
    scale_y_discrete(labels = function(x) sub(".*::", "", x)) +
    scale_x_continuous(breaks = scales::extended_breaks(n = 3),
                       labels = function(v) sub("^0\\.", ".", formatC(v, format = "g", digits = 2)),
                       guide = guide_axis(check.overlap = TRUE)) +
    labs(x = expression(paste("bootstrap ", eta^2)), y = NULL,
         title = "Bootstrap robustness of region-driver effect (54 programs)",
         subtitle = "top 6 programs per representative driver subclass; diamond = full-data estimate") +
    theme_fig4(10) + theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10),
      strip.text = element_text(size = 10.5, face = "bold"),
      legend.position = "top", legend.justification = "left", panel.spacing = unit(2.4, "mm"))
}

## ===================================================================
## panel h : driver network subclass -> program -> region (3 layers)
##   legibility fix: node text 1.5 -> 2.0
## ===================================================================
build_h <- function() {
  eta <- read.delim(file.path(RES, "within_subclass_region_eta2.tsv")); eta$program <- as.integer(eta$program)
  m <- read.delim(file.path(RES, "region_subclass_program_mean.tsv")); m$program <- as.integer(m$program)
  # RENUMBER 2026-06-20: filter excluded
  eta <- eta[!eta$program %in% .EXCLUDED_OLD, ]
  m   <- m[!m$program %in% .EXCLUDED_OLD, ]
  sp <- eta %>% filter(fdr < 0.05) %>% arrange(desc(eta2)) %>% head(22)
  prog_keep <- sort(unique(sp$program)); sc_keep <- unique(sp$subclass)
  pr <- m %>% filter(program %in% prog_keep) %>% group_by(program, region) %>%
    summarise(mean = mean(mean), .groups = "drop") %>% group_by(program) %>%
    mutate(spec = (mean - min(mean)) / (max(mean) - min(mean) + 1e-9)) %>%
    arrange(program, desc(mean)) %>% slice_head(n = 2) %>% ungroup()
  reg_keep <- unique(pr$region)
  nodes <- bind_rows(
    data.frame(name = sc_keep, type = "Subclass", layer = 1),
    data.frame(name = paste0("p", prog_keep), type = "Program", layer = 2),
    data.frame(name = reg_keep, type = "Region", layer = 3)) %>% distinct(name, .keep_all = TRUE)
  nodes$class <- ifelse(nodes$type == "Subclass", sc_class(nodes$name), nodes$type)
  nodes$disp <- nodes$name
  .is_prog <- nodes$type == "Program"
  nodes$disp[.is_prog] <- prog_label(as.integer(sub("^p", "", nodes$name[.is_prog])))
  e1 <- sp %>% transmute(from = subclass, to = paste0("p", program), weight = eta2, etype = "sc_prog")
  e2 <- pr %>% transmute(from = paste0("p", program), to = region, weight = spec, etype = "prog_reg")
  edges <- bind_rows(e1, e2)
  g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  lay <- nodes %>% group_by(layer) %>%
    mutate(y = scales::rescale(rank(name, ties.method = "first"), to = c(0, 1))) %>% ungroup()
  pal_node <- c(pal_class, Program = "#7F8C8D", Region = "#34495E")
  ggraph(g, layout = "manual", x = lay$layer, y = lay$y) +
    geom_edge_diagonal(aes(edge_width = weight, edge_colour = etype, edge_alpha = weight)) +
    scale_edge_width(range = c(0.15, 1.4), guide = "none") +
    scale_edge_alpha(range = c(0.22, 0.85), guide = "none") +
    scale_edge_colour_manual(values = c(sc_prog = "#C0392B", prog_reg = "#2471A3"),
      labels = c(sc_prog = "subclass→program", prog_reg = "program→region"), name = NULL) +
    geom_node_point(aes(colour = class, shape = type), size = 1.8) +
    scale_colour_manual(values = pal_node, name = NULL, breaks = c("Excitatory", "Inhibitory", "Non-neuronal")) +
    scale_shape_manual(values = c(Subclass = 16, Program = 15, Region = 17), name = NULL) +
    geom_node_text(aes(label = disp, hjust = ifelse(layer == 1, 1.10, ifelse(layer == 3, -0.10, 0.5)),
      vjust = ifelse(layer == 2, -0.9, 0.5)), size = 2.0, colour = "grey15") +
    scale_x_continuous(expand = expansion(c(0.20, 0.20))) +
    scale_y_continuous(expand = expansion(c(0.05, 0.09))) +
    annotate("text", x = 1, y = 1.09, label = "Subclass", size = 2.4, fontface = "bold") +
    annotate("text", x = 2, y = 1.09, label = "Program", size = 2.4, fontface = "bold") +
    annotate("text", x = 3, y = 1.09, label = "Region", size = 2.4, fontface = "bold") +
    labs(title = "Driver network", subtitle = "subclass → program → region") +
    theme_void(base_size = 7) + theme(plot.title = element_text(face = "bold", size = 9, hjust = 0),
      plot.subtitle = element_text(size = 6.5, colour = "grey35"), legend.position = "bottom",
      legend.text = element_text(size = 6), legend.key.size = unit(3, "mm"),
      legend.margin = margin(2, 2, 2, 2), legend.box.spacing = unit(2, "mm"),
      plot.margin = margin(0, 0, 0, 0))
}

## ===================================================================
## render each panel to ZERO-margin svglite SVG
## ===================================================================
cat("rendering fig4 panels to zero-margin SVG...\n")

# a : heatmap grob -> measure native then save (square cells fix height)
gA <- build_a_grob()
svglite(svgf("a"), width = 9.5, height = 3.6, bg = "white")
grid.draw(gA); invisible(dev.off()); cat("  a done\n")

save_panel <- function(id, p, w, h) {
  svglite(svgf(id), width = w, height = h, bg = "white")
  print(nomar(p)); invisible(dev.off()); cat(sprintf("  %s done\n", id))
}

save_panel("b", build_b(), 3.0, 4.6)
save_panel("c", build_c(), 3.2, 4.2)
save_panel("d", build_d(), 4.2, 4.2)
save_panel("e", build_e(), 3.4, 4.6)

# f : patchwork twin (two coord_fixed UMAPs) -> wide
pf <- build_f()
svglite(svgf("f"), width = 9.0, height = 3.0, bg = "white")
print(pf & theme(plot.margin = margin(0, 0, 0, 0))); invisible(dev.off()); cat("  f done\n")

save_panel("g", build_g(), 9.0, 4.2)  # 6-facet 3x2 wide banner; bigger fonts (legibility rebuild 2026-06-10)
save_panel("h", build_h(), 4.6, 4.4)

cat("ALL fig4 svg panels written to", SVGD, "\n")
