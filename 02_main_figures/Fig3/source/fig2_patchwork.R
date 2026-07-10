#!/usr/bin/env Rscript
# =====================================================================
# Fig.2  Spatial program backbone of human cortex (bin50, 50 µm)
# 8 panels (a-h). Submission-grade, Nature-family style.
# SINGLE NATIVE PATCHWORK script (replaces python+R+PIL-montage pipeline).
# Panels a/b/f were python(matplotlib) -> here PORTED to ggplot.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
  library(circlize); library(ggridges); library(scatterpie)
  library(ggrastr); library(data.table)
  library(dplyr); library(tidyr); library(grid)
  library(scales); library(RColorBrewer)
})

set.seed(42)

## ---- paths ----------------------------------------------------------
WORK <- "CORTEX_PROGRAM_ROOT/scripts/fig2"
RES  <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUT  <- "CORTEX_PROGRAM_ROOT/figures/fig2"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- palettes (unified family across panels) ------------------------
LAYER_ORDER  <- c("ARACHNOID","L1","L2","L3","L4","L5","L6","WM")
LAYER_COLORS <- c(ARACHNOID="#9E9E9E", L1="#3B4CC0", L2="#5A78D6", L3="#7DA0E0",
                  L4="#36A66B", L5="#E8A93B", L6="#D65A3B", WM="#7B3294")
CELLTYPES <- c("AST","CHANDELIER","ENDO","ET","L2-L3 IT LINC00507","L3-L4 IT RORB",
  "L4-L5 IT RORB","L6 CAR3","L6 CT","L6 IT","L6B","LAMP5","MICRO","NDNF","NP",
  "OLIGO","OPC","PAX6","PVALB","SST","VIP","VLMC")
CT_COLORS <- c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df8a",
  "#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b","#c49c94","#e377c2",
  "#f7b6d2","#7f7f7f","#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5",
  "#393b79","#637939")
names(CT_COLORS) <- CELLTYPES
# z-activity diverging palette (CNS blue-white-red), matches old ZCMAP intent
div_pal <- colorRampPalette(c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7",
                              "#FDDBC7","#EF8A62","#B2182B"))(256)

## ---- Nature theme (mirror fig3 theme_nat) ---------------------------
theme_nat <- function(base = 6.6) {
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      axis.line   = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.35, colour = "black"),
      axis.title  = element_text(size = base),
      axis.text   = element_text(size = base - 0.6, colour = "black"),
      legend.title= element_text(size = base - 0.4),
      legend.text = element_text(size = base - 1.0),
      legend.key.size = unit(3.0, "mm"),
      strip.text  = element_text(size = base - 0.2, face = "bold"),
      strip.background = element_blank(),
      plot.title  = element_text(size = base + 1, face = "bold"),
      plot.subtitle = element_text(size = base - 0.6, colour = "grey30"),
      panel.grid  = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}
theme_set(theme_nat())

## ---- choices (rep chip + exemplar programs) -------------------------
ch <- readLines(file.path(WORK, "_choices.txt"))
REP <- sub("REP=", "", grep("^REP=", ch, value = TRUE)[1])
ex_lines <- grep("^EX\t", ch, value = TRUE)
EX <- setNames(sub(".*\t.*\t", "", ex_lines), sub("^EX\t", "", sub("\t[^\t]*$", "", ex_lines)))

## ---- program -> functional GO:BP label lookup ----------------------
## label = "P{n} {name_short}" (strip trailing " P{n}"; "*" if brain-weak)
NM <- fread(file.path(RES, "program_names.tsv"))
.mk_label <- function(n) {
  r <- NM[program == as.integer(n)]
  if (nrow(r) == 0) return(paste0("P", n))
  ns   <- trimws(sub("\\s+P[0-9]+$", "", as.character(r$name_short[1])))
  star <- if (as.character(r$confidence[1]) == "brain-weak") "*" else ""
  paste0("P", n, " ", ns, star)
}
PLABEL <- setNames(vapply(NM$program, .mk_label, character(1)), as.character(NM$program))
plabel <- function(prog_key) {        # accepts "program_37" or "37"
  n <- sub("program_", "", as.character(prog_key))
  out <- PLABEL[n]; out[is.na(out)] <- paste0("P", n[is.na(out)]); unname(out)
}
## anatomical/cell-type context tag per exemplar (spatial anchor)
EX_CONTEXT <- c(program_5="L2/3", program_35="L4", program_53="L6",
                program_37="OLIGO / WM", program_56="ENDO / vascular",
                program_29="Inhibitory")
ex_title <- function(prog_key, sep = " · ") {
  ctx <- EX_CONTEXT[as.character(prog_key)]
  lab <- plabel(prog_key)
  ifelse(is.na(ctx), lab, paste0(ctx, sep, lab))
}

## =====================================================================
## load rep-chip spatial data (shared by panels a, b)
## =====================================================================
meta <- fread(file.path(WORK, "repchip_meta.tsv"))
prog <- fread(file.path(WORK, "repchip_progscores.tsv"))
df <- merge(meta, prog, by = "bin")
df[, yy := max(y) - y]                # flip y so pia/top reads naturally
cat("rep chip bins for plot:", nrow(df), "\n")

## =====================================================================
## PANEL a -- rep-chip tissue scatter colored by majorDomain (layers)
## PORTED from python matplotlib. ggplot geom_point, rasterized via ggrastr.
## =====================================================================
df_a <- df %>% mutate(majorDomain = factor(majorDomain, levels = LAYER_ORDER))
# scale bar 1 mm = 1000 um, bottom-left
xr <- range(df_a$x); yr <- range(df_a$yy)
sb_x0 <- xr[1] + diff(xr) * 0.05
sb_y0 <- yr[1] + diff(yr) * 0.03
p_a <- ggplot(df_a, aes(x, yy, colour = majorDomain)) +
  ggrastr::rasterise(geom_point(shape = 15, size = 0.18, stroke = 0), dpi = 450) +
  geom_segment(aes(x = sb_x0, xend = sb_x0 + 1000, y = sb_y0, yend = sb_y0),
               inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
  annotate("text", x = sb_x0 + 500, y = sb_y0 + diff(yr) * 0.03,
           label = "1 mm", size = 1.9, vjust = 0) +
  scale_colour_manual(values = LAYER_COLORS, name = "majorDomain",
                      drop = FALSE,
                      guide = guide_legend(override.aes = list(size = 1.8))) +
  coord_equal(clip = "off") +
  labs(title = sprintf("Cortical layers (chip %s)", REP)) +
  theme_void(base_size = 6.6) +
  theme(plot.title = element_text(size = 7.0, face = "bold", hjust = 0,
                                  margin = margin(b = 3)),
        legend.title = element_text(size = 6.2),
        legend.text  = element_text(size = 5.6),
        legend.key.size = unit(2.6, "mm"),
        legend.position = "right",
        plot.background = element_rect(fill = "white", colour = NA))
cat("panel a done\n")

## =====================================================================
## PANEL b -- 6-program spatial small-multiples (z activity) on rep chip
## PORTED from python matplotlib gridspec. ggplot facet_wrap.
## =====================================================================
exprogs <- names(EX)
# per-program symmetric vlim (98th pct of |z|), then squish for shared scale display
b_long <- df %>% select(x, yy, all_of(exprogs)) %>%
  pivot_longer(all_of(exprogs), names_to = "program", values_to = "z")
# clamp each program to its own 98th-pct |z| so facets are individually legible,
# then map onto a shared symmetric diverging scale in [-1,1] (z / vlim)
b_vlim <- b_long %>% group_by(program) %>%
  summarise(vlim = quantile(abs(z), 0.98, na.rm = TRUE), .groups = "drop")
b_long <- b_long %>% left_join(b_vlim, by = "program") %>%
  mutate(zc = pmax(pmin(z / vlim, 1), -1),
         prog_label = factor(ex_title(program, sep = "\n"),
                             levels = ex_title(exprogs, sep = "\n")))
p_b <- ggplot(b_long, aes(x, yy, colour = zc)) +
  ggrastr::rasterise(geom_point(shape = 15, size = 0.14, stroke = 0), dpi = 400) +
  facet_wrap(~ prog_label, nrow = 2) +
  scale_colour_gradientn(colours = div_pal, limits = c(-1, 1),
                         breaks = c(-1, 0, 1), labels = c("low", "0", "high"),
                         name = "z-activity\n(per-program)") +
  coord_equal(clip = "off") +
  labs(title = sprintf("Spatial program activity (chip %s)", REP)) +
  theme_void(base_size = 6.6) +
  theme(plot.title = element_text(size = 7.4, face = "bold", hjust = 0,
                                  margin = margin(b = 3)),
        strip.text = element_text(size = 5.8, face = "bold", lineheight = 0.95,
                                  margin = margin(b = 1, t = 1)),
        legend.title = element_text(size = 5.8),
        legend.text  = element_text(size = 5.2),
        legend.key.height = unit(4.0, "mm"), legend.key.width = unit(2.4, "mm"),
        legend.position = "right",
        panel.spacing = unit(1.2, "mm"),
        plot.background = element_rect(fill = "white", colour = NA))
cat("panel b done\n")

## =====================================================================
## PANEL c -- scatterpie of RCTD composition on coarse grid (already R)
## =====================================================================
rctd <- fread(file.path(WORK, "repchip_rctd.tsv"))
m <- merge(meta[, .(bin, x, y)], rctd[rctd_pass_mask == TRUE], by = "bin")
gx <- 700
m[, gx2 := floor(x / gx) * gx]
m[, gy2 := floor(y / gx) * gx]
agg <- m[, c(lapply(.SD, mean), .(n = .N)), by = .(gx2, gy2), .SDcols = CELLTYPES]
agg <- agg[n >= 4]
comp <- as.matrix(agg[, ..CELLTYPES]); comp <- comp / rowSums(comp)
agg[, (CELLTYPES) := as.data.table(comp)]
agg[, r := gx * 0.42]
agg[, yy := max(gy2) - gy2]
cat("panel c pies:", nrow(agg), "\n")
p_c <- ggplot() +
  geom_scatterpie(aes(x = gx2, y = yy, r = r), data = agg, cols = CELLTYPES,
                  color = NA, linewidth = 0) +
  scale_fill_manual(values = CT_COLORS, name = "cell type") +
  coord_equal(clip = "off") +
  labs(title = sprintf("RCTD composition\n(chip %s, %d-µm grid)", REP, gx)) +
  theme_void(base_size = 6.6) +
  theme(plot.title = element_text(size = 7.0, face = "bold", hjust = 0,
                                  lineheight = 1.0, margin = margin(b = 3)),
        plot.margin = margin(t = 4, r = 2, b = 2, l = 2),
        legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 5.0),
        legend.title = element_text(size = 6.0),
        plot.background = element_rect(fill = "white", colour = NA)) +
  guides(fill = guide_legend(ncol = 2))
cat("panel c done\n")

## =====================================================================
## PANEL d -- program activity along cortical depth (ridgelines)
## functional names, pia->WM (already R / ggridges)
## =====================================================================
pcl <- fread(file.path(WORK, "prog_x_layer_per_chip.tsv"))
selprog <- names(EX)
pdd <- pcl[program %in% selprog & n >= 50]
pdd[, majorDomain := factor(majorDomain, levels = rev(LAYER_ORDER))]
pdd[, prog_label := factor(ex_title(program), levels = ex_title(selprog))]
p_d <- ggplot(pdd, aes(x = mean_z, y = majorDomain, fill = after_stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 1.6, rel_min_height = 0.01,
      linewidth = 0.2, colour = "grey30") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, name = "mean z") +
  facet_wrap(~ prog_label, nrow = 2) +
  labs(title = "Program activity along cortical depth (pia→WM, across 44 chips)",
       x = "mean program z-activity", y = "cortical depth") +
  theme_nat() +
  theme(strip.text = element_text(size = 5.8, face = "bold"),
        panel.spacing = unit(1.5, "mm"),
        legend.key.size = unit(3.0, "mm"))
cat("panel d done\n")

## =====================================================================
## PANEL e -- program x layer ComplexHeatmap (~60 programs, 2-COLUMN)
## spotlight full functional names, rest "P{n} . {peakct}".
## captured via grid.grabExpr (mirror fig3 panel a).
## IMPORTANT: built as a SINGLE Heatmap with column_split into 2 blocks
## (programs 1-30 | 31-60 of the clustered order), NOT a `ht1 + ht2`
## HeatmapList. A HeatmapList draw() corrupts the grid viewport stack on
## grid.grabExpr (viewport-overwrite -> empty grob -> blanks every panel
## after e in the patchwork). A single Heatmap with 2 column blocks grabs
## cleanly while still giving the same wide 2-column 30-row layout.
## =====================================================================
G <- fread(file.path(WORK, "prog_x_layer_global.tsv"))
G <- G[match(LAYER_ORDER, majorDomain)]
progcols <- paste0("program_", 1:60)
Me <- t(as.matrix(G[, ..progcols]))        # 60 programs x 8 layers
rownames(Me) <- sub("program_", "P", progcols)
colnames(Me) <- G$majorDomain
col_fun <- colorRamp2(c(-1.2, -0.4, 0, 0.4, 1.2),
   c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"))
cc <- fread(file.path(WORK, "program_celltype_corr.tsv"))
setnames(cc, 1, "program")
peakct <- apply(as.matrix(cc[, ..CELLTYPES]), 1, function(r) CELLTYPES[which.max(r)])
names(peakct) <- cc$program
SPOTLIGHT <- c(5, 29, 35, 37, 53, 56, 26)
pnum <- as.integer(sub("program_", "", progcols))
rowlab_all <- ifelse(pnum %in% SPOTLIGHT,
                 paste0(plabel(progcols), " · ", peakct[progcols]),
                 paste0(rownames(Me), " · ", peakct[progcols]))
names(rowlab_all) <- rownames(Me)
face_all <- ifelse(pnum %in% SPOTLIGHT, "bold", "plain"); names(face_all) <- rownames(Me)
size_all <- ifelse(pnum %in% SPOTLIGHT, 5.6, 5.2);        names(size_all) <- rownames(Me)
# cluster once, split clustered order into two equal halves -> side-by-side blocks
hc   <- hclust(dist(Me, method = "euclidean"), method = "ward.D2")
ord  <- hc$order
half <- ceiling(length(ord) / 2)           # 30
idx1 <- ord[1:half]; idx2 <- ord[(half + 1):length(ord)]
rn1  <- rownames(Me)[idx1]; rn2 <- rownames(Me)[idx2]
# 30 x 16 matrix: cols 1-8 = layers of block1 progs, cols 9-16 = block2 progs
W <- cbind(Me[idx1, , drop = FALSE], Me[idx2, , drop = FALSE])
rownames(W) <- NULL
col_block <- factor(rep(c("b1", "b2"), each = length(LAYER_ORDER)),
                    levels = c("b1", "b2"))
# top layer-color strip (8 layers repeated per block)
top_e <- HeatmapAnnotation(layer = rep(LAYER_ORDER, 2),
   col = list(layer = LAYER_COLORS), show_legend = FALSE,
   annotation_name_gp = gpar(fontsize = 5.5), simple_anno_size = unit(2.2, "mm"))
# row labels: block1 names on LEFT (face/size per spotlight), block2 names on RIGHT
gp1 <- gpar(fontsize = size_all[rn1], fontface = face_all[rn1])
gp2 <- gpar(fontsize = size_all[rn2], fontface = face_all[rn2])
ra_left  <- rowAnnotation(b1 = anno_text(rowlab_all[rn1], gp = gp1,
                                         location = 1, just = "right"))
ra_right <- rowAnnotation(b2 = anno_text(rowlab_all[rn2], gp = gp2,
                                         location = 0, just = "left"))
ht_e <- Heatmap(W, name = "mean z", col = col_fun,
   cluster_columns = FALSE, cluster_rows = FALSE,
   row_order = seq_len(nrow(W)), column_order = seq_len(ncol(W)),
   column_split = col_block, column_gap = unit(4, "mm"),
   column_labels = rep(LAYER_ORDER, 2),
   column_names_gp = gpar(fontsize = 5.6), column_names_rot = 45,
   top_annotation = top_e, left_annotation = ra_left, right_annotation = ra_right,
   show_heatmap_legend = TRUE,
   heatmap_legend_param = list(title_gp = gpar(fontsize = 6.5),
       labels_gp = gpar(fontsize = 5.5), legend_height = unit(15, "mm"),
       grid_width = unit(2.5, "mm")),
   width = unit(74, "mm"), height = unit(64, "mm"),
   column_title = c("Program × cortical layer  (block 1)", "(block 2)"),
   column_title_gp = gpar(fontsize = 6.4))
# capture as grob for patchwork assembly (vector-preserving, single Heatmap)
PANEL_E_W <- 11.0; PANEL_E_H <- 4.0
gb_e <- grid.grabExpr(
  draw(ht_e, heatmap_legend_side = "right", padding = unit(c(2, 2, 2, 2), "mm")),
  width = PANEL_E_W, height = PANEL_E_H)
p_e <- patchwork::wrap_elements(full = gb_e)
cat("panel e done (single-Heatmap 2-block, grabExpr)\n")

## =====================================================================
## PANEL f -- program-z vs RCTD-weight hexbin, exemplar pairs
## PORTED from python hexbin. ggplot geom_hex, facet by pair.
## =====================================================================
fdat <- fread(file.path(WORK, "panelf_pairs.tsv"))
PAIRS <- list(
  list(pp = "program_37", cc = "OLIGO",               xl = "OLIGO weight"),
  list(pp = "program_56", cc = "ENDO",                 xl = "ENDO weight"),
  list(pp = "program_5",  cc = "L2-L3 IT LINC00507",   xl = "L2/3 IT weight"))
f_list <- lapply(PAIRS, function(P) {
  x <- fdat[[P$cc]]; y <- fdat[[P$pp]]
  ok <- is.finite(x) & is.finite(y)
  r  <- cor(x[ok], y[ok])
  data.frame(wt = x[ok], z = y[ok],
             facet = sprintf("%s  (r=%.2f)", plabel(P$pp), r),
             xlab  = P$xl, prog = plabel(P$pp))
})
# facets need their own x-axis title; encode via free panels + per-facet xlab.
# Build 3 separate ggplots so each gets its correct x-axis label, then combine.
mk_hex <- function(d, P) {
  r <- cor(d$wt, d$z)
  ggplot(d, aes(wt, z)) +
    geom_hex(bins = 42) +
    scale_fill_gradient(low = "#fcdbb8", high = "#5b1a6b",
                        trans = "log10", name = expression(log[10]~"bins")) +
    annotate("text", x = -Inf, y = Inf, label = sprintf("r = %.2f", r),
             hjust = -0.15, vjust = 1.4, size = 1.9) +
    labs(x = P$xl, y = sprintf("%s z", plabel(P$pp)),
         title = sprintf("%s", plabel(P$pp))) +
    theme_nat() +
    theme(plot.title = element_text(size = 6.0, face = "bold"),
          axis.title = element_text(size = 5.8),
          legend.key.size = unit(2.6, "mm"),
          legend.position = "right")
}
pf_list <- mapply(function(d, P) mk_hex(d, P), f_list, PAIRS, SIMPLIFY = FALSE)
# collect with one shared legend on the right
p_f <- (pf_list[[1]] + pf_list[[2]] + pf_list[[3]]) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(title = "Program activity vs RCTD cell-type weight",
                  theme = theme(plot.title = element_text(size = 7.4, face = "bold",
                                                          family = "Helvetica")))
p_f <- patchwork::wrap_elements(full = p_f)
cat("panel f done\n")

## =====================================================================
## PANEL g -- cross-chip reproducibility boxplot (already R)
## =====================================================================
gsum  <- fread(file.path(WORK, "panelg_summary.tsv"))
grep_ <- fread(file.path(WORK, "panelg_reproducibility.tsv"))
gord  <- gsum[order(median)]$program
grep_[, program := factor(program, levels = gord)]
medmap <- setNames(gsum$median, gsum$program)
grep_[, medv := medmap[as.character(program)]]
p_g <- ggplot(grep_, aes(x = program, y = corr_to_mean)) +
  geom_hline(yintercept = c(0, 0.5, 1), colour = "grey85", linewidth = 0.25) +
  geom_boxplot(aes(fill = medv), outlier.size = 0.2, outlier.colour = "grey60",
               linewidth = 0.15, width = 0.72) +
  scale_fill_gradientn(colours = c("#D73027", "#FEE08B", "#1A9850"),
      limits = c(min(grep_$medv), 1), name = "median r") +
  scale_x_discrete(labels = function(x) sub("program_", "P", x)) +
  labs(title = "Cross-chip reproducibility of program layer profiles (n=44 chips)",
       x = "program (sorted by median)", y = "correlation to mean layer profile") +
  coord_cartesian(ylim = c(min(0, quantile(grep_$corr_to_mean, 0.01)), 1)) +
  theme_nat() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3.8),
        legend.key.size = unit(3.0, "mm"))
cat("panel g done\n")

## =====================================================================
## PANEL h -- WM-distance radial/line profile (already R)
## =====================================================================
rings <- fread(file.path(WORK, "panelh_rings.tsv"))
rings[, ring := as.numeric(as.character(ring))]
wmp <- c("program_37", "program_26", "program_5", "program_56")
.role <- c(program_37 = "OLIGO/WM", program_26 = "OLIGO/WM",
           program_5  = "L2/3 ctrl", program_56 = "vascular")
plab <- setNames(paste0(plabel(wmp), " (", .role[wmp], ")"), wmp)
longh <- melt(rings, id.vars = c("ring", "chip"), measure.vars = wmp,
              variable.name = "program", value.name = "mean_z")
longh[, prog_label := factor(plab[as.character(program)], levels = plab[wmp])]
prog_pal <- setNames(c("#7B3294", "#C2A5CF", "#5A78D6", "#D6604D"), plab[wmp])
p_h <- ggplot(longh, aes(x = ring, y = mean_z, colour = prog_label,
                         group = interaction(prog_label, chip))) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_line(aes(linetype = chip), linewidth = 0.5) +
  geom_point(size = 0.6) +
  scale_color_manual(values = prog_pal, name = "program") +
  scale_linetype_manual(values = c("solid", "31"), name = "chip") +
  labs(title = "Program activity vs distance from white-matter seeds (OLIGO≥0.5)",
       x = "distance from WM (µm)", y = "mean program z-activity") +
  theme_nat() +
  theme(legend.position = "right", legend.box = "vertical",
        legend.key.size = unit(3.0, "mm"))
cat("panel h done\n")

## =====================================================================
## ASSEMBLE -- native patchwork design (mirror fig3 layout strategy)
## balanced grid; e gets a wide row (2-col heatmap); f/g/h wide rows.
## =====================================================================
design <- "
AAABBBB
AAABBBB
CCCDDDD
CCCDDDD
EEEEEEE
EEEEEEE
FFFFFFF
GGGGGGG
HHHHHHH
"
fig <- p_a + p_b + p_c + p_d + p_e + p_f + p_g + p_h +
  plot_layout(design = design,
              heights = c(1, 1, 1, 1, 0.85, 0.85, 1.6, 0.86, 1.05)) +
  plot_annotation(tag_levels = "a",
    title = "Fig. 2  Spatial program backbone of human cortex (bin50, 50 µm)",
    theme = theme(plot.title = element_text(size = 11, face = "bold",
                                            family = "Helvetica"))) &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(face = "bold", size = 10, family = "Helvetica",
                                hjust = 0, vjust = 1))

FIG_W <- 13.0; FIG_H <- 17.5
ggsave(file.path(OUT, "fig2_full.pdf"), fig,
       width = FIG_W, height = FIG_H, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "fig2_full.png"), fig,
       width = FIG_W, height = FIG_H, dpi = 150, bg = "white", limitsize = FALSE)
cat("assembly done\n")
cat("\nALL PANELS COMPLETE -> figures/fig2/fig2_full.pdf\n")
