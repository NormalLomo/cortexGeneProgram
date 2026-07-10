#!/usr/bin/env Rscript
# =====================================================================
# Fig.2  -> 8 STANDALONE zero-margin SVG panels (a-h) for svgutils compose.
# REUSES the per-panel construction code from fig2_patchwork.R, but each
# panel is rendered as its OWN tight zero-margin SVG to:
#   figures/fig2/svg_panels/fig2_{a..h}.svg
# (does NOT modify fig2_patchwork.R; this is a sibling script.)
#
# Rules honoured:
#   * each panel device sized close to natural content (tight crop)
#   * theme(plot.margin = margin(0,0,0,0)) for ggplot panels; for the
#     ComplexHeatmap panel e: draw(ht, padding = unit(c(0,0,0,0),"mm"))
#   * spatial / x-y panels (a,b,c,f) use coord_fixed(1)/coord_equal so the
#     content is square (no stretch)
#   * titles SHORTENED to concise <=~22 char single-line forms (avoids the
#     title-overflow-into-neighbour bug)
#   * program labels from CURRENT program_names.tsv name_short (+ "*")
#   * fonts >= 5pt; dense point layers rasterised (ggrastr) dpi ~400-450
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
  library(circlize); library(ggridges); library(scatterpie)
  library(ggrastr); library(data.table)
  library(dplyr); library(tidyr); library(grid)
  library(scales); library(RColorBrewer); library(svglite)
})

set.seed(42)
mm2in <- function(x) as.numeric(x) / 25.4

## ---- paths ----------------------------------------------------------
WORK <- "CORTEX_PROGRAM_ROOT/scripts/fig2"
RES  <- "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUT  <- "CORTEX_PROGRAM_ROOT/figures/fig2/svg_panels"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
svgf <- function(id) file.path(OUT, sprintf("fig2_%s.svg", id))

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
div_pal <- colorRampPalette(c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7",
                              "#FDDBC7","#EF8A62","#B2182B"))(256)

## ---- Nature theme ---------------------------------------------------
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
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
}
theme_set(theme_nat())
ZM <- theme(plot.margin = margin(0, 0, 0, 0))   # zero-margin overlay

## ---- choices (rep chip + exemplar programs) -------------------------
ch <- readLines(file.path(WORK, "_choices.txt"))
REP <- sub("REP=", "", grep("^REP=", ch, value = TRUE)[1])
ex_lines <- grep("^EX\t", ch, value = TRUE)
EX <- setNames(sub(".*\t.*\t", "", ex_lines), sub("^EX\t", "", sub("\t[^\t]*$", "", ex_lines)))

## ---- old->new cNMF component renumber map (used throughout) ----------
## EXCLUDED old components: 9, 18, 19, 35, 52, 57
OLD2NEW <- c(`1`=1,`2`=2,`3`=3,`4`=4,`5`=5,`6`=6,`7`=7,`8`=8,
             `10`=9,`11`=10,`12`=11,`13`=12,`14`=13,`15`=14,`16`=15,`17`=16,
             `20`=17,`21`=18,`22`=19,`23`=20,`24`=21,`25`=22,`26`=23,`27`=24,
             `28`=25,`29`=26,`30`=27,`31`=28,`32`=29,`33`=30,`34`=31,
             `36`=32,`37`=33,`38`=34,`39`=35,`40`=36,`41`=37,`42`=38,`43`=39,
             `44`=40,`45`=41,`46`=42,`47`=43,`48`=44,`49`=45,`50`=46,`51`=47,
             `53`=48,`54`=49,`55`=50,`56`=51,`58`=52,`59`=53,`60`=54)

## ---- program -> functional GO:BP label lookup (CURRENT program_names.tsv) ---
## new_P col = "P1".."P54"; label = "P{n} {name_short}" + "*" if brain-weak
NM <- fread(file.path(RES, "program_names.tsv"))
NM[, program_int := as.integer(sub("^P", "", new_P))]  # "P1"->1, "P54"->54
.mk_label <- function(n) {
  r <- NM[program_int == as.integer(n)]
  if (nrow(r) == 0) return(paste0("P", n))
  ns   <- trimws(sub("\\s+P[0-9]+$", "", as.character(r$name_short[1])))
  star <- if (as.character(r$confidence[1]) == "brain-weak") "*" else ""
  paste0("P", n, " ", ns, star)
}
PLABEL <- setNames(vapply(NM$program_int, .mk_label, character(1)), as.character(NM$program_int))
plabel <- function(prog_key) {
  n <- sub("program_", "", as.character(prog_key))
  out <- PLABEL[n]; out[is.na(out)] <- paste0("P", n[is.na(out)]); unname(out)
}
## FLAG: old program_35 -> EXCLUDED; removed from EX_CONTEXT. Keys are OLD cNMF nums (match _choices.txt / data).
## Context descriptions unchanged; new_P numbers will appear via d_short/ex_title renaming.
EX_CONTEXT <- c(program_5="L2/3", program_53="L6",
                program_37="OLIGO / WM", program_56="ENDO / vascular",
                program_29="Inhibitory")
## plabel_by_old: convert old cNMF prog_key -> new_P -> plabel (since PLABEL is keyed by new_P)
.plabel_by_old <- function(prog_key) {
  old_n <- as.character(as.integer(sub("program_", "", as.character(prog_key))))
  new_n <- OLD2NEW[old_n]
  plabel(paste0("program_", new_n))
}
ex_title <- function(prog_key, sep = " · ") {
  ctx <- EX_CONTEXT[as.character(prog_key)]
  lab <- .plabel_by_old(prog_key)
  ifelse(is.na(ctx), lab, paste0(ctx, sep, lab))
}

## ---- SHORT panel titles (<=~22 chars, single line) ------------------
TT_A <- "Cortical layers"
TT_B <- "Spatial program activity"
TT_C <- "RCTD composition"
TT_D <- "Depth profile"
TT_E <- "Program × layer"
TT_F <- "z vs RCTD weight"
TT_G <- "Cross-chip reproducibility"
TT_H <- "Distance from WM"

## =====================================================================
## load rep-chip spatial data (shared by a, b)
## =====================================================================
meta <- fread(file.path(WORK, "repchip_meta.tsv"))
prog <- fread(file.path(WORK, "repchip_progscores.tsv"))
df <- merge(meta, prog, by = "bin")
df[, yy := max(y) - y]
cat("rep chip bins for plot:", nrow(df), "\n")

## =====================================================================
## PANEL a -- rep-chip tissue scatter colored by layer (square)
## =====================================================================
df_a <- df %>% mutate(majorDomain = factor(majorDomain, levels = LAYER_ORDER))
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
  coord_fixed(ratio = 1, clip = "off") +
  labs(title = TT_A) +
  theme_void(base_size = 6.6) +
  theme(plot.title = element_text(size = 7.0, face = "bold", hjust = 0,
                                  margin = margin(b = 3)),
        legend.title = element_text(size = 6.2),
        legend.text  = element_text(size = 5.6),
        legend.key.size = unit(2.6, "mm"),
        legend.position = "bottom",
        legend.box.margin = margin(t = 1),
        plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.8)))
svglite(svgf("a"), width = mm2in(42), height = mm2in(41), bg = "white")
print(p_a); invisible(dev.off())
cat("panel a SVG written\n")

## =====================================================================
## PANEL b -- 6-program spatial small-multiples (z), one per cell class.
## SPATIALLY-VALIDATED exemplars (renumbered: old->new):
##   P31 EX [old34] | P21 INH(VIP)† [old24] | P14 AST [old15] | P33 OLIGO [old37] | P51 ENDO [old56] | P36 MICRO [old40]
## P21 caveat: inhibitory programs sparsely resolved at bin50 -> "†" + footnote.
## Arranged 3 rows x 2 cols (square facets via coord_fixed) -> panel ~1:1.
## NOTE: repchip_progscores_smooth.tsv columns use OLD cNMF nums. B_PROGS uses old nums for data lookup.
## =====================================================================
## B_PROGS = old nums (data column names); B_CLASS keyed by old nums for B_PROGS lookup.
B_PROGS   <- c("program_34", "program_24", "program_15",
               "program_37", "program_56", "program_40")
B_CLASS   <- c(program_34 = "Excitatory", program_24 = "Inhibitory (VIP)†",
               program_15 = "Astrocyte",  program_37 = "Oligodendrocyte",
               program_56 = "Endothelial", program_40 = "Microglia")
## plabel_new: lookup by new_P number for display labels on panel b
.plabel_new_b <- function(old_prog_key) {
  old_n <- as.character(as.integer(sub("program_", "", old_prog_key)))
  new_n <- OLD2NEW[old_n]
  plabel(paste0("program_", new_n))
}
b_title <- function(prog_key) {           # "Class\nP{n} name_short" (wrapped); prog_key = old num
  cls <- B_CLASS[as.character(prog_key)]
  ## wrap the (sometimes long) program label so it never clips the strip;
  ## width ~20 chars keeps "Astrocyte glutamate transport" / "Microglial
  ## immune activation" on 2 short lines instead of overflowing.
  lab <- vapply(.plabel_new_b(prog_key),  # use new_P label lookup
                function(s) paste(strwrap(s, width = 20), collapse = "\n"),
                character(1))
  paste0(cls, "\n", lab)
}
exprogs <- B_PROGS
## DISPLAY-LAYER smoothing: use per-chip KNN-MEDIAN-smoothed program-z (k=25
## nearest bins in x,y) from repchip_progscores_smooth.py, NOT the raw per-bin
## scores. Suppresses single-bin speckle while preserving anatomy. Raw parquet
## untouched. Panel a (layer scatter) still uses the un-smoothed shared df.
prog_sm <- fread(file.path(WORK, "repchip_progscores_smooth.tsv"))
df_b <- merge(meta[, .(bin, x, y)], prog_sm, by = "bin")
df_b[, yy := max(y) - y]
b_long <- df_b %>% select(x, yy, all_of(exprogs)) %>%
  pivot_longer(all_of(exprogs), names_to = "program", values_to = "z")
b_vlim <- b_long %>% group_by(program) %>%
  summarise(vlim = quantile(abs(z), 0.98, na.rm = TRUE), .groups = "drop")
b_long <- b_long %>% left_join(b_vlim, by = "program") %>%
  mutate(zc = pmax(pmin(z / vlim, 1), -1),
         prog_label = factor(b_title(program),
                             levels = b_title(exprogs)))
p_b <- ggplot(b_long, aes(x, yy, colour = zc)) +
  ggrastr::rasterise(geom_point(shape = 15, size = 0.14, stroke = 0), dpi = 400) +
  facet_wrap(~ prog_label, nrow = 1, ncol = 6) +
  scale_colour_gradientn(colours = div_pal, limits = c(-1, 1),
                         breaks = c(-1, 0, 1), labels = c("low", "0", "high"),
                         name = "z-activity (per-program)") +
  coord_fixed(ratio = 1, clip = "off") +
  labs(title = TT_B,
       caption = "† inhibitory programs sparsely resolved at bin50") +
  theme_void(base_size = 6.6) +
  theme(plot.title = element_text(size = 7.4, face = "bold", hjust = 0,
                                  margin = margin(b = 3)),
        plot.caption = element_text(size = 5.4, colour = "grey35", hjust = 0,
                                    margin = margin(t = 2)),
        strip.text = element_text(size = 5.3, face = "bold", lineheight = 0.92,
                                  margin = margin(b = 1.4, t = 1)),
        legend.title = element_text(size = 5.8),
        legend.text  = element_text(size = 5.2),
        legend.key.height = unit(2.4, "mm"), legend.key.width = unit(4.0, "mm"),
        legend.position = "bottom",
        legend.box.margin = margin(t = 1),
        panel.spacing = unit(1.2, "mm"),
        plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA)) +
  guides(colour = guide_colourbar(title.position = "top"))
## 1 row x 6 cols (USER FB v6): wide banner; six square tiles side-by-side
## with bottom colourbar legend; height includes 2-line strip + legend.
svglite(svgf("b"), width = mm2in(172), height = mm2in(58), bg = "white")
print(p_b); invisible(dev.off())
cat("panel b SVG written 1x6 banner (renumbered exemplars P31/P21/P14/P33/P51/P36)\n")

## =====================================================================
## PANEL c -- scatterpie of RCTD composition on coarse grid (square)
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
  coord_fixed(ratio = 1, clip = "off") +
  labs(title = TT_C) +
  theme_void(base_size = 6.6) +
  theme(plot.title = element_text(size = 7.0, face = "bold", hjust = 0,
                                  lineheight = 1.0, margin = margin(b = 3)),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "center",
        legend.box.margin = margin(t = 1),
        legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 5.4),
        legend.title = element_text(size = 6.0),
        plot.background = element_rect(fill = "white", colour = NA)) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE,
                             title.position = "top"))
## square tissue + bottom multi-row cell-type legend -> overall panel ~1:1
svglite(svgf("c"), width = mm2in(84), height = mm2in(93), bg = "white")
print(p_c); invisible(dev.off())
cat("panel c SVG written\n")

## =====================================================================
## PANEL d -- program activity along cortical depth (ridgelines, 2x3 wide)
## =====================================================================
pcl <- fread(file.path(WORK, "prog_x_layer_per_chip.tsv"))
## Remove EXCLUDED old programs (35=EXCLUDED) from EX programs before filtering
EXCLUDED_OLD_PROGS <- paste0("program_", c(9, 18, 19, 35, 52, 57))
selprog <- names(EX)[!names(EX) %in% EXCLUDED_OLD_PROGS]
pdd <- pcl[program %in% selprog & n >= 50]
pdd[, majorDomain := factor(majorDomain, levels = rev(LAYER_ORDER))]
## panel d short strip label: context + new P-number (long name clips at 41mm)
d_short <- function(prog_key) {
  ctx <- EX_CONTEXT[as.character(prog_key)]
  old_n <- as.character(as.integer(sub("program_", "", as.character(prog_key))))
  new_n <- OLD2NEW[old_n]
  pn  <- paste0("P", ifelse(is.na(new_n), old_n, new_n))
  ifelse(is.na(ctx), pn, paste0(ctx, " \u00b7 ", pn))
}
pdd[, prog_label := factor(d_short(program), levels = d_short(selprog))]
p_d <- ggplot(pdd, aes(x = mean_z, y = majorDomain, fill = after_stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 1.6, rel_min_height = 0.01,
      linewidth = 0.2, colour = "grey30") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, name = "median z") +
  facet_wrap(~ prog_label, nrow = 2, ncol = 3) +
  labs(title = TT_D,
       x = "median program z-activity", y = "cortical depth") +
  theme_nat() +
  theme(strip.text = element_text(size = 5.8, face = "bold"),
        panel.spacing = unit(1.5, "mm"),
        legend.key.size = unit(3.0, "mm"),
        plot.margin = margin(0, 0, 0, 0))
svglite(svgf("d"), width = mm2in(56), height = mm2in(42), bg = "white")
print(p_d); invisible(dev.off())
cat("panel d SVG written\n")

## =====================================================================
## PANEL e -- program x layer ComplexHeatmap (2-COLUMN, wide banner)
## rendered DIRECTLY to svglite with zero padding (no grabExpr/patchwork).
## =====================================================================
G <- fread(file.path(WORK, "prog_x_layer_global.tsv"))
G <- G[match(LAYER_ORDER, majorDomain)]
progcols_all <- paste0("program_", 1:60)
EXCLUDED_OLD <- c(9, 18, 19, 35, 52, 57)   # old cNMF components, excluded
progcols <- progcols_all[!as.integer(sub("program_","",progcols_all)) %in% EXCLUDED_OLD]
Me <- t(as.matrix(G[, ..progcols]))
## OLD2NEW already defined at top of script
old_nums <- as.integer(sub("program_","",progcols))
rownames(Me) <- paste0("P", OLD2NEW[as.character(old_nums)])
colnames(Me) <- G$majorDomain
col_fun <- colorRamp2(c(-1.2, -0.4, 0, 0.4, 1.2),
   c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"))
cc <- fread(file.path(WORK, "program_celltype_corr.tsv"))
setnames(cc, 1, "program")
peakct <- apply(as.matrix(cc[, ..CELLTYPES]), 1, function(r) CELLTYPES[which.max(r)])
names(peakct) <- cc$program
## FLAG: old P35 -> EXCLUDED, removed. New numbers: 5->5, 29->26, 37->33, 53->48, 56->51, 26->23
SPOTLIGHT <- c(5, 26, 33, 48, 51, 23)   # new_P (old: 5,29,37,53,56,26); old35=EXCLUDED removed
pnum <- OLD2NEW[as.character(as.integer(sub("program_","",progcols)))]
## routeD (rule 8 dense-axis): SPOTLIGHT rows keep full functional + peak-cell
## label (bold); the other ~53 rows keep ONLY their P-number (drop the long
## "· celltype" suffix) so labels are short and legible at >=5pt. No row dropped.
## plabel uses new_P numbers (PLABEL keyed by new_P int), so pass paste0("program_", pnum).
rowlab_all <- ifelse(pnum %in% SPOTLIGHT,
                 paste0(plabel(paste0("program_", pnum)), " · ", peakct[progcols]),
                 rownames(Me))
names(rowlab_all) <- rownames(Me)
face_all <- ifelse(pnum %in% SPOTLIGHT, "bold", "plain"); names(face_all) <- rownames(Me)
size_all <- ifelse(pnum %in% SPOTLIGHT, 5.8, 5.2);        names(size_all) <- rownames(Me)
hc   <- hclust(dist(Me, method = "euclidean"), method = "ward.D2")
ord  <- hc$order
half <- ceiling(length(ord) / 2)
idx1 <- ord[1:half]; idx2 <- ord[(half + 1):length(ord)]
rn1  <- rownames(Me)[idx1]; rn2 <- rownames(Me)[idx2]
W <- cbind(Me[idx1, , drop = FALSE], Me[idx2, , drop = FALSE])
rownames(W) <- NULL
col_block <- factor(rep(c("b1", "b2"), each = length(LAYER_ORDER)),
                    levels = c("b1", "b2"))
top_e <- HeatmapAnnotation(layer = rep(LAYER_ORDER, 2),
   col = list(layer = LAYER_COLORS), show_legend = FALSE,
   annotation_name_gp = gpar(fontsize = 5.6), simple_anno_size = unit(2.2, "mm"))
gp1 <- gpar(fontsize = size_all[rn1], fontface = face_all[rn1])
gp2 <- gpar(fontsize = size_all[rn2], fontface = face_all[rn2])
ra_left  <- rowAnnotation(b1 = anno_text(rowlab_all[rn1], gp = gp1,
                                         location = 1, just = "right"))
ra_right <- rowAnnotation(b2 = anno_text(rowlab_all[rn2], gp = gp2,
                                         location = 0, just = "left"))
ht_e <- Heatmap(W, name = "median z", col = col_fun,
   cluster_columns = FALSE, cluster_rows = FALSE,
   row_order = seq_len(nrow(W)), column_order = seq_len(ncol(W)),
   column_split = col_block, column_gap = unit(4, "mm"),
   column_labels = rep(LAYER_ORDER, 2),
   column_names_gp = gpar(fontsize = 5.6), column_names_rot = 45,
   top_annotation = top_e, left_annotation = ra_left, right_annotation = ra_right,
   show_heatmap_legend = TRUE,
   heatmap_legend_param = list(title_gp = gpar(fontsize = 6.5),
       labels_gp = gpar(fontsize = 5.6), legend_height = unit(15, "mm"),
       grid_width = unit(2.5, "mm")),
   ## routeD (skill rule 1): SQUARE cells. 16 cols x 30 rows at 2.4mm each ->
   ## wide-short banner of square cells (was 1.93x2.47 non-square stretched-tall).
   width = unit(ncol(W) * 2.4, "mm"), height = unit(nrow(W) * 2.4, "mm"),
   column_title = c(paste0(TT_E, "  (block 1)"), "(block 2)"),
   column_title_gp = gpar(fontsize = 6.4))
## ===== USER FB v6: MERGE h (distance-from-WM) INTO e =================
## Build the WM-distance line subplot HERE (formerly panel h) and stack it
## BELOW the heatmap in a single composite SVG -> e_combined. h dropped.
rings <- fread(file.path(WORK, "panelh_rings.tsv"))
rings[, ring := as.numeric(as.character(ring))]
## panelh_rings.tsv columns use OLD cNMF nums; wmp uses old nums for data lookup.
## old->new: 37->33, 26->23, 5->5, 56->51
wmp <- c("program_37", "program_26", "program_5", "program_56")  # old nums (data columns)
wmp_new <- paste0("program_", c(OLD2NEW["37"], OLD2NEW["26"], OLD2NEW["5"], OLD2NEW["56"]))  # new nums for labels
.role <- c(program_37 = "OLIGO/WM", program_26 = "OLIGO/WM",
           program_5  = "L2/3 ctrl", program_56 = "vascular")
plab_h <- setNames(plabel(wmp_new), wmp)  # keyed by old prog name, value = new-P label
longh <- melt(rings, id.vars = c("ring", "chip"), measure.vars = wmp,
              variable.name = "program", value.name = "mean_z")
longh[, prog_label := factor(plab_h[as.character(program)], levels = plab_h[wmp])]
prog_pal <- setNames(c("#7B3294", "#C2A5CF", "#5A78D6", "#D6604D"), plab_h[wmp])
p_wm <- ggplot(longh, aes(x = ring, y = mean_z, colour = prog_label,
                         group = interaction(prog_label, chip))) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_line(aes(linetype = chip), linewidth = 0.45) +
  geom_point(size = 0.4) +
  scale_color_manual(values = prog_pal, name = "program") +
  scale_linetype_manual(values = c("solid", "31"), name = "chip") +
  labs(title = "Distance from WM (annotation)",
       x = "distance from WM (µm)", y = "median z") +
  theme_nat() +
  theme(plot.title = element_text(size = 6.4, face = "bold"),
        axis.title = element_text(size = 5.6),
        legend.position = "right", legend.box = "vertical",
        legend.key.size = unit(2.2, "mm"),
        legend.text = element_text(size = 5.4),
        legend.title = element_text(size = 5.6),
        plot.margin = margin(2, 0, 0, 0))

## render heatmap to a grid grob via grid.grabExpr so it can stack with ggplot.
ht_grob <- grid::grid.grabExpr(
  draw(ht_e, heatmap_legend_side = "right",
       padding = unit(c(0, 0, 0, 0), "mm"))
)
## stack heatmap (taller) + WM-distance line plot (slim) = combined panel e.
p_e_combined <- patchwork::wrap_plots(
  patchwork::wrap_elements(full = ht_grob),
  p_wm,
  ncol = 1, heights = c(3.0, 1)
) + plot_annotation(theme = theme(plot.margin = margin(0,0,0,0)))
svglite(svgf("e"), width = mm2in(87), height = mm2in(74), bg = "white")
print(p_e_combined); invisible(dev.off())
cat("panel e SVG written (HEATMAP + WM-distance merged; h dropped)\n")

## =====================================================================
## PANEL f -- program-z vs RCTD-weight hexbin, 3 exemplar facets (wide row)
## square panels via coord_fixed, rasterised hex; standalone patchwork row.
## =====================================================================
fdat <- fread(file.path(WORK, "panelf_pairs.tsv"))
PAIRS <- list(
  ## panelf_pairs.tsv columns use OLD cNMF nums; pp = old num for data lookup.
  ## st uses plabel with new_P (via OLD2NEW): old37->new33, old56->new51, old5->5
  list(pp = "program_37", cc = "OLIGO",               xl = "OLIGO weight",  st = plabel(paste0("program_", OLD2NEW["37"]))),
  list(pp = "program_56", cc = "ENDO",                 xl = "ENDO weight",  st = plabel(paste0("program_", OLD2NEW["56"]))),
  list(pp = "program_5",  cc = "L2-L3 IT LINC00507",   xl = "L2/3 IT weight", st = plabel(paste0("program_", OLD2NEW["5"]))))
f_list <- lapply(PAIRS, function(P) {
  x <- fdat[[P$cc]]; y <- fdat[[P$pp]]
  ok <- is.finite(x) & is.finite(y)
  data.frame(wt = x[ok], z = y[ok], prog = P$st)  # P$st already holds new-P label
})
mk_hex <- function(d, P) {
  rp <- cor(d$wt, d$z, method = "pearson")          # Pearson r
  rs <- cor(d$wt, d$z, method = "spearman")         # Spearman rho
  rng_x <- range(d$wt); rng_y <- range(d$z)
  asp <- diff(rng_x) / diff(rng_y)
  lab <- sprintf("Pearson r = %.2f\nSpearman rs = %.2f", rp, rs)
  ggplot(d, aes(wt, z)) +
    ggrastr::rasterise(geom_hex(bins = 42), dpi = 400) +
    scale_fill_gradient(low = "#fcdbb8", high = "#5b1a6b",
                        trans = "log10", name = "log10 (bins)") +
    ## trend line (NOT identity): linear fit of program-z on RCTD weight.
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "#111111", fill = "grey75", linewidth = 0.5,
                alpha = 0.35) +
    ## routeD: enlarge correlation annotation so it renders >=5pt after layout
    ## scale (geom text size is in mm; 2.0mm~=5.7pt body). Keep both r values.
    annotate("text", x = -Inf, y = Inf, label = lab,
             hjust = -0.06, vjust = 1.2, size = 2.0, lineheight = 0.95) +
    coord_fixed(ratio = asp, clip = "off") +
    labs(x = P$xl, y = sprintf("%s z", P$st),  # P$st holds new-P label
         title = P$st) +
    theme_nat() +
    theme(plot.title = element_text(size = 6.4, face = "bold"),
          axis.title = element_text(size = 6.0),
          axis.text  = element_text(size = 5.4),
          legend.key.size = unit(2.6, "mm"),
          legend.position = "bottom",
          legend.title = element_text(size = 5.6),
          legend.text  = element_text(size = 5.2),
          plot.margin = margin(0, 0, 0, 0))
}
pf_list <- mapply(function(d, P) mk_hex(d, P), f_list, PAIRS, SIMPLIFY = FALSE)
## routeD: 1x3 wide row (was 2x2). Wider, shorter panel (AR up ~2.9) frees page
## height for the dense panels; shared bottom legend; each hexbin still square.
p_f <- (pf_list[[1]] + pf_list[[2]] + pf_list[[3]]) +
  plot_layout(ncol = 3, nrow = 1, guides = "collect") +
  plot_annotation(title = TT_F,
                  theme = theme(plot.title = element_text(size = 7.4, face = "bold",
                                                          family = "Helvetica"),
                                plot.margin = margin(0, 0, 0, 0))) &
  theme(plot.margin = margin(0, 0, 0, 0), legend.position = "bottom")
svglite(svgf("f"), width = mm2in(91), height = mm2in(40), bg = "white")
print(p_f); invisible(dev.off())
cat("panel f SVG written (routeD 1x3 wide)\n")

## =====================================================================
## PANEL g -- cross-chip reproducibility boxplot (very wide)
## =====================================================================
gsum  <- fread(file.path(WORK, "panelg_summary.tsv"))
grep_ <- fread(file.path(WORK, "panelg_reproducibility.tsv"))
gord  <- gsum[order(median)]$program
grep_[, program := factor(program, levels = gord)]
medmap <- setNames(gsum$median, gsum$program)
grep_[, medv := medmap[as.character(program)]]
## VERTICAL orientation: programs on y-axis (~56 rows -> tall), r on x.
## routeD (rule 8 dense-axis): with ~56 boxes a per-program y label is illegible
## (<3pt). Label ONLY the SPOTLIGHT/highlighted programs + the top/bottom few by
## median r; the rest keep their P-number tick so the axis stays traceable but
## the labels enlarge to >=5pt. No data dropped — every program still plotted.
## FLAG: old P35 -> EXCLUDED, removed. SPOTLIGHT_G_OLD for data matching (old cNMF nums in panelg_summary.tsv).
## old->new: 5->5, 29->26, 37->33, 53->48, 56->51, 26->23, 34->31, 15->14, 40->36
SPOTLIGHT_G_OLD <- c(5, 29, 37, 53, 56, 26, 34, 15, 40)  # old_P (35=EXCLUDED removed); match panelg_summary.tsv
SPOTLIGHT_G <- sapply(as.character(SPOTLIGHT_G_OLD), function(x) OLD2NEW[x])  # new_P for display
gmed <- gsum[order(median)]
n_g  <- nrow(gmed)
extremes <- unique(c(head(gmed$program, 3), tail(gmed$program, 3)))  # 3 worst + 3 best
hl_prog  <- union(paste0("program_", SPOTLIGHT_G_OLD), extremes)
g_lab <- function(x) {
  old_n <- as.character(as.integer(sub("program_", "", x)))
  new_n <- OLD2NEW[old_n]
  pn <- ifelse(is.na(new_n), sub("program_", "P", x), paste0("P", new_n))
  ifelse(x %in% hl_prog, pn, "")          # blank non-highlighted (tick kept)
}
g_face <- ifelse(gord %in% hl_prog, "bold", "plain")
p_g <- ggplot(grep_, aes(x = program, y = corr_to_mean)) +
  geom_hline(yintercept = c(0, 0.5, 1), colour = "grey85", linewidth = 0.25) +
  geom_boxplot(aes(fill = medv), outlier.size = 0.15, outlier.colour = "grey60",
               linewidth = 0.13, width = 0.72) +
  scale_fill_gradientn(colours = c("#D73027", "#FEE08B", "#1A9850"),
      limits = c(min(grep_$medv), 1), name = "median r") +
  scale_x_discrete(labels = g_lab) +
  labs(title = TT_G,
       x = "program (sorted by median; key programs labelled)",
       y = "correlation to mean layer profile") +
  coord_flip(ylim = c(min(0, quantile(grep_$corr_to_mean, 0.01)), 1)) +
  theme_nat() +
  theme(axis.text.y = element_text(size = 5.4, face = g_face),
        axis.ticks.y = element_line(linewidth = 0.2),
        legend.key.size = unit(3.0, "mm"),
        legend.position = "right",
        plot.margin = margin(0, 0, 0, 0))
svglite(svgf("g"), width = mm2in(78), height = mm2in(89), bg = "white")
print(p_g); invisible(dev.off())
cat("panel g SVG written (routeD: highlight-only y labels)\n")

## =====================================================================
## PANEL h -- DROPPED (USER FB v6): WM-distance now embedded inside panel e
## as a stacked annotation row. h.svg deliberately NOT written.
## =====================================================================
cat("panel h SKIPPED (merged into e per user FB v6)\n")

## =====================================================================
## PANEL i -- spatial-bin UMAP in EXPRESSION space + Harmony (REPLACES the
## old cNMF program-score UMAP). Each point = one bin50 spatial bin embedded
## by its RAW SCT counts -> Pearson residuals -> HVG -> PCA -> Harmony
## integrate by chip -> scanpy neighbors+umap (seed=0).
## Main map coloured by majorDomain (anatomy: GM layers / WM / arachnoid);
## THREE program-feature insets: P33 OLIGO/WM (old P37, expect WM), P31 excitatory (old P34),
## P7 inhibitory. Embedding precomputed by stage_spatial_umap_expr.py.
## Square via coord_fixed; fonts >= 5pt; hexbin-aggregated (no speckle).
## =====================================================================
UM <- fread(file.path("CORTEX_PROGRAM_ROOT/figures/fig2/_intermediate",
                      "spatial_umap_expr_harmony.csv"))
UM[, majorDomain := factor(majorDomain, levels = LAYER_ORDER)]
asp_um <- diff(range(UM$UMAP1)) / diff(range(UM$UMAP2))

## insets need the program-z columns; the embedding CSV lacks them, so join
## from a small companion CSV (bin + program_37/34/7 for the embedding bins),
## precomputed by scripts/fig2/make_umap_feat.py (R 'arrow' pkg unavailable).
## routeC (option C, 羅老師裁決): RESTORE all 3 feature insets (no content loss).
## P33 OLIGO/WM (old P37) + P31 excitatory (old P34) + P7 inhibitory. Height is
## relaxed (<=290mm page) so every inset still keeps a >=5pt box at 3-up.
## FEAT_I_OLD = old cNMF component numbers (column names in spatial_umap_feat.csv are old nums).
FEAT_I_OLD <- c(37, 34, 7)    # old cNMF component numbers (data columns in spatial_umap_feat.csv)
FEAT_I_NEW <- sapply(as.character(FEAT_I_OLD), function(x) OLD2NEW[x])  # new_P: 33,31,7
have_feat <- all(paste0("program_", FEAT_I_OLD) %in% names(UM))
if (!have_feat) {
  featf <- file.path("CORTEX_PROGRAM_ROOT/figures/fig2/_intermediate",
                     "spatial_umap_feat.csv")
  pj <- fread(featf)
  UM <- merge(UM, pj, by = "bin", all.x = TRUE, sort = FALSE)
}
## alias FEAT_I for legacy compatibility (used in lapply below)
FEAT_I <- FEAT_I_OLD

## ---- HEXBIN aggregation (replaces speckly per-bin scatter) ----------
## Anatomy map: each hex coloured by the MAJORITY majorDomain of bins in it.
## Feature insets: each hex coloured by the MEDIAN program-z of bins in it.
## hexbin -> per-point cell id; aggregate mode (domain) / median (z) per cell;
## draw filled hexagons via geom_hex(stat="identity"). coord_fixed kept.
NBINS_HEX <- 60                       # ~60 hexes across the longer UMAP axis
xr_um <- range(UM$UMAP1); yr_um <- range(UM$UMAP2)
hb <- hexbin::hexbin(UM$UMAP1, UM$UMAP2, xbins = NBINS_HEX,
                     xbnds = xr_um, ybnds = yr_um, IDs = TRUE)
hcell  <- hb@cID                      # per-point hex cell id
hcent  <- hexbin::hcell2xy(hb)            # cell centroids (x,y) by @cell order
cellxy <- data.table(cell = hb@cell, hx = hcent$x, hy = hcent$y)
UM[, hex := hcell]

.mode_domain <- function(v) {
  v <- v[!is.na(v)]; if (!length(v)) return(NA_character_)
  names(sort(table(v), decreasing = TRUE))[1]
}
dom_hex <- UM[, .(majorDomain = .mode_domain(as.character(majorDomain)), n = .N),
              by = .(cell = hex)]
dom_hex <- merge(dom_hex, cellxy, by = "cell")
dom_hex[, majorDomain := factor(majorDomain, levels = LAYER_ORDER)]

p_i_main <- ggplot(dom_hex, aes(hx, hy, fill = majorDomain)) +
  geom_hex(stat = "identity", colour = NA) +
  scale_fill_manual(values = LAYER_COLORS, name = "majorDomain", drop = FALSE,
                    guide = guide_legend(override.aes = list(size = 1.6),
                                         nrow = 2, byrow = TRUE)) +
  coord_fixed(ratio = asp_um, clip = "off") +
  labs(title = "Spatial-bin UMAP (expression + Harmony)", x = "UMAP1", y = "UMAP2") +
  theme_nat() +
  theme(plot.title = element_text(size = 7.0, face = "bold", hjust = 0),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom", legend.box.margin = margin(t = 1),
        legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 5.2),
        legend.title = element_text(size = 5.8),
        plot.margin = margin(0, 0, 0, 0))

## mk_feat: pnum_old = old cNMF component number (column in spatial_umap_feat.csv)
## pnum_new = new renumbered P number (for label display via PLABEL)
mk_feat <- function(pnum_old) {
  pnum_new <- OLD2NEW[as.character(pnum_old)]
  col_name <- paste0("program_", pnum_old)
  z <- UM[[col_name]]
  if (is.null(z)) return(NULL)
  vl <- quantile(abs(z), 0.98, na.rm = TRUE)
  zhex <- UM[, .(z = median(get(col_name), na.rm = TRUE)),
             by = .(cell = hex)]
  zhex <- merge(zhex, cellxy, by = "cell")
  zhex[, zc := pmax(pmin(z / vl, 1), -1)]
  ggplot(zhex, aes(hx, hy, fill = zc)) +
    geom_hex(stat = "identity", colour = NA) +
    scale_fill_gradientn(colours = div_pal, limits = c(-1, 1),
                         breaks = c(-1, 0, 1), labels = c("lo", "0", "hi"),
                         name = "median z") +
    coord_fixed(ratio = asp_um, clip = "off") +
    ## single-source title (name_short via plabel using new_P number); wrap to >=2 lines so the
    ## 3 side-by-side insets don't collide (htitlefix: legibility only, no
    ## content/name change). width chosen so each line <= one inset width.
    labs(title = paste(strwrap(plabel(paste0("program_", pnum_new)), width = 22),
                       collapse = "\n")) +
    theme_nat() +
    theme(plot.title = element_text(size = 5.2, face = "bold", hjust = 0.5,
                                    lineheight = 0.85, margin = margin(b = 0.6)),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(2.0, "mm"), legend.key.width = unit(3.0, "mm"),
          legend.position = "bottom", legend.title = element_text(size = 5.4),
          legend.text = element_text(size = 5.4),
          plot.margin = margin(0, 0, 0, 0))
}
feat_list <- Filter(Negate(is.null), lapply(FEAT_I_OLD, mk_feat))

## layout: big anatomy map on top, THREE feature insets in a row below
## (P33 OLIGO/WM [old P37] | P31 excitatory [old P34] | P7 inhibitory)  -- routeC: 3 insets restored
p_i <- p_i_main / (feat_list[[1]] | feat_list[[2]] | feat_list[[3]]) +
  patchwork::plot_layout(heights = c(1.7, 1.06)) +   # htitlefix: +space for 2-line inset titles
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
svglite(svgf("i"), width = mm2in(88), height = mm2in(96), bg = "white")  # htitlefix: 93->96 for wrapped titles
print(p_i); invisible(dev.off())
cat("panel i SVG written (routeC: spatial-bin UMAP + P33/P31/P7 insets [renumbered from P37/P34/P7], 3 restored)\n")

cat("\nALL 8 STANDALONE SVG PANELS COMPLETE (h dropped, merged into e) ->", OUT, "\n")
