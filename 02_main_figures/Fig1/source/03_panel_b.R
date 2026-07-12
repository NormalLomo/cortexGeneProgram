#!/usr/bin/env Rscript
# =====================================================================
# Fig.1 panel b -- MERGED "program backbone" panel.
# Fuses old panels b (program x subclass identity heatmap) + d (top
# genes) + e (GO:BP NES) into ONE ComplexHeatmap keyed on the 60 cNMF
# programs, rows ordered by ward.D2 hierarchical clustering with the row
# dendrogram kept. Right rowAnnotations carry the d/e content.
# Standalone -- does NOT source fig1_svg_panels.R. Only writes fig1_b.svg.
#
# Map 60 cNMF components to the retained 54-program identifiers.
# =====================================================================
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid)
  library(svglite); library(dplyr)
})
set.seed(42)

PROJ <- "CORTEX_PROGRAM_ROOT"
RES  <- file.path(PROJ, "results/crossregion_v1")
INT  <- file.path(PROJ, "figures/fig1/_intermediate")
SVGD <- file.path(PROJ, "figures/fig1/svg_panels")
mm2in <- function(x) x / 25.4

# ---- owner-specified palette ---------------------------------------
CLASS_LV  <- c("excitatory", "inhibitory", "non-neuronal")
CLASS_COL <- c(excitatory = "#2166AC", inhibitory = "#B2182B",
               `non-neuronal` = "#1B7837")

# ---- data: subclass -> class map -----------------------------------
scc <- read.csv(file.path(INT, "subclass_class.csv"), stringsAsFactors = FALSE)
class_map <- setNames(scc$class, scc$subclass)

# ---- data: program specificity (dominant_class) --------------------
spec <- read.csv(file.path(INT, "program_specificity.csv"))
spec$program <- as.integer(spec$program)
dom_cls <- setNames(spec$dominant_class, spec$program)

# ---- OLD2NEW mapping (cNMF component old -> new_P number) ----------
# 6 EXCLUDED: old 9,18,19,35,52,57
EXCLUDED_OLD <- c(9L, 18L, 19L, 35L, 52L, 57L)
OLD2NEW <- c(
  `1`=1, `2`=2, `3`=3, `4`=4, `5`=5, `6`=6, `7`=7, `8`=8,
  `10`=9, `11`=10, `12`=11, `13`=12, `14`=13, `15`=14, `16`=15, `17`=16,
  `20`=17, `21`=18, `22`=19, `23`=20, `24`=21, `25`=22, `26`=23, `27`=24,
  `28`=25, `29`=26, `30`=27, `31`=28, `32`=29, `33`=30, `34`=31,
  `36`=32, `37`=33, `38`=34, `39`=35, `40`=36, `41`=37, `42`=38, `43`=39,
  `44`=40, `45`=41, `46`=42, `47`=43, `48`=44, `49`=45, `50`=46, `51`=47,
  `53`=48, `54`=49, `55`=50, `56`=51, `58`=52, `59`=53, `60`=54
)

# ---- data: program functional names + NES + fdr --------------------
# program_names.tsv cols: new_P, cnmf_component, name_full, name_short, confidence, ...
pn <- read.delim(file.path(RES, "program_names.tsv"), quote = "",
                 comment.char = "", stringsAsFactors = FALSE)
# Filter EXCLUDED (new_P == "EXCLUDED")
pn <- pn[pn$new_P != "EXCLUDED", ]
# cnmf_component = old cNMF integer (key for data lookups)
pn$program <- as.integer(pn$cnmf_component)
# new_P integer for display
pn$new_p_int <- as.integer(sub("^P", "", pn$new_P))
# row label = "P{new_n} {name_short}"  (+ '*' for brain-weak programs)
.BRAIN_WEAK_NEW <- as.integer(pn$new_p_int[pn$confidence == "brain-weak"])
trim_prog_name <- function(x, max_n = 24L) {
  ifelse(nchar(x) <= max_n, x, paste0(substr(x, 1L, max_n - 3L), "..."))
}
mk_prog_id <- function(p_old) {
  new_n <- OLD2NEW[as.character(p_old)]
  star <- ifelse(new_n %in% .BRAIN_WEAK_NEW, "*", "")
  sprintf("P%d%s", new_n, star)
}
mk_prog_name <- function(p_old) {
  i <- match(p_old, pn$program)
  nm <- pn$name_short[i]
  trim_prog_name(nm)
}
nes_of <- function(p) pn$brain_term_NES[match(p, pn$program)]

# ---- data: top-3 gene symbols per program --------------------------
top <- read.csv(file.path(INT, "top100_genes_per_program.csv"),
                stringsAsFactors = FALSE)
top$program <- as.integer(top$program)
top3 <- top %>% filter(rank <= 3) %>% arrange(program, rank) %>%
  group_by(program) %>%
  summarise(genes = paste(gene, collapse = ", "),
            mean_top3_loading = mean(loading), .groups = "drop")
gene_of <- function(p) top3$genes[match(p, top3$program)]
top3load_of <- function(p) top3$mean_top3_loading[match(p, top3$program)]

# ---- MAIN MATRIX: program x subclass identity ----------------------
m <- read.csv(file.path(INT, "program_x_subclass_mean.csv"),
              check.names = FALSE)
progs <- as.integer(m$program)
mat <- as.matrix(m[, -1])
rownames(mat) <- as.character(progs)          # temp rownames = program id
# ---- DROP all 6 EXCLUDED programs (old P9,18,19,35,52,57) ----------
.keep <- !(progs %in% EXCLUDED_OLD)
progs <- progs[.keep]
mat   <- mat[.keep, , drop = FALSE]            # row-subset matrix consistently
# row-z-score (matches old panel b), clamp +-3
z <- t(scale(t(mat))); z[is.na(z)] <- 0
z[z > 3] <- 3; z[z < -3] <- -3

# column order: by class block then subclass name (keep major-class grouping)
col_cls <- class_map[colnames(z)]
col_ord <- order(factor(col_cls, levels = CLASS_LV), colnames(z))
z <- z[, col_ord]; col_cls <- col_cls[col_ord]

# functional row names (replace the program-id rownames with new_P labels)
row_ids  <- vapply(progs, mk_prog_id, character(1))
row_prog <- vapply(progs, mk_prog_name, character(1))
rownames(z) <- paste(row_ids, row_prog)

# per-program (row-aligned) vectors for annotations
row_nes  <- nes_of(progs)
row_gene <- gene_of(progs)
row_t3ld <- top3load_of(progs)     # mean cNMF loading of each program's top-3 genes
row_dom  <- dom_cls[as.character(progs)]
row_dom[is.na(row_dom)] <- "non-neuronal"
row_dom <- factor(row_dom, levels = CLASS_LV)

# ---- top class color strip -----------------------------------------
ha_col <- HeatmapAnnotation(
  class = col_cls,
  col = list(class = CLASS_COL),
  show_annotation_name = FALSE,
  annotation_legend_param = list(class = list(
      title = "class", title_gp = gpar(fontsize = 5.6),
      labels_gp = gpar(fontsize = 5.2))),
  simple_anno_size = unit(2.4, "mm"))

ha_row_left <- rowAnnotation(
  class = row_dom,
  col = list(class = CLASS_COL),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = unit(2.2, "mm"),
  border = TRUE)

# ---- RIGHT row annotations -----------------------------------------
ha_row <- rowAnnotation(
  # (a) GO:BP NES barplot  (old panel e)
  NES = anno_barplot(row_nes, gp = gpar(fill = CLASS_COL[as.character(row_dom)], col = NA),
                     border = FALSE, baseline = 0, bar_width = 0.84,
                     width = unit(11.5, "mm"),
                     axis_param = list(at = c(0, 1, 2), labels = c("0", "1", "2"),
                                       gp = gpar(fontsize = 4.2))),
  # (d) mean cNMF loading of the program's top-3 genes (mako/magma ramp)
  top3load = row_t3ld,
  col = list(top3load = circlize::colorRamp2(
               quantile(row_t3ld, c(0, .5, 1), na.rm = TRUE),
               c("#FCFDBF", "#B5367A", "#1D1147"))),
  # (b) top-3 gene symbols  (old panel d content)
  genes = anno_text(row_gene, gp = gpar(fontsize = 4.6, fontface = "italic"),
                    just = "left", location = unit(0, "npc")),
  annotation_label = c("NES", "load", "genes"),
  annotation_name_gp = gpar(fontsize = 4.4, fontface = "bold"),
  annotation_name_rot = 0,
  annotation_name_side = "top",
  simple_anno_size = unit(1.8, "mm"),
  gap = unit(1.3, "mm"),
  annotation_legend_param = list(
    top3load = list(
      title = "mean top-3 loading", title_gp = gpar(fontsize = 5.6),
      labels_gp = gpar(fontsize = 5.2))),
  show_legend = c(NES = FALSE, top3load = TRUE, genes = FALSE))

# ---- HEATMAP -------------------------------------------------------
B_CELL <- 2.0   # mm per column cell
R_CELL <- 2.6   # mm per program row (taller so 54 rows + anno_text legible)
ht_b <- Heatmap(
  z, name = "row z",
  col = colorRamp2(c(-3, 0, 3), c("#3B4CC0", "#F7F7F7", "#B40426")),
  width  = unit(ncol(z) * B_CELL, "mm"),
  height = unit(nrow(z) * R_CELL, "mm"),
  cluster_rows = TRUE, clustering_method_rows = "ward.D2",
  row_split = row_dom,
  row_gap = unit(1.0, "mm"),
  cluster_row_slices = FALSE,
  show_row_dend = TRUE, row_dend_width = unit(7, "mm"),
  cluster_columns = FALSE,
  column_split = factor(col_cls, levels = CLASS_LV),
  column_gap = unit(1.4, "mm"),
  row_names_gp = gpar(fontsize = 4.8),
  row_names_side = "left",
  row_title_gp = gpar(fontsize = 5.2, fontface = "bold"),
  row_title_rot = 0,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 5.2),
  column_names_rot = 90, column_names_side = "bottom",
  column_title_gp = gpar(fontsize = 6, fontface = "bold"),
  top_annotation = ha_col,
  left_annotation = ha_row_left,
  right_annotation = ha_row,
  heatmap_legend_param = list(title = "row z", title_gp = gpar(fontsize = 5.6),
      labels_gp = gpar(fontsize = 5.2), legend_height = unit(12, "mm")))

# ---- RENDER zero-margin to SVG -------------------------------------
tmp <- file.path(SVGD, "_fig1_b_merged_tmp.svg")
# device generously sized; ink-crop trims to true drawn extent afterwards
svglite(tmp, width = mm2in(190), height = mm2in(185), bg = "white")
draw(ht_b,
     heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legend = TRUE,
     column_title = "Program backbone",
     column_title_gp = gpar(fontsize = 7.2, fontface = "bold"),
     padding = unit(c(0, 0, 0, 0), "mm"))
dev.off()

# overwrite ONLY fig1_b.svg
fb <- file.path(SVGD, "fig1_b.svg")
file.copy(tmp, fb, overwrite = TRUE)
file.remove(tmp)
cat(sprintf("merged panel b SVG -> %s (%d programs)\n", fb, nrow(z)))
