.annotation_group <- function(value) ifelse(endsWith(as.character(value), "sig"), "class_a", "class_b")
PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
    library(svglite)
    library(dplyr)
})
set.seed(42)
PROJ <- PROJECT_ROOT
RES <- file.path(PROJ, "results/crossregion_v1")
INT <- file.path(PROJ, "figures/fig1/_intermediate")
SVGD <- file.path(PROJ, "figures/fig1/svg_panels")
dir.create(SVGD, recursive = TRUE, showWarnings = FALSE)
mm2in <- function(x) x/25.4
CLASS_LV <- c("excitatory", "inhibitory", "non-neuronal")
CLASS_COL <- c(excitatory = "#2166AC", inhibitory = "#B2182B", `non-neuronal` = "#1B7837")
scc <- read.csv(file.path(INT, "subclass_class.csv"), stringsAsFactors = FALSE)
class_map <- setNames(scc$class, scc$subclass)
spec <- read.csv(file.path(INT, "program_specificity.csv"))
spec$program <- as.integer(spec$program)
dom_cls <- setNames(spec$dominant_class, spec$program)
pn <- read.delim(file.path(RES, "program_names.tsv"), quote = "", comment.char = "", stringsAsFactors = FALSE)
pn <- pn[pn$new_P != "EXCLUDED", ]
pn$program <- as.integer(pn$cnmf_component)
pn$new_p_int <- as.integer(sub("^P", "", pn$new_P))
old_to_new <- setNames(pn$new_p_int, as.character(pn$program))
.CLASS_B_PROGRAMS <- as.integer(pn$new_p_int[.annotation_group(pn$confidence) == "class_b"])
trim_prog_name <- function(x, max_n = 24L) {
    ifelse(nchar(x) <= max_n, x, paste0(substr(x, 1L, max_n - 3L), "..."))
}
mk_prog_id <- function(p_old) {
    new_n <- old_to_new[as.character(p_old)]
    star <- ifelse(new_n %in% .CLASS_B_PROGRAMS, "*", "")
    sprintf("P%d%s", new_n, star)
}
mk_prog_name <- function(p_old) {
    i <- match(p_old, pn$program)
    nm <- pn$name_short[i]
    trim_prog_name(nm)
}
nes_of <- function(p) pn$brain_term_NES[match(p, pn$program)]
top <- read.csv(file.path(INT, "top100_genes_per_program.csv"), stringsAsFactors = FALSE)
top$program <- as.integer(top$program)
top3 <- top %>% filter(rank <= 3) %>% arrange(program, rank) %>% group_by(program) %>% summarise(genes = paste(gene, collapse = ", "), mean_top3_loading = mean(loading), .groups = "drop")
gene_of <- function(p) top3$genes[match(p, top3$program)]
top3load_of <- function(p) top3$mean_top3_loading[match(p, top3$program)]
m <- read.csv(file.path(INT, "program_x_subclass_mean.csv"), check.names = FALSE)
progs <- as.integer(m$program)
mat <- as.matrix(m[, -1])
rownames(mat) <- as.character(progs)
.keep <- progs %in% pn$program
progs <- progs[.keep]
mat <- mat[.keep, , drop = FALSE]
z <- t(scale(t(mat)))
z[is.na(z)] <- 0
z[z > 3] <- 3
z[z < -3] <- -3
col_cls <- class_map[colnames(z)]
col_ord <- order(factor(col_cls, levels = CLASS_LV), colnames(z))
z <- z[, col_ord]
col_cls <- col_cls[col_ord]
row_ids <- vapply(progs, mk_prog_id, character(1))
row_prog <- vapply(progs, mk_prog_name, character(1))
rownames(z) <- paste(row_ids, row_prog)
row_nes <- nes_of(progs)
row_gene <- gene_of(progs)
row_t3ld <- top3load_of(progs)
row_dom <- dom_cls[as.character(progs)]
row_dom[is.na(row_dom)] <- "non-neuronal"
row_dom <- factor(row_dom, levels = CLASS_LV)
ha_col <- HeatmapAnnotation(class = col_cls, col = list(class = CLASS_COL), show_annotation_name = FALSE, annotation_legend_param = list(class = list(title = "class", title_gp = gpar(fontsize = 5.6), labels_gp = gpar(fontsize = 5.2))), simple_anno_size = unit(2.4, "mm"))
ha_row_left <- rowAnnotation(class = row_dom, col = list(class = CLASS_COL), show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(2.2, "mm"), border = TRUE)
ha_row <- rowAnnotation(NES = anno_barplot(row_nes, gp = gpar(fill = CLASS_COL[as.character(row_dom)], col = NA), border = FALSE, baseline = 0, bar_width = 0.84, width = unit(11.5, "mm"), axis_param = list(at = c(0, 1, 2), labels = c("0", "1", "2"), gp = gpar(fontsize = 4.2))), top3load = row_t3ld, col = list(top3load = circlize::colorRamp2(quantile(row_t3ld, c(0, 0.5, 1), na.rm = TRUE), c("#FCFDBF", "#B5367A", "#1D1147"))), genes = anno_text(row_gene, gp = gpar(fontsize = 4.6, fontface = "italic"), 
    just = "left", location = unit(0, "npc")), annotation_label = c("NES", "load", "genes"), annotation_name_gp = gpar(fontsize = 4.4, fontface = "bold"), annotation_name_rot = 0, annotation_name_side = "top", simple_anno_size = unit(1.8, "mm"), gap = unit(1.3, "mm"), annotation_legend_param = list(top3load = list(title = "mean top-3 loading", title_gp = gpar(fontsize = 5.6), labels_gp = gpar(fontsize = 5.2))), show_legend = c(NES = FALSE, top3load = TRUE, genes = FALSE))
B_CELL <- 2
R_CELL <- 2.6
ht_b <- Heatmap(z, name = "row z", col = colorRamp2(c(-3, 0, 3), c("#3B4CC0", "#F7F7F7", "#B40426")), width = unit(ncol(z) * B_CELL, "mm"), height = unit(nrow(z) * R_CELL, "mm"), cluster_rows = TRUE, clustering_method_rows = "ward.D2", row_split = row_dom, row_gap = unit(1, "mm"), cluster_row_slices = FALSE, show_row_dend = TRUE, row_dend_width = unit(7, "mm"), cluster_columns = FALSE, column_split = factor(col_cls, levels = CLASS_LV), column_gap = unit(1.4, "mm"), row_names_gp = gpar(fontsize = 4.8), 
    row_names_side = "left", row_title_gp = gpar(fontsize = 5.2, fontface = "bold"), row_title_rot = 0, show_column_names = TRUE, column_names_gp = gpar(fontsize = 5.2), column_names_rot = 90, column_names_side = "bottom", column_title_gp = gpar(fontsize = 6, fontface = "bold"), top_annotation = ha_col, left_annotation = ha_row_left, right_annotation = ha_row, heatmap_legend_param = list(title = "row z", title_gp = gpar(fontsize = 5.6), labels_gp = gpar(fontsize = 5.2), legend_height = unit(12, "mm")))
assembled <- file.path(SVGD, "fig1_b_assembled.svg")
svglite(assembled, width = mm2in(190), height = mm2in(185), bg = "white")
draw(ht_b, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE, column_title = "Program backbone", column_title_gp = gpar(fontsize = 7.2, fontface = "bold"), padding = unit(c(0, 0, 0, 0), "mm"))
dev.off()
fb <- file.path(SVGD, "fig1_b.svg")
file.copy(assembled, fb, overwrite = TRUE)
file.remove(assembled)
