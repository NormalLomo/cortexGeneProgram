# Panel A: subclass x program eta2 heatmap (ComplexHeatmap)
source("CORTEX_PROGRAM_ROOT/figures/fig4/fig4_theme.R")
suppressMessages({library(ComplexHeatmap); library(circlize); library(grid)})

eta <- read.delim(file.path(RES, "within_subclass_region_eta2.tsv"))
eta$program <- as.character(eta$program)
M <- eta %>% select(subclass, program, eta2) %>%
  pivot_wider(names_from = program, values_from = eta2) %>%
  as.data.frame()
rownames(M) <- M$subclass; M$subclass <- NULL
# order columns numerically
M <- M[, order(as.integer(colnames(M)))]
M <- as.matrix(M)

cls <- sc_class(rownames(M))
# order rows by driver rank (median eta2)
drv <- read.delim(file.path(RES, "subclass_driver_rank.tsv"))
M <- M[drv$subclass, ]
cls <- sc_class(rownames(M))

col_fun <- colorRamp2(
  quantile(M, c(0, 0.5, 0.8, 0.95, 1)),
  c("#FCFDBF","#FEC287","#F1605D","#8C2981","#2D1160"))

ra <- rowAnnotation(
  Class = cls,
  col = list(Class = pal_class),
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(Class = list(title = "Class",
      labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 8))),
  simple_anno_size = unit(3, "mm"))

# mark top cells: top 12 (subclass,program) by eta2
top <- eta %>% arrange(desc(eta2)) %>% head(12)
mark_mat <- matrix("", nrow = nrow(M), ncol = ncol(M),
                   dimnames = dimnames(M))
for (i in seq_len(nrow(top))) {
  s <- top$subclass[i]; p <- top$program[i]
  if (s %in% rownames(mark_mat) && p %in% colnames(mark_mat))
    mark_mat[s, p] <- "*"
}

# Dense 60-wide axis: full functional name only for the highlighted (starred)
# programs; all others stay as bare "P{n}" to avoid overflow/overlap.
hi_progs <- sort(unique(as.integer(top$program)))
col_labs <- prog_label_selective(colnames(M), hi_progs)

ht <- Heatmap(M, name = "eta2", col = col_fun,
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_row_dend = TRUE, show_column_dend = TRUE,
  row_dend_width = unit(8, "mm"), column_dend_height = unit(8, "mm"),
  row_names_gp = gpar(fontsize = 6.5),
  column_labels = col_labs,
  column_names_gp = gpar(fontsize = 5),
  column_title = "Programs (K=60)", column_title_gp = gpar(fontsize = 8),
  row_title = "Subclass", row_title_gp = gpar(fontsize = 8),
  left_annotation = ra,
  heatmap_legend_param = list(title = expression(eta^2),
    labels_gp = gpar(fontsize = 7.5), title_gp = gpar(fontsize = 9),
    legend_height = unit(28, "mm"), legend_width = unit(6, "mm"),
    grid_width = unit(5, "mm")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (mark_mat[i, j] == "*")
      grid.text("*", x, y, gp = gpar(fontsize = 8, col = "white", fontface = "bold"))
  },
  rect_gp = gpar(col = NA))

pdf(file.path(FIG, "fig4_a.pdf"), width = 8.2, height = 4.0)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
cat("panel a done\n")
