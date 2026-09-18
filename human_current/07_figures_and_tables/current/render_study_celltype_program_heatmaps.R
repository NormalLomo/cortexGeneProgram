if (!nzchar(Sys.getenv("NMF_WORK_ROOT"))) stop("Set NMF_WORK_ROOT; original source directories must remain read-only.")
#!/usr/bin/env Rscript
options(future.globals.maxSize = 50 * 1024^3)

# Native ComplexHeatmap entry adapted from the actual Fig.1b core:
# Historical source: scripts/fig1/fig1_panel_b_merged.R
# Reused here: class strips, Ward.D2 row slices, Heatmap geometry, palette, and draw().
# The matrix below is the retained column-zscore field; no row-zscore, NA->0, or clamping
# is performed on the source values.
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

SOURCE <- paste0(Sys.getenv("NMF_WORK_ROOT"), "/analysis/external_single_cell_program_scores/study_celltype_program_score_heatmap.tsv")
OUTDIR <- paste0(Sys.getenv("NMF_WORK_ROOT"), "/analysis/external_single_cell_program_scores/program_score_heatmaps_by_dataset")
SUBCLASS_CLASS <- paste0(Sys.getenv("NMF_WORK_ROOT"), "/inputs/cortex_nmf_program/archived/figures/fig1/_intermediate/subclass_class.csv")

STUDIES <- c(
  "Allen",
  "SEAAD",
  "Grubman_2019_GSE138852",
  "Lake_2018_GSE97930",
  "Schirmer_2019_MS_UCSC",
  "SingleSoma_AD_PFC_CELLxGENE",
  "Tran_2021_reward_cortex_CELLxGENE",
  "GSE144136",
  "GSE174367",
  "GSE291605"
)
STUDY_LABELS <- c(
  Reference_1M = "1M reference",
  Allen = "Allen",
  SEAAD = "SEA-AD",
  Grubman_2019_GSE138852 = "Grubman",
  Lake_2018_GSE97930 = "Lake",
  Schirmer_2019_MS_UCSC = "Schirmer",
  SingleSoma_AD_PFC_CELLxGENE = "SingleSoma",
  Tran_2021_reward_cortex_CELLxGENE = "Tran",
  GSE144136 = "Nagy",
  GSE174367 = "Morabito",
  GSE291605 = "Mural"
)
CLASS_LEVELS <- c("excitatory", "inhibitory", "non-neuronal", "Mixed", "Unknown")
CLASS_COL <- c(
  excitatory = "#2166AC",
  inhibitory = "#B2182B",
  `non-neuronal` = "#1B7837",
  Mixed = "#7F7F7F",
  Unknown = "#BDBDBD"
)
PROGRAMS <- paste0("P", seq_len(54))

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
src <- read.delim(
  SOURCE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = "NA"
)
src$program <- as.character(src$program)
src$column_key <- as.character(src$column_key)
src$column_order <- as.integer(src$column_order)
src$program_order <- as.integer(src$program_order)

column_meta <- src[!duplicated(src$column_key), c(
  "column_key", "column_order", "column_type", "column_label",
  "study", "study_display", "author_label", "author_label_display",
  "reference_subclass", "column_major_class", "column_annotation_status",
  "column_annotation_raw", "column_annotation_source"
), drop = FALSE]
column_meta <- column_meta[order(column_meta$column_order), , drop = FALSE]
column_meta$column_major_class <- as.character(column_meta$column_major_class)
column_meta$column_major_class[
  !column_meta$column_major_class %in% CLASS_LEVELS
] <- "Unknown"
column_meta$column_label <- as.character(column_meta$column_label)
column_meta$column_label[is.na(column_meta$column_label) | column_meta$column_label == ""] <-
  "Missing original annotation"

program_meta <- src[match(PROGRAMS, src$program), c(
  "program", "program_dominant_class", "program_dominant_subclass",
  "program_confidence", "program_annotation_source"
), drop = FALSE]
program_class <- as.character(program_meta$program_dominant_class)
program_class[is.na(program_class) | !program_class %in% CLASS_LEVELS] <- "Unknown"
program_class <- factor(program_class, levels = CLASS_LEVELS)

z_matrix <- matrix(
  NA_real_,
  nrow = length(PROGRAMS),
  ncol = nrow(column_meta),
  dimnames = list(PROGRAMS, column_meta$column_key)
)
for (j in seq_len(nrow(column_meta))) {
  rows <- src[src$column_key == column_meta$column_key[j], , drop = FALSE]
  rows <- rows[match(PROGRAMS, rows$program), , drop = FALSE]
  z_matrix[, j] <- as.numeric(rows$column_zscore)
}

subclass_class <- read.csv(
  SUBCLASS_CLASS,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
reference_class <- setNames(subclass_class$class, subclass_class$subclass)

render_one <- function(target_name, indices) {
  z <- z_matrix[, indices, drop = FALSE]
  cols <- column_meta[indices, , drop = FALSE]
  col_class <- factor(cols$column_major_class, levels = CLASS_LEVELS)
  row_cluster <- any(is.finite(z))
  row_split <- if (row_cluster) program_class else NULL

  ha_col <- HeatmapAnnotation(
    class = col_class,
    col = list(class = CLASS_COL),
    show_annotation_name = FALSE,
    annotation_legend_param = list(class = list(
      title = "class",
      title_gp = gpar(fontsize = 5.6),
      labels_gp = gpar(fontsize = 5.2)
    )),
    simple_anno_size = unit(2.4, "mm")
  )
  ha_row_left <- rowAnnotation(
    class = program_class,
    col = list(class = CLASS_COL),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(2.2, "mm"),
    border = TRUE
  )

  ht <- Heatmap(
    z,
    name = "column z",
    col = colorRamp2(c(-3, 0, 3), c("#3B4CC0", "#F7F7F7", "#B40426")),
    na_col = "#E6E6E6",
    width = unit(ncol(z) * 2.0, "mm"),
    height = unit(nrow(z) * 2.6, "mm"),
    cluster_rows = row_cluster,
    clustering_method_rows = "ward.D2",
    row_split = row_split,
    cluster_row_slices = FALSE,
    show_row_dend = row_cluster,
    row_dend_width = unit(7, "mm"),
    cluster_columns = FALSE,
    column_split = col_class,
    column_gap = unit(1.4, "mm"),
    column_title_gp = gpar(fontsize = 5.2, fontface = "bold"),
    row_names_gp = gpar(fontsize = 4.8),
    row_names_side = "left",
    row_title_gp = gpar(fontsize = 5.2, fontface = "bold"),
    row_title_rot = 0,
    show_column_names = TRUE,
    column_labels = cols$column_label,
    column_names_gp = gpar(fontsize = 5.2),
    column_names_rot = 90,
    column_names_side = "bottom",
    top_annotation = ha_col,
    left_annotation = ha_row_left,
    heatmap_legend_param = list(
      title = "column z (±3 saturated)",
      title_gp = gpar(fontsize = 5.6),
      labels_gp = gpar(fontsize = 5.2),
      at = c(-3, 0, 3),
      labels = c("<= -3", "0", ">= 3")
    )
  )

  width_mm <- max(120, ncol(z) * 2.0 + 70)
  height_mm <- 54 * 2.6 + 60
  pdf_path <- file.path(OUTDIR, paste0(target_name, ".pdf"))
  png_path <- file.path(OUTDIR, paste0(target_name, ".png"))
  pdf(pdf_path, width = width_mm / 25.4, height = height_mm / 25.4, bg = "white", useDingbats = FALSE)
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE,
    column_title = STUDY_LABELS[[target_name]],
    column_title_gp = gpar(fontsize = 7.2, fontface = "bold"),
    padding = unit(c(0, 0, 0, 0), "mm")
  )
  dev.off()
  png(
    png_path,
    width = round(width_mm / 25.4 * 300),
    height = round(height_mm / 25.4 * 300),
    res = 300,
    type = "cairo-png",
    bg = "white"
  )
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE,
    column_title = STUDY_LABELS[[target_name]],
    column_title_gp = gpar(fontsize = 7.2, fontface = "bold"),
    padding = unit(c(0, 0, 0, 0), "mm")
  )
  dev.off()
  finite_cells <- sum(is.finite(z))
  class_count <- table(factor(as.character(col_class), levels = CLASS_LEVELS))
  cat(sprintf(
    "study_complete\t%s\tcolumns\t%d\tfinite_cells\t%d\tna_cells\t%d\n",
    target_name, ncol(z), finite_cells, length(z) - finite_cells
  ))
  for (class_name in names(class_count)) {
    if (class_count[[class_name]] > 0) {
      cat(sprintf("column_class\t%s\t%s\t%d\n", target_name, class_name, class_count[[class_name]]))
    }
  }
  cat(sprintf("output\t%s\noutput\t%s\n", pdf_path, png_path))
}

reference_indices <- which(column_meta$column_type == "reference_subclass")
render_one("Reference_1M", reference_indices)
for (study in STUDIES) {
  indices <- which(column_meta$study == study & column_meta$column_type != "reference_subclass")
  render_one(study, indices)
}
cat(sprintf("program_annotation_source\t%s\n", unique(program_meta$program_annotation_source)[1]))
cat(sprintf("program_rows\t%d\nexternal_columns\t%d\nreference_columns\t%d\n", length(PROGRAMS), length(STUDIES) * 0 + sum(column_meta$column_type != "reference_subclass"), length(reference_indices)))
