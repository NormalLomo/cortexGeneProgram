#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[[2]] %in% c("reference", "spatial")) {
  stop("Usage: 02_prepare_rctd_objects.R INPUT_H5SEURAT reference|spatial OUTPUT_RDS", call. = FALSE)
}

input_h5seurat <- args[[1]]
object_kind <- args[[2]]
output_rds <- args[[3]]
object <- LoadH5Seurat(input_h5seurat)
counts <- GetAssayData(object, assay = DefaultAssay(object), layer = "counts")
metadata <- object[[]]

if (object_kind == "reference") {
  if (!"cell_type" %in% colnames(metadata)) stop("Reference metadata requires cell_type.", call. = FALSE)
  output <- list(counts = counts, cell_type = metadata$cell_type)
} else {
  if (!all(c("x", "y") %in% colnames(metadata))) {
    stop("Spatial H5Seurat metadata requires x and y coordinates from 00_bin_spatial_counts.py.", call. = FALSE)
  }
  output <- list(counts = counts, coordinates = data.frame(x = metadata$x, y = metadata$y, row.names = rownames(metadata)))
}
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(output, output_rds)
message("Wrote RCTD ", object_kind, " object: ", output_rds)
