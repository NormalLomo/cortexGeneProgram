#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 02_run_rctd.R REFERENCE_RDS SPATIAL_RDS OUTPUT_RDS", call. = FALSE)
}
reference <- readRDS(args[[1]])
spatial <- readRDS(args[[2]])
if (is.list(reference) && all(c("counts", "cell_type") %in% names(reference))) {
  reference_counts <- reference$counts
  cell_type <- reference$cell_type
} else {
  reference_counts <- GetAssayData(reference, assay = DefaultAssay(reference), layer = "counts")
  if (!"cell_type" %in% colnames(reference[[]])) stop("Reference metadata requires cell_type.", call. = FALSE)
  cell_type <- reference$cell_type
}
if (is.list(spatial) && all(c("counts", "coordinates") %in% names(spatial))) {
  spatial_counts <- spatial$counts
  coordinates <- spatial$coordinates
} else {
  spatial_counts <- GetAssayData(spatial, assay = DefaultAssay(spatial), layer = "counts")
  coordinates <- GetTissueCoordinates(spatial)
}
rctd <- create.RCTD(
  puck = SpatialRNA(coords = coordinates, counts = spatial_counts),
  reference = Reference(counts = reference_counts, cell_types = cell_type),
  max_cores = 1
)
rctd <- run.RCTD(rctd, doublet_mode = "full")
dir.create(dirname(args[[3]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(rctd, args[[3]])
message("Wrote RCTD result: ", args[[3]])
