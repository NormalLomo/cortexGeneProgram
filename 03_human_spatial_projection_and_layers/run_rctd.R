PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(Sys.getenv("CORTEX_PROGRAM_ROOT", "/DATA/cortex_nmf_program"), mustWork = FALSE)
suppressPackageStartupMessages({
    library(Seurat)
    library(spacexr)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop(call. = FALSE)
}
reference <- readRDS(args[[1]])
spatial <- readRDS(args[[2]])
if (is.list(reference) && all(c("counts", "cell_type") %in% names(reference))) {
    reference_counts <- reference$counts
    cell_type <- reference$cell_type
} else {
    reference_counts <- GetAssayData(reference, assay = DefaultAssay(reference), layer = "counts")
    if (!"cell_type" %in% colnames(reference[[]])) 
        stop(call. = FALSE)
    cell_type <- reference$cell_type
}
if (is.list(spatial) && all(c("counts", "coordinates") %in% names(spatial))) {
    spatial_counts <- spatial$counts
    coordinates <- spatial$coordinates
} else {
    spatial_counts <- GetAssayData(spatial, assay = DefaultAssay(spatial), layer = "counts")
    coordinates <- GetTissueCoordinates(spatial)
}
rctd <- create.RCTD(puck = SpatialRNA(coords = coordinates, counts = spatial_counts), reference = Reference(counts = reference_counts, cell_types = cell_type), max_cores = 1)
rctd <- run.RCTD(rctd, doublet_mode = "full")
dir.create(dirname(args[[3]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(rctd, args[[3]])
