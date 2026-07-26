PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(Sys.getenv("CORTEX_PROGRAM_ROOT", "/DATA/cortex_nmf_program"), mustWork = FALSE)
suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratDisk)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[[2]] %in% c("reference", "spatial")) {
    stop(call. = FALSE)
}
input_h5seurat <- args[[1]]
object_kind <- args[[2]]
output_rds <- args[[3]]
object <- LoadH5Seurat(input_h5seurat)
counts <- GetAssayData(object, assay = DefaultAssay(object), layer = "counts")
metadata <- object[[]]
if (object_kind == "reference") {
    if (!"cell_type" %in% colnames(metadata)) 
        stop(call. = FALSE)
    output <- list(counts = counts, cell_type = metadata$cell_type)
} else {
    if (!all(c("x", "y") %in% colnames(metadata))) {
        stop(call. = FALSE)
    }
    output <- list(counts = counts, coordinates = data.frame(x = metadata$x, y = metadata$y, row.names = rownames(metadata)))
}
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(output, output_rds)
