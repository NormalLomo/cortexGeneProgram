PROJECT_ROOT <- Sys.getenv("CORTEX_PROGRAM_ROOT", unset = "/DATA/cortex_nmf_program")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(Sys.getenv("CORTEX_PROGRAM_ROOT", "/DATA/cortex_nmf_program"), mustWork = FALSE)
suppressPackageStartupMessages(library(SeuratDisk))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop(call. = FALSE)
}
input_h5ad <- args[[1]]
output_h5seurat <- args[[2]]
if (!file.exists(input_h5ad)) stop(call. = FALSE)
dir.create(dirname(output_h5seurat), recursive = TRUE, showWarnings = FALSE)
Convert(input_h5ad, dest = "h5seurat", filename = output_h5seurat, overwrite = TRUE)
if (!file.exists(output_h5seurat)) stop(call. = FALSE)
