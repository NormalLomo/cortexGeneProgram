#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(SeuratDisk))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 01_convert_h5ad_to_h5seurat.R INPUT_H5AD OUTPUT_H5SEURAT", call. = FALSE)
}

input_h5ad <- args[[1]]
output_h5seurat <- args[[2]]
if (!file.exists(input_h5ad)) stop("Input H5AD does not exist: ", input_h5ad, call. = FALSE)
dir.create(dirname(output_h5seurat), recursive = TRUE, showWarnings = FALSE)
Convert(input_h5ad, dest = "h5seurat", filename = output_h5seurat, overwrite = TRUE)
if (!file.exists(output_h5seurat)) stop("SeuratDisk did not create: ", output_h5seurat, call. = FALSE)
message("Wrote H5Seurat: ", output_h5seurat)
