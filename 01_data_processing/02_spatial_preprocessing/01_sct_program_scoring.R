#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 01_sct_program_scoring.R INPUT_H5SEURAT PROGRAM_SPECTRA_TSV OUTPUT_TSV", call. = FALSE)
}
input_h5seurat <- args[[1]]
spectra_tsv <- args[[2]]
output_tsv <- args[[3]]

object <- LoadH5Seurat(input_h5seurat)
object <- SCTransform(object, assay = DefaultAssay(object), verbose = FALSE)
expression <- GetAssayData(object, assay = "SCT", layer = "data")
spectra <- as.matrix(read.delim(spectra_tsv, row.names = 1, check.names = FALSE))
genes <- intersect(colnames(spectra), rownames(expression))
if (length(genes) == 0L) stop("No overlapping program-spectrum and SCT genes.", call. = FALSE)
scores <- t(spectra[, genes, drop = FALSE] %*% expression[genes, , drop = FALSE])
scores <- as.data.frame(scores)
scores$bin_id <- rownames(scores)
dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)
write.table(scores, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote SCT program scores: ", output_tsv)
