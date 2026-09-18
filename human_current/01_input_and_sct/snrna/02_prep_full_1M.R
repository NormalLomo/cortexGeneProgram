options(future.globals.maxSize = 50 * 1024^3)
# =============================================================================
# STEP 02 / 03  —  extract per-donor count chunks from the SCT Seurat object
# IN : snRNA.seu.sct.qs        (this dir, produced by 01_prepareData_snRNA.R)
# DOES: JoinLayers; GetAssayData(assay="RNA", layer="counts") per donor -> .mtx chunks
# OUT: chunks_1M/ , snRNA_1M_obs.csv , snRNA_1M_var.csv   (this dir)
# PREV: 01_prepareData_snRNA.R     NEXT: 03_build_h5ad.py
# (full chain documented in 01_prepareData_snRNA.R header)
# =============================================================================
# P10 v2: write per-donor chunks then concatenate
suppressPackageStartupMessages({
  library(qs); library(Seurat); library(SeuratObject); library(Matrix)
})
out_dir <- Sys.getenv("SNRNA_WORK_DIR")
if (!nzchar(out_dir)) stop("Set SNRNA_WORK_DIR.")

message(sprintf("[%s] Loading 160 GB .qs", Sys.time()))
seu <- qread(file.path(out_dir, "snRNA.seu.sct.qs"), nthreads=8)
message(sprintf("[%s] Loaded: %d cells x %d genes", Sys.time(), ncol(seu), nrow(seu)))

message("JoinLayers...")
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

# Save full metadata once
write.csv(seu@meta.data, file.path(out_dir, "snRNA_1M_obs.csv"), row.names=TRUE)
write.csv(data.frame(gene=rownames(seu)), file.path(out_dir, "snRNA_1M_var.csv"), row.names=FALSE)
message("Metadata + var written")

donors <- unique(seu$donor)
message(sprintf("Donors: %d", length(donors)))

chunk_dir <- file.path(out_dir, "chunks_1M")
dir.create(chunk_dir, showWarnings=FALSE)

for (d in donors) {
  message(sprintf("[%s] Donor %s", Sys.time(), d))
  sub <- seu[, seu$donor == d]
  message(sprintf("  cells: %d", ncol(sub)))
  cnt <- GetAssayData(sub, assay="RNA", layer="counts")
  cnt <- as(cnt, "dgCMatrix")
  writeMM(cnt, file.path(chunk_dir, sprintf("counts_%s.mtx", d)))
  # also save cell barcodes for this chunk to verify ordering
  writeLines(colnames(cnt), file.path(chunk_dir, sprintf("barcodes_%s.txt", d)))
  message(sprintf("  Done: nnz=%d", length(cnt@x)))
  rm(sub, cnt); gc()
}
message("DONE.")
