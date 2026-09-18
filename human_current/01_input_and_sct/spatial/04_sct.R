#!/usr/bin/env Rscript
# 04_sct.R — per-chip SCTransform on bin50 raw counts, matching the bin200 recipe
# (sctransform::vst -> sctransform::correct at a GLOBAL scale_factor, do_round + do_pos).
# Pure sctransform (no Seurat needed). Outputs SCT-corrected integer counts (genes x cells).
#
# Args: <chip> <counts_mtx(genes x cells)> <genes.txt> <bcs.txt> <scale_factor> <out_corrected_mtx>
suppressMessages({library(Matrix); library(sctransform)})

# sctransform 0.4.3 dispatches its per-gene NB regression through `future`, which caps
# exported globals at 500 MiB by default. At bin50 scale the closure is ~1.7 GiB -> raise it.
# Pin sequential so future spawns no workers (cross-chip parallelism is handled by xargs).
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(if (requireNamespace("future", quietly = TRUE)) try(future::plan("sequential"), silent = TRUE))

args <- commandArgs(trailingOnly = TRUE)
chip <- args[1]; mtx <- args[2]; gp <- args[3]; bp <- args[4]
scale_factor <- as.numeric(args[5]); out <- args[6]

umi <- as(Matrix::readMM(mtx), "CsparseMatrix")        # genes x cells, integer counts
rownames(umi) <- readLines(gp)
colnames(umi) <- readLines(bp)
cat(sprintf("[04_sct] %s vst: %d genes x %d cells (scale_factor=%g)\n",
            chip, nrow(umi), ncol(umi), scale_factor))

vst_out <- sctransform::vst(umi = umi,
                            return_corrected_umi = FALSE,
                            return_gene_attr = TRUE,
                            return_cell_attr = TRUE,
                            residual_type = "pearson",
                            verbosity = 1)

# correct() reverses the regularized NB GLM to the common scale_factor depth.
corrected <- sctransform::correct(vst_out,
                                  scale_factor = scale_factor,
                                  do_round = TRUE, do_pos = TRUE,
                                  verbosity = 1)
corrected <- as(corrected, "CsparseMatrix")

Matrix::writeMM(corrected, out)
writeLines(rownames(corrected), paste0(out, ".genes"))
writeLines(colnames(corrected), paste0(out, ".bcs"))
cat(sprintf("[04_sct] %s corrected %d genes x %d cells -> %s\n",
            chip, nrow(corrected), ncol(corrected), out))

# --- Pearson residuals (dense matrix from vst) -> HDF5 for the Python side ---
# bin200 SCT discarded these; we now store them in the completo h5ad layers['sct_residuals'].
res_out <- sub("\\.mtx$", "_residuals.h5", out)
y <- as.matrix(vst_out$y)   # modeled genes x cells, real-valued
ok <- FALSE
if (requireNamespace("hdf5r", quietly = TRUE)) {
  h <- hdf5r::H5File$new(res_out, mode = "w")
  # write float32 (halves disk + downstream RAM; full precision not needed for residuals)
  h$create_dataset("residuals", robj = y, dtype = hdf5r::h5types$H5T_NATIVE_FLOAT)
  h[["genes"]] <- rownames(y)
  h[["bcs"]]   <- colnames(y)
  h$close_all(); ok <- TRUE
} else if (requireNamespace("rhdf5", quietly = TRUE)) {
  if (file.exists(res_out)) stop("Output already exists: ", res_out)
  rhdf5::h5createFile(res_out)
  rhdf5::h5createDataset(res_out, "residuals", dim(y),
                         storage.mode = "double", H5type = "H5T_IEEE_F32LE")
  rhdf5::h5write(y, res_out, "residuals")
  rhdf5::h5write(rownames(y), res_out, "genes")
  rhdf5::h5write(colnames(y), res_out, "bcs")
  rhdf5::H5close(); ok <- TRUE
}
if (!ok) stop("[04_sct] need either 'hdf5r' or 'rhdf5' R package to write residuals")
cat(sprintf("[04_sct] %s DONE: residuals %d x %d -> %s\n",
            chip, nrow(y), ncol(y), res_out))
