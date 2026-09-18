# =============================================================================
# STEP 01 / 03  —  snRNA preprocessing (SCTransform)
# Provenance chain for the data in this directory (neuropeptide_cortex/data/human/snrna/):
#   [upstream, NOT in this dir] cortex/SnRNA/3_SnRNA_seurat_merged_1m_Cells.RDS
#        merged human cortex snRNA: edlein 814,034 + us 222,005 = 1,036,039 cells
#   --01_prepareData_snRNA.R (THIS)--> snRNA.seu.sct.qs
#   --02_prep_full_1M.R-->             chunks_1M/ , snRNA_1M_obs.csv , snRNA_1M_var.csv
#   --03_build_h5ad.py-->         snRNA_1M.h5ad
# THIS STEP: per-library_prep SCTransform vst + correct (+ 2nd ScaleData on residuals),
#   then merge -> snRNA.seu.sct.qs.
#   NOTE: the RNA 'counts' layer written here is SCT-CORRECTED, depth-equalized, rounded
#   counts (NOT raw UMI) — this is intentional.
# IN : cortex/SnRNA/3_SnRNA_seurat_merged_1m_Cells.RDS  (external upstream input)
#      SeuratObject and Seurat source trees supplied by environment variables
# OUT: snRNA.seu.sct.qs , vst.list_snRNA.qs  (this dir)
# NEXT: 02_prep_full_1M.R
# =============================================================================
work_dir <- Sys.getenv("SNRNA_WORK_DIR")
if (!nzchar(work_dir)) stop("Set SNRNA_WORK_DIR to a writable analysis directory.")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)
devtools::load_all(Sys.getenv("SEURAT_OBJECT_SOURCE"))
devtools::load_all(Sys.getenv("SEURAT_SOURCE"))
library(magrittr)
library(tidyverse)
library(flock)
options(future.globals.maxSize = 50 * 1024^3)

# === read upstream merged snRNA (raw counts), join split layers ===
snRNA <- readRDS(Sys.getenv("SNRNA_MERGED_RDS"))
snRNA %<>% JoinLayers()

# === per-library_prep SCTransform (vst + correct + 2nd ScaleData) ===
snRNA.seu <- CreateSeuratObject(counts = snRNA[["RNA"]]$counts, meta.data = snRNA@meta.data)
snRNA.seu.list <- snRNA.seu %>% SplitObject(split.by = "library_prep")

sample_median_umis <- lapply(snRNA.seu.list, FUN = function(x) median(x$nCount_RNA)) %>% unlist %>% median(na.rm = TRUE) %>% round()

vst.list <- parallel::mclapply(snRNA.seu.list, mc.cores = 10, mc.preschedule = F,
                               FUN = function(seu) {

                                 lock <- flock::lock("tmp.lock", exclusive = TRUE)
                                 seu[[DefaultAssay(seu)]]$counts %<>% as("dgCMatrix")
                                 flock::unlock(lock)

                                 vst_out <- sctransform::vst(
                                   umi = seu[[DefaultAssay(seu)]]$counts,
                                   return_corrected_umi = FALSE,
                                   return_gene_attr = TRUE,
                                   return_cell_attr = TRUE,
                                   residual_type = "pearson", verbosity = 2)

                                 umi_corrected <- sctransform::correct(
                                   vst_out,
                                   scale_factor = sample_median_umis, # 取所有样本里最小/中位的，不能往大了取
                                   do_round = TRUE, do_pos = T, verbosity = 2
                                 )

                                 scale.data <- vst_out$y
                                 # clip the residuals
                                 clip.range = c(-sqrt(x = ncol(x = seu) / 30), sqrt(x = ncol(x = seu) / 30))
                                 scale.data[scale.data < clip.range[1]] <- clip.range[1]
                                 scale.data[scale.data > clip.range[2]] <- clip.range[2]

                                 # 2nd regression (center only)
                                 scale.data %<>% ScaleData(
                                   features = NULL,
                                   vars.to.regress = NULL,
                                   latent.data = NULL,
                                   model.use = 'linear',
                                   use.umi = FALSE,
                                   do.scale = F,
                                   do.center = T,
                                   scale.max = Inf )

                                 vst_out$y <- scale.data
                                 vst_out$y %<>% as("dgCMatrix") %>% BPCells::write_matrix_memory()
                                 vst_out$umi_corrected <- umi_corrected %>% as("dgCMatrix") %>% BPCells::write_matrix_memory()
                                 vst_out$library_prep <- seu$library_prep[[1]]

                                 vst_out
                               })

# === assemble merged SCT object ===
corrected_counts.list <- lapply(vst.list, function(x){ x$umi_corrected })
corrected_data.list   <- lapply(corrected_counts.list, log1p)
scale.data.list <- lapply(vst.list, function(x){ x$y })
scaled.features <- lapply(scale.data.list, rownames) %>% purrr::reduce(intersect)
scale.data.list <- lapply(scale.data.list, function(x) x[scaled.features, , drop = FALSE])
scale.data <- do.call(cbind, scale.data.list)

snRNA.seu.sct <- CreateSeuratObject(
  counts = corrected_counts.list, data = corrected_data.list,
  assay = "RNA",
  meta.data = snRNA@meta.data
)

snRNA.seu.sct[["RNA"]] %<>% SetAssayData(
  slot = 'scale.data',
  new.data = scale.data
)

snRNA.seu.sct %>% qs::qsave("snRNA.seu.sct.qs")
