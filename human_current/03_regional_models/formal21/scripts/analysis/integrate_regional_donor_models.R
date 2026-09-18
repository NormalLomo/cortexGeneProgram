options(future.globals.maxSize = 50 * 1024^3)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript integrate_regional_donor_models.R <project_root>")
root <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- file.path(root, "results", "crossregion_v2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

programs <- paste0("P", seq_len(54L))
overall <- list(
  nature_lmm = read_tsv(file.path(root, "celltype_regional/overall_nature_lmm/all_54_overall_nature_lmm.tsv")),
  limma_dupcor = read_tsv(file.path(root, "limma_dupcor_validation/all_54_limma_dupcor.tsv")),
  dream = read_tsv(file.path(root, "dream_validation/all_54_dream.tsv"))
)
subclass <- list(
  nature_lmm = read_tsv(file.path(root, "celltype_regional/nature_lmm/all_54_by_22_nature_lmm.tsv")),
  limma_dupcor = read_tsv(file.path(root, "regional_donor_validation_v2/methods/limma_dupcor/all_54_by_22_limma_dupcor.tsv")),
  dream = read_tsv(file.path(root, "regional_donor_validation_v2/methods/dream/all_54_by_22_dream.tsv"))
)

# Export the unchanged six model tables under the names consumed by Table S4.
for (method in names(overall)) {
  write_tsv(overall[[method]], file.path(output_dir, paste0("overall_", method, ".tsv")))
  write_tsv(subclass[[method]], file.path(output_dir, paste0("subclass_", method, ".tsv")))
}

region_p <- function(x) {
  if ("p_region" %in% names(x)) as.numeric(x$p_region) else as.numeric(x[["P.Value"]])
}

overall_flag <- function(x) {
  x$BH_54_recomputed <- p.adjust(region_p(x), method = "BH", n = 54L)
  flag <- x$BH_54_recomputed < 0.05
  names(flag) <- x$program
  as.integer(flag[programs])
}

subclass_count <- function(x) {
  x$BH_global_22x54_recomputed <- p.adjust(region_p(x), method = "BH", n = 54L * 22L)
  flag <- x$BH_global_22x54_recomputed < 0.05
  count <- tapply(flag, factor(x$program, levels = programs), sum)
  as.integer(count[programs])
}

matrix <- data.frame(
  program = programs,
  overall_nature_lmm = overall_flag(overall$nature_lmm),
  overall_limma_dupcor = overall_flag(overall$limma_dupcor),
  overall_dream = overall_flag(overall$dream),
  n_subclasses_nature_lmm = subclass_count(subclass$nature_lmm),
  n_subclasses_limma_dupcor = subclass_count(subclass$limma_dupcor),
  n_subclasses_dream = subclass_count(subclass$dream),
  check.names = FALSE
)
consensus <- matrix[rowSums(matrix[, c("overall_nature_lmm", "overall_limma_dupcor", "overall_dream"), drop = FALSE]) == 3L, , drop = FALSE]
write_tsv(matrix, file.path(output_dir, "program_method_matrix.tsv"))
write_tsv(consensus, file.path(output_dir, "regional_consensus_21.tsv"))
