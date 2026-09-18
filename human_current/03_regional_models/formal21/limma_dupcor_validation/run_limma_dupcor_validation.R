if (!nzchar(Sys.getenv("CORTEX_PROGRAM_ROOT"))) stop("Set CORTEX_PROGRAM_ROOT; original source directories must remain read-only.")
options(future.globals.maxSize = 50 * 1024^3)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript run_limma_dupcor_validation.R <audit_root> <output_dir>")

root <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(root, "Rlib"), .libPaths()))

suppressPackageStartupMessages(library(limma))

programs <- paste0("P", seq_len(54L))
old8 <- c("P1", "P3", "P4", "P6", "P8", "P9", "P13", "P33")
input_path <- file.path(root, "output", "donor_region_means.tsv")
metadata_path <- paste0(Sys.getenv("CORTEX_PROGRAM_ROOT"), "/inputs/snRNA_1M_obs.csv")
name_path <- file.path(root, "scripts", "program_names.tsv")
formula_text <- "~ region + age + sex"
command_text <- paste("Rscript run_limma_dupcor_validation.R", root, output_dir)

write_notes <- function(lines) writeLines(lines, file.path(output_dir, "METHOD_NOTES.md"))

means <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
required_metadata <- c("donor", "age", "sex")
if (!all(required_metadata %in% names(metadata))) stop("Source metadata lacks donor, age, or sex")

donor_covariates <- unique(metadata[, required_metadata])
if (anyDuplicated(donor_covariates$donor)) {
  donor_check <- aggregate(cbind(age, sex) ~ donor, metadata[, required_metadata], function(x) length(unique(x)))
  if (any(donor_check$age > 1L | donor_check$sex > 1L)) stop("Age or sex is inconsistent within a donor")
  donor_covariates <- donor_covariates[!duplicated(donor_covariates$donor), , drop = FALSE]
}

dat <- merge(means, donor_covariates, by = "donor", all.x = TRUE, sort = FALSE)
if (nrow(dat) != nrow(means)) stop("Metadata merge changed the donor-region observation count")
dat$region <- droplevels(factor(dat$region))
dat$sex <- droplevels(factor(dat$sex))
dat$donor <- droplevels(factor(dat$donor))
dat$age <- suppressWarnings(as.numeric(dat$age))
if (anyNA(dat[, c("region", "sex", "donor", "age")])) stop("Missing covariate after donor metadata merge")

design <- model.matrix(~ region + age + sex, data = dat)
region_columns <- grep("^region", colnames(design))
if (length(region_columns) < 1L || qr(design)$rank != ncol(design)) stop("Region-adjusted design matrix is not full rank")

y <- t(as.matrix(dat[, programs, drop = FALSE]))
storage.mode(y) <- "double"
colnames(y) <- paste(dat$donor, dat$region, sep = "__")

fit_status <- "ok"
failure_message <- ""
corfit <- tryCatch(
  limma::duplicateCorrelation(y, design = design, block = dat$donor),
  error = function(e) e
)

if (inherits(corfit, "error")) {
  fit_status <- "duplicateCorrelation_error"
  failure_message <- conditionMessage(corfit)
  results <- data.frame(
    program = programs,
    F_region = NA_real_,
    p_region = NA_real_,
    BH_54 = NA_real_,
    significant_BH_54 = FALSE,
    consensus_correlation = NA_real_,
    n_donor_region_observations = nrow(dat),
    n_donors = nlevels(dat$donor),
    n_regions = nlevels(dat$region),
    fixed_design_rank = qr(design)$rank,
    fixed_design_columns = ncol(design),
    formula = formula_text,
    status = fit_status,
    error = failure_message,
    stringsAsFactors = FALSE
  )
} else {
  fit <- tryCatch(
    limma::lmFit(y, design = design, block = dat$donor, correlation = corfit$consensus),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    fit_status <- "lmFit_error"
    failure_message <- conditionMessage(fit)
    results <- data.frame(
      program = programs,
      F_region = NA_real_,
      p_region = NA_real_,
      BH_54 = NA_real_,
      significant_BH_54 = FALSE,
      consensus_correlation = corfit$consensus,
      n_donor_region_observations = nrow(dat),
      n_donors = nlevels(dat$donor),
      n_regions = nlevels(dat$region),
      fixed_design_rank = qr(design)$rank,
      fixed_design_columns = ncol(design),
      formula = formula_text,
      status = fit_status,
      error = failure_message,
      stringsAsFactors = FALSE
    )
  } else {
    fit <- limma::eBayes(fit)
    joint <- limma::topTable(fit, coef = region_columns, number = Inf, sort.by = "none")
    joint <- joint[programs, , drop = FALSE]
    adj <- p.adjust(joint$P.Value, method = "BH", n = length(programs))
    results <- data.frame(
      program = rownames(joint),
      F_region = joint$F,
      p_region = joint$P.Value,
      BH_54 = adj,
      significant_BH_54 = adj < 0.05,
      consensus_correlation = corfit$consensus,
      n_donor_region_observations = nrow(dat),
      n_donors = nlevels(dat$donor),
      n_regions = nlevels(dat$region),
      fixed_design_rank = qr(design)$rank,
      fixed_design_columns = ncol(design),
      formula = formula_text,
      status = "ok",
      error = "",
      stringsAsFactors = FALSE
    )
  }
}

if (file.exists(name_path)) {
  names_tbl <- read.delim(name_path, stringsAsFactors = FALSE)
  keep <- intersect(c("new_P", "name_short", "name_full", "confidence"), names(names_tbl))
  names_tbl <- names_tbl[, keep, drop = FALSE]
  names(names_tbl)[names(names_tbl) == "new_P"] <- "program"
  results <- merge(results, names_tbl, by = "program", all.x = TRUE, sort = FALSE)
}
results <- results[order(as.integer(sub("P", "", results$program))), , drop = FALSE]

write.table(results, file.path(output_dir, "all_54_limma_dupcor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(results[results$program %in% old8, , drop = FALSE], file.path(output_dir, "old8_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
