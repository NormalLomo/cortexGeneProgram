if (!nzchar(Sys.getenv("CORTEX_PROGRAM_ROOT"))) stop("Set CORTEX_PROGRAM_ROOT; original source directories must remain read-only.")
options(future.globals.maxSize = 50 * 1024^3)
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else stop("Supply the audit root path")
remote_rlib <- file.path(root, "Rlib")
.libPaths(c(remote_rlib, .libPaths()))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lmerTest))
options(contrasts = c("contr.sum", "contr.poly"))

input_path <- file.path(root, "output", "donor_region_means.tsv")
metadata_path <- paste0(Sys.getenv("CORTEX_PROGRAM_ROOT"), "/inputs/snRNA_1M_obs.csv")
output_dir <- file.path(root, "celltype_regional", "overall_nature_lmm")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

programs <- paste0("P", seq_len(54L))
formula_text <- "score ~ region + age + sex + (1 | donor)"

means <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
required_metadata <- c("donor", "age", "sex")
if (!all(required_metadata %in% names(metadata))) stop("Source metadata lacks donor, age, or sex")
if ("PMI" %in% names(metadata) || "pmi" %in% names(metadata) || "post_mortem_interval" %in% names(metadata)) stop("PMI detected unexpectedly; update the audit specification before fitting")

donor_covariates <- unique(metadata[, required_metadata])
if (anyDuplicated(donor_covariates$donor)) {
  donor_check <- aggregate(cbind(age, sex) ~ donor, metadata[, required_metadata], function(x) length(unique(x)))
  if (any(donor_check$age > 1L | donor_check$sex > 1L)) stop("Age or sex is inconsistent within at least one donor")
  donor_covariates <- donor_covariates[!duplicated(donor_covariates$donor), , drop = FALSE]
}

dat <- merge(means, donor_covariates, by = "donor", all.x = TRUE, sort = FALSE)
if (nrow(dat) != nrow(means)) stop("Metadata merge changed the number of donor-region rows")
dat$age <- suppressWarnings(as.numeric(dat$age))
dat$sex <- droplevels(factor(dat$sex))
dat$region <- droplevels(factor(dat$region))
dat$donor <- droplevels(factor(dat$donor))
if (anyNA(dat$age) || anyNA(dat$sex) || anyNA(dat$region) || anyNA(dat$donor)) stop("Missing age, sex, region, or donor after merge")

fixed_matrix <- model.matrix(~ region + age + sex, data = dat)
fixed_rank <- qr(fixed_matrix)$rank
fixed_columns <- ncol(fixed_matrix)
full_rank <- fixed_rank == fixed_columns
if (!full_rank) stop("Overall fixed-effect design matrix is rank deficient")
if (nlevels(dat$region) < 2L || nlevels(dat$donor) < 2L || nrow(dat) <= fixed_rank) stop("Overall model is not estimable")

result_rows <- vector("list", length(programs))
for (i in seq_along(programs)) {
  program_value <- programs[[i]]
  fit_data <- dat[, c("donor", "region", "age", "sex", program_value), drop = FALSE]
  names(fit_data)[ncol(fit_data)] <- "score"
  result <- data.frame(
    program = program_value,
    formula = formula_text,
    n_combinations = nrow(dat),
    n_cells = sum(dat$n_cells),
    n_donors = nlevels(dat$donor),
    n_regions = nlevels(dat$region),
    n_age_values = length(unique(dat$age)),
    n_sex_values = length(unique(dat$sex)),
    PMI_available = FALSE,
    fixed_design_rank = fixed_rank,
    fixed_design_columns = fixed_columns,
    fixed_design_full_rank = full_rank,
    F_region = NA_real_,
    df_num_region = NA_real_,
    df_den_region = NA_real_,
    p_region = NA_real_,
    singular = NA,
    estimable = FALSE,
    status = "not_run",
    error = NA_character_,
    stringsAsFactors = FALSE
  )
  fit <- tryCatch(
    suppressWarnings(lmerTest::lmer(score ~ region + age + sex + (1 | donor), data = fit_data, REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000L)))),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    result$status <- "fit_error"
    result$error <- substr(conditionMessage(fit), 1L, 300L)
  } else {
    region_test <- tryCatch(
      suppressWarnings(as.data.frame(anova(fit, type = 3, ddf = "Satterthwaite"))),
      error = function(e) e
    )
    if (inherits(region_test, "error") || !("region" %in% rownames(region_test))) {
      result$status <- "region_test_error"
      result$error <- if (inherits(region_test, "error")) substr(conditionMessage(region_test), 1L, 300L) else "Region term absent from Satterthwaite table"
    } else {
      region_row <- region_test["region", , drop = FALSE]
      result$F_region <- as.numeric(region_row[["F value"]])
      result$df_num_region <- as.numeric(region_row[["NumDF"]])
      result$df_den_region <- as.numeric(region_row[["DenDF"]])
      result$p_region <- as.numeric(region_row[["Pr(>F)"]])
      result$singular <- lme4::isSingular(fit, tol = 1e-4)
      result$estimable <- TRUE
      result$status <- if (isTRUE(result$singular)) "ok_singular" else "ok"
    }
  }
  result_rows[[i]] <- result
}

results <- do.call(rbind, result_rows)
results$BH_54 <- p.adjust(results$p_region, method = "BH", n = 54L)
results$significant_BH_54 <- !is.na(results$BH_54) & results$BH_54 < 0.05

name_path <- file.path(root, "scripts", "program_names.tsv")
if (file.exists(name_path)) {
  program_names <- read.delim(name_path, stringsAsFactors = FALSE)
  program_names <- program_names[, intersect(c("new_P", "name_short", "name_full", "confidence"), names(program_names)), drop = FALSE]
  names(program_names)[names(program_names) == "new_P"] <- "program"
  results <- merge(results, program_names, by = "program", all.x = TRUE, sort = FALSE)
  results <- results[order(as.integer(sub("P", "", results$program))), , drop = FALSE]
}

coverage <- data.frame(
  n_combinations = nrow(dat),
  n_cells = sum(dat$n_cells),
  n_donors = nlevels(dat$donor),
  n_regions = nlevels(dat$region),
  n_age_values = length(unique(dat$age)),
  n_sex_values = length(unique(dat$sex)),
  fixed_design_rank = fixed_rank,
  fixed_design_columns = fixed_columns,
  fixed_design_full_rank = full_rank,
  PMI_available = FALSE,
  covariates = "region + age + sex; donor random intercept; PMI unavailable",
  formula = formula_text,
  stringsAsFactors = FALSE
)

write.table(results, file.path(output_dir, "all_54_overall_nature_lmm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(coverage, file.path(output_dir, "coverage_and_estimability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
