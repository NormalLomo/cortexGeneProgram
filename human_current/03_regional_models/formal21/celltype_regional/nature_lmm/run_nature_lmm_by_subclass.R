if (!nzchar(Sys.getenv("CORTEX_PROGRAM_ROOT"))) stop("Set CORTEX_PROGRAM_ROOT; original source directories must remain read-only.")
options(future.globals.maxSize = 50 * 1024^3)
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else stop("Supply the audit root path")
remote_rlib <- file.path(root, "Rlib")
.libPaths(c(remote_rlib, .libPaths()))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lmerTest))
options(contrasts = c("contr.sum", "contr.poly"))

input_path <- file.path(root, "celltype_regional", "aggregated", "donor_region_subclass_means.tsv")
metadata_path <- paste0(Sys.getenv("CORTEX_PROGRAM_ROOT"), "/inputs/snRNA_1M_obs.csv")
output_dir <- file.path(root, "celltype_regional", "nature_lmm")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

programs <- paste0("P", seq_len(54L))
original_eight <- c("P1", "P3", "P4", "P6", "P8", "P9", "P13", "P33")
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
if (nrow(dat) != nrow(means)) stop("Metadata merge changed the number of donor-region-subclass rows")
dat$age <- suppressWarnings(as.numeric(dat$age))
dat$sex <- factor(dat$sex)
dat$region <- factor(dat$region)
dat$donor <- factor(dat$donor)
if (anyNA(dat$age) || anyNA(dat$sex) || anyNA(dat$region) || anyNA(dat$donor)) stop("Missing age, sex, region, or donor after merge")

subclasses <- sort(unique(dat$subclass))
result_template <- expand.grid(subclass = subclasses, program = programs, stringsAsFactors = FALSE)
result_template <- result_template[order(result_template$subclass, as.integer(sub("P", "", result_template$program))), , drop = FALSE]

coverage_rows <- vector("list", length(subclasses))
result_rows <- vector("list", nrow(result_template))
row_index <- 1L

for (subclass_value in subclasses) {
  all_subclass <- dat[dat$subclass == subclass_value, , drop = FALSE]
  eligible <- all_subclass[all_subclass$n_cells >= 10L, , drop = FALSE]
  eligible$region <- droplevels(eligible$region)
  eligible$donor <- droplevels(eligible$donor)
  eligible$sex <- droplevels(eligible$sex)
  n_total <- nrow(all_subclass)
  n_eligible <- nrow(eligible)
  n_excluded <- n_total - n_eligible
  n_donors <- if (n_eligible > 0L) nlevels(droplevels(eligible$donor)) else 0L
  n_regions <- if (n_eligible > 0L) nlevels(droplevels(eligible$region)) else 0L
  n_age_values <- if (n_eligible > 0L) length(unique(eligible$age)) else 0L
  n_sex_values <- if (n_eligible > 0L) length(unique(eligible$sex)) else 0L
  fixed_rank <- NA_integer_
  fixed_columns <- NA_integer_
  full_rank <- FALSE
  design_status <- "not_evaluated"
  if (n_eligible > 0L) {
    fixed_matrix <- model.matrix(~ region + age + sex, data = eligible)
    fixed_rank <- qr(fixed_matrix)$rank
    fixed_columns <- ncol(fixed_matrix)
    full_rank <- fixed_rank == fixed_columns
    design_status <- if (!full_rank) "rank_deficient" else if (n_regions < 2L || n_donors < 2L || n_eligible <= fixed_rank) "insufficient_observations" else "estimable"
  }
  coverage_rows[[which(subclasses == subclass_value)]] <- data.frame(
    subclass = subclass_value,
    n_combinations_total = n_total,
    n_combinations_min_cells_10 = n_eligible,
    n_combinations_excluded_min_cells_10 = n_excluded,
    n_cells_total = sum(all_subclass$n_cells),
    n_cells_min_cells_10 = sum(eligible$n_cells),
    n_donors = n_donors,
    n_regions = n_regions,
    n_age_values = n_age_values,
    n_sex_values = n_sex_values,
    fixed_design_rank = fixed_rank,
    fixed_design_columns = fixed_columns,
    fixed_design_full_rank = full_rank,
    PMI_available = FALSE,
    covariates = "region + age + sex; donor random intercept; PMI unavailable",
    design_status = design_status,
    stringsAsFactors = FALSE
  )
  for (program_value in programs) {
    result <- data.frame(
      subclass = subclass_value,
      program = program_value,
      formula = formula_text,
      n_combinations_total = n_total,
      n_combinations_min_cells_10 = n_eligible,
      n_combinations_excluded_min_cells_10 = n_excluded,
      n_donors = n_donors,
      n_regions = n_regions,
      n_age_values = n_age_values,
      n_sex_values = n_sex_values,
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
      status = design_status,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
    if (identical(design_status, "estimable")) {
      fit_data <- eligible[, c("donor", "region", "age", "sex", program_value), drop = FALSE]
      names(fit_data)[ncol(fit_data)] <- "score"
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
    }
    result_rows[[row_index]] <- result
    row_index <- row_index + 1L
  }
}

results <- do.call(rbind, result_rows)
coverage <- do.call(rbind, coverage_rows)
results$BH_within_subclass_54 <- NA_real_
for (subclass_value in subclasses) {
  idx <- results$subclass == subclass_value
  results$BH_within_subclass_54[idx] <- p.adjust(results$p_region[idx], method = "BH", n = 54L)
}
results$BH_global_22x54 <- p.adjust(results$p_region, method = "BH", n = length(subclasses) * length(programs))
results$significant_BH_within_subclass_54 <- !is.na(results$BH_within_subclass_54) & results$BH_within_subclass_54 < 0.05
results$significant_BH_global_22x54 <- !is.na(results$BH_global_22x54) & results$BH_global_22x54 < 0.05

name_path <- file.path(root, "scripts", "program_names.tsv")
if (file.exists(name_path)) {
  program_names <- read.delim(name_path, stringsAsFactors = FALSE)
  program_names <- program_names[, intersect(c("new_P", "name_short", "name_full", "confidence"), names(program_names)), drop = FALSE]
  names(program_names)[names(program_names) == "new_P"] <- "program"
  results <- merge(results, program_names, by = "program", all.x = TRUE, sort = FALSE)
  results <- results[order(results$subclass, as.integer(sub("P", "", results$program))), , drop = FALSE]
}

collapse_subclasses <- function(x) if (length(x) == 0L) "" else paste(sort(unique(x)), collapse = ";")
program_summary <- do.call(rbind, lapply(programs, function(program_value) {
  x <- results[results$program == program_value, , drop = FALSE]
  data.frame(
    program = program_value,
    program_name = if ("name_short" %in% names(x)) x$name_short[[1L]] else NA_character_,
    n_subclasses_total = length(subclasses),
    n_subclasses_estimable = sum(x$estimable),
    n_subclasses_BH_within_subclass_54 = sum(x$significant_BH_within_subclass_54),
    n_subclasses_BH_global_22x54 = sum(x$significant_BH_global_22x54),
    subclasses_BH_within_subclass_54 = collapse_subclasses(x$subclass[x$significant_BH_within_subclass_54]),
    subclasses_BH_global_22x54 = collapse_subclasses(x$subclass[x$significant_BH_global_22x54]),
    stringsAsFactors = FALSE
  )
}))

subclass_summary <- do.call(rbind, lapply(subclasses, function(subclass_value) {
  x <- results[results$subclass == subclass_value, , drop = FALSE]
  data.frame(
    subclass = subclass_value,
    n_programs_total = length(programs),
    n_programs_estimable = sum(x$estimable),
    n_programs_BH_within_subclass_54 = sum(x$significant_BH_within_subclass_54),
    n_programs_BH_global_22x54 = sum(x$significant_BH_global_22x54),
    programs_BH_within_subclass_54 = collapse_subclasses(x$program[x$significant_BH_within_subclass_54]),
    programs_BH_global_22x54 = collapse_subclasses(x$program[x$significant_BH_global_22x54]),
    stringsAsFactors = FALSE
  )
}))

original_eight_detail <- results[results$program %in% original_eight, , drop = FALSE]
write.table(results, file.path(output_dir, "all_54_by_22_nature_lmm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(coverage, file.path(output_dir, "coverage_and_estimability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(original_eight_detail, file.path(output_dir, "original_eight_detail.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(program_summary, file.path(output_dir, "program_significance_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(subclass_summary, file.path(output_dir, "subclass_significance_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
