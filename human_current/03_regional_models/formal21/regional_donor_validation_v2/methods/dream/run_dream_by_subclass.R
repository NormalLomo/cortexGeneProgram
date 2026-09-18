if (!nzchar(Sys.getenv("CORTEX_PROGRAM_ROOT"))) stop("Set CORTEX_PROGRAM_ROOT; original source directories must remain read-only.")
options(future.globals.maxSize = 50 * 1024^3)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = FALSE) else normalizePath(getwd(), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_root <- normalizePath(file.path(out_dir, "..", ".."), mustWork = TRUE)
source_root <- normalizePath(file.path(analysis_root, ".."), mustWork = TRUE)
remote_rlib <- file.path(source_root, "Rlib")
configured_rlibs <- strsplit(Sys.getenv("FORMAL21_R_LIBS"), .Platform$path.sep, fixed = TRUE)[[1L]]
configured_rlibs <- configured_rlibs[nzchar(configured_rlibs) & dir.exists(configured_rlibs)]
.libPaths(unique(c(remote_rlib, configured_rlibs, .libPaths())))

suppressPackageStartupMessages(library(variancePartition))
options(digits = 17)

input_path <- file.path(source_root, "celltype_regional", "aggregated", "donor_region_subclass_means.tsv")
metadata_path <- paste0(Sys.getenv("CORTEX_PROGRAM_ROOT"), "/inputs/snRNA_1M_obs.csv")
name_path <- file.path(source_root, "scripts", "program_names.tsv")

programs <- paste0("P", seq_len(54L))
formula_text <- "score ~ region + age + sex + (1 | donor)"
model_formula <- stats::as.formula("~ region + age + sex + (1 | donor)")
response_scale_factor <- 1000
scale_tolerance <- 1e-6

if (!file.exists(input_path)) stop("Missing donor-region-subclass program mean input: ", input_path)
if (!file.exists(metadata_path)) stop("Missing metadata input: ", metadata_path)

means <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
required_mean_columns <- c("donor", "region", "subclass", "n_cells", programs)
required_metadata <- c("donor", "age", "sex")
missing_mean_columns <- setdiff(required_mean_columns, names(means))
missing_metadata_columns <- setdiff(required_metadata, names(metadata))
if (length(missing_mean_columns) > 0L) stop("Mean input is missing columns: ", paste(missing_mean_columns, collapse = ", "))
if (length(missing_metadata_columns) > 0L) stop("Metadata input is missing columns: ", paste(missing_metadata_columns, collapse = ", "))
if (anyNA(means[, c("donor", "region", "subclass", "n_cells")]) || any(means$donor == "" | means$region == "" | means$subclass == "")) stop("Mean input contains empty donor, region, subclass, or n_cells values")
if (anyDuplicated(means[, c("donor", "region", "subclass")])) stop("Mean input contains duplicate donor-region-subclass keys")
if (any(!is.finite(as.matrix(means[, programs, drop = FALSE])))) stop("Mean input contains non-finite program scores")
if (any(!is.finite(means$n_cells)) || any(means$n_cells < 0)) stop("Mean input contains invalid n_cells values")

for (covariate in c("age", "sex")) {
  values_per_donor <- tapply(metadata[[covariate]], metadata$donor, function(x) length(unique(x[!is.na(x)])))
  if (any(values_per_donor > 1L)) stop("Metadata has donor-inconsistent ", covariate, " values")
}
donor_covariates <- unique(metadata[, required_metadata, drop = FALSE])
donor_covariates <- donor_covariates[!duplicated(donor_covariates$donor), , drop = FALSE]
dat <- merge(means, donor_covariates, by = "donor", all.x = TRUE, sort = FALSE)
if (nrow(dat) != nrow(means)) stop("Metadata merge changed the donor-region-subclass row count")
dat$age <- suppressWarnings(as.numeric(dat$age))
dat$sex <- factor(dat$sex)
dat$region <- factor(dat$region)
dat$donor <- factor(dat$donor)
if (anyNA(dat[, c("age", "sex", "region", "donor")])) stop("Missing age, sex, region, or donor after metadata merge")

subclasses <- sort(unique(as.character(dat$subclass)))
if (length(subclasses) != 22L) stop("Expected 22 cortical subclasses, found ", length(subclasses))
if (!all(programs %in% names(dat))) stop("Program columns are not exactly P1 through P54")

collapse_values <- function(x) {
  if (length(x) == 0L) return("")
  paste(sort(unique(as.character(x))), collapse = ";")
}

error_message <- function(x) {
  substr(conditionMessage(x), 1L, 1000L)
}

backend_error <- function(errors, program) {
  if (is.null(errors) || is.null(names(errors)) || !(program %in% names(errors))) return("Program was not returned by dream")
  message <- paste(as.character(errors[[program]]), collapse = " ")
  if (is.na(message) || trimws(message) == "") "Program was not returned by dream without a backend message" else message
}

run_dream <- function(expression, fit_data) {
  tryCatch(
    variancePartition::dream(
      exprObj = expression,
      formula = model_formula,
      data = fit_data,
      ddf = "Satterthwaite",
      useWeights = FALSE,
      REML = TRUE,
      hideErrorsInBackend = TRUE
    ),
    error = function(e) e
  )
}

moderate_dream <- function(fit) {
  tryCatch(variancePartition::eBayes(fit), error = function(e) e)
}

extract_joint_region <- function(fit) {
  coefficient_names <- colnames(fit[["coefficients"]])
  region_coefficients <- grep("^region", coefficient_names, value = TRUE)
  if (length(region_coefficients) < 1L) return(list(error = "No region coefficients are available for the joint test"))
  joint <- tryCatch(
    variancePartition::topTable(fit, coef = region_coefficients, number = Inf, sort.by = "none"),
    error = function(e) e
  )
  if (inherits(joint, "error")) return(list(error = error_message(joint)))
  joint <- as.data.frame(joint, stringsAsFactors = FALSE)
  if (!all(c("F", "P.Value") %in% names(joint))) return(list(error = "Joint region table lacks F or P.Value"))
  list(table = joint, region_coefficients = region_coefficients)
}


fit_subclass <- function(subclass_value) {
  all_subclass <- dat[dat$subclass == subclass_value, , drop = FALSE]
  eligible <- all_subclass[all_subclass$n_cells >= 10L, , drop = FALSE]
  n_total <- nrow(all_subclass)
  n_eligible <- nrow(eligible)
  n_cells_total <- sum(all_subclass$n_cells)
  n_cells_eligible <- sum(eligible$n_cells)
  n_donors <- 0L
  n_regions <- 0L
  fixed_rank <- NA_integer_
  fixed_columns <- NA_integer_
  design_status <- "not_evaluated"
  design_error <- ""

  if (n_eligible == 0L) {
    design_status <- "no_rows_after_n_cells_filter"
  } else {
    eligible$region <- droplevels(eligible$region)
    eligible$donor <- droplevels(eligible$donor)
    eligible$sex <- droplevels(eligible$sex)
    n_donors <- nlevels(eligible$donor)
    n_regions <- nlevels(eligible$region)
    fixed_matrix <- tryCatch(model.matrix(~ region + age + sex, data = eligible), error = function(e) e)
    if (inherits(fixed_matrix, "error")) {
      design_status <- "fixed_design_error"
      design_error <- error_message(fixed_matrix)
    } else {
      fixed_rank <- qr(fixed_matrix)$rank
      fixed_columns <- ncol(fixed_matrix)
      if (fixed_rank != fixed_columns) {
        design_status <- "rank_deficient"
        design_error <- "Fixed-effect design is rank deficient"
      } else if (n_regions < 2L || n_donors < 2L || n_eligible <= fixed_rank) {
        design_status <- "insufficient_observations"
        design_error <- "Fewer than two regions or donors, or insufficient rows after the n_cells filter"
      } else {
        design_status <- "estimable"
      }
    }
  }

  rows <- data.frame(
    subclass = rep(subclass_value, length(programs)),
    program = programs,
    formula = rep(formula_text, length(programs)),
    n_combinations_total = rep(n_total, length(programs)),
    n_combinations_min_cells_10 = rep(n_eligible, length(programs)),
    n_combinations_excluded_min_cells_10 = rep(n_total - n_eligible, length(programs)),
    n_cells_total = rep(n_cells_total, length(programs)),
    n_cells_min_cells_10 = rep(n_cells_eligible, length(programs)),
    n_donors = rep(n_donors, length(programs)),
    n_regions = rep(n_regions, length(programs)),
    fixed_design_rank = rep(fixed_rank, length(programs)),
    fixed_design_columns = rep(fixed_columns, length(programs)),
    fixed_design_full_rank = rep(identical(design_status, "estimable"), length(programs)),
    minimum_cell_rule = rep("n_cells >= 10", length(programs)),
    response_scale_factor = rep(response_scale_factor, length(programs)),
    region_coefficients_tested = rep(NA_character_, length(programs)),
    F_region = rep(NA_real_, length(programs)),
    P.Value = rep(NA_real_, length(programs)),
    status = rep(design_status, length(programs)),
    error = rep(if (design_status == "estimable") "" else design_error, length(programs)),
    backend_message = rep("", length(programs)),
    stringsAsFactors = FALSE
  )

  if (design_status != "estimable") return(list(rows = rows, summary = list(
    subclass = subclass_value,
    n_combinations_total = n_total,
    n_combinations_min_cells_10 = n_eligible,
    n_combinations_excluded_min_cells_10 = n_total - n_eligible,
    n_cells_total = n_cells_total,
    n_cells_min_cells_10 = n_cells_eligible,
    n_donors = n_donors,
    n_regions = n_regions,
    fixed_design_rank = fixed_rank,
    fixed_design_columns = fixed_columns,
    design_status = design_status
  )))

  expression <- t(as.matrix(eligible[, programs, drop = FALSE]))
  storage.mode(expression) <- "double"
  colnames(expression) <- paste0("obs", seq_len(nrow(eligible)))
  rownames(eligible) <- colnames(expression)

  fit <- run_dream(expression * response_scale_factor, eligible)
  if (inherits(fit, "error")) {
    rows$status <- "dream_error"
    rows$error <- error_message(fit)
  } else {
    returned_programs <- rownames(fit[["coefficients"]])
    fit_errors <- attr(fit, "errors")
    fit <- moderate_dream(fit)
    if (inherits(fit, "error")) {
      rows$status <- "ebayes_error"
      rows$error <- error_message(fit)
    } else {
      joint <- extract_joint_region(fit)
      if (!is.null(joint$error)) {
        rows$status <- "joint_region_error"
        rows$error <- joint$error
      } else {
        rows$region_coefficients_tested <- paste(joint$region_coefficients, collapse = ";")
        for (program in programs) {
          idx <- rows$program == program
          if (!(program %in% returned_programs)) {
            rows$status[idx] <- "not_returned_by_dream"
            rows$error[idx] <- backend_error(fit_errors, program)
          } else if (!(program %in% rownames(joint$table))) {
            rows$status[idx] <- "not_returned_by_joint_region_test"
            rows$error[idx] <- "Program was returned by dream but absent from the joint-region table"
          } else {
            statistic <- joint$table[program, , drop = FALSE]
            rows$F_region[idx] <- statistic$F[[1L]]
            rows$P.Value[idx] <- statistic$P.Value[[1L]]
            if (!is.finite(rows$F_region[idx]) || !is.finite(rows$P.Value[idx])) {
              rows$status[idx] <- "invalid_joint_statistics"
              rows$error[idx] <- "Joint region F or P.Value is non-finite"
            } else {
              rows$status[idx] <- "ok"
              rows$error[idx] <- ""
              if (!is.null(fit_errors) && !is.null(names(fit_errors)) && program %in% names(fit_errors)) rows$backend_message[idx] <- paste(as.character(fit_errors[[program]]), collapse = " ")
            }
          }
        }
      }
    }
  }
  list(rows = rows, summary = list(
    subclass = subclass_value,
    n_combinations_total = n_total,
    n_combinations_min_cells_10 = n_eligible,
    n_combinations_excluded_min_cells_10 = n_total - n_eligible,
    n_cells_total = n_cells_total,
    n_cells_min_cells_10 = n_cells_eligible,
    n_donors = n_donors,
    n_regions = n_regions,
    fixed_design_rank = fixed_rank,
    fixed_design_columns = fixed_columns,
    design_status = design_status
  ))
}

subclass_fits <- lapply(subclasses, fit_subclass)
results <- do.call(rbind, lapply(subclass_fits, function(x) x$rows))
summary_input <- lapply(subclass_fits, function(x) x$summary)

subclass_summary <- do.call(rbind, lapply(summary_input, function(x) {
  data.frame(
    subclass = x$subclass,
    n_combinations_total = x$n_combinations_total,
    n_combinations_min_cells_10 = x$n_combinations_min_cells_10,
    n_combinations_excluded_min_cells_10 = x$n_combinations_excluded_min_cells_10,
    n_cells_total = x$n_cells_total,
    n_cells_min_cells_10 = x$n_cells_min_cells_10,
    n_donors = x$n_donors,
    n_regions = x$n_regions,
    fixed_design_rank = x$fixed_design_rank,
    fixed_design_columns = x$fixed_design_columns,
    design_status = x$design_status,
    stringsAsFactors = FALSE
  )
}))


if (nrow(results) != length(subclasses) * length(programs)) stop("Result grid does not contain 54 x 22 rows")
if (anyNA(results$subclass) || anyNA(results$program) || any(results$subclass == "" | results$program == "")) stop("Result grid contains an empty program or subclass identifier")
if (anyDuplicated(results[, c("subclass", "program")])) stop("Result grid contains duplicate program-subclass keys")

results$BH_within_subclass_54 <- NA_real_
for (subclass_value in subclasses) {
  idx <- results$subclass == subclass_value
  results$BH_within_subclass_54[idx] <- stats::p.adjust(results$P.Value[idx], method = "BH", n = 54L)
}
results$BH_global_22x54 <- stats::p.adjust(results$P.Value, method = "BH", n = 54L * 22L)
results$significant_BH_within_subclass_54 <- !is.na(results$BH_within_subclass_54) & results$BH_within_subclass_54 < 0.05
results$significant_BH_global_22x54 <- !is.na(results$BH_global_22x54) & results$BH_global_22x54 < 0.05

key_statistics <- c("F_region", "P.Value", "BH_within_subclass_54", "BH_global_22x54")
failed_rows <- results$status != "ok"
if (any(rowSums(is.na(results[, key_statistics, drop = FALSE])) > 0L & !failed_rows)) stop("An ok result row has missing key statistics")
if (any(failed_rows & (is.na(results$error) | trimws(results$error) == ""))) stop("A failed result row lacks an explicit error")

program_annotations <- data.frame(program = programs, program_name = NA_character_, program_name_full = NA_character_, annotation_confidence = NA_character_, stringsAsFactors = FALSE)
if (file.exists(name_path)) {
  annotations <- read.delim(name_path, check.names = FALSE, stringsAsFactors = FALSE)
  annotation_columns <- intersect(c("new_P", "name_short", "name_full", "confidence"), names(annotations))
  if ("new_P" %in% annotation_columns) {
    annotations <- annotations[, annotation_columns, drop = FALSE]
    names(annotations)[names(annotations) == "new_P"] <- "program"
    annotations <- annotations[!duplicated(annotations$program), , drop = FALSE]
    annotations <- annotations[match(programs, annotations$program), , drop = FALSE]
    if ("name_short" %in% names(annotations)) program_annotations$program_name <- annotations$name_short
    if ("name_full" %in% names(annotations)) program_annotations$program_name_full <- annotations$name_full
    if ("confidence" %in% names(annotations)) program_annotations$annotation_confidence <- annotations$confidence
  }
}
results <- merge(results, program_annotations, by = "program", all.x = TRUE, sort = FALSE)
results <- results[order(results$subclass, as.integer(sub("P", "", results$program))), , drop = FALSE]

program_subclass_counts <- do.call(rbind, lapply(programs, function(program_value) {
  x <- results[results$program == program_value, , drop = FALSE]
  data.frame(
    program = program_value,
    program_name = x$program_name[[1L]],
    n_subclasses_total = length(subclasses),
    n_subclasses_estimable = sum(x$status == "ok"),
    n_subclasses_significant_BH_within_subclass_54 = sum(x$significant_BH_within_subclass_54),
    n_subclasses_significant_BH_global_22x54 = sum(x$significant_BH_global_22x54),
    subclasses_significant_BH_within_subclass_54 = collapse_values(x$subclass[x$significant_BH_within_subclass_54]),
    subclasses_significant_BH_global_22x54 = collapse_values(x$subclass[x$significant_BH_global_22x54]),
    integrated_heatmap_primary = "BH_global_22x54",
    stringsAsFactors = FALSE
  )
}))

subclass_summary$n_programs_total <- 54L
subclass_summary$n_programs_estimable <- vapply(subclass_summary$subclass, function(subclass_value) sum(results$subclass == subclass_value & results$status == "ok"), integer(1L))
subclass_summary$n_programs_failed <- 54L - subclass_summary$n_programs_estimable
subclass_summary$n_programs_significant_BH_within_subclass_54 <- vapply(subclass_summary$subclass, function(subclass_value) sum(results$subclass == subclass_value & results$significant_BH_within_subclass_54), integer(1L))
subclass_summary$n_programs_significant_BH_global_22x54 <- vapply(subclass_summary$subclass, function(subclass_value) sum(results$subclass == subclass_value & results$significant_BH_global_22x54), integer(1L))
subclass_summary$programs_significant_BH_within_subclass_54 <- vapply(subclass_summary$subclass, function(subclass_value) collapse_values(results$program[results$subclass == subclass_value & results$significant_BH_within_subclass_54]), character(1L))
subclass_summary$programs_significant_BH_global_22x54 <- vapply(subclass_summary$subclass, function(subclass_value) collapse_values(results$program[results$subclass == subclass_value & results$significant_BH_global_22x54]), character(1L))
subclass_summary <- subclass_summary[order(subclass_summary$subclass), , drop = FALSE]

write.table(results, file.path(out_dir, "all_54_by_22_dream.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(program_subclass_counts, file.path(out_dir, "program_subclass_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(subclass_summary, file.path(out_dir, "subclass_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")


cat("Completed donor-aware dream analysis for", nrow(results), "program-subclass rows.\n")
