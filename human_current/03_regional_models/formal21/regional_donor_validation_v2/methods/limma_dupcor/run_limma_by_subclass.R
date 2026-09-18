if (!nzchar(Sys.getenv("CORTEX_PROGRAM_ROOT"))) stop("Set CORTEX_PROGRAM_ROOT; original source directories must remain read-only.")
options(future.globals.maxSize = 50 * 1024^3)
#!/usr/bin/env Rscript

default_source_root <- Sys.getenv("FORMAL21_MODEL_ROOT")
default_metadata_path <- paste0(Sys.getenv("CORTEX_PROGRAM_ROOT"), "/inputs/snRNA_1M_obs.csv")
expected_programs <- paste0("P", seq_len(54L))
expected_subclass_count <- 22L
design_formula_text <- "~ region + age + sex"
block_text <- "donor"

abort <- function(...) stop(..., call. = FALSE)

script_directory <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) abort("Could not determine the script path")
  dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
}

short_error <- function(x) substr(as.character(x), 1L, 500L)

empty_fit <- function(programs, status, error, consensus_correlation = NA_real_,
                      fixed_design_rank = NA_integer_, fixed_design_columns = NA_integer_,
                      fixed_design_full_rank = NA, n_region_coefficients = NA_integer_) {
  results <- data.frame(
    F_region = rep(NA_real_, length(programs)),
    df_num_region = rep(NA_real_, length(programs)),
    df_den_region = rep(NA_real_, length(programs)),
    p_region = rep(NA_real_, length(programs)),
    consensus_correlation = rep(consensus_correlation, length(programs)),
    estimable = rep(FALSE, length(programs)),
    status = rep(status, length(programs)),
    error = rep(error, length(programs)),
    row.names = programs,
    stringsAsFactors = FALSE
  )
  list(
    results = results,
    consensus_correlation = consensus_correlation,
    fixed_design_rank = fixed_design_rank,
    fixed_design_columns = fixed_design_columns,
    fixed_design_full_rank = fixed_design_full_rank,
    n_region_coefficients = n_region_coefficients,
    fit_status = status
  )
}

fit_subclass_limma <- function(eligible, programs) {
  required_columns <- c("donor", "region", "age", "sex", programs)
  missing_columns <- setdiff(required_columns, names(eligible))
  if (length(missing_columns) > 0L) {
    return(empty_fit(
      programs, "missing_required_columns",
      paste("Missing columns:", paste(missing_columns, collapse = ", "))
    ))
  }
  if (nrow(eligible) == 0L) {
    return(empty_fit(programs, "no_eligible_observations", "No donor-region observations passed n_cells >= 10"))
  }

  fit_data <- eligible[, required_columns, drop = FALSE]
  id_columns <- c("donor", "region", "sex")
  if (any(vapply(fit_data[id_columns], function(x) anyNA(x) || any(trimws(as.character(x)) == ""), logical(1L)))) {
    return(empty_fit(programs, "missing_covariate", "At least one donor, region, or sex value is missing or empty"))
  }
  fit_data$age <- suppressWarnings(as.numeric(fit_data$age))
  if (anyNA(fit_data$age) || any(!is.finite(fit_data$age))) {
    return(empty_fit(programs, "missing_or_nonfinite_age", "At least one age value is missing or non-finite"))
  }
  fit_data$donor <- droplevels(factor(as.character(fit_data$donor)))
  fit_data$region <- droplevels(factor(as.character(fit_data$region)))
  fit_data$sex <- droplevels(factor(as.character(fit_data$sex)))

  n_donors <- nlevels(fit_data$donor)
  n_regions <- nlevels(fit_data$region)
  n_sex_values <- nlevels(fit_data$sex)
  if (n_regions < 2L) {
    return(empty_fit(programs, "insufficient_regions", "Fewer than two retained regions"))
  }
  if (n_donors < 2L) {
    return(empty_fit(programs, "insufficient_donors", "Fewer than two retained donors"))
  }
  if (n_sex_values < 2L) {
    return(empty_fit(programs, "single_sex_level", "Sex has fewer than two levels after filtering"))
  }

  design <- tryCatch(
    stats::model.matrix(~ region + age + sex, data = fit_data),
    error = function(e) e
  )
  if (inherits(design, "error")) {
    return(empty_fit(programs, "design_matrix_error", short_error(conditionMessage(design))))
  }
  fixed_design_rank <- qr(design)$rank
  fixed_design_columns <- ncol(design)
  fixed_design_full_rank <- fixed_design_rank == fixed_design_columns
  region_columns <- grep("^region", colnames(design))
  n_region_coefficients <- length(region_columns)
  if (n_region_coefficients < 1L) {
    return(empty_fit(
      programs, "region_term_absent", "The design matrix has no estimable region coefficient",
      fixed_design_rank = fixed_design_rank, fixed_design_columns = fixed_design_columns,
      fixed_design_full_rank = fixed_design_full_rank, n_region_coefficients = n_region_coefficients
    ))
  }
  if (!fixed_design_full_rank) {
    return(empty_fit(
      programs, "rank_deficient_design", "The region + age + sex design matrix is rank deficient",
      fixed_design_rank = fixed_design_rank, fixed_design_columns = fixed_design_columns,
      fixed_design_full_rank = fixed_design_full_rank, n_region_coefficients = n_region_coefficients
    ))
  }
  if (nrow(fit_data) <= fixed_design_rank) {
    return(empty_fit(
      programs, "insufficient_residual_df", "Retained observations do not exceed the fixed-effect design rank",
      fixed_design_rank = fixed_design_rank, fixed_design_columns = fixed_design_columns,
      fixed_design_full_rank = fixed_design_full_rank, n_region_coefficients = n_region_coefficients
    ))
  }

  y <- t(as.matrix(fit_data[, programs, drop = FALSE]))
  storage.mode(y) <- "double"
  rownames(y) <- programs
  colnames(y) <- paste(as.character(fit_data$donor), as.character(fit_data$region), seq_len(nrow(fit_data)), sep = "__")
  nonfinite_programs <- programs[rowSums(!is.finite(y)) > 0L]
  if (length(nonfinite_programs) > 0L) {
    return(empty_fit(
      programs, "nonfinite_program_scores",
      paste("Non-finite scores prevent a complete 54-program matrix:", paste(nonfinite_programs, collapse = ", ")),
      fixed_design_rank = fixed_design_rank, fixed_design_columns = fixed_design_columns,
      fixed_design_full_rank = fixed_design_full_rank, n_region_coefficients = n_region_coefficients
    ))
  }

  corfit <- tryCatch(
    limma::duplicateCorrelation(y, design = design, block = fit_data$donor),
    error = function(e) e
  )
  if (inherits(corfit, "error")) {
    return(empty_fit(
      programs, "duplicateCorrelation_error", short_error(conditionMessage(corfit)),
      fixed_design_rank = fixed_design_rank, fixed_design_columns = fixed_design_columns,
      fixed_design_full_rank = fixed_design_full_rank, n_region_coefficients = n_region_coefficients
    ))
  }
  consensus_correlation <- as.numeric(corfit$consensus)
  if (length(consensus_correlation) != 1L || !is.finite(consensus_correlation)) {
    return(empty_fit(
      programs, "invalid_consensus_correlation", "duplicateCorrelation returned a missing or non-finite consensus correlation",
      fixed_design_rank = fixed_design_rank, fixed_design_columns = fixed_design_columns,
      fixed_design_full_rank = fixed_design_full_rank, n_region_coefficients = n_region_coefficients
    ))
  }

  fit <- tryCatch(
    limma::lmFit(y, design = design, block = fit_data$donor, correlation = consensus_correlation),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(empty_fit(
      programs, "lmFit_error", short_error(conditionMessage(fit)), consensus_correlation,
      fixed_design_rank, fixed_design_columns, fixed_design_full_rank, n_region_coefficients
    ))
  }
  fit <- tryCatch(limma::eBayes(fit), error = function(e) e)
  if (inherits(fit, "error")) {
    return(empty_fit(
      programs, "eBayes_error", short_error(conditionMessage(fit)), consensus_correlation,
      fixed_design_rank, fixed_design_columns, fixed_design_full_rank, n_region_coefficients
    ))
  }

  joint <- tryCatch({
    if (n_region_coefficients == 1L) {
      single_term <- limma::topTable(fit, coef = region_columns, number = Inf, sort.by = "none")
      data.frame(F = single_term$t^2, P.Value = single_term$P.Value, row.names = rownames(single_term))
    } else {
      multi_term <- limma::topTable(fit, coef = region_columns, number = Inf, sort.by = "none")
      if (!all(c("F", "P.Value") %in% names(multi_term))) abort("limma did not return moderated F and P.Value for the region coefficient block")
      multi_term[, c("F", "P.Value"), drop = FALSE]
    }
  }, error = function(e) e)
  if (inherits(joint, "error")) {
    return(empty_fit(
      programs, "joint_region_test_error", short_error(conditionMessage(joint)), consensus_correlation,
      fixed_design_rank, fixed_design_columns, fixed_design_full_rank, n_region_coefficients
    ))
  }

  joint_index <- match(programs, rownames(joint))
  fit_index <- match(programs, rownames(fit$coefficients))
  F_region <- rep(NA_real_, length(programs))
  p_region <- rep(NA_real_, length(programs))
  df_den_region <- rep(NA_real_, length(programs))
  present <- !is.na(joint_index) & !is.na(fit_index)
  F_region[present] <- as.numeric(joint$F[joint_index[present]])
  p_region[present] <- as.numeric(joint$P.Value[joint_index[present]])
  df_den_region[present] <- as.numeric(fit$df.total[fit_index[present]])
  valid <- present & is.finite(F_region) & is.finite(p_region) & p_region >= 0 & p_region <= 1 & is.finite(df_den_region) & df_den_region > 0
  F_region[!is.finite(F_region)] <- NA_real_
  p_region[!is.finite(p_region)] <- NA_real_
  df_den_region[!is.finite(df_den_region)] <- NA_real_
  status <- ifelse(valid, "ok", ifelse(!present, "missing_joint_result", "nonfinite_or_invalid_joint_statistic"))
  error <- ifelse(
    valid,
    "none",
    ifelse(!present, "limma did not return this program in the region test table", "limma returned a missing, non-finite, or invalid moderated joint statistic")
  )
  results <- data.frame(
    F_region = F_region,
    df_num_region = rep(as.numeric(n_region_coefficients), length(programs)),
    df_den_region = df_den_region,
    p_region = p_region,
    consensus_correlation = rep(consensus_correlation, length(programs)),
    estimable = valid,
    status = status,
    error = error,
    row.names = programs,
    stringsAsFactors = FALSE
  )
  list(
    results = results,
    consensus_correlation = consensus_correlation,
    fixed_design_rank = fixed_design_rank,
    fixed_design_columns = fixed_design_columns,
    fixed_design_full_rank = fixed_design_full_rank,
    n_region_coefficients = n_region_coefficients,
    fit_status = if (all(valid)) "ok" else "partial_nonfinite_or_missing_statistics"
  )
}

read_donor_covariates <- function(metadata_path) {
  metadata <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  required_metadata <- c("donor", "age", "sex")
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata) > 0L) abort("Source metadata lacks: ", paste(missing_metadata, collapse = ", "))
  metadata <- metadata[, required_metadata, drop = FALSE]
  metadata$donor <- trimws(as.character(metadata$donor))
  metadata$sex <- trimws(as.character(metadata$sex))
  metadata$age <- suppressWarnings(as.numeric(metadata$age))
  if (anyNA(metadata$donor) || any(metadata$donor == "")) abort("Metadata contains a missing or empty donor ID")
  if (anyNA(metadata$sex) || any(metadata$sex == "")) abort("Metadata contains a missing or empty sex value")
  if (anyNA(metadata$age) || any(!is.finite(metadata$age))) abort("Metadata contains a missing or non-finite age value")
  age_consistent <- vapply(split(metadata$age, metadata$donor), function(x) length(unique(x)) == 1L, logical(1L))
  sex_consistent <- vapply(split(metadata$sex, metadata$donor), function(x) length(unique(x)) == 1L, logical(1L))
  inconsistent <- union(names(age_consistent)[!age_consistent], names(sex_consistent)[!sex_consistent])
  if (length(inconsistent) > 0L) abort("Age or sex is inconsistent within donor(s): ", paste(inconsistent, collapse = ", "))
  metadata[!duplicated(metadata$donor), , drop = FALSE]
}

read_program_annotations <- function(name_path, programs) {
  annotations <- data.frame(
    program = programs,
    program_name_short = NA_character_,
    program_name_full = NA_character_,
    program_confidence = NA_character_,
    stringsAsFactors = FALSE
  )
  if (!file.exists(name_path)) return(annotations)
  name_table <- utils::read.delim(name_path, check.names = FALSE, stringsAsFactors = FALSE)
  id_column <- intersect(c("new_P", "program"), names(name_table))
  if (length(id_column) != 1L) abort("Program-name table lacks exactly one supported identifier column (new_P or program)")
  name_table$program <- as.character(name_table[[id_column]])
  name_table <- name_table[name_table$program %in% programs, , drop = FALSE]
  if (anyDuplicated(name_table$program)) abort("Program-name table has duplicate program identifiers")
  mapping <- c(name_short = "program_name_short", name_full = "program_name_full", confidence = "program_confidence")
  name_index <- match(programs, name_table$program)
  for (source_name in names(mapping)) {
    if (source_name %in% names(name_table)) annotations[[mapping[[source_name]]]] <- as.character(name_table[[source_name]][name_index])
  }
  annotations
}

collapse_ids <- function(x) {
  if (length(x) == 0L) "none" else paste(sort(unique(as.character(x))), collapse = ";")
}

assert_results <- function(results, programs, subclasses) {
  expected_rows <- length(programs) * length(subclasses)
  if (nrow(results) != expected_rows) abort("Result row count is ", nrow(results), ", expected ", expected_rows)
  if (length(unique(results$program)) != length(programs) || !setequal(unique(results$program), programs)) abort("Result table does not contain exactly the expected 54 programs")
  if (length(unique(results$subclass)) != length(subclasses) || !setequal(unique(results$subclass), subclasses)) abort("Result table does not contain exactly the expected subclasses")
  key <- paste(results$program, results$subclass, sep = "\r")
  if (anyDuplicated(key)) abort("Result table has duplicate program-subclass keys")
  if (anyNA(results$program) || anyNA(results$subclass) || any(trimws(results$program) == "") || any(trimws(results$subclass) == "")) abort("Result table has a missing or empty program/subclass ID")
  if (anyNA(results$status) || any(results$status == "") || anyNA(results$error) || any(results$error == "")) abort("Every result row must state status and error")
  if (!all(results$estimable == (results$status == "ok"))) abort("Estimability and status are inconsistent")
  success_columns <- c("F_region", "df_num_region", "df_den_region", "p_region", "BH_within_subclass_54", "BH_global_22x54")
  if (any(results$estimable & !stats::complete.cases(results[, success_columns, drop = FALSE]))) abort("An estimable row has a missing key statistic")
  if (any(!results$estimable & results$error == "none")) abort("A failed row does not explain its failure")

}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

main <- function() {
  source_root <- Sys.getenv("LIMMA_DUCPOR_SOURCE_ROOT", unset = default_source_root)
  if (!nzchar(source_root)) source_root <- default_source_root
  if (!nzchar(source_root)) abort("Set FORMAL21_MODEL_ROOT or LIMMA_DUCPOR_SOURCE_ROOT.")
  output_dir <- Sys.getenv("LIMMA_DUCPOR_OUTPUT_DIR", unset = file.path(source_root, "regional_donor_validation_v2/methods/limma_dupcor"))
  if (!nzchar(output_dir)) output_dir <- file.path(source_root, "regional_donor_validation_v2/methods/limma_dupcor")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metadata_path <- Sys.getenv("CORTEX_METADATA_PATH", unset = default_metadata_path)
  source_root <- normalizePath(source_root, mustWork = TRUE)
  metadata_path <- normalizePath(metadata_path, mustWork = TRUE)
  rlib <- file.path(source_root, "Rlib")
  if (dir.exists(rlib)) .libPaths(c(rlib, .libPaths()))
  suppressPackageStartupMessages(library(limma))

  input_path <- file.path(source_root, "celltype_regional", "aggregated", "donor_region_subclass_means.tsv")
  name_path <- file.path(source_root, "scripts", "program_names.tsv")
  if (!file.exists(input_path)) abort("Input donor-region-subclass means file does not exist: ", input_path)

  means <- utils::read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
  required_input <- c("donor", "region", "subclass", "n_cells", expected_programs)
  missing_input <- setdiff(required_input, names(means))
  if (length(missing_input) > 0L) abort("Input means lacks columns: ", paste(missing_input, collapse = ", "))
  for (id_column in c("donor", "region", "subclass")) {
    means[[id_column]] <- trimws(as.character(means[[id_column]]))
    if (anyNA(means[[id_column]]) || any(means[[id_column]] == "")) abort("Input means has a missing or empty ", id_column, " ID")
  }
  means$n_cells <- suppressWarnings(as.numeric(means$n_cells))
  if (anyNA(means$n_cells) || any(!is.finite(means$n_cells)) || any(means$n_cells < 0)) abort("Input means has an invalid n_cells value")
  input_key <- paste(means$donor, means$region, means$subclass, sep = "\r")
  if (anyDuplicated(input_key)) abort("Input means has duplicate donor-region-subclass rows")

  donor_covariates <- read_donor_covariates(metadata_path)
  covariate_index <- match(means$donor, donor_covariates$donor)
  if (anyNA(covariate_index)) abort("At least one donor in the means table is absent from source metadata")
  dat <- means
  dat$age <- donor_covariates$age[covariate_index]
  dat$sex <- donor_covariates$sex[covariate_index]
  subclasses <- sort(unique(dat$subclass))
  if (length(subclasses) != expected_subclass_count) abort("Observed ", length(subclasses), " subclasses, expected ", expected_subclass_count)

  annotation <- read_program_annotations(name_path, expected_programs)
  result_rows <- vector("list", length(subclasses))
  coverage_rows <- vector("list", length(subclasses))
  for (i in seq_along(subclasses)) {
    subclass_value <- subclasses[[i]]
    all_subclass <- dat[dat$subclass == subclass_value, , drop = FALSE]
    eligible <- all_subclass[all_subclass$n_cells >= 10L, , drop = FALSE]
    fit <- fit_subclass_limma(eligible, expected_programs)
    subclass_results <- fit$results
    subclass_results$program <- rownames(subclass_results)
    rownames(subclass_results) <- NULL
    subclass_results$subclass <- subclass_value
    subclass_results$design_formula <- design_formula_text
    subclass_results$block <- block_text
    subclass_results$n_combinations_total <- nrow(all_subclass)
    subclass_results$n_combinations_min_cells_10 <- nrow(eligible)
    subclass_results$n_combinations_excluded_min_cells_10 <- nrow(all_subclass) - nrow(eligible)
    subclass_results$n_donors <- if (nrow(eligible) == 0L) 0L else length(unique(eligible$donor))
    subclass_results$n_regions <- if (nrow(eligible) == 0L) 0L else length(unique(eligible$region))
    subclass_results$n_age_values <- if (nrow(eligible) == 0L) 0L else length(unique(eligible$age))
    subclass_results$n_sex_values <- if (nrow(eligible) == 0L) 0L else length(unique(eligible$sex))
    subclass_results$fixed_design_rank <- fit$fixed_design_rank
    subclass_results$fixed_design_columns <- fit$fixed_design_columns
    subclass_results$fixed_design_full_rank <- fit$fixed_design_full_rank
    subclass_results$n_region_coefficients <- fit$n_region_coefficients
    result_rows[[i]] <- subclass_results[, c(
      "subclass", "program", "design_formula", "block",
      "n_combinations_total", "n_combinations_min_cells_10", "n_combinations_excluded_min_cells_10",
      "n_donors", "n_regions", "n_age_values", "n_sex_values",
      "fixed_design_rank", "fixed_design_columns", "fixed_design_full_rank", "n_region_coefficients",
      "F_region", "df_num_region", "df_den_region", "p_region", "consensus_correlation",
      "estimable", "status", "error"
    )]
    coverage_rows[[i]] <- data.frame(
      subclass = subclass_value,
      n_combinations_total = nrow(all_subclass),
      n_combinations_min_cells_10 = nrow(eligible),
      n_combinations_excluded_min_cells_10 = nrow(all_subclass) - nrow(eligible),
      n_cells_total = sum(all_subclass$n_cells),
      n_cells_min_cells_10 = sum(eligible$n_cells),
      n_donors = if (nrow(eligible) == 0L) 0L else length(unique(eligible$donor)),
      n_regions = if (nrow(eligible) == 0L) 0L else length(unique(eligible$region)),
      n_age_values = if (nrow(eligible) == 0L) 0L else length(unique(eligible$age)),
      n_sex_values = if (nrow(eligible) == 0L) 0L else length(unique(eligible$sex)),
      fixed_design_rank = fit$fixed_design_rank,
      fixed_design_columns = fit$fixed_design_columns,
      fixed_design_full_rank = fit$fixed_design_full_rank,
      n_region_coefficients = fit$n_region_coefficients,
      analysis_status = fit$fit_status,
      consensus_correlation = fit$consensus_correlation,
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, result_rows)
  rownames(results) <- NULL
  results$BH_within_subclass_54 <- NA_real_
  for (subclass_value in subclasses) {
    idx <- results$subclass == subclass_value
    results$BH_within_subclass_54[idx] <- stats::p.adjust(results$p_region[idx], method = "BH", n = length(expected_programs))
  }
  results$BH_global_22x54 <- stats::p.adjust(results$p_region, method = "BH", n = length(expected_programs) * length(subclasses))
  results$significant_BH_within_subclass_54 <- !is.na(results$BH_within_subclass_54) & results$BH_within_subclass_54 < 0.05
  results$significant_BH_global_22x54 <- !is.na(results$BH_global_22x54) & results$BH_global_22x54 < 0.05
  results$heatmap_significant_global_BH <- results$significant_BH_global_22x54
  annotation_index <- match(results$program, annotation$program)
  results$program_name_short <- annotation$program_name_short[annotation_index]
  results$program_name_full <- annotation$program_name_full[annotation_index]
  results$program_confidence <- annotation$program_confidence[annotation_index]
  results <- results[order(results$subclass, as.integer(sub("^P", "", results$program))), , drop = FALSE]
  rownames(results) <- NULL
  assert_results(results, expected_programs, subclasses)

  coverage <- do.call(rbind, coverage_rows)
  rownames(coverage) <- NULL
  subclass_summary <- do.call(rbind, lapply(subclasses, function(subclass_value) {
    x <- results[results$subclass == subclass_value, , drop = FALSE]
    y <- coverage[coverage$subclass == subclass_value, , drop = FALSE]
    data.frame(
      subclass = subclass_value,
      n_combinations_total = y$n_combinations_total,
      n_combinations_min_cells_10 = y$n_combinations_min_cells_10,
      n_combinations_excluded_min_cells_10 = y$n_combinations_excluded_min_cells_10,
      n_cells_total = y$n_cells_total,
      n_cells_min_cells_10 = y$n_cells_min_cells_10,
      n_donors = y$n_donors,
      n_regions = y$n_regions,
      n_age_values = y$n_age_values,
      n_sex_values = y$n_sex_values,
      fixed_design_rank = y$fixed_design_rank,
      fixed_design_columns = y$fixed_design_columns,
      fixed_design_full_rank = y$fixed_design_full_rank,
      n_region_coefficients = y$n_region_coefficients,
      consensus_correlation = y$consensus_correlation,
      analysis_status = y$analysis_status,
      n_programs_total = length(expected_programs),
      n_programs_estimable = sum(x$estimable),
      n_programs_failed = sum(!x$estimable),
      n_programs_significant_BH_within_subclass_54 = sum(x$significant_BH_within_subclass_54),
      n_programs_significant_BH_global_22x54 = sum(x$significant_BH_global_22x54),
      programs_significant_BH_within_subclass_54 = collapse_ids(x$program[x$significant_BH_within_subclass_54]),
      programs_significant_BH_global_22x54 = collapse_ids(x$program[x$significant_BH_global_22x54]),
      integrated_heatmap_primary = "BH_global_22x54",
      stringsAsFactors = FALSE
    )
  }))
  program_subclass_counts <- do.call(rbind, lapply(expected_programs, function(program_value) {
    x <- results[results$program == program_value, , drop = FALSE]
    annotation_row <- annotation[annotation$program == program_value, , drop = FALSE]
    data.frame(
      program = program_value,
      program_name_short = annotation_row$program_name_short,
      program_name_full = annotation_row$program_name_full,
      program_confidence = annotation_row$program_confidence,
      n_subclasses_total = length(subclasses),
      n_subclasses_estimable = sum(x$estimable),
      n_subclasses_failed = sum(!x$estimable),
      n_subclasses_significant_BH_within_subclass_54 = sum(x$significant_BH_within_subclass_54),
      n_subclasses_significant_BH_global_22x54 = sum(x$significant_BH_global_22x54),
      subclasses_significant_BH_within_subclass_54 = collapse_ids(x$subclass[x$significant_BH_within_subclass_54]),
      subclasses_significant_BH_global_22x54 = collapse_ids(x$subclass[x$significant_BH_global_22x54]),
      integrated_heatmap_primary = "BH_global_22x54",
      stringsAsFactors = FALSE
    )
  }))

  write_tsv(results, file.path(output_dir, "all_54_by_22_limma_dupcor.tsv"))
  write_tsv(program_subclass_counts, file.path(output_dir, "program_subclass_counts.tsv"))
  write_tsv(subclass_summary, file.path(output_dir, "subclass_summary.tsv"))

}

if (sys.nframe() == 0L) main()
