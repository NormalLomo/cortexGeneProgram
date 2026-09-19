if (!nzchar(Sys.getenv("CORTEX_PROGRAM_ROOT"))) stop("Set CORTEX_PROGRAM_ROOT; original source directories must remain read-only.")
options(future.globals.maxSize = 50 * 1024^3)
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else stop("Supply audit root")
out_dir <- file.path(root, "dream_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

remote_rlib <- file.path(root, "Rlib")
configured_rlibs <- strsplit(Sys.getenv("FORMAL21_R_LIBS"), .Platform$path.sep, fixed = TRUE)[[1L]]
configured_rlibs <- configured_rlibs[nzchar(configured_rlibs) & dir.exists(configured_rlibs)]
.libPaths(unique(c(remote_rlib, configured_rlibs, .libPaths())))

suppressPackageStartupMessages(library(variancePartition))

input_path <- file.path(root, "output", "donor_region_means.tsv")
metadata_path <- paste0(Sys.getenv("CORTEX_PROGRAM_ROOT"), "/inputs/snRNA_1M_obs.csv")
name_path <- file.path(root, "scripts", "program_names.tsv")
programs <- paste0("P", seq_len(54L))
formula_text <- "~ region + age + sex + (1 | donor)"
response_scale_factor <- 1000

means <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
required_metadata <- c("donor", "age", "sex")
if (!all(required_metadata %in% names(metadata))) stop("Source metadata lacks donor, age, or sex")

donor_covariates <- unique(metadata[, required_metadata])
if (anyDuplicated(donor_covariates$donor)) {
  consistency <- aggregate(cbind(age, sex) ~ donor, metadata[, required_metadata], function(x) length(unique(x)))
  if (any(consistency$age > 1L | consistency$sex > 1L)) stop("Age or sex varies within donor")
  donor_covariates <- donor_covariates[!duplicated(donor_covariates$donor), , drop = FALSE]
}

dat <- merge(means, donor_covariates, by = "donor", all.x = TRUE, sort = FALSE)
if (nrow(dat) != nrow(means)) stop("Metadata merge changed donor-region observation count")
dat$age <- suppressWarnings(as.numeric(dat$age))
dat$sex <- droplevels(factor(dat$sex))
dat$region <- droplevels(factor(dat$region))
dat$donor <- droplevels(factor(dat$donor))
if (anyNA(dat[, c("age", "sex", "region", "donor")])) stop("Missing covariate after merge")

fixed_design <- model.matrix(~ region + age + sex, data = dat)
if (qr(fixed_design)$rank != ncol(fixed_design)) stop("Fixed-effect design is rank deficient")

expr <- t(as.matrix(dat[, programs, drop = FALSE]))
storage.mode(expr) <- "double"
rownames(expr) <- programs
colnames(expr) <- paste0("obs", seq_len(nrow(dat)))
rownames(dat) <- colnames(expr)
original_variances <- apply(expr, 1L, stats::var)
expr <- expr * response_scale_factor

fit <- variancePartition::dream(
  exprObj = expr,
  formula = stats::as.formula(formula_text),
  data = dat,
  ddf = "Satterthwaite",
  useWeights = FALSE,
  REML = TRUE,
  hideErrorsInBackend = TRUE
)
fit_errors <- attr(fit, "errors")
returned_programs <- rownames(fit$coefficients)
missing_programs <- setdiff(programs, returned_programs)
if (length(missing_programs) > 0L) {
  missing_messages <- vapply(missing_programs, function(program) {
    message <- fit_errors[[program]]
    if (is.null(message)) "missing without a recorded backend error" else as.character(message)
  }, character(1L))
  stop(paste("dream did not return all 54 programs:", paste(paste(missing_programs, missing_messages, sep = "="), collapse = "; ")))
}
if (length(unique(returned_programs)) != 54L) stop("dream returned duplicate program rows")
fit <- variancePartition::eBayes(fit)
region_coefficients <- grep("^region", colnames(fit$coefficients), value = TRUE)
if (length(region_coefficients) < 1L) stop("No region coefficients available for joint test")
joint <- variancePartition::topTable(fit, coef = region_coefficients, number = Inf, sort.by = "none")
joint$program <- rownames(joint)
joint <- joint[match(programs, joint$program), , drop = FALSE]
joint$BH_54 <- p.adjust(joint$P.Value, method = "BH", n = 54L)
joint$significant_BH_54 <- !is.na(joint$BH_54) & joint$BH_54 < 0.05
joint$formula <- "score ~ region + age + sex + (1 | donor)"
joint$n_combinations <- nrow(dat)
joint$n_cells <- sum(dat$n_cells)
joint$n_donors <- nlevels(dat$donor)
joint$n_regions <- nlevels(dat$region)
joint$region_coefficients_tested <- paste(region_coefficients, collapse = ";")
joint$weights <- "none; continuous program-score matrix supplied directly to dream(useWeights=FALSE)"
joint$ddf <- "Satterthwaite"
joint$response_scale_factor <- response_scale_factor
joint$original_response_variance <- original_variances[joint$program]

joint$status <- ifelse(joint$program %in% names(fit_errors), "fit_error", "ok")
joint$error <- ifelse(joint$program %in% names(fit_errors), unname(fit_errors[joint$program]), "")
n_model_errors <- sum(joint$status != "ok", na.rm = TRUE)

if (file.exists(name_path)) {
  names <- read.delim(name_path, check.names = FALSE, stringsAsFactors = FALSE)
  names <- names[, intersect(c("new_P", "name_short", "name_full", "confidence"), names(names)), drop = FALSE]
  names(names)[names(names) == "new_P"] <- "program"
  joint <- merge(joint, names, by = "program", all.x = TRUE, sort = FALSE)
  joint <- joint[match(programs, joint$program), , drop = FALSE]
}

if (nrow(joint) != 54L) stop("Final result does not contain 54 rows")
if (anyNA(joint$program) || any(joint$program == "")) stop("Final result contains an empty program identifier")
if (!setequal(joint$program, programs) || anyDuplicated(joint$program)) stop("Final program identifiers are incomplete or duplicated")
key_columns <- c("F", "P.Value", "BH_54", "significant_BH_54")
if (anyNA(joint[, key_columns, drop = FALSE])) stop("Final result contains NA in a key statistical column")
if (n_model_errors != 0L) stop("At least one dream model retained an error")

write.table(joint, file.path(out_dir, "all_54_dream.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
