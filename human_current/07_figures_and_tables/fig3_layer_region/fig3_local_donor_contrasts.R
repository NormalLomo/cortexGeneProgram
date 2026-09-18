#!/usr/local/bin/Rscript
options(future.globals.maxSize = 50 * 1024^3)
options(stringsAsFactors = FALSE, digits = 17)
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(pbkrtest))

ROOT <- Sys.getenv("NMF_WORK_ROOT")
OLD <- Sys.getenv("NMF_ARCHIVE_SOURCE_ROOT")
if (!nzchar(ROOT) || !nzchar(OLD)) {
  stop("Set NMF_WORK_ROOT and NMF_ARCHIVE_SOURCE_ROOT; source directories remain read-only.")
}
OUT <- file.path(ROOT, "analysis/fig3_layer_region/fig3_local_donor_contrasts.tsv")
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
ann <- read_tsv(file.path(ROOT, "tables/TableS3_program_annotation.tsv"))
ann <- ann[order(as.integer(sub("^P", "", ann$new_P))), ]
ids <- unique(read_tsv(file.path(ROOT, "inputs/cortex_nmf_program/archived/revisions/2026-08-02_v49_professor_review/human_validation/HUMAN_VALIDATION_spatial_section_by_layer_aggregates_all54.tsv"))[, c("chip", "region", "donor")])
chip <- read_tsv(file.path(OLD, "scripts/fig2/prog_x_layer_per_chip.tsv"))
profiles <- read_tsv(file.path(ROOT, "analysis/fig3_layer_region/layer_region_profiles.tsv"))
donors <- c("DonorB", "DonorF", "DonorG", "DonorH", "DonorI")
layers <- paste0("L", 1:6)
domains <- c(layers, "WM")
areas <- c("ACC", "AG", "DLPFC", "FPPFC", "ITG", "M1", "PoCG", "S1",
           "S1E", "SMG", "SPL", "STG", "V1", "VLPFC")
infer_areas <- c("AG", "DLPFC", "FPPFC", "M1", "S1", "SMG", "SPL", "V1", "VLPFC")
chip <- merge(chip, ids, by = "chip", all.x = TRUE, sort = FALSE)
finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
clean_text <- function(x) gsub("[\t\r\n]", " ", paste(x, collapse = " | "))
observed <- tolower(as.character(profiles$observed)) %in% c("true", "1")
donor_profiles <- profiles[profiles$record_type == "donor" & observed, ]
records <- list()
append_record <- function(x) records[[length(records) + 1L]] <<- x

make_record <- function(panel, a, region = "", domain = "", eligible = TRUE) {
  data.frame(panel = panel, program = a$new_P, raw_component = a$cnmf_component,
    functional_name = a$functional_name, region = region, domain = domain,
    comparison = "", in_inference_family = eligible,
    family_size = if (eligible) if (panel == "b") 378L else 2916L else NA_integer_,
    method = if (panel == "b") "two-sided one-sample donor-contrast t"
             else "Gaussian LMM REML; Kenward-Roger two-sided linear contrast",
    effect = NA_real_, standard_error = NA_real_, df = NA_real_, statistic = NA_real_,
    ci95_lower = NA_real_, ci95_upper = NA_real_, p = NA_real_, q_BY = NA_real_,
    direction = "", significant_BY_0_05 = FALSE, status = "",
    n_donors = 0L, n_donors_model = 0L, n_donor_region_units = 0L,
    n_sections = 0L, donor_ids = "", regions_by_donor = "", donor_contrasts = "",
    model_singular = NA, detail = "", stringsAsFactors = FALSE)
}

# b: each biological donor contributes one mean of its observable within-area
# domain contrasts. No region, domain or donor values are imputed.
for (ai in seq_len(nrow(ann))) {
  a <- ann[ai, ]
  raw <- paste0("program_", a$cnmf_component)
  s <- chip[chip$program == raw & chip$majorDomain %in% domains, ]
  ag <- aggregate(s$mean_z, s[, c("donor", "region", "majorDomain")], finite_mean)
  names(ag)[4] <- "score"
  for (target in domains) {
    row <- make_record("b", a, domain = target)
    others <- setdiff(domains, target)
    row$comparison <- paste0(target, " minus equal mean of ", paste(others, collapse = ","),
      "; equal areas within each observed donor; equal five donors")
    dvalues <- setNames(rep(NA_real_, length(donors)), donors)
    coverage <- setNames(vector("list", length(donors)), donors)
    for (person in donors) {
      d <- ag[ag$donor == person, ]
      area_values <- c()
      kept <- character()
      for (area in sort(unique(d$region))) {
        z <- d[d$region == area, ]
        v <- setNames(z$score, z$majorDomain)[domains]
        if (length(v) == 7L && all(is.finite(v))) {
          area_values <- c(area_values, v[target] - mean(v[others]))
          kept <- c(kept, area)
        }
      }
      if (length(area_values)) dvalues[person] <- mean(area_values)
      coverage[[person]] <- kept
    }
    available <- is.finite(dvalues)
    row$n_donors <- sum(available)
    row$n_donors_model <- sum(available)
    row$donor_ids <- paste(names(dvalues)[available], collapse = ";")
    row$regions_by_donor <- paste(vapply(donors, function(d)
      paste0(d, "=", paste(coverage[[d]], collapse = ",")), character(1)), collapse = ";")
    row$donor_contrasts <- paste(vapply(donors, function(d)
      paste0(d, "=", if (is.finite(dvalues[d])) format(dvalues[d], digits = 17) else "NA"),
      character(1)), collapse = ";")
    row$n_donor_region_units <- sum(lengths(coverage))
    kept_keys <- unlist(lapply(donors, function(d) paste(d, coverage[[d]], sep = "|")),
                        use.names = FALSE)
    row$n_sections <- length(unique(s$chip[paste(s$donor, s$region, sep = "|") %in% kept_keys]))
    if (sum(available) != length(donors)) {
      row$status <- "NA_required_five_donor_contrasts_unavailable"
      row$detail <- "A complete seven-domain within-area contrast is required; no donor imputation."
    } else {
      row$effect <- mean(dvalues)
      ans <- tryCatch(stats::t.test(unname(dvalues), mu = 0, alternative = "two.sided",
                                   conf.level = 0.95), error = identity)
      if (inherits(ans, "error")) {
        row$status <- "NA_t_contrast_not_estimable"
        row$detail <- clean_text(conditionMessage(ans))
      } else if (!all(is.finite(c(ans$stderr, ans$parameter, ans$statistic, ans$p.value)))
                 || ans$stderr <= 0 || ans$parameter <= 0) {
        row$status <- "NA_invalid_standard_error_or_df"
      } else {
        row$standard_error <- unname(ans$stderr)
        row$df <- unname(ans$parameter)
        row$statistic <- unname(ans$statistic)
        row$ci95_lower <- ans$conf.int[1]
        row$ci95_upper <- ans$conf.int[2]
        row$p <- ans$p.value
        row$status <- "ok"
      }
    }
    append_record(row)
  }
}

# c: a separate same-layer Gaussian random-intercept model per program.
# Single-donor regions remain in the complete display grid, outside inference.
for (ai in seq_len(nrow(ann))) {
  a <- ann[ai, ]
  for (layer in layers) {
    all_d <- donor_profiles[donor_profiles$program == a$new_P &
                             donor_profiles$layer == layer, ]
    dat <- all_d[all_d$region %in% infer_areas & is.finite(all_d$score), ]
    dat$region <- factor(dat$region, levels = infer_areas)
    dat$donor <- factor(dat$donor, levels = donors)
    warnings <- character()
    fit <- withCallingHandlers(
      tryCatch(lme4::lmer(score ~ region + (1 | donor), data = dat, REML = TRUE),
               error = identity),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
      message = function(m) { warnings <<- c(warnings, conditionMessage(m)); invokeRestart("muffleMessage") })
    fit_status <- "ok"
    singular <- NA
    detail <- clean_text(warnings)
    if (inherits(fit, "error")) {
      fit_status <- "NA_model_fit_error"
      detail <- clean_text(c(detail, conditionMessage(fit)))
    } else {
      singular <- lme4::isSingular(fit)
      conv <- fit@optinfo$conv
      conv_messages <- conv$lme4$messages
      substantive <- conv_messages[!grepl("boundary.*singular", conv_messages, ignore.case = TRUE)]
      if ((length(conv$opt) && any(conv$opt != 0)) || length(substantive)) {
        fit_status <- "NA_model_nonconvergence"
        detail <- clean_text(c(detail, substantive))
      } else if (singular) fit_status <- "ok_boundary"
    }
    grid <- data.frame(region = factor(infer_areas, levels = infer_areas))
    design <- model.matrix(~ region, data = grid)
    for (area in areas) {
      eligible <- area %in% infer_areas
      row <- make_record("c", a, region = area, domain = layer, eligible = eligible)
      area_dat <- all_d[all_d$region == area & is.finite(all_d$score), ]
      row$n_donors <- length(unique(area_dat$donor))
      row$donor_ids <- paste(sort(unique(area_dat$donor)), collapse = ";")
      row$n_sections <- sum(area_dat$n_sections)
      if (!eligible) {
        row$status <- "NA_single_donor_region_outside_inference_scope"
        row$comparison <- "Display only; not a member of the nine-region contrast family"
        row$detail <- "No population-level regional inference; not a nonsignificant result."
      } else {
        row$n_donors_model <- length(unique(dat$donor))
        row$n_donor_region_units <- nrow(dat)
        row$regions_by_donor <- paste(vapply(donors, function(person)
          paste0(person, "=", paste(as.character(dat$region[dat$donor == person]), collapse = ",")),
          character(1)), collapse = ";")
        row$comparison <- paste0(area, " minus equal adjusted mean of ",
          paste(setdiff(infer_areas, area), collapse = ","), "; same layer ", layer)
        row$model_singular <- singular
        row$status <- fit_status
        row$detail <- detail
        if (row$n_donors < 2L) {
          row$status <- "NA_required_region_donor_coverage_unavailable"
        } else if (startsWith(fit_status, "ok")) {
          fixed_names <- names(lme4::fixef(fit))
          if (!setequal(colnames(design), fixed_names)) {
            row$status <- "NA_nonestimable_fixed_effect_contrast"
          } else {
            w <- rep(-1/8, length(infer_areas))
            w[match(area, infer_areas)] <- 1
            L <- as.numeric(w %*% design[, fixed_names, drop = FALSE])
            cwarnings <- character()
            ans <- withCallingHandlers(
              tryCatch({
                # Direct installed pbkrtest API used by lmerTest's KR1D path.
                # No automatic Satterthwaite fallback is permitted.
                adjusted <- pbkrtest::vcovAdj(fit)
                variance <- as.numeric(t(L) %*% adjusted %*% L)
                ddf <- pbkrtest::Lb_ddf(L = L, V0 = vcov(fit), Vadj = adjusted)
                estimate <- sum(L * lme4::fixef(fit))
                se <- sqrt(variance)
                tvalue <- estimate / se
                data.frame("Estimate" = estimate, "Std. Error" = se, "df" = ddf,
                  "t value" = tvalue,
                  "Pr(>|t|)" = 2 * stats::pt(abs(tvalue), df = ddf, lower.tail = FALSE),
                  check.names = FALSE)
              }, error = identity),
              warning = function(w) { cwarnings <<- c(cwarnings, conditionMessage(w)); invokeRestart("muffleWarning") },
              message = function(m) { cwarnings <<- c(cwarnings, conditionMessage(m)); invokeRestart("muffleMessage") })
            row$detail <- clean_text(c(detail, cwarnings))
            if (inherits(ans, "error")) {
              row$status <- "NA_Kenward_Roger_contrast_failed"
              row$detail <- clean_text(c(row$detail, conditionMessage(ans)))
            } else {
              vals <- unlist(ans[1, c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")],
                             use.names = FALSE)
              if (length(vals) != 5L || !all(is.finite(vals)) || vals[2] <= 0 ||
                  vals[3] <= 0 || vals[5] < 0 || vals[5] > 1) {
                row$status <- "NA_invalid_Kenward_Roger_standard_error_or_df"
              } else {
                row$effect <- vals[1]
                row$standard_error <- vals[2]
                row$df <- vals[3]
                row$statistic <- vals[4]
                row$p <- vals[5]
                radius <- stats::qt(0.975, vals[3]) * vals[2]
                row$ci95_lower <- vals[1] - radius
                row$ci95_upper <- vals[1] + radius
              }
            }
          }
        }
      }
      append_record(row)
    }
  }
}

result <- do.call(rbind, records)
for (panel in c("b", "c")) {
  family <- if (panel == "b") 378L else 2916L
  ix <- which(result$panel == panel & result$in_inference_family)
  result$q_BY[ix] <- p.adjust(result$p[ix], method = "BY", n = family)
}
result$direction <- ifelse(is.na(result$effect), "", ifelse(result$effect > 0, "positive",
                              ifelse(result$effect < 0, "negative", "zero")))
result$significant_BY_0_05 <- is.finite(result$q_BY) & result$q_BY < 0.05
write.table(result, OUT, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA",
            fileEncoding = "UTF-8")
cat("RETAINED", OUT, "rows", nrow(result), "\n")
for (panel in c("b", "c")) {
  s <- result[result$panel == panel, ]
  cat("PANEL", panel, "rows", nrow(s), "family", sum(s$in_inference_family),
      "valid", sum(is.finite(s$q_BY)), "NA", sum(!is.finite(s$q_BY)),
      "significant", sum(s$significant_BY_0_05),
      "positive", sum(s$significant_BY_0_05 & s$effect > 0, na.rm = TRUE),
      "negative", sum(s$significant_BY_0_05 & s$effect < 0, na.rm = TRUE), "\n")
  print(table(s$status, useNA = "ifany"))
}
