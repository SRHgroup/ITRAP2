#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

step <- function(msg) message(sprintf("[stress-test] %s", msg))

load_itrap2 <- function() {
  if (requireNamespace("ITRAP2", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ITRAP2))
    return("package")
  }
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  for (f in sort(r_files)) sys.source(f, envir = globalenv())
  "source"
}

run_with_capture <- function(expr) {
  warns <- character()
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  t1 <- proc.time()[["elapsed"]]

  if (inherits(out, "error")) {
    list(ok = FALSE, value = NULL, error = conditionMessage(out), warnings = warns, elapsed = t1 - t0)
  } else {
    list(ok = TRUE, value = out, error = NA_character_, warnings = warns, elapsed = t1 - t0)
  }
}

validate_assignment <- function(obj) {
  md <- obj@meta.data
  req <- c("pMHC_classification", "pMHC_pvalues")
  missing <- setdiff(req, colnames(md))
  if (length(missing) > 0) {
    return(list(ok = FALSE, detail = paste("missing columns:", paste(missing, collapse = ", "))))
  }

  non_na_class <- sum(!is.na(md$pMHC_classification))
  non_na_pvals <- sum(!is.na(md$pMHC_pvalues))

  if (non_na_class == 0 || non_na_pvals == 0) {
    return(list(ok = FALSE, detail = "pMHC output columns are all NA"))
  }

  list(ok = TRUE, detail = sprintf("non-NA class=%d; non-NA pvals=%d", non_na_class, non_na_pvals))
}

step("Loading ITRAP2")
mode <- load_itrap2()
step(sprintf("Load mode: %s", mode))

step("Loading data")
load("data/viral_screen.RData")

step("Preparing baseline object")
viral_screen <- score_pmhc_noise(viral_screen, how = "per_clone", downsample_rate = 0.1, verbose = FALSE)
viral_screen <- ScaleDataNoOutliers(viral_screen, verbose = FALSE)
obj_sub <- subset(viral_screen, BC == "BC390")

if (ncol(obj_sub) == 0) {
  stop("subset BC == 'BC390' returned zero cells")
}

smooth_grid <- list(
  list(id = "S1_baseline", params = list(assay = "pMHC", cap_upper_quantiles = TRUE, replacement_q = 0.85, replace_ones = TRUE, cl_size_thresh = 3, span_val = 1, degree_val = 1, family_val = "symmetric", verbose = FALSE)),
  list(id = "S2_gaussian_nocap", params = list(assay = "pMHC", cap_upper_quantiles = FALSE, replacement_q = 0.85, replace_ones = FALSE, cl_size_thresh = 3, span_val = 0.75, degree_val = 1, family_val = "gaussian", verbose = FALSE)),
  list(id = "S3_strict_clone", params = list(assay = "pMHC", cap_upper_quantiles = TRUE, replacement_q = 0.95, replace_ones = TRUE, cl_size_thresh = 5, span_val = 1, degree_val = 1, family_val = "symmetric", verbose = FALSE))
)

assign_grid <- list(
  list(id = "A1_rosner_default", params = list(assay = "pMHC", slot = "scale.data", assignment = "rosner", rosner_alpha = 0.001, rosner_pval_dist = "normal", remove_for_params = 2, params = "remove_topx", pseudobulk_fun = median, adjust_permutation = FALSE, adjust_test_within_clone = TRUE, padj_method = "BH", assign_small_clones = TRUE, verbose = FALSE)),
  list(id = "A2_rosner_adjust_both", params = list(assay = "pMHC", slot = "scale.data", assignment = "rosner", rosner_alpha = 0.001, rosner_pval_dist = "normal", remove_for_params = 2, params = "remove_topx", pseudobulk_fun = median, adjust_permutation = TRUE, adjust_test_within_clone = TRUE, padj_method = "BH", assign_small_clones = TRUE, verbose = FALSE)),
  list(id = "A3_extreme", params = list(assay = "pMHC", slot = "scale.data", assignment = "extreme_distribution", extreme_alpha = 0.001, extreme_type = "regular", pseudobulk_fun = median, adjust_permutation = FALSE, adjust_test_within_clone = FALSE, assign_small_clones = TRUE, verbose = FALSE))
)

results <- list()

for (s in smooth_grid) {
  step(sprintf("Running smooth case: %s", s$id))
  smooth_res <- run_with_capture(do.call(smooth_pmhc, c(list(object = obj_sub), s$params)))

  if (!smooth_res$ok) {
    results[[length(results) + 1]] <- tibble(
      smooth_case = s$id,
      assign_case = NA_character_,
      status = "smooth_error",
      elapsed_sec = smooth_res$elapsed,
      n_warnings = length(smooth_res$warnings),
      error = smooth_res$error,
      validation = NA_character_,
      n_pairs = NA_integer_
    )
    next
  }

  smoothed <- smooth_res$value

  for (a in assign_grid) {
    case_id <- paste(s$id, a$id, sep = "__")
    step(sprintf("Running assign case: %s", case_id))

    assign_res <- run_with_capture(do.call(assign_pmhc, c(list(object = smoothed), a$params)))

    if (!assign_res$ok) {
      results[[length(results) + 1]] <- tibble(
        smooth_case = s$id,
        assign_case = a$id,
        status = "assign_error",
        elapsed_sec = assign_res$elapsed,
        n_warnings = length(smooth_res$warnings) + length(assign_res$warnings),
        error = assign_res$error,
        validation = NA_character_,
        n_pairs = NA_integer_
      )
      next
    }

    assigned <- assign_res$value
    v <- validate_assignment(assigned)

    fp_res <- run_with_capture({
      filtered_obj <- filter_pmhc(assigned, condition = "pvalues < 0.1")
      pairs <- extract_pairs(filtered_obj, custom_columns = c("BC", "cdr3_beta", "cdr3_alpha"))
      nrow(pairs)
    })

    status <- if (v$ok && fp_res$ok) "ok" else "validation_error"
    err <- NA_character_
    if (!v$ok) err <- v$detail
    if (!fp_res$ok) err <- paste(c(err, fp_res$error), collapse = " | ")

    results[[length(results) + 1]] <- tibble(
      smooth_case = s$id,
      assign_case = a$id,
      status = status,
      elapsed_sec = smooth_res$elapsed + assign_res$elapsed + fp_res$elapsed,
      n_warnings = length(smooth_res$warnings) + length(assign_res$warnings) + length(fp_res$warnings),
      error = err,
      validation = v$detail,
      n_pairs = if (fp_res$ok) as.integer(fp_res$value) else NA_integer_
    )
  }
}

report <- bind_rows(results) %>%
  arrange(factor(status, levels = c("assign_error", "smooth_error", "validation_error", "ok")), desc(n_warnings), elapsed_sec)

out_dir <- "scripts/output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "smooth_assign_stress_report.tsv")
readr::write_tsv(report, out_file)

cat("\n=== Stress Test Summary ===\n")
summary_lines <- apply(as.data.frame(report), 1, function(x) paste(x, collapse = " | "))
writeLines(summary_lines)
cat("\nReport saved to:", out_file, "\n")

if (any(report$status != "ok")) {
  q(status = 1)
}
