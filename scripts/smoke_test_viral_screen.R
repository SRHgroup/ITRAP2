#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(EnvStats)
})

step <- function(msg) message(sprintf("[smoke-test] %s", msg))

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

assert_col_has_non_na <- function(df, col) {
  assert_true(col %in% colnames(df), sprintf("Missing expected column: %s", col))
  assert_true(sum(!is.na(df[[col]])) > 0, sprintf("Column %s is all NA", col))
}

load_itrap2 <- function() {
  if (requireNamespace("ITRAP2", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ITRAP2))
    return("package")
  }

  step("ITRAP2 not installed; sourcing local R/*.R files")
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  assert_true(length(r_files) > 0, "No R source files found under ./R")
  for (f in sort(r_files)) {
    sys.source(f, envir = globalenv())
  }
  "source"
}

step("Loading ITRAP2")
mode <- load_itrap2()
step(sprintf("ITRAP2 load mode: %s", mode))

step("Loading data/viral_screen.RData")
assert_true(file.exists("data/viral_screen.RData"), "Missing data/viral_screen.RData")
load("data/viral_screen.RData")
assert_true(exists("viral_screen"), "Object 'viral_screen' not found after loading RData")

step("Running score_pmhc_noise")
viral_screen <- score_pmhc_noise(viral_screen, how = "per_clone", downsample_rate = 0.1, verbose = FALSE)

step("Running ScaleDataNoOutliers")
viral_screen <- ScaleDataNoOutliers(viral_screen, verbose = FALSE)

step("Subsetting BC == 'BC390'")
obj_sub <- subset(viral_screen, BC == "BC390")
assert_true(ncol(obj_sub) > 0, "Subset BC == 'BC390' returned zero cells")

step("Running smooth_pmhc")
obj_sub <- smooth_pmhc(
  object = obj_sub,
  assay = "pMHC",
  cap_upper_quantiles = TRUE,
  replacement_q = 0.85,
  replace_ones = TRUE
)

step("Running assign_pmhc")
obj_sub <- assign_pmhc(
  object = obj_sub,
  assay = "pMHC",
  slot = "scale.data",
  assignment = "rosner",
  rosner_alpha = 0.001,
  verbose = FALSE,
  rosner_pval_dist = "normal",
  remove_for_params = 2,
  params = "remove_topx",
  pseudobulk_fun = median,
  adjust_permutation = FALSE,
  adjust_test_within_clone = TRUE,
  padj_method = "BH",
  assign_small_clones = TRUE
)

step("Validating pMHC output columns")
meta <- obj_sub@meta.data
pmhc_cols <- grep("^pMHC", colnames(meta), value = TRUE)
assert_true(length(pmhc_cols) > 0, "No metadata columns starting with 'pMHC' were created")
assert_col_has_non_na(meta, "pMHC_classification")
assert_col_has_non_na(meta, "pMHC_pvalues")

step("Running filter_pmhc")
filtered_obj <- filter_pmhc(obj_sub, condition = "pvalues < 0.1")
assert_true(ncol(filtered_obj) > 0, "filter_pmhc returned object with zero cells")

step("Running extract_pairs")
pairs <- extract_pairs(
  filtered_obj,
  custom_columns = c("BC", "cdr3_beta", "cdr3_alpha")
)
assert_true(is.data.frame(pairs), "extract_pairs did not return a data.frame")

step("Smoke test passed")
message("\nOK: workflow executed without errors and expected pMHC outputs are present.")
message(sprintf("Rows in pairs table: %d", nrow(pairs)))
message(sprintf("pMHC columns detected: %s", paste(pmhc_cols, collapse = ", ")))
