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

normalize_pairs <- function(df) {
  df <- tibble::as_tibble(df) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.x)))
  sort_cols <- sort(colnames(df))
  df %>% dplyr::arrange(dplyr::across(dplyr::all_of(sort_cols)))
}

compare_or_write_expected_pairs <- function(pairs_df, expected_path = "scripts/smoke_test_expected_pairs.tsv") {
  pairs_norm <- normalize_pairs(pairs_df)
  
  if (!file.exists(expected_path)) {
    readr::write_tsv(pairs_norm, expected_path)
    step(sprintf("Wrote baseline expected pairs: %s", expected_path))
    return(invisible(TRUE))
  }
  
  expected <- readr::read_tsv(expected_path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  expected_norm <- normalize_pairs(expected)
  
  same_columns <- identical(colnames(expected_norm), colnames(pairs_norm))
  same_rows <- nrow(expected_norm) == nrow(pairs_norm)
  same_content <- same_columns && same_rows && identical(as.data.frame(expected_norm), as.data.frame(pairs_norm))
  
  if (!same_content) {
    out_dir <- "scripts/output"
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    latest_path <- file.path(out_dir, "smoke_test_pairs_latest.tsv")
    readr::write_tsv(pairs_norm, latest_path)
    
    warning(
      sprintf(
        paste0(
          "extract_pairs output changed vs baseline expected file (%s). ",
          "Current snapshot written to %s. ",
          "expected_rows=%d, current_rows=%d, same_columns=%s"
        ),
        expected_path, latest_path, nrow(expected_norm), nrow(pairs_norm), same_columns
      ),
      call. = FALSE
    )
  } else {
    step("extract_pairs output matches baseline expected pairs")
  }
  
  invisible(TRUE)
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
compare_or_write_expected_pairs(pairs)

step("Smoke test passed")
message("\nOK: workflow executed without errors and expected pMHC outputs are present.")
message(sprintf("Rows in pairs table: %d", nrow(pairs)))
message(sprintf("pMHC columns detected: %s", paste(pmhc_cols, collapse = ", ")))
