#' Compute Average Expression for Pseudobulks Based on Metadata Grouping
#'
#' This function calculates the average expression for each feature across pseudobulk samples,
#' where pseudobulks are defined by groups in a specified metadata column of a Seurat object.
#' It preserves the original column names and ensures that the averages are computed exactly
#' according to the specified assay and slot without any unwanted alterations to column names.
#'
#' @param object A Seurat object containing the assay data to be analyzed.
#' @param assay The name of the assay to use for expression data (default: "RNA").
#' @param slot The slot containing the expression data ("data", "counts", or "scale.data").
#' @param clone_col The name of the metadata column used to define groups (default: "clone_id").
#'
#' @return A matrix containing the average expression values for each feature across the defined
#' pseudobulk groups. Columns are named after the unique identifiers found in the `clone_col` metadata column.
#'
#' @examples
#' # Assuming 'seurat_obj' is a Seurat object with a metadata column named 'clone_id'
#' pseudobulk_matrix <- ClonePseudobulk(seurat_obj, assay = "RNA", slot = "data", clone_col = "clone_id")
#' # This will return a matrix with average expressions for each 'clone_id' group.
#'
#' @export
#'
#' @importFrom Seurat GetAssayData
ClonePseudobulk <- function(object, assay="pMHC", slot="scale.data", 
                            clone_col="clone_id", z_score_threshold=2,
                            filter_by_zscore=FALSE, FUN=median) {
  
  if (filter_by_zscore){
    warning("filter_by_zscore should only be used if slot='scale.data'")
    cells <- Cells(object)
    cells_to_subset <- apply(GetAssayData(object = object, assay = assay, layer = slot), 
                             2, function(x) any(x > threshold))
    cells_subset <- cells[cells_to_subset]
    object <- subset(object, cells = cells_subset)
  }
  
  data <- GetAssayData(object = object, assay = assay, layer = slot)
  
  if (!clone_col %in% colnames(object@meta.data)) {
    stop(paste("Metadata column", clone_col, "not found in the Seurat object."))
  }
  
  group_ids <- setNames(object = object@meta.data[[clone_col]], nm = Cells(object)) 
  unique_groups <- unique(group_ids) %>% na.omit()  # Changed from drop.na() to na.omit() which is correct for omitting NA values.
  
  pseudobulk_means <- lapply(unique_groups, function(group) {
    cells_in_group <- names(group_ids[group_ids == group]) %>% drop.na()
    
    group_data <- data[, cells_in_group, drop = FALSE]  # Ensure it is always treated as a matrix.
    
    # Calculate row means, ensuring output is a vector if only one row.
    if (is.matrix(group_data)) {
      return(apply(group_data, 1, FUN, na.rm = TRUE))
      # Directly return the column as is if only one column present.
    } else {
      return(group_data[,1])
    }
  })
  
  # Convert list of vectors to a matrix
  pseudobulk_means <- do.call(cbind, pseudobulk_means)
  
  # Assign the group names to the columns of the result
  colnames(pseudobulk_means) <- unique_groups
  
  return(pseudobulk_means)
}

#' Drop NA Values from a Vector
#'
#' This function removes `NA` values from a given vector and returns the cleaned vector.
#'
#' @param vec A vector of any type (numeric, character, etc.) from which `NA` values are to be removed.
#'
#' @return A vector of the same type as the input `vec`, but with `NA` values removed.
#'
#' @examples
#' # Remove NA values from a numeric vector
#' drop.na(c(1, 2, NA, 4, 5, NA))
#' # Output: c(1, 2, 4, 5)
#'
#' # Remove NA values from a character vector
#' drop.na(c("a", "b", NA, "d"))
#' # Output: c("a", "b", "d")
#'
#' @export
drop.na <- function(vec, drop.neg=F){
  
  if (drop.neg){
    vec <- vec[vec!='Negative']
  }
  
  vec[!is.na(vec)]
}

#' Normalize a Numeric Vector
#'
#' This function normalizes a numeric vector to a range between 0 and 1.
#'
#' @param vec A numeric vector that you want to normalize.
#'
#' @return A numeric vector with values scaled to the range [0, 1].
#'
#' @examples
#' # Normalize a numeric vector
#' normalize_vector(c(1, 2, 3, 4, 5))
#' # Output: c(0, 0.25, 0.5, 0.75, 1)
#'
#' @export
normalize_vector <- function(vec) {
  if (!is.numeric(vec)) {
    stop("Input vector must be numeric")
  }
  
  min_val <- min(vec)
  max_val <- max(vec)
  
  normalized_vec <- (vec - min_val) / (max_val - min_val)
  
  return(normalized_vec)
}


#' Extract pMHC-TCR Pairs from Seurat Object
#'
#' This function extracts pMHC-TCR pairs from a Seurat object, ensuring each pair is written with each pair in a single row.
#'
#' @param object A Seurat object containing TCR and pMHC data.
#' @param custom_columns a character vector of columns from your `objects meta.data you want to include 
#' @param check_notation check and pring possible errors in storing pMHC related columns, with incinsistent : signs, separating multiple pMHC 
#' @return A data frame with each row representing a unique pMHC-TCR pair. The data frame includes columns for pMHC classification, confidence, p-values, clonotype details, and TCR details.
#'
#' @examples
#' # Assuming `seurat_obj` is a Seurat object with relevant data
#' pairs_df <- extract_pairs(seurat_obj)
#'
#' @export
extract_pairs <- function(object, include_negative = FALSE, custom_columns = NULL, check_notation = FALSE) {
  md <- object@meta.data
  existing_cols <- colnames(md)
  
  maybe_etm <- if ("epitope_type_multiple" %in% existing_cols) "epitope_type_multiple" else NULL
  
  default_columns <- c(
    "pMHC_classification", "pMHC_confidence", "pMHC_pvalues",
    "pMHC_wclone_pvalues", "pMHC_deltas", "pMHC_scaled_umis",
    "clone_id", "junction_beta", "junction_alpha", "clone_size",
    "v_call_beta", "c_call_beta", "j_call_beta", "d_call_beta",
    "v_call_alpha", "c_call_alpha", "j_call_alpha", "d_call_alpha"
  )
  
  all_columns <- c(default_columns, custom_columns, maybe_etm)
  
  split_targets <- intersect(
    c("pMHC_classification", "pMHC_confidence", "pMHC_pvalues",
      "pMHC_wclone_pvalues", "pMHC_deltas", "pMHC_scaled_umis",
      "epitope_type_multiple"),
    existing_cols
  )
  
  # ---------------------------
  # ⭐ PRE-CHECK FOR BAD COLS ⭐
  # ---------------------------
  if (check_notation){
    if (length(split_targets) > 1) {
      
      # count number of ":" items per column
      sep_counts <- md %>%
        mutate(row_id_original = row_number()) %>%
        mutate(across(all_of(split_targets), ~ ifelse(is.na(.x) | .x == "", 1L,
                                                      stringr::str_count(.x, ":") + 1L)))
      
      # for each row, check if all split-target columns have same count
      inconsistent <- sep_counts %>%
        rowwise() %>%
        mutate(is_inconsistent = length(unique(c_across(all_of(split_targets)))) > 1) %>%
        ungroup() %>%
        filter(is_inconsistent)
      
      if (nrow(inconsistent) > 0) {
        message("\n❌ ERROR: Inconsistent ':' counts detected before separate_rows().\n")
        
        purrr::walk(1:nrow(inconsistent), function(i) {
          row <- inconsistent[i, ]
          row_id <- row$row_id_original
          clone <- row$clone_id %||% "NA"
          
          cat("-----\n")
          cat("Original meta.data row:", row_id, "\n")
          cat("Clone ID:", clone, "\n")
          cat("Column item counts:\n")
          
          for (col in split_targets) {
            cat("  •", col, "→", row[[col]], "\n")
          }
        })
        
        stop("\nFix the inconsistent ':'-separated columns before splitting.\n")
      }
    }
  }
  
  md %>%
    { if (!include_negative) filter(., pMHC_classification != "Negative" & !is.na(pMHC_classification)) else . } %>%
    filter(productive_beta & productive_alpha) %>%
    select(all_of(intersect(all_columns, existing_cols))) %>%
    tidyr::separate_rows(all_of(split_targets), sep = ":") %>%
    mutate(
      HLA = gsub("_.+", "", pMHC_classification),
      peptide = gsub(".+_", "", pMHC_classification),
      tcr_pmhc = paste0(junction_beta, "_", junction_alpha, "_", peptide)
    ) %>%
    distinct()
}


#' Extract pMHC-TCR Pairs from Seurat Object
#'
#' This function extracts pMHC-TCR pairs from a Seurat object, ensuring each pair is written with each pair in a single row.
#'
#' @param object A Seurat object containing TCR and pMHC data.
#' @param custom_columns a character vector of columns from your `objects meta.data you want to include 
#' @return A data frame with each row representing a unique pMHC-TCR pair. The data frame includes columns for pMHC classification, confidence, p-values, clonotype details, and TCR details.
#'
#' @examples
#' # Assuming `seurat_obj` is a Seurat object with relevant data
#' pairs_df <- extract_pairs(seurat_obj)
#'
#' @export
extract_pairs_old <- function(object, include_negative = FALSE, custom_columns = NULL) {
  md <- object@meta.data
  existing_cols <- colnames(md)
  
  # include epitope_type_multiple only if present
  maybe_etm <- if ("epitope_type_multiple" %in% existing_cols) "epitope_type_multiple" else NULL
  
  default_columns <- c(
    "pMHC_classification", "pMHC_confidence", "pMHC_pvalues",
    "pMHC_wclone_pvalues", "pMHC_deltas", "pMHC_scaled_umis",
    "clone_id", "junction_beta", "junction_alpha", "clone_id", "clone_size",
    "v_call_beta", "c_call_beta", "j_call_beta", "d_call_beta",
    "v_call_alpha", "c_call_alpha", "j_call_alpha", "d_call_alpha"
  )
  
  all_columns <- c(default_columns, custom_columns, maybe_etm)
  
  # columns to split by ":" (only those that actually exist)
  split_targets <- intersect(
    c("pMHC_classification", "pMHC_confidence", "pMHC_pvalues",
      "pMHC_wclone_pvalues", "pMHC_deltas", "pMHC_scaled_umis",
      "epitope_type_multiple"),
    existing_cols
  )
  
  md %>%
    { if (!include_negative) dplyr::filter(., pMHC_classification != "Negative" & !is.na(pMHC_classification)) else . } %>%
    dplyr::filter(productive_beta & productive_alpha) %>%
    dplyr::select(dplyr::all_of(intersect(all_columns, existing_cols))) %>%
    tidyr::separate_rows(dplyr::all_of(split_targets), sep = ":") %>%
    dplyr::mutate(HLA = gsub("_.+", "", pMHC_classification)) %>%
    dplyr::mutate(peptide = gsub(".+_", "", pMHC_classification)) %>%
    dplyr::mutate(tcr_pmhc = paste0(junction_beta, "_", junction_alpha, "_", peptide)) %>%
    dplyr::distinct()
}

#' Compute per-clone pMHC scaled UMI values and delta-scores
#'
#' This function aggregates pMHC signals per T-cell clone in a Seurat object.
#' For every clone, it computes:
#'
#' * **scaled UMI expression**: mean per-feature values from the assay `@data`
#'   collapsed into a ":"-separated string in the order of pMHC tokens.
#'
#' * **delta scores**: feature-wise positive-part scaled expression
#'   from `@scale.data`, normalized and background-subtracted using the mean
#'   expression of all *non-target* pMHC features.
#'
#' The pMHC features belonging to each clone are extracted from a metadata
#' column such as `"pMHC_classification"`. If the metadata encodes multiple
#' pMHCs for a clone in a single string (e.g. `"A0201_GVLDAVWRV:A0201_KVDDTFYYV"`),
#' the function splits and maps each token to the true assay feature (Barcode)
#' using `object@misc$pmhc`.
#'
#' The resulting per-clone collapsed values are written back into
#' `object@meta.data` under new columns defined by `out_scaled_col`
#' and `out_delta_col`. Any previous versions of these columns are removed.
#'
#' @param object A Seurat object containing a pMHC assay. The metadata must
#'   include `clone_id` and a column with pMHC tokens (default:
#'   `"pMHC_classification"`). If present, `object@misc$pmhc` must contain
#'   a mapping between human-readable pMHC names and feature `Barcode`s.
#'
#' @param assay Character. Name of the assay storing pMHC data (default: `"pMHC"`).
#'
#' @param pmhc_class_col Character. Metadata column specifying pMHC tokens
#'   collapsed per clone (e.g. `"A0201_GVLDAWVRV:A0201_KVDDTFYYV"`).
#'
#' @param out_scaled_col Character. Name of the metadata column to store
#'   per-clone collapsed scaled-UMI values (default `"pMHC_scaled_umis"`).
#'
#' @param out_delta_col Character. Name of the metadata column to store
#'   per-clone delta-scores (default `"pMHC_deltas"`).
#'
#' @param clone_col Character. Metadata column that uniquely identifies a clone
#'   (default `"clone_id"`).
#'
#' @param verbose Logical. Whether to display a progress bar (default `TRUE`).
#'
#' @details
#'
#'
#' **Returned columns**
#' * `out_scaled_col`: e.g. `"0.124:0.522:0.991"`
#' * `out_delta_col`:  e.g. `"0.001:0.551:-0.110"`
#'
#' If a clone has no valid pMHC tokens or no mapped features,
#' `"Negative"` is returned for both fields.
#'
#' @return The input Seurat object with two updated metadata columns:
#'   * `out_scaled_col` (character): ":"-collapsed per-clone UMI means  
#'   * `out_delta_col`  (character): ":"-collapsed per-clone delta-scores
#'
#' @examples
#' \dontrun{
#'   seurat_obj <- calculate_pmhc_scaled_umis_and_deltas(seurat_obj)
#'   head(seurat_obj@meta.data[c("pMHC_scaled_umis", "pMHC_deltas")])
#' }
#'
#' @export
calculate_pmhc_scaled_umis_and_deltas <- function(
    object,
    assay = "pMHC",
    pmhc_class_col = "pMHC_classification",
    out_scaled_col = "pMHC_scaled_umis",
    out_delta_col  = "pMHC_deltas",
    clone_col = "clone_id",
    verbose = TRUE
){
  stopifnot(!is.null(object@meta.data[[clone_col]]))
  if (is.null(object@meta.data[[pmhc_class_col]]))
    stop(sprintf("meta.data missing '%s'", pmhc_class_col))
  
  # pull matrices once
  pm_data  <- tryCatch(GetAssayData(object, assay = assay, layer = "data"),
                       error = function(e) GetAssayData(object, assay = assay, slot = "data"))
  pm_scale <- tryCatch(GetAssayData(object, assay = assay, layer = "scale.data"),
                       error = function(e) GetAssayData(object, assay = assay, slot = "scale.data"))
  
  # map human-readable -> feature id (Barcode)
  pmhc_map <- NULL
  if (!is.null(object@misc$pmhc) && all(c("pmhc","Barcode") %in% colnames(object@misc$pmhc)))
    pmhc_map <- setNames(object@misc$pmhc$Barcode, object@misc$pmhc$pmhc)
  
  map_token <- function(tok){
    if (tok %in% rownames(pm_data))  return(tok)
    if (tok %in% rownames(pm_scale)) return(tok)
    if (!is.null(pmhc_map) && tok %in% names(pmhc_map)) {
      bc <- pmhc_map[[tok]]
      if (bc %in% rownames(pm_data) || bc %in% rownames(pm_scale)) return(bc)
    }
    NA_character_
  }
  
  make_row <- function(cl, scaled, delta){
    df <- data.frame(clone_id = cl, check.names = FALSE)
    df[[out_scaled_col]] <- scaled
    df[[out_delta_col]]  <- delta
    df
  }
  
  clones <- unique(stats::na.omit(object@meta.data[[clone_col]]))
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = length(clones), style = 3, width = 50, char = "+")
  
  res <- lapply(seq_along(clones), function(i){
    cl <- clones[i]
    if (verbose) utils::setTxtProgressBar(pb, i)
    
    cl_idx <- which(object@meta.data[[clone_col]] == cl)
    if (!length(cl_idx)) return(make_row(cl, "Negative", "Negative"))
    
    # tokens = your already-found outliers
    cls <- unique(stats::na.omit(object@meta.data[[pmhc_class_col]][cl_idx]))
    if (!length(cls) || cls[1] %in% c("Negative","")) return(make_row(cl, "Negative", "Negative"))
    
    toks <- strsplit(cls[1], ":", fixed = TRUE)[[1]]
    toks <- toks[nzchar(toks) & toks != "Negative"]
    feats <- vapply(toks, map_token, character(1))
    feats <- feats[!is.na(feats)]
    if (!length(feats)) return(make_row(cl, "Negative", "Negative"))
    
    cl_cells <- Cells(object)[cl_idx]
    m_data   <- Matrix::rowMeans(pm_data[,  cl_cells, drop = FALSE])
    m_scale  <- Matrix::rowMeans(pm_scale[, cl_cells, drop = FALSE])
    
    # ensure names exist (some combos drop them)
    if (is.null(names(m_data)))  names(m_data)  <- rownames(pm_data)
    if (is.null(names(m_scale))) names(m_scale) <- rownames(pm_scale)
    
    # order features to tokens order and keep only present
    feats <- feats[feats %in% names(m_data)]
    if (!length(feats)) return(make_row(cl, "Negative", "Negative"))
    
    # scaled_umis: raw means from @data
    scaled_str <- paste(round(m_data[feats], 3), collapse = ":")
    
    # deltas: clamp negatives, normalize, pos - mean(others)
    v <- m_scale
    v[!is.finite(v)] <- 0
    v <- pmax(v, 0)
    norm_all <- normalize_vector(v)
    names(norm_all) <- names(v)  # keep names after normalization
    
    others <- setdiff(names(norm_all), feats)
    bg_mean <- if (length(others)) mean(norm_all[others], na.rm = TRUE) else 0
    
    deltas <- norm_all[feats] - bg_mean
    deltas[!is.finite(deltas)] <- 0
    delta_str <- paste(round(deltas, 3), collapse = ":")
    
    make_row(cl, scaled_str, delta_str)
  })
  
  if (verbose) close(pb)
  
  df <- do.call(rbind, res)
  # remove any old versions of these columns, then join fresh
  pmhc_cols_re <- sprintf("^(%s|%s)$", out_scaled_col, out_delta_col)
  
  object@meta.data <- object@meta.data %>%
    tibble::rownames_to_column("row_id") %>%
    dplyr::select(-tidyselect::matches(pmhc_cols_re)) %>%
    dplyr::left_join(df, by = setNames("clone_id", clone_col)) %>%
    tibble::column_to_rownames("row_id")
  
  object
}


#' Filter pMHC–TCR pairs by donor HLA matching and collapse matched pMHC features per clone
#'
#' This function takes a Seurat object containing TCR–pMHC assignments and donor
#' HLA genotypes (stored in `object@misc$hla_genotype`). For each cell, it checks
#' whether the presented pMHC allele matches one of the donor's typed HLA alleles.
#' Only HLA-matched pMHC–TCR pairs are retained, and pMHC-related metadata
#' (`pMHC_classification`, confidence scores, deltas, UMIs, etc.) are collapsed
#' per clone into a single ":"-separated string.
#'
#' Optionally, the function preserves the unfiltered pMHC metadata columns before
#' filtering and collapsing, allowing downstream comparison between filtered and
#' unfiltered data.
#'
#' @param object Seurat object containing pMHC–TCR metadata in `@meta.data` and
#'   HLA genotypes in `@misc$hla_genotype`. The metadata must contain
#'   `clone_id`, `HLA`, and pMHC fields.
#' @param keep_unfiltered Logical; if `TRUE` (default), the function stores the
#'   original pMHC metadata columns in new columns ending with `_unf` before
#'   filtering.
#' @param donor_col Name of the metadata column identifying the donor for each
#'   cell. Must match the first column name of the HLA genotype table stored in
#'   `object@misc$hla_genotype`.
#'
#' @return A Seurat object with updated `@meta.data`, where:
#'   * HLA-matching pMHC–TCR pairs are retained,
#'   * pMHC metadata are collapsed per `clone_id`,
#'   * optional unfiltered metadata columns are added with `_unf` suffix.
#'
#' @examples
#' \dontrun{
#'   seurat_filtered <- filter_mismatch_responses(seurat_obj, donor_col = "donor")
#'   head(seurat_filtered@meta.data)
#' }
#'
#' @export
filter_mismatch_responses <- function(object, keep_unfiltered=T, donor_col='donor'){
  if (is.null(object@misc$hla_genotype)) {
    stop('no hla genotyping provided in object@misc$gla_genotype')
  } else {
    hlas <- object@misc$hla_genotype
  }
  
  if (colnames(hlas)[1] != donor_col){
    stop('you have to name 1st column in the hla_genotype, the same way as your donor column in the object, and donor_col param in this function')
  }
  
  if (keep_unfiltered){
    object@meta.data <- object@meta.data %>%
      mutate(
        pMHC_classification_unf = pMHC_classification,
        pMHC_confidence_unf     = pMHC_confidence,
        pMHC_pvalues_unf        = pMHC_pvalues,
        pMHC_wclone_pvalues_unf = pMHC_wclone_pvalues,  
        pMHC_deltas_unf         = pMHC_deltas,
        pMHC_scaled_umis_unf    = pMHC_scaled_umis
      )
  }
  
  object_pairs <- extract_pairs(object, custom_columns = donor_col)
  
  object_pairs <- object_pairs %>% 
    left_join(hlas, by = donor_col) %>%
    rowwise() %>%
    mutate(
      hla_match = case_when(
        all(is.na(c(`HLA-A1`, `HLA-A2`,
                    `HLA-B1`, `HLA-B2`,
                    `HLA-C1`, `HLA-C2`))) ~ NA,
        
        TRUE ~ HLA %in% c(`HLA-A1`, `HLA-A2`,
                          `HLA-B1`, `HLA-B2`,
                          `HLA-C1`, `HLA-C2`)
      )
    ) %>%
    ungroup()
  
  pmhc_cols <- c("pMHC_classification", "pMHC_classification", "pMHC_confidence",
                 "pMHC_pvalues", "pMHC_wclone_pvalues", "pMHC_deltas", "pMHC_scaled_umis")
  
  object_pairs <- object_pairs %>% 
    filter(hla_match) 
  
  object_pairs_collapsed <- object_pairs %>%
    group_by(clone_id) %>%
    summarise(
      across(
        all_of(pmhc_cols),
        ~ paste(unique(trimws(as.character(.x))), collapse = ":")
      ),
      .groups = "drop"
    )
  
  object_meta <- object@meta.data %>% 
    dplyr::select(-all_of(pmhc_cols)) %>%
    rownames_to_column('cell_id') %>%
    left_join(object_pairs_collapsed) %>%
    column_to_rownames('cell_id')
  
  object@meta.data <- object_meta[Cells(object),]
  
  return(object)
}




mask_pmhc_to_screen <- function(
    object,
    screened_hla,
    assay = "pMHC",
    donor_col = "donor",
    barcode_col = "Barcode",
    dict_hla_col = "HLA"
) {
  stopifnot(assay %in% names(object@assays))
  if (is.null(object@misc$pmhc)) {
    stop("object@misc$pmhc is missing; expected multimer annotation table there.")
  }
  
  # Pull counts (sparse) and convert to dense so we can write NAs
  X <- SeuratObject::GetAssayData(object, assay = assay, slot = "counts")
  barcodes <- rownames(X)
  cells    <- colnames(X)
  
  # ---- Get multimer map from misc
  multimer_map <- object@misc$pmhc
  mm <- multimer_map %>%
    dplyr::select(!!barcode_col, !!dict_hla_col) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      !!barcode_col := as.character(!!rlang::sym(barcode_col)),
      !!dict_hla_col := as.character(!!rlang::sym(dict_hla_col))
    )
  
  # Align rows to assay barcodes
  mm_rows <- dplyr::right_join(
    mm,
    tibble::tibble(!!barcode_col := barcodes),
    by = barcode_col
  ) %>%
    dplyr::arrange(factor(!!rlang::sym(barcode_col), levels = barcodes))
  
  if (!all(mm_rows[[barcode_col]] == barcodes)) {
    stop("Row alignment failure between multimer_map and assay rows.")
  }
  hla_by_row <- mm_rows[[dict_hla_col]]
  
  # ---- Build donor per column
  md <- object@meta.data
  if (!donor_col %in% colnames(md)) stop("meta.data has no column: ", donor_col)
  donors_by_cell <- md[[donor_col]][match(cells, rownames(md))]
  
  # ---- donor → allowed HLA
  keyed <- screened_hla %>%
    dplyr::mutate(
      donor_key = as.character(.data$donor_id),
      screened_hla = as.character(.data$screened_hla)
    ) %>%
    dplyr::select(donor_key, screened_hla) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(donor_key) & donor_key != "")
  
  allowed_map <- stats::setNames(
    strsplit(keyed$screened_hla, "\\s*,\\s*"),
    keyed$donor_key
  )
  
  # ---- Mask disallowed rows per donor
  M <- as.matrix(X)
  rows_by_hla <- split(seq_along(hla_by_row), hla_by_row)
  
  for (d in unique(donors_by_cell)) {
    cols_d <- which(donors_by_cell == d)
    if (!length(cols_d)) next
    if (!(d %in% names(allowed_map))) next
    
    allowed <- allowed_map[[d]]
    allowed_rows <- unlist(rows_by_hla[intersect(names(rows_by_hla), allowed)], use.names = FALSE)
    disallowed_rows <- setdiff(seq_along(hla_by_row), allowed_rows)
    if (length(disallowed_rows) == 0) next
    
    M[disallowed_rows, cols_d] <- NA_real_
  }
  
  # ---- add as new assay
  M <-  Matrix(M, sparse = TRUE)
  new_assay <- SeuratObject::CreateAssayObject(counts = M)
  object[["pMHC_masked"]] <- new_assay
  attr(object, "pMHC_masked_matrix") <- M
  return(object)
}

