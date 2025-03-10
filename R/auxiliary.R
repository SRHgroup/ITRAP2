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
drop.na <- function(vec){
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
#' @return A data frame with each row representing a unique pMHC-TCR pair. The data frame includes columns for pMHC classification, confidence, p-values, clonotype details, and TCR details.
#'
#' @examples
#' # Assuming `seurat_obj` is a Seurat object with relevant data
#' pairs_df <- extract_pairs(seurat_obj)
#'
#' @export
extract_pairs <- function(object, custom_columns = NULL){
  
  default_columns <- c("pMHC_classification", "pMHC_confidence", "pMHC_pvalues", "pMHC_wclone_pvalues",
                       "clone_id", "junction_beta", "junction_alpha", "clone_id", 'clone_size',
                       "v_call_beta", "c_call_beta", "j_call_beta", "d_call_beta",
                       "v_call_alpha", "c_call_alpha", "j_call_alpha", "d_call_alpha")
  
  all_columns <- c(default_columns, custom_columns)
  
  return(
    object@meta.data %>%
      filter(pMHC_classification != 'Negative' & !is.na(pMHC_classification)) %>%
      filter(productive_beta & productive_alpha) %>%
      select(!!!syms(all_columns)) %>%
      separate_rows(c(pMHC_classification, pMHC_confidence, pMHC_pvalues, pMHC_wclone_pvalues), sep = ':') %>%
      mutate(HLA = gsub('_.+', '', pMHC_classification)) %>%
      mutate(peptide = gsub('.+_', '', pMHC_classification)) %>%
      mutate(tcr_pmhc = paste0(junction_beta, '_', junction_alpha, '_', peptide)) %>%
      distinct()
  )
}

