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

#' Plot TCR-pMHC Permutation Tests
#'
#' This function creates histograms for permutation tests and highlights, within simulated values per clone, where the observed clone value and stability lie.
#'
#' @param object A Seurat object containing TCR and pMHC data.
#' @param clonal_id The ID of the clonotype to analyze.
#' @param assay The assay to use for retrieving data. Default is 'pMHC'.
#' @param slot The slot to use for retrieving data. Default is 'data'.
#' @param stab_z The stability threshold for determining non-zero UMI counts. Default is 1.
#' @param n_permutations The number of permutations to perform for the test. Default is 1000.
#'
#' @return A combined ggplot object showing histograms of simulated means and stabilities for each pMHC, with observed values highlighted.
#'
#' @examples
#' # Assuming `seurat_obj` is a Seurat object with relevant data
#' plot_tcr_pmhc_permutation(seurat_obj, clonal_id = 'clone_1')
#'
#' @export
plot_tcr_pmhc_permutation <- function(object, clonal_id, assay = 'pMHC', 
                          slot='data', stab_z=1, n_permutations=1000){
  
  library(patchwork)
  
  clone_coords <- which(object$clone_id == clonal_id)
  clone_size <- object$clone_size[object$clone_id == clonal_id] %>% drop.na() %>% unique() 
  
  meta <- object@meta.data %>%
    filter(clone_id == clonal_id) %>% 
    dplyr::select(pMHC_classification, pMHC_confidence, pMHC_pvalues) %>%
    distinct() %>%
    separate_rows(c(pMHC_classification, pMHC_confidence, pMHC_pvalues), sep = ':')
  
  if (sum(meta$pMHC_classification != 'Negative' | !is.na(meta$pMHC_classification)) == 0){
    stop(paste0('clone ', clonal_id, ' doesnt have an assign pmhc specificity'))
  }
  
  pmhcs <- meta %>% pull(pMHC_classification)
  pmhcs <- object@misc$pmhc %>% filter(pmhc %in% pmhcs) %>% pull(Barcode)
  p_values <- meta %>% pull(pMHC_pvalues)
  names(p_values) <- pmhcs
  
  pmhc_mat <- GetAssayData(object, assay = assay, layer=slot)[pmhcs,]
  observed_means <- rowMeans(pmhc_mat[, clone_coords, drop=FALSE])
  observed_stabs <- apply(pmhc_mat[, clone_coords, drop=FALSE], 1, function(x) mean(x > stab_z))
  
  simulated_means <- matrix(NA, nrow = nrow(pmhc_mat), ncol = n_permutations)
  simulated_stabs <- matrix(NA, nrow = nrow(pmhc_mat), ncol = n_permutations)
  
  rownames(simulated_means) <- rownames(pmhc_mat)
  rownames(simulated_stabs) <- rownames(pmhc_mat)
  
  for (i in 1:n_permutations) {
    sampled_columns <- sample(ncol(pmhc_mat), clone_size, replace = FALSE)
    simulated_means[, i] <- rowMeans(pmhc_mat[, sampled_columns, drop=FALSE])
    simulated_stabs[, i] <- apply(pmhc_mat[, sampled_columns, drop=FALSE], 1, function(x) mean(x > stab_z))
  }
  
  plot_list <- list()
  
  for (pmhc in rownames(simulated_means)) {
    print(pmhc)
    current_means <- data.frame(value = simulated_means[pmhc, ])
    current_obs_mean <- observed_means[pmhc]
    current_stabs <- data.frame(value = simulated_stabs[pmhc, ])
    current_obs_stab <- observed_stabs[pmhc]
    
    means_plot <- ggplot(current_means, aes(x = value)) + 
      geom_histogram(binwidth = 0.1, color = "black", fill = "grey") +
      geom_vline(xintercept = current_obs_mean, linetype = "dashed", color = "red", size = 2) +
      labs(title = paste0('simulated means for ', pmhc, 
                          ',\nclone_size=', clone_size,
                          ' p_value = ', p_values[pmhc])) +
      theme_minimal()
    
    stabs_plot <- ggplot(current_stabs, aes(x = value)) + 
      geom_histogram(binwidth = 0.01, color = "black", fill = "blue") +
      geom_vline(xintercept = current_obs_stab, linetype = "dashed", color = "green", size = 2) +
      labs(title =  paste0('simulated stabilities for ', pmhc, ',\nclone_size=', clone_size)) +
      theme_minimal()
    
    combined_plot <- means_plot + stabs_plot + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
    plot_list[[pmhc]] <- combined_plot
    
  }
  final_plot <- patchwork::wrap_plots(plot_list, ncol = 1)
  final_plot 
}


#' Volcano Plot for pMHC-TCR Pairs
#'
#' This function creates a volcano plot using permutation p-values and confidence of pMHC-TCR pairs as an effect size metric.
#'
#' @param object A Seurat object containing TCR and pMHC data.
#' @param conf_threshold The confidence threshold for highlighting significant points. Default is 0.5.
#' @param pval_threshold The p-value threshold for highlighting significant points. Default is -log10(0.05).
#'
#' @return A ggplot object representing the volcano plot.
#'
#' @examples
#' # Assuming `seurat_obj` is a Seurat object with relevant data
#' volcano_plot <- pmhc_volcano_plot(seurat_obj)
#' print(volcano_plot)
#'
#' @export
pmhc_volcano_plot <- function(meta, conf_threshold=.5, pval_threshold=-log10(0.05)) {
  # Convert confidence to log scale if needed
  meta <- extract_pairs(object)
  meta <- meta %>%
    mutate(pMHC_confidence = as.numeric(pMHC_confidence)) %>%
    mutate(pMHC_pvalues = as.numeric(pMHC_pvalues)) %>%
    mutate(log_conf = log10(pMHC_confidence)) %>%
    mutate(min10logp = -log10(pMHC_pvalues))         
  
  # Create the plot
  plot <- ggplot(meta, aes(x = log_conf, y = min10logp)) +
    geom_point(aes(color = (log_conf >= conf_threshold & pMHC_pvalues <= pval_threshold)), size = 2) +
    scale_color_manual(values = c("grey", "red")) +
    geom_vline(xintercept = conf_threshold, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "blue") +
    labs(
      title = "Volcano Plot",
      x = "Log10(Confidence)",
      y = "-Log10(p-value)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Add text annotations for significant points
  significant_points <- subset(meta, log_conf >= conf_threshold & pMHC_pvalues <= pval_threshold)
  plot <- plot + 
    geom_text_repel(
      data = significant_points,
      aes(label = pMHC_classification),
      size = 3,
      box.padding = 0.3,
      point.padding = 0.3,
      max.overlaps = 10
    )
  
  return(plot)
}






