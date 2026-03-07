#' Centers and scales a rows in numeric matrix
#'
#' Centers and scales a rows in numeric matrix, you have an option to replace 0 sd
#' with 0.001 not to cause a division by zero during the scaling
#'
#' @param mat A numeric matrix (`matrix`).
#' @param sd_fun, function (`function`) that calculates sd, functions are expected to
#' have na.rm parameter to remove NA values
#' @param replace_zero (`logical`) if TRUE, replaces zero SD values with 1e-4 before scaling.
#'
#' @return a mean ceantered and median scaled matrix (`matrix`)
#'
#' @examples
#' scale2(mat)
#' scale2(mat, sd_fun=sd)
#'
#' @export
scale2 <- function(mat, sd_fun = sd, replace_zero = FALSE) {
  
  if(is.vector(mat)) {
    mean <- mean(mat, na.rm = T)
    sd <- sd_fun(mat, na.rm = T)
    
    if (replace_zero && sd == 0) {
      sd <- .0001
    }
    return( (mat - mean) / sd )
  }
  else {
    means <- apply(mat, 1, mean, na.rm = T)
    sds <- apply(mat, 1, function(x) sd_fun(x, na.rm = T))
    if (replace_zero) {
      sds[sds == 0] <- .0001
    }
    
    return( (mat - means) / sds )
  }
}

#' Scale without outliers
#'
#' Performs scaling of the matrix iteratively, 
#' removing outliers from sd estimation
#'
#'
#' @param umi_matrix A numeric matrix (`matrix`).
#' @param outlier (`numeric`) from what z-score a value is considered an outlier
#' @param sd_fun (`function`) what function to use to calculate standard deviation
#' @return a mean ceatered and median scaled matrix (`matrix`)
#'
#' @examples
#' super_scale(umi_matrix=mat)
#'
#' @export

super_scale <- function(umi_matrix, outlier=3, sd_fun=sd) {
  scaled_umi_matrix <- scale2(umi_matrix, sd_fun)  
  umi_matrix[abs(scaled_umi_matrix) > outlier] <- NA
  return(umi_matrix)
}

#' Perform the scaling with no outliers on your Seurat object
#'
#' Perform the scaling with no outliers on your Seurat object. Your Seurat object
#' is expected to have a certain structure. 
#' pMHC count matrix should be stored as separate assay called "pMHC"
#' slot "misc" is expected to have an annotation for your pMHC barcodes with following 
#' columns: "Barcode" - barcode ID (same as rownames in pMHC matrix) 
#'          "Sequence" - peptide sequence
#'          "HLA" - HLA allele
#'
#' @param object (`Seurat`) A Seurat object that must be of the above described 
#' structure 
#' @param outlier (`numeric`) from what z-score a value is considered an outlier
#' @param return_params (`logical`) whether you want to return a lust sds and 
#' means along with the object
#' @param verbose whether to show a progress bar
#'
#' @return (`Seurat`) with the ScaleData slot in pMHC assay filled with scaled 
#' pMGC matrix
#'
#' @examples
#' ScaleDataNoOutliers(object)
#'
#' @export
ScaleDataNoOutliers <- function(object, slot = 'counts', assay = 'pMHC', 
                                outlier = 3, return_params = FALSE, verbose = TRUE) {
  
  if (!is(object, "Seurat")) {
    stop("The provided object is not a Seurat object.")
  }
  
  input_matrix <- GetAssayData(object, assay = assay, layer = slot)
  umi_matrix <- input_matrix
  iter0 <- apply(umi_matrix, 1, sd, na.rm = TRUE)
  mean_iter0 <- apply(umi_matrix, 1, mean, na.rm = TRUE)
  iter <- 0
  sds_iter <- list()
  sds_stats_iter <- list()
  means_iter <- list()
  
  initial_outliers <- sum(abs(scale2(umi_matrix, sd_fun = sd)) > outlier, na.rm = TRUE)
  current_outliers <- initial_outliers
  
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = initial_outliers, style = 3)
  }
  
  while (any(abs(scale2(umi_matrix, sd_fun = sd)) > outlier, na.rm = TRUE)) {
    iter <- iter + 1
    if (verbose) {
      cat(sprintf("\nIteration %d\n", iter))
    }
    umi_matrix <- super_scale(umi_matrix, sd_fun = sd)
    nzeroes <- apply(umi_matrix, 1, function(x) sum(x[x != 0], na.rm = TRUE))
    
    sds_stats_iter[[paste0('iteration=', iter)]] <- nzeroes
    sds_iter[[paste0('iteration=', iter)]] <- apply(umi_matrix, 1, sd, na.rm = TRUE)
    means_iter[[paste0('iteration=', iter)]] <- apply(umi_matrix, 1, mean, na.rm = TRUE)
    
    if (verbose) {
      current_outliers <- sum(abs(scale2(umi_matrix, sd_fun = sd)) > outlier, na.rm = TRUE)
      setTxtProgressBar(pb, initial_outliers - current_outliers)
    }
  }
  
  if (verbose) {
    close(pb)
  }
  
  sds_iter[['iteration=0']] <- iter0
  means_iter[['iteration=0']] <- mean_iter0
  sds_iter <- do.call(cbind, sds_iter)
  means_iter <- do.call(cbind, means_iter)
  sds_iter <- sds_iter[, paste('iteration=', 0:(ncol(sds_iter)-1), sep = '')]
  means_iter <- means_iter[, paste('iteration=', 0:(ncol(means_iter)-1), sep = '')]
  
  get_last_nonzero <- function(row) {
    nonzero_values <- row[row != 0]
    if (length(nonzero_values) > 0) {
      return(tail(nonzero_values, n = 1))
    } else {
      return(0)
    }
  }
  
  sds <- apply(sds_iter, 1, get_last_nonzero)
  
  centered <- as.matrix(
    input_matrix - apply(umi_matrix, 1, mean, na.rm = TRUE)
  )
  
  object@assays[[assay]]@data <- centered / sds
  
  object@commands$ScaleDataNoOutliers <- TRUE
  
  if (return_params) {
    return(list('object' = object, 'sds' = sds_iter, 'means' = means_iter))
  } else {
    return(object)
  }
}

cap_high_outliers <- function(values, upper_quantile = 0.95, replacement = "quantile") {
  
  upper_threshold <- quantile(values, upper_quantile, na.rm = TRUE)
  
  if (replacement == "quantile") {
    values[values > upper_threshold] <- upper_threshold
  } else if (replacement == "mean") {
    non_outliers <- values[values <= upper_threshold]
    values[values > upper_threshold] <- mean(non_outliers, na.rm = TRUE)
  } else {
    stop("Invalid replacement method. Choose 'quantile', 'mean', or 'scaled'.")
  }
  
  return(values)
}


#' Perform smoothing of scaled pMHC counts for every big clone in your Seurat object
#'
#' Perform the smoothing of scaled pMHC counts for every big clone in your 
#' Seurat object. Your Seurat object is expected "clone_id" column in your meta.data
#' with clonotype assigments.
#'
#' @param object (`Seurat`) A Seurat object that must be of the above described 
#' structure 
#' @param best_params (`list`) ignore it for now
#' smooth scaled by ScaleDataNoOutliers counts
#' @param cap_upper_quantiles (`bool`) whether to trim higher than .9 quantile values with .9 quantile
#' @param assay (`character`) name of the slot with pMHC counts, default "pMHC"
#' @param normalise (`bool`) parameter passed to LOESS
#' @param cl_size_thresh (`numeric`) from what clone size a clone will be smoothed
#' @param span_val (`numeric`) value for span parameter in loess
#' @param replace_ones (`bool`) whether you want to replace 1 UMI count with 0
#' @param degree_val (`numeric`) value for degree parameter in loess
#' @param family_val (`character`) value for family parameter in loess
#' @param verbose whether to show a progress bar
#' 
#' @return (`Seurat`) with the ScaleData slot with scaled pMHC counts smoothed 
#'
#' @examples
#' smooth_pmhc(object)
#'
#' @export

smooth_pmhc <- function(object, assay = 'pMHC', cap_upper_quantiles=T,
                        cl_size_thresh = 3, normalise = FALSE, span_val = 1, replace_ones=F,
                        degree_val = 1, family_val = "symmetric", replacement_q=.85, verbose = TRUE) {
  
  if (cl_size_thresh<2){
    stop('cl_size_thresh must be minimum 2')
  }
  
  if (is.null(object@commands$ScaleDataNoOutliers)){
    stop('you have to run ScaleDataNoOutliers() before running smooth_pmhc')
  }
  
  if (!"clone_size" %in% names(object@meta.data)) {
    object@meta.data <- object@meta.data %>%
      tibble::rownames_to_column('row_id') %>% 
      group_by(clone_id) %>%
      mutate(clone_size = n()) %>%
      ungroup() %>% 
      tibble::column_to_rownames('row_id')
    
    object$clone_size[is.na(object$clone_id)] <- NA
  }
  
  bigger_thanXclones <- object$clone_id[object$clone_size > cl_size_thresh] %>% drop.na() %>% unique()
  
  if (length(bigger_thanXclones) == 0){
    stop('there are no expanded clones in the dataset, hence the pMHC counts are 
          unsmoothable, if you study naive repertoirs, you can try to skip this step
          doing assign_mhc(..., force=T)')
  }
  
  smoothable_cells <- Cells(object)[object$clone_id %in% bigger_thanXclones]
  unsmoothable_cells <- Cells(object)[!Cells(object) %in% smoothable_cells]
  
  scaled_counts <- GetAssayData(object, assay = assay, layer = 'data')
  raw_countss <- GetAssayData(object, assay = assay, layer = 'counts')
  pmhcs <- rownames(scaled_counts)
  smoothed_counts <- matrix(0, nrow = length(pmhcs), ncol = length(smoothable_cells), dimnames = list(pmhcs, smoothable_cells))
  
  unsmoothable_counts <- scaled_counts[, unsmoothable_cells]
  scaled_counts <- scaled_counts[, smoothable_cells]
  raw_countss <- raw_countss[, smoothable_cells]
  
  clone_assignments <- object@meta.data[smoothable_cells,]$clone_id
  
  if (verbose) {
    total_iterations <- length(bigger_thanXclones)
    pb <- txtProgressBar(min = 0, max = total_iterations, style = 3, width = 50, char = "+")
    cat(sprintf('\nsmoothing each pMHC within each clone bigger than %d \n', cl_size_thresh))
  }
  
  index <- 0
  object@commands$smooth_pmhc <- list()
  
  for (clone_id in bigger_thanXclones) {
    clone_idx <- which(clone_assignments == clone_id)
    clone_cells <- smoothable_cells[clone_idx]
    
    for (pmhc in pmhcs) {
      counts <- scaled_counts[pmhc, clone_idx]
      
      if (all(is.nan(counts) | is.na(counts))){
        next
      }
      
      smoothed_counts[pmhc, clone_idx] <- counts
      
      
      
      ratio_zero <- sum(counts == 0)/length(counts)
      one_outlier <- sum(counts < 0) == 1 | ratio_zero > 0.9
      if (replace_ones & one_outlier){
        smoothed_counts[pmhc, clone_idx] <- rep(0, length(counts))
      }
      
      error_occurred <- tryCatch({
        if (cap_upper_quantiles){
          counts <- cap_high_outliers(counts, upper_quantile = replacement_q)
        }
        t <- try({
          loess_fit <- loess(counts ~ seq_along(counts), normalize = normalise, 
                             span = span_val, degree = degree_val, family = family_val)
          smoothed_counts[pmhc, clone_idx] <- predict(loess_fit)
          object@commands$smooth_pmhc[[clone_id]] <- 'done'
        }, silent = TRUE)
        
        if("try-error" %in% class(t)) {
          if (family_val == "symmetric") {
            t_gaussian <- try({
              loess_fit <- loess(counts ~ seq_along(counts), normalize = normalise, 
                                 span = span_val, degree = degree_val, family = "gaussian")
              smoothed_counts[pmhc, clone_idx] <- predict(loess_fit)
              object@commands$smooth_pmhc[[clone_id]] <- 'done_with_gaussian'
            }, silent = TRUE)
            
            if("try-error" %in% class(t_gaussian)) {
              cat(sprintf('\nError in fitting loess for pmhc %s in clone %s with both symmetric and gaussian families', pmhc, clone_id))
              object@commands$smooth_pmhc[[clone_id]] <- 'error'
            }
          } else {
            cat(sprintf('\nError in fitting loess for pmhc %s in clone %s with gaussian family', pmhc, clone_id))
            object@commands$smooth_pmhc[[clone_id]] <- 'error'
          }
        }
        ####
        
      }, error = function(e) {
        error_both_families <- ("try-error" %in% class(t)) & (exists("t_gaussian") && "try-error" %in% class(t_gaussian))
        if (verbose & error_both_families) {
          cat(sprintf('\nerror in fitting loess for pmhc %s in clone %s', pmhc, clone_id))
        }
        object@commands$smooth_pmhc[[clone_id]] <- 'error'
      })
      
      if (verbose) {
        setTxtProgressBar(pb, index)
      }
    }
    index <- index + 1
  }
  
  if (verbose) {
    close(pb)
  }
  
  smoothed_counts <- cbind(smoothed_counts, unsmoothable_counts)
  smoothed_counts[is.nan(smoothed_counts)] <- 0
  object@assays[[assay]]@scale.data <- as.matrix(smoothed_counts[, Cells(object)])
  
  return(object)
}


filter_by_z_score <- function(scaled_mat, z_score_threshold){
  valid_columns <- apply(scaled_mat, 2, function(x) any(x > z_score_threshold))
  valid_columns[is.na(valid_columns)] <- FALSE
  
  filtered_mat <- scaled_mat[, valid_columns]
  
  return(filtered_mat)
}

#' Perform permutation test for pairs within assign_pmhc
#' 
#' @param pmhc_mat (`matrix`) pMHC UMI scaled matrix
#' @param barcodes_to_test (`character` vector) List of barcodes for specific multimers
#' @param clone_coords (`integer` vector) Indices of the clone in the matrix
#' @param stab_z (`numeric`) Z-score threshold to call the gem positive for a given pMHC
#' @param df (`integer`) Degrees of freedom for Fisher method, uniting p-values, must be 4 for 2 tests
#' @param sample_permutation (`character`) Whether to downsample values within the same pMHC 
#' (`"within_pmhc"`) or across any random pMHC (`"random"`)
#' @param n_permutations (`integer`) Number of generated pseudoclones
#' @param p_adj_method (`character`) Method to use for multiple testing adjustment (e.g., `"BH"`)
#' @return (`numeric` vector) Adjusted p-values for each barcode
#'
#' @examples
#' tcr_pmhc_permutation.test()
#'
#' @export
tcr_pmhc_permutation.test <- function(pmhc_mat, barcodes_to_test, 
                                      clone_coords, stab_z=1, df=4,
                                      sample_permutation='within_pmhc',
                                      n_permutations=1000, p_adj_method='none') {
  
  cl_size <- length(clone_coords)
  
  sub_obs <- pmhc_mat[, clone_coords, drop = FALSE]
  
  observed_means <- setNames(rowMeans(sub_obs), rownames(sub_obs))[barcodes_to_test]
  observed_stabs <- setNames(rowMeans(sub_obs > stab_z), rownames(sub_obs))[barcodes_to_test]
  
  simulated_means <- matrix(NA, nrow = nrow(pmhc_mat), ncol = n_permutations)
  simulated_stabs <- matrix(NA, nrow = nrow(pmhc_mat), ncol = n_permutations)
  
  rownames(simulated_means) <- rownames(pmhc_mat)
  rownames(simulated_stabs) <- rownames(pmhc_mat)
  
  if (sample_permutation == 'random'){
    barcodes <- rownames(pmhc_mat)
    pmhc_mat <- apply(pmhc_mat, 2, sample)
    rownames(pmhc_mat) <- barcodes
  }
  
  for (i in 1:n_permutations) {
    sampled_columns <- sample(ncol(pmhc_mat), cl_size, replace = FALSE)
    sub_sim <- pmhc_mat[, sampled_columns, drop = FALSE]
    simulated_means[, i] <- rowMeans(sub_sim, na.rm = T)
    simulated_stabs[, i] <- rowMeans(sub_sim > stab_z, na.rm = T)
  }
  
  means_p_values <- sapply(barcodes_to_test, 
                           function(pmhc) mean(simulated_means[pmhc,] >= observed_means[pmhc], na.rm=T))
  stabs_p_values <- sapply(barcodes_to_test, 
                           function(pmhc) mean(simulated_stabs[pmhc,] >= observed_stabs[pmhc], na.rm=T))
  
  combined_p_values <- mapply(function(p1, p2) {
    chi_sq_statistic <- -2 * (log(p1) + log(p2))
    pchisq(chi_sq_statistic, df = df, lower.tail = FALSE)
  }, means_p_values, stabs_p_values) 
  
  return(
    p.adjust(combined_p_values, method = p_adj_method) %>% round(., digits = 3)
  )
}


#' Perform permutation test for a specific pMHC-TCR pair
#' 
#' Perform permutation test for a specific pMHC-TCR pair
#' @param object (`Seurat`) Seurat object 
#' @param bc_or_pmhc (`character` vector) vector of barcodes or pmhc names
#' @use_pmhc (`bool`) whether bc_pr_pmhc params supplied with bc or pmhc names
#' @clone (`character`) name of the clone within which to test, keep as NULL if you want to identify clone with full_tcr
#' @full_tcr (`full_tcr`) name of the full_tcr of the clone totest within, keep as NULL if you want to identify clone with clone_id
#'
#' @return (`numeric`)  
#'
#' @examples
#' permutation_test_specific_pair()
#'
#' @export
permutation_test_specific_pair <- function(object, bc_or_pmhc, use_pmhc=T, drop.na=F,
                                           clone=NULL, full_tcr=NULL, perm_test_adj_method='BH', ...){
  
  pmhc_mat <- GetAssayData(object, layer = 'data', assay = 'pMHC')
  
  if (use_pmhc){
    bc_to_pmhc <- setNames(object@misc$pmhc$pmhc, object@misc$pmhc$Barcode)
    rownames(pmhc_mat) <- recode(rownames(pmhc_mat), !!!bc_to_pmhc)
  }
  
  if (drop.na){
    pmhc_mat <- pmhc_mat[bc_or_pmhc, ,drop=FALSE]
    pmhc_mat <- pmhc_mat[, !is.na(pmhc_mat[1, ]), drop = FALSE] 
    
    clone_coords <- which(colnames(pmhc_mat) %in% Cells(object)[object$clone_id == clone])
  } else{
    if (is.null(full_tcr)){
      clone_coords <- which(object$clone_id == clone)
    } else if (is.null(clone)){
      clone_coords <- which(object$full_tcr == full_tcr)
    }
  }
  
  p_values <- tcr_pmhc_permutation.test(pmhc_mat=pmhc_mat, 
                                        barcodes_to_test=bc_or_pmhc,
                                        sample_permutation='within_pmhc',
                                        clone_coords=clone_coords, 
                                        p_adj_method = perm_test_adj_method,
                                        ...)
  
  return(p_values)
}

#' Adjust pMHC P-values Within Clones
#'
#' This function processes a data frame containing pMHC classifications, confidence scores, and p-values,
#' separating multiple p-values per row and adjusting them using `p.adjust()`. The results are re-collapsed
#' into the original format, keeping clone-level grouping.
#'
#' @param df A data frame containing `pMHC_pvalues` and `pMHC_wclone_pvalues` as `:`-separated strings.
#' @param adjust_permutation Logical. If `TRUE`, adjust `pMHC_pvalues` using `p.adjust()`.
#' @param adjust_test_within_clone Logical. If `TRUE`, adjust `pMHC_wclone_pvalues` using `p.adjust()`.
#' @param method Character string specifying the multiple testing correction method for `p.adjust()`. Default is `"BH"`.
#'
#' @return A modified data frame with adjusted p-values, where each `pMHC_pvalues` and `pMHC_wclone_pvalues`
#'         are re-collapsed back into `:`-separated format after adjustment.
#'
#' @details The function first splits the `pMHC_pvalues` and `pMHC_wclone_pvalues` columns into separate rows,
#'          applies `p.adjust()` per `clone_id`, and then re-aggregates them back into their original structure.
#'
#' @examples
#' df_adjusted <- adjust_pmhc_pvalues(df, adjust_permutation = TRUE, adjust_test_within_clone = TRUE)
#'
#' @export
adjust_pmhc_pvalues <- function(df, adjust_permutation, adjust_test_within_clone, method = "BH") {
  
  library(purrr)
  
  df_long <- df %>%
    separate_rows(c(pMHC_classification, pMHC_pvalues, 
                    pMHC_wclone_pvalues, pMHC_confidence, 
                    pMHC_deltas, pMHC_scaled_umis), sep=":")
  
  if (adjust_permutation) {
    df_long <- df_long %>%
      group_by(clone_id) %>%
      mutate(pMHC_pvalues = p.adjust(pMHC_pvalues, method = method)) %>%
      ungroup()
  }
  
  if (adjust_test_within_clone) {
    df_long <- df_long %>%
      group_by(clone_id) %>%
      mutate(pMHC_wclone_pvalues = p.adjust(pMHC_wclone_pvalues, method = method)) %>%
      ungroup()
  }
  
  df_wide <- df_long %>%
    group_by(clone_id) %>%
    summarise(across(starts_with("pMHC_"), ~paste(., collapse = ":"), .names = "collapsed_{col}")) %>%
    rename_with(~ gsub("^collapsed_", "", .), starts_with("collapsed_"))
  
  return(df_wide)
}


#' Perform pMHC Assignment for Each Clone
#'
#' This function performs pMHC assignment for each "big" clone in a Seurat object, with the option to 
#' include small clones as well. The Seurat object must contain a pMHC count matrix as a separate assay 
#' named "pMHC", and metadata in the "misc" slot should include columns for pMHC barcode annotation:
#' - "Barcode": barcode ID (matching row names in the pMHC matrix)
#' - "Sequence": peptide sequence
#' - "HLA": HLA allele
#'
#' @param object (`Seurat`) A Seurat object with the required structure described above.
#' @param slot (`character`) The data slot to use within the specified assay, default is "scale.data".
#' @param assay (`character`) The name of the assay with pMHC counts, default is "pMHC".
#' @param assignment (`character`) The method for identifying outliers, can be one of 
#'   "z-score", "delta_gap", "rosner", or "extreme_distribution".
#' @param entropy_thresh (`numeric`) Threshold above which pMHCs are excluded from TCR-pMHC pairing for 
#'   clones smaller than `cl_size_thresh` due to noise.
#' @param cl_size_thresh (`numeric`) Minimum clone size for using smoothed z-score values, default is 3.
#' @param delta_threshold (`numeric`) Threshold to separate background from specificity signal, 
#'   recommended values between 0.7 and 1. Used in delta_gap assigment.
#' @param assign_small_clones (`logical`) Whether to assign pMHCs to small clones (smaller than `cl_size_thresh`), 
#'   which may reduce precision. Defaults to `FALSE`.
#' @param kmax (`numeric`) The maximum number of potential outliers, used in `rosner` assignment.
#' @param rosner_alpha (`numeric`) Significance level for Rosner's test, used in `rosner` assignment.
#' @param extreme_alpha (`numeric`) Significance level for `extreme_distribution` assignment.
#' @param predefined_extreme_params (`list` or `NULL`) Predefined extreme distribution parameters for `extreme_distribution` assignment.
#' @param filter_by_zscore (`logical`) Whether to filter small clones by z-score threshold, 
#' assuming that the absense of any extreme values would result in no specificities, thus saving
#' computational time. 
#' @param filter_by_zscore_big_clones (`logical`) Whether to filter large clones by z-score threshold,
#' #' assuming that the absense of any extreme values would result in no specificities, thus saving
#' computational time.
#' @param z_score_threshold (`numeric`) z-score threshold for filtering clones.
#' @param top_down (`logical`) If `TRUE`, uses a top-down approach in delta gap assignment, if `FALSE` then uses bottom-up approach,
#' that leads to lessstrict assigments, default is `FALSE`.
#' @param calculate_pvalue (`logical`) Whether to calculate p-values for TCR-pMHC assignments.
#' @param calculate_confidence_score (`logical`) Whether to calculate confidence scores for TCR-pMHC assignments.
#' @param alpha, beta, gamma, delta (`numeric`) Parameters used in confidence score calculations.
#' @param force (`logical`) If `TRUE`, forces the function to run even if preconditions are not met.
#' @param verbose (`logical`) If `TRUE`, displays progress information.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return (`Seurat`) The modified Seurat object, with TCR-pMHC classification and scores added to `meta.data`.
#'
#' @examples
#' assign_pmhc(object)
#'
#' @export

assign_pmhc <- function(object, slot='scale.data', assay='pMHC', clones_to_analyse=NULL,
                        assignment='rosner', params = 'remove_topx', remove_for_params = 2, assign_small_clones=F, 
                        cl_size_thresh=3, entropy_thresh=1, kmax=10, rosner_alpha=0.01, rosner_pval_dist='t',
                        extreme_alpha=0.001, delta_threshold = 1, pseudobulk_fun=median, n_tests=10, extreme_type = 'regular',
                        double_loc_scale=F, filter_by_zscore=FALSE, filter_by_zscore_big_clones=FALSE,
                        adjust_permutation = FALSE, adjust_test_within_clone = FALSE, 
                        perm_test_adj_method = 'none', padj_method = "BH",
                        z_score_threshold=2, calculate_pvalue=TRUE, calculate_confidence_score=TRUE,
                        alpha=1, beta=1, gamma=1, delta=1, force=F, print_clone_id=FALSE, verbose=TRUE, ...) {
  library(shades)
  library(EnvStats)
  
  if (cl_size_thresh<2) {
    stop('cl_size_thresh must be minimum 2')
  }
  
  if (!params %in% c('iteratively', 'remove_topx', 'extreme_params')) {
    stop('params  must me either iteratively, remove_topx or extreme_params')
  }
  
  if (is.null(object@commands$smooth_pmhc) & !force) {
    stop('you have to run smooth_pmhc() before running assign_pmhc')
  }
  if (is.null(object@misc$noise_score) & assign_small_clones & !force) {
    stop('you have to run score_pmhc_noise() before running assign_pmhc with assign_small_clones=T')
  }
  if (is.null(object@misc$pmhc)) {
    object$pMHC_classification <- 'Negative'
    if (verbose) cat('No pmhc barcode annotation is found in object@misc$pmhc\n')
    return(object)
  }
  if (!assignment %in% c('delta_gap', 'rosner', 'extreme_distribution')) {
    stop('assigment must be either z-score, or delta_gap or extreme_distribution, or check your spelling')
  }
  
  if (!"clone_size" %in% names(object@meta.data)) {
    object@meta.data <- object@meta.data %>%
      tibble::rownames_to_column('row_id') %>% 
      group_by(clone_id) %>%
      mutate(clone_size = n()) %>%
      ungroup() %>% 
      tibble::column_to_rownames('row_id')
    
    object$clone_size[is.na(object$clone_id)] <- NA
  }
  
  clones <- object$clone_id[!is.na(object$clone_id)]
  
  if (!is.null(clones_to_analyse)){
    submited_clones <- clones[clones %in% clones_to_analyse] %>% unique()
    clone_bulk <- ClonePseudobulk(object = object %>% subset(clone_id %in% submited_clones), 
                                  assay = assay, clone_col = 'clone_id', 
                                  slot = slot, filter_by_zscore=FALSE)
  } else{
    
    big_clones <- names(table(clones)[table(clones) >= cl_size_thresh])
    
    clone_bulk <- NULL
    
    if (length(big_clones) > 0) {
      if (verbose) {
        cat('\nAggregating clonotype aggregated pseudobulk pMHC matrix\n')
      }
      clone_bulk <- ClonePseudobulk(object %>% subset(clone_id %in% big_clones), 
                                    assay = assay, clone_col = 'clone_id', 
                                    slot = slot, 
                                    filter_by_zscore=filter_by_zscore_big_clones,
                                    z_score_threshold=z_score_threshold)
    } else if (verbose) {
      cat('\nNo big clones found, skipping pseudobulk aggregation\n')
    }
    
    if (assign_small_clones){
      small_clones <- object$clone_id[!object$clone_id %in% big_clones] %>% drop.na() %>% unique()
      
      if (!is.null(clones_to_analyse)){
        small_clones <- small_clones[small_clones %in% clones_to_analyse]
      }
      
      single_gems <- ClonePseudobulk(object %>% subset(clone_id %in% small_clones), 
                                     assay = assay, clone_col = 'clone_id', 
                                     slot = slot, 
                                     filter_by_zscore=filter_by_zscore_big_clones,
                                     z_score_threshold=z_score_threshold)
      
      if (filter_by_zscore){
        single_gems <- filter_by_z_score(single_gems, 
                                         z_score_threshold = z_score_threshold)
      }
      
      if (!is.null(clone_bulk)) {
        clone_bulk <- cbind(clone_bulk, single_gems)
      } else {
        clone_bulk <- single_gems
      }
  }
    
    if (verbose) {
      cat(sprintf("\nRemoving features with entropy smaller than %g for clones > than %g gems\n", entropy_thresh, cl_size_thresh))
    }
    high_entropy_pmhc <- rownames(object@misc$noise_score)[object@misc$noise_score$entropy > entropy_thresh]
    small_clones <- colnames(clone_bulk)[colnames(clone_bulk) %in% object$clone_id[object$clone_size < cl_size_thresh]]
    clone_bulk[,small_clones][high_entropy_pmhc,] <- -10
  }
  
  n_iter <- ncol(clone_bulk)
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "+")
    cat('\nassigning pMHC-TCR pairs\n')
  }
  clone_outliers <- list()
  data_layer <- GetAssayData(object, assay = assay, layer = "data")
  ##################
  ####  L O O P ####
  ##################
  for (i in 1:n_iter) {
    if (print_clone_id){
      clone_print <- object$clone_id[object$clone_id == colnames(clone_bulk)[i]] %>% unique() %>% drop.na()
      print(clone_print)
    }
    cl_size <- object$clone_size[object$clone_id == colnames(clone_bulk)[i]] %>% unique() %>% drop.na()
    scores <- sort(clone_bulk[,i], decreasing = F)
    scores[scores < 0] <- 0
    
    scores_interpolated <- normalize_vector(scores[names(scores)])
    
    clone_coords <- which(object$clone_id == colnames(clone_bulk)[i]) 
    
    if (assignment == 'rosner') {
      
      if (length(scores[scores != 0]) < 25 & length(scores)>25){
        scores <- scores[order(scores)[(length(scores)-26):length(scores)]]
      } else if (length(scores)>25) { # remove else if (length(scores)>25) to just else if big changes
        scores <- scores[scores != 0]
      }
      
      rosner <- tryCatch({
        library(shades)
        library(EnvStats)
        rosnerTest2(scores, k = kmax, alpha = rosner_alpha, pval_dist=rosner_pval_dist,
                    params = params, remove_for_params = remove_for_params)
      }, error = function(e) {
        rosner <- 'failed'
      })
      
      if (is.character(rosner) && rosner == "failed") {
        sprintf('rosner test failed for clone %s', colnames(clone_bulk)[i])
        clone_outliers[[colnames(clone_bulk)[i]]] <- list('confidence_scores' = 'Negative', 
                                                          'p_values' = 'Negative',
                                                          'within_clone_p'='Negative')
        next
      }
      
      if (all(!rosner$all.stats$Outlier)) {
        clone_outliers[[colnames(clone_bulk)[i]]] <- list('confidence_scores'='Negative', 
                                                          'p_values'='Negative',
                                                          'within_clone_pvalues'='Negative')
        next
      }
      
      gap_pos <- rosner$all.stats %>% filter(Outlier) %>% slice_min(order_by = Value, n = 1) %>% pull(Value)
      within_clone_pvalues <- rosner$all.stats %>% filter(Outlier) %>% pull()
      
      outliers <- scores[scores >= gap_pos]
      
      background <- scores_interpolated[!names(scores_interpolated) %in% names(outliers)]
      positives <- scores_interpolated[names(scores_interpolated) %in% names(outliers)]  %>% round(., digits=3)
      
      deltas <- positives - mean(background)%>% round(., digits=3)
      
      within_clone_pvalues <- rosner$all.stats %>% filter(Outlier) %>% pull(P.Values) %>% rev() %>% round(., digits=3)
      
    } else if (assignment == 'extreme_distribution') {
      
      extreme_p <- tryCatch({
        extreme_outlier_test(x = scores, type = extreme_type, double_loc_scale = double_loc_scale)
      }, error = function(e) {
        extreme_p <- 'failed'
      })
      
      if (is.character(extreme_p) && extreme_p == "failed") {
        sprintf('extreme test failed for clone %s', colnames(clone_bulk)[i])
        clone_outliers[[colnames(clone_bulk)[i]]] <- list('confidence_scores' = 'Negative', 
                                                          'p_values' = 'Negative',
                                                          'within_clone_p'='Negative')
        next
      }
      
      if (sum((extreme_p)<extreme_alpha)==0) {
        clone_outliers[[colnames(clone_bulk)[i]]] <- list('confidence_scores'='Negative', 
                                                          'p_values'='Negative',
                                                          'within_clone_pvalues'='Negative')
        next
      }
      
      outliers <- scores[which(extreme_p<extreme_alpha)]
      
      background <- scores_interpolated[!names(scores_interpolated) %in% names(outliers)]
      positives <- scores_interpolated[names(scores_interpolated) %in% names(outliers)] %>% round(., digits=3)
      
      deltas <- (positives - mean(background)) %>% round(., digits=3)
      
      within_clone_pvalues <- round(extreme_p[which(extreme_p<extreme_alpha)], digits = 3)
    }
    
    pmhc_mat <- data_layer[names(outliers), , drop = FALSE]
    
    if (is.null(dim(pmhc_mat))) {
      pmhc_mat <- matrix(pmhc_mat, nrow = 1, byrow = TRUE,
                         dimnames = list(names(outliers), names(pmhc_mat)))  
    }
    
    p_values <- tryCatch(
      {
        tcr_pmhc_permutation.test(
          pmhc_mat = pmhc_mat,
          barcodes_to_test = names(outliers),
          sample_permutation = "within_pmhc",
          p_adj_method = perm_test_adj_method,
          clone_coords = clone_coords
        )
      },
      error = function(e) {
        message("Skipping iteration due to error: ", conditionMessage(e))
      }
    )
    if (is.null(p_values)){
      next
    }
    
    entropy <- object@misc$noise_score[names(outliers),]$entropy
    clone_data <- data_layer[, clone_coords, drop = FALSE]
    concordance <- calculate_pmhc_concordance(clone_data, assay = assay, slot = "data")
    
    concordance <- concordance[names(outliers)]
    
    confidence_scores <- calculate_confidence(entropy=entropy, concordances=concordance, 
                                              deltas=deltas, clone_size = cl_size, 
                                              alpha = alpha, beta = beta, gamma = gamma, delta = delta)
    
    clone_outliers[[colnames(clone_bulk)[i]]] <- list('confidence_scores'=confidence_scores[names(outliers)], 
                                                      'p_values'=p_values[names(outliers)],
                                                      'within_clone_pvalues'=within_clone_pvalues,
                                                      'deltas'=deltas,
                                                      'scaled_umis'=scores[names(positives)] %>% round(., digits = 3))
    
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose) {
    close(pb)
    cat('\nmerging TCR-pMHC information into objects @meta.data\n')
  }
  
  pmhc_bc <- setNames(object@misc$pmhc$pmhc, object@misc$pmhc$Barcode)
  df <- do.call(rbind, lapply(names(clone_outliers), function(x) {
    confidences <- clone_outliers[[x]]$confidence_scores
    p_values <- clone_outliers[[x]]$p_values
    wclone_pvalues <- clone_outliers[[x]]$within_clone_pvalues
    deltas <- clone_outliers[[x]]$deltas
    scaled_umis <- clone_outliers[[x]]$scaled_umis
    
    confidence_is_negative <- length(confidences) == 1 && confidences == "Negative"
    pvalue_is_negative <- length(p_values) == 1 && p_values == "Negative"
    wpvalue_is_negative <- length(wclone_pvalues) == 1 && wclone_pvalues == "Negative"
    deltas_is_negative <- length(deltas) == 1 && deltas == "Negative"
    scaled_umis_is_negative <- length(scaled_umis) == 1 && scaled_umis == "Negative"
    
    data.frame(clone_id = x,
               pMHC_classification = if (confidence_is_negative) "Negative" else 
                 paste(names(clone_outliers[[x]]$confidence_scores) %>% recode(!!!pmhc_bc), collapse=":"),
               pMHC_confidence = if (confidence_is_negative) 'Negative' else 
                 paste(clone_outliers[[x]]$confidence_scores, collapse=":"),
               pMHC_pvalues = if (pvalue_is_negative) 'Negative' else 
                 paste(clone_outliers[[x]]$p_values, collapse=":"),
               pMHC_wclone_pvalues = if (wpvalue_is_negative) 'Negative' else 
                 paste(clone_outliers[[x]]$within_clone_pvalues, collapse=":"),
               pMHC_deltas = if (deltas_is_negative) 'Negative' else 
                 paste(clone_outliers[[x]]$deltas, collapse=":"),
               pMHC_scaled_umis = if (scaled_umis_is_negative) 'Negative' else 
                 paste(clone_outliers[[x]]$scaled_umis, collapse=":")   
    )
  }))
  df$clone_id <- gsub('-', '_', df$clone_id)
  
  if (adjust_permutation | adjust_test_within_clone){
    df <- adjust_pmhc_pvalues(df, 
                              adjust_permutation=adjust_permutation,
                              adjust_test_within_clone=adjust_test_within_clone,
                              method=padj_method)
  }
  
  pmhc_cols <- c(
    "pMHC_classification",
    "pMHC_confidence",
    "pMHC_pvalues",
    "pMHC_wclone_pvalues",
    "pMHC_deltas",
    "pMHC_scaled_umis"
  )
  
  if (is.null(clones_to_analyse)) {
    object@meta.data <- object@meta.data %>%
      tibble::rownames_to_column("row_id") %>%
      dplyr::select(-any_of(pmhc_cols)) %>%
      dplyr::left_join(df, by = "clone_id") %>%
      tibble::column_to_rownames("row_id")
  } else {
    old_metadata <- object@meta.data %>%
      dplyr::filter(!clone_id %in% clones_to_analyse)
    
    new_metadata <- object@meta.data %>%
      dplyr::filter(clone_id %in% clones_to_analyse) %>%
      tibble::rownames_to_column("row_id") %>%
      dplyr::select(-any_of(pmhc_cols)) %>%
      dplyr::left_join(df, by = "clone_id") %>%
      tibble::column_to_rownames("row_id")
    
    metadata <- rbind(old_metadata, new_metadata)
    object@meta.data <- metadata[Cells(object), , drop = FALSE]
  }
  
  return(object)
}


check_pmhc_consistency <- function(object, sep = ":", stop_on_error = FALSE) {
  md <- object@meta.data
  rn <- rownames(md)
  
  cols_base <- c(
    "pMHC_classification", "pMHC_confidence", "pMHC_pvalues",
    "pMHC_wclone_pvalues", "pMHC_deltas", "pMHC_scaled_umis"
  )
  cols <- cols_base[cols_base %in% colnames(md)]
  if ("epitope_type_multiple" %in% colnames(md))
    cols <- c(cols, "epitope_type_multiple")
  
  if (!length(cols)) {
    warning("No pMHC columns found to check.")
    return(tibble::tibble())
  }
  
  cnt <- function(x) {
    sx <- as.character(x)
    if (length(sx) == 0L || is.na(sx) || sx == "" || sx == "Negative" || sx == "NA") return(0L)
    parts <- strsplit(sx, sep, fixed = TRUE)[[1]]
    sum(nzchar(parts))
  }
  
  counts <- as.data.frame(lapply(md[, cols, drop = FALSE], function(v) {
    vapply(v, cnt, integer(1))
  }))
  rownames(counts) <- rn
  
  all_equal <- apply(counts, 1, function(r) length(unique(r)) == 1)
  issues <- counts[!all_equal, , drop = FALSE]
  
  if (nrow(issues)) {
    issues <- tibble::as_tibble(issues, rownames = "cell_id") %>%
      dplyr::mutate(clone_id = md[cell_id, "clone_id", drop = TRUE]) %>%
      dplyr::relocate(clone_id, .after = cell_id)
    
    attr(issues, "summary") <- list(
      total_checked = nrow(md),
      problems = nrow(issues),
      ok = nrow(md) - nrow(issues)
    )
    
    if (stop_on_error) {
      stop(sprintf("Consistency check failed for %d rows. Inspect returned tibble.", nrow(issues)), call. = FALSE)
    }
  } else {
    issues <- tibble::tibble()
    attr(issues, "summary") <- list(total_checked = nrow(md), problems = 0, ok = nrow(md))
  }
  
  issues
}
