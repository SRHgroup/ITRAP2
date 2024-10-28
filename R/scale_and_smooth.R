#' Centers and scales a rows in numeric matrix
#'
#' Centers and scales a rows in numeric matrix, you have an option to replace 0 sd
#' with 0.001 not to cause a division by zero during the scaling
#'
#' @param mat A numeric matrix (`matrix`).
#' @param sd_fun, function (`function`) that calculates sd, functions are expected to
#' have na.rm parameter to remove NA values
#'
#' @return a mean ceantered and median scaled matrix (`matrix`)
#'
#' @examples
#' scale2(mat)
#' scale2(mat, sd_fun=sd)
#'
#' @export
scale2 <- function(mat, sd_fun = sd) {
  
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
ScaleDataNoOutliers <- function(object, outlier = 3, return_params = FALSE, verbose = TRUE) {
  
  if (!is(object, "Seurat")) {
    stop("The provided object is not a Seurat object.")
  }
  
  umi_matrix <- object@assays$pMHC@counts
  iter0 <- apply(umi_matrix, 1, sd, na.rm = TRUE)
  mean_iter0 <- apply(umi_matrix, 1, mean, na.rm = TRUE)
  iter <- 0
  sds_iter <- list()
  sds_stats_iter <- list()
  means_iter <- list()
  
  initial_outliers <- sum(abs(scale2(umi_matrix, sd_fun = sd)) > outlier, na.rm = TRUE)
  current_outliers <- initial_outliers
  
  # Initialize progress bar only if verbose is TRUE
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
    
    # Update for progress bar only if verbose is TRUE
    if (verbose) {
      current_outliers <- sum(abs(scale2(umi_matrix, sd_fun = sd)) > outlier, na.rm = TRUE)
      setTxtProgressBar(pb, initial_outliers - current_outliers)
    }
  }
  
  # Close the progress bar only if verbose is TRUE
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
    object@assays$pMHC@counts - apply(umi_matrix, 1, mean, na.rm = TRUE)
                        )

  object@assays$pMHC@data <- centered / sds
  
  object@commands$ScaleDataNoOutliers <- TRUE
  
  if (return_params) {
    return(list('object' = object, 'sds' = sds_iter, 'means' = means_iter))
  } else {
    return(object)
  }
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
#' @param slot (`character`) within what clot to perform smoothing, we recommend to
#' smooth scaled by ScaleDataNoOutliers counts
#' @param assay (`character`) name of the slot with pMHC counts, default "pMHC"
#' @param cl_size_thresh (`numeric`) from what clone size a clone will be smoothed
#' @param span_val (`numeric`) value for span parameter in loess
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

# Define the function as an S3 method for Seurat class
smooth_pmhc <- function(object, best_params = NULL, assay = 'pMHC', 
                        cl_size_thresh = 3, normalise = FALSE, span_val = 1, 
                        degree_val = 1, family_val = "gaussian", verbose = TRUE) {
  
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
  
  bigger_thanXclones <- object$clone_id[object$clone_size > cl_size_thresh & !is.na(object$clone_size)] %>% unique()
  
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
  
  # Initialize the progress bar if verbose is TRUE
  if (verbose) {
    total_iterations <- length(bigger_thanXclones)
    pb <- txtProgressBar(min = 0, max = total_iterations, style = 3, width = 50, char = "+")
    cat(sprintf('\nsmoothing each pMHC within each clone bigger than %d \n', cl_size_thresh))
  }
  
  index <- 0
  object@commands$smooth_pmhc <- list()
  
  for (clone_id in bigger_thanXclones) {
    clone_cells <- rownames(object@meta.data[smoothable_cells,])[which(clone_assignments == clone_id)]
    
    for (pmhc in pmhcs) {
      counts <- scaled_counts[pmhc,][clone_cells]
      raw_counts <- raw_countss[pmhc,][clone_cells]
      smoothed_counts[pmhc,][clone_cells] <- counts
      
      error_occurred <- tryCatch({
        # Apply custom parameters if available
        if (!is.null(best_params[[pmhc]])) {
          span_val <- as.numeric(best_params[[pmhc]][1])
          degree_val <- as.numeric(best_params[[pmhc]][2])
          family_val <- best_params[[pmhc]][3]
        }
        
        ######
        t <- try({
          loess_fit <- loess(counts ~ seq_along(counts), normalize = normalise, 
                             span = span_val, degree = degree_val, family = family_val)
          smoothed_counts[pmhc,][clone_cells] <- predict(loess_fit)
          object@commands$smooth_pmhc[[clone_id]] <- 'done'
        }, silent = TRUE)
        
        # If there is an error, retry with the "gaussian" family
        if("try-error" %in% class(t)) {
          if (family_val == "symmetric") {
            t_gaussian <- try({
              loess_fit <- loess(counts ~ seq_along(counts), normalize = normalise, 
                                 span = span_val, degree = degree_val, family = "gaussian")
              smoothed_counts[pmhc,][clone_cells] <- predict(loess_fit)
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
        if (verbose) {
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

# store scaled in data and use it here

#' Perform permutation test for a specific pMHC-TCR pair
#' 
#' Perform permutation test for a specific pMHC-TCR pair
#' @param sample_permutation whether to downsample values within the same pMHC, 
#' ('within_pmhc') or across any random pMHC ('random')
#' @return (`numeric`)  
#'
#' @examples
#' tcr_pmhc_permutation.test()
#'
#' @export

tcr_pmhc_permutation.test <- function(pmhc_mat, barcodes_to_test, 
                                      clone_coords, stab_z=1, 
                                      sample_permutation='within_pmhc',
                                      n_permutations=1000, p_adj_method='BH') {
  
  cl_size <- length(clone_coords)
  
  observed_means <- rowMeans(pmhc_mat[, clone_coords, drop=FALSE])[barcodes_to_test]
  observed_stabs <- apply(pmhc_mat, 1, function(x) mean(x > stab_z))[barcodes_to_test]
  
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
      simulated_means[, i] <- rowMeans(pmhc_mat[, sampled_columns, drop=FALSE])
      simulated_stabs[, i] <- apply(pmhc_mat[, sampled_columns, drop=FALSE], 1, function(x) mean(x > stab_z))
    }
  
  means_p_values <- sapply(barcodes_to_test, 
                           function(pmhc) mean(simulated_means[pmhc,] >= observed_means[pmhc]))
  stabs_p_values <- sapply(barcodes_to_test, 
                           function(pmhc) mean(simulated_stabs[pmhc,] >= observed_stabs[pmhc]))

  combined_p_values <- mapply(function(p1, p2) {
    chi_sq_statistic <- -2 * (log(p1) + log(p2))
    pchisq(chi_sq_statistic, df = 4, lower.tail = FALSE)
  }, means_p_values, stabs_p_values) 
  
  return(
    p.adjust(combined_p_values, method = p_adj_method) %>% round(., digits = 3)
    )
}


#' Perform pMHC assigment /each big clone
#'
#' Perform pMHC assigment each big clone, your Seurat object is expected to
#' have a certain structure
#' pMHC count matrix should be stored as separate assay called "pMHC"
#' slot "misc" is expected to have an annotation for your pMHC barcodes with following 
#' columns: "Barcode" - barcode ID (same as rownames in pMHC matrix) 
#'          "Sequence" - peptide sequence
#'          "HLA" - HLA allele
#'
#' @param object (`Seurat`) A Seurat object that must be of the above described 
#' structure 
#' @param assay (`character`) name of the slot with pMHC counts, default "pMHC"
#' @param entropy_thresh (`numeric`) pMHCs with higher than this threshold will be
#' removed from TCR-pMHC pairing for clones smaller than cl_size_thresh. Meant for
#' extremely noisy pMHCs.
#' @param cl_size_thresh (`numeric`) from what clone size a clone will be smoothed
#' @param delta_threshold what should be a threshold that differentiates 
#' background from specificity signal, recommended from 0.7 to 1
#' @param assign_small_clones (`logical`) whether you want to assign pMHC to small
#' clones (<cl_size_thresh), smaller clones result in less confidence, so if 
#' you aim for high precision you may want to skip them, this will also make the 
#' assignment run way faster
#' @param sample_permutation whether to downsample values within the same pMHC, 
#' ('within_pmhc') or across any random pMHC ('random')
#' @param verbose whether to show a progress bar
#' 
#' @return (`Seurat`) with the ScaleData slot with scaled pMHC counts smoothed 
#'
#' @examples
#' assign_pmhc(object)
#'
#' @export

assign_pmhc <- function(object, slot='scale.data', assay='pMHC', 
                        assignment='rosner', params = 'iteratively', remove_for_params = 3, assign_small_clones=F, 
                        trim_outliers_for_zassignment=NULL, assigment_z_score_threshold=2,
                        cl_size_thresh=3, entropy_thresh=1, kmax=10, rosner_alpha=0.05, extreme_alpha=0.05, delta_threshold = 1, 
                        n_tests=10, p_alpha=0.1, sample_permutation='within_pmhc', predefined_extreme_params=NULL,
                        filter_by_zscore=FALSE, filter_by_zscore_big_clones=FALSE, top_down=F, 
                        z_score_threshold=2, calculate_pvalue=TRUE, calculate_confidence_score=TRUE,
                        alpha=1, beta=1, gamma=1, delta=1, force=F, verbose=TRUE, ...){
  
  if (cl_size_thresh<2){
    stop('cl_size_thresh must be minimum 2')
  }
  
  if (is.null(object@commands$smooth_pmhc) & !force){
    stop('you have to run smooth_pmhc() before running assign_pmhc')
  }
  if (is.null(object@misc$noise_score) & assign_small_clones & !force){
    stop('you have to run score_pmhc_noise() before running assign_pmhc with assign_small_clones=T')
  }
  if (is.null(object@misc$pmhc)){
    object$pMHC_classification <- 'Negative'
    if (verbose) cat('No pmhc barcode annotation is found in object@misc$pmhc\n')
    return(object)
  }
  if (!assignment %in% c('z-score', 'permutation', 'delta_gap', 'rosner', 'extreme_distribution')){
    stop('assigment must be either z-score, permutation or delta_gap or extreme_distribution, or check your spelling')
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
  big_clones <- names(table(clones)[table(clones) > cl_size_thresh])
  
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
    single_gem_obj <- subset(object, clone_size == 1)
    single_gems <- GetAssayData(single_gem_obj, layer='scale.data', assay='pMHC')
    colnames(single_gems) <- single_gem_obj@meta.data$clone_id
    rm(single_gem_obj)
    
    if (filter_by_zscore){
      single_gems <- filter_by_z_score(single_gems, 
                                       z_score_threshold = z_score_threshold)
    }
    
    # Combine single-gem data with clone_bulk if it exists
    if (!is.null(clone_bulk)) {
      clone_bulk <- cbind(clone_bulk, single_gems)
    } else {
      clone_bulk <- single_gems
    }
    
    if (verbose) {
      cat(sprintf("\nRemoving features with entropy smaller than %g for clones > than %g gems\n", entropy_thresh, cl_size_thresh))
    }
    high_entropy_pmhc <- rownames(object@misc$noise_score)[object@misc$noise_score$entropy > entropy_thresh]
    small_clones <- colnames(clone_bulk)[colnames(clone_bulk) %in% object$clone_id[object$clone_size < cl_size_thresh]]
    clone_bulk[,small_clones][high_entropy_pmhc,] <- -10
  }
  
  n_iter <- ncol(clone_bulk)
  if (verbose) pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "+")
  if (verbose) cat('\nassigning pMHC-TCR pairs\n')
  
  clone_outliers <- list()
  
  for (i in 1:n_iter) {
    #print(i)
    cl_size <- object$clone_size[object$clone_id == colnames(clone_bulk)[i]] %>% unique() %>% drop.na()
    scores <- sort(clone_bulk[,i], decreasing = F)
    scores[scores < 0] <- 0
    diffs <- diff(scores)
    
    if (top_down){
      gap_pos <- rev(which(diffs > delta_threshold))[1]
    } else {
      gap_pos <- which(diffs > delta_threshold)[1]
    }
    
    scores_interpolated <- normalize_vector(scores[names(diffs)])
    
    if (!is.na(gap_pos) | assignment == 'permutation') {
      
      clone_coords <- which(object$clone_id == colnames(clone_bulk)[i]) 
      
      if (assignment == 'delta_gap'){
        outliers <- scores[(gap_pos+1):length(scores)]
        
        background <- scores_interpolated[1:(gap_pos-1)]
        positives <- scores_interpolated[names(diffs)][(gap_pos):length(diffs)]
        
        deltas <- positives - mean(background)
      } else if (assignment == 'permutation'){
        outliers <- diffs[(length(diffs)-n_tests):length(diffs)]
      } else if (assignment == 'z-score'){
        
        if (is.numeric(trim_outliers_for_zassignment)){
          zscores <- sort((scores-mean(scores))/sd(scores))
          mean_score <- mean(scores[abs(zscores) < 3])
          sd_score <- sd(scores[abs(zscores) < 3])
        } else if (is.character(trim_outliers_for_zassignment)){
          trim_outliers_for_zassignment <- as.numeric(trim_outliers_for_zassignment)
          
          trimmed_zscores <- zscores[!(zscores %in% 
                              sort(zscores, decreasing = TRUE)[1:trim_outliers_for_zassignment])]
          
          mean_score <- mean(trimmed_zscores)
          sd_score <- sd(trimmed_zscores)
          
        } else {
          mean_score <- mean(scores)
          sd_score <- sd(scores)
        }
        
        zscores <- sort((scores-mean_score)/sd_score)
        
        outliers <- scores[zscores > assigment_z_score_threshold]
        
        background <- scores_interpolated[!names(scores_interpolated) %in% names(outliers)]
        positives <- scores_interpolated[names(scores_interpolated) %in% names(outliers)]
        
        deltas <- positives - mean(background)
      } else if (assignment == 'rosner'){
        
        if (length(scores[scores != 0]) < 25){
          scores <- scores[order(scores)[(length(scores)-26):length(scores)]]
        } else {
          scores <- scores[scores != 0]
        }
        
        rosner <- rosnerTest(scores, k=kmax, alpha = rosner_alpha, 
                             params = params, remove_for_params = remove_for_params)
        
        if (all(!rosner$all.stats$Outlier)){
          clone_outliers[[colnames(clone_bulk)[i]]] <- list('confidence_scores'='Negative', 
                                                            'p_values'='Negative')
          next
        }
        
        gap_pos <- rosner$all.stats %>% filter(Outlier) %>% slice_min(order_by = Value, n = 1) %>% pull(Value)
        
        outliers <- scores[scores >= gap_pos]
        
        background <- scores_interpolated[!names(scores_interpolated) %in% names(outliers)]
        positives <- scores_interpolated[names(scores_interpolated) %in% names(outliers)]
        
        deltas <- positives - mean(background)
      } else if (assignment == 'extreme_distribution'){
        ##################
        ##################
        
        if (is.null(predefined_extreme_params)){
          extr <- extRemes::fevd(scores)
          
          location <- extr$results$par["location"]
          scale <- extr$results$par["scale"]
          shape <- extr$results$par["shape"]
        } else {
          location <- predefined_extreme_params["location"]
          scale <- predefined_extreme_params["scale"]
          shape <- predefined_extreme_params["shape"]
        }
        
        
        pvalues <- 1-extRemes::pevd(q = scores, loc = location, scale = scale, shape = shape)
        
        if (sum((pvalues)<extreme_alpha)==0){
          clone_outliers[[colnames(clone_bulk)[i]]] <- list('confidence_scores'='Negative', 
                                                            'p_values'='Negative')
          next
        }
        
        outliers <- scores[which(pvalue<extreme_alpha)]
        
        background <- scores_interpolated[!names(scores_interpolated) %in% names(outliers)]
        positives <- scores_interpolated[names(scores_interpolated) %in% names(outliers)]
        
        deltas <- positives - mean(background)
        
        ##################
        ##################
      }
    
      if (sample_permutation == 'within_pmhc'){
        pmhc_mat <- GetAssayData(object, assay = 'pMHC', layer = 'data')[names(outliers),]
        
        if (is.null(dim(pmhc_mat))) {
          pmhc_mat <- matrix(pmhc_mat, nrow = 1, byrow = TRUE,
                             dimnames = list(names(outliers), names(pmhc_mat)))  
        }
        
        p_values <- tcr_pmhc_permutation.test(pmhc_mat=pmhc_mat, 
                                              barcodes_to_test=names(outliers),
                                              sample_permutation=sample_permutation,
                                              clone_coords=clone_coords)
      } else if (sample_permutation == 'random'){
        pmhc_mat <- GetAssayData(object, assay = 'pMHC', layer = 'data')
        
        if (is.null(dim(pmhc_mat))) {
          pmhc_mat <- matrix(pmhc_mat, nrow = 1, byrow = TRUE,
                             dimnames = list(names(outliers), names(pmhc_mat)))  
        }
        
        p_values <- tcr_pmhc_permutation.test(pmhc_mat=pmhc_mat, 
                                              barcodes_to_test=names(outliers),
                                              sample_permutation=sample_permutation,
                                              clone_coords=clone_coords)
      } else {
        stop('sample_permutation must be either "random" or "within_pmhc"')
      }
      
      
      if (assignment == 'permutation') {
        if (sum(p_values < p_alpha) == 0){
          clone_outliers[[colnames(clone_bulk)[i]]] <- list('confidence_scores'='Negative', 
                                                            'p_values'='Negative')
          next
        } else {
          outliers <- p_values[p_values < p_alpha]
          positives <- scores_interpolated[names(outliers)]
          background <- scores_interpolated[!names(scores_interpolated) %in% names(outliers)]
          deltas <- positives - mean(background) 
          }
        }
      
      entropy <- object@misc$noise_score[names(outliers),]$entropy
      concordance <- calculate_pmhc_concordance(clone_obj = subset(object, clone_id == colnames(clone_bulk)[i]), 
                                                assay='pMHC', slot='scale.data')
      
      concordance <- concordance[names(outliers)]
      
      ########
      ########
      
      confidence_scores <- calculate_confidence(entropy=entropy, concordances=concordance, 
                                                deltas=deltas, clone_size = cl_size, 
                                                alpha = alpha, beta = beta, gamma = gamma, delta = delta)
    
      clone_outliers[[colnames(clone_bulk)[i]]] <- list('confidence_scores'=confidence_scores[names(outliers)], 
                                                        'p_values'=p_values[names(outliers)])
    } else {
      clone_outliers[[colnames(clone_bulk)[i]]] <- list('confidence_scores'='Negative', 
                                                        'p_values'='Negative')
    }
    if (verbose) setTxtProgressBar(pb, i)
  }
  
  if (verbose) close(pb)
  if (verbose) cat('\nmerging TCR-pMHC information into objects @meta.data\n')
  
  pmhc_bc <- setNames(object@misc$pmhc$pmhc, object@misc$pmhc$Barcode)
  df <- do.call(rbind, lapply(names(clone_outliers), function(x) {
    #print(x)
    confidences <- clone_outliers[[x]]$confidence_scores
    p_values <- clone_outliers[[x]]$p_values
    confidence_is_negative <- length(confidences) == 1 && confidences == "Negative"
    pvalue_is_negative <- length(p_values) == 1 && p_values == "Negative"
    
    data.frame(clone_id = x,
               pMHC_classification = if (confidence_is_negative) "Negative" else 
                 paste(names(clone_outliers[[x]]$confidence_scores) %>% recode(!!!pmhc_bc), collapse=":"),
               pMHC_confidence = if (confidence_is_negative) 'Negative' else 
                 paste(clone_outliers[[x]]$confidence_scores, collapse=":"),
               pMHC_pvalues = if (pvalue_is_negative) 'Negative' else 
                 paste(clone_outliers[[x]]$p_values, collapse=":")
               )
  }))
  df$clone_id <- gsub('-', '_', df$clone_id)
  
  object@meta.data <- object@meta.data %>%
    tibble::rownames_to_column('row_id') %>%
    dplyr::select(-matches('^(pMHC_classification|pMHC_confidence|pMHC_pvalues)$')) %>%
    left_join(df, by = 'clone_id') %>%
    tibble::column_to_rownames('row_id')
  
  return(object)
}









