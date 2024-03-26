#' Centers and scales a rows in numeric matrix
#'
#' Centers and scales a rows in numeric matrix, you have an option to replace 0 sd
#' with 0.001 not to cause a division by zero during the scaling
#'
#' @param mat A numeric matrix (`matrix`).
#' @param sd_sun, function (`function`) that calculates sd, functions are expected to
#' have na.rm parameter to remove NA values
#'
#' @return a mean ceantered and median scaled matrix (`matrix`)
#'
#' @examples
#' scale2(mat)
#' scale2(mat, sd_fun=sd, replace_zero=T)
#'
#' @export
scale2 <- function(mat, sd_fun = sd, replace_zero=F) {
  
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
    
    if (replace_zero){
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
#' @param mat A numeric matrix (`matrix`).
#' @param outlier (`numeric`) from what z-score a value is considered an outlier
#'
#' @return a mean ceatered and median scaled matrix (`matrix`)
#'
#' @examples
#' super_scale(mat)
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
#'
#' @return (`Seurat`) with the ScaleData slot in pMHC assay filled with scaled 
#' pMGC matrix
#'
#' @examples
#' ScaleDataNoOutliers(object)
#'
#' @export

ScaleDataNoOutliers <- function(object, outlier=3, return_params=F) {
  
   if(!is(object, "Seurat")) {
    stop("The provided object is not a Seurat object.")
  }
  
  umi_matrix <- object@assays$pMHC@counts
  iter0 <- apply(umi_matrix, 1, sd, na.rm = TRUE)
  iter <- 0
  sds_iter <- list()
  sds_stats_iter <- list()
  
  initial_outliers <- sum(abs(scale2(umi_matrix, sd_fun = sd)) > outlier, na.rm = TRUE)
  current_outliers <- initial_outliers
  pb <- txtProgressBar(min = 0, max = initial_outliers, style = 3)
  
  while(any(abs(scale2(umi_matrix, sd_fun = sd)) > outlier, na.rm = TRUE)) {
    iter <- iter + 1
    cat(sprintf("\nIteration %d\n", iter))
    umi_matrix <- super_scale(umi_matrix, sd_fun = sd)
    nzeroes <- apply(umi_matrix, 1, function(x) sum(x[x != 0], na.rm = T))
    sds_stats_iter[[paste0('iteration=', iter)]] <- nzeroes
    
    sds_iter[[paste0('iteration=', iter)]] <- apply(umi_matrix, 1, sd, na.rm=T)
    
    # Update for progress bar
    current_outliers <- sum(abs(scale2(umi_matrix, sd_fun = sd)) > outlier, na.rm = TRUE)
    setTxtProgressBar(pb, initial_outliers - current_outliers)
  }
  
  close(pb)
  
  sds_iter[['iteration=0']] <- iter0
  sds_iter <- do.call(cbind, sds_iter)
  sds_iter <- sds_iter[, paste('iteration=', 0:(ncol(sds_iter)-1), sep = '')]
  
  get_last_nonzero <- function(row) {
    nonzero_values <- row[row != 0]
    if(length(nonzero_values) > 0) {
      return(tail(nonzero_values, n = 1))
    } else {
      return(0)
    }
  }
  
  sds <- apply(sds_iter, 1, get_last_nonzero)
  #rownames(sds_iter) <- recode(rownames(sds_iter), !!!bc_pmhc)
  
  object@assays$pMHC@scale.data <- as.matrix((object@assays$pMHC@data - apply(umi_matrix, 1, mean, na.rm = T)) / sds)
  
  object@commands$ScaleDataNoOutliers <- T
  
  if (return_params){
    return(list('object'=object,
                'params'=sds_iter))
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
#' 
#' @return (`Seurat`) with the ScaleData slot with scaled pMHC counts smoothed 
#'
#' @examples
#' smooth_pmhc(object)
#'
#' @export

# Define the function as an S3 method for Seurat class
smooth_pmhc <- function(object, best_params = NULL, slot='scale.data', assay = 'pMHC', cl_size_thresh = 3, 
                        span_val = 1, degree_val = 1, family_val = "gaussian") {
  
  if (is.null(object@commands$ScaleDataNoOutliers)){
    stop('you have to run ScaleDataNoOutluers() before running smooth_pmhc')
  }
  
  scaled_counts <- GetAssayData(object, assay = assay, slot = slot)
  
  clone_assignments <- object$clone_id#[!is.na(object$clone_id)]
  
  clone_ids <- unique(clone_assignments)
  
  if (!"clone_size" %in% names(object@meta.data)) {
    object@meta.data <- object@meta.data %>%
      tibble::rownames_to_column('row_id') %>% 
      group_by(clone_id) %>%
      mutate(clone_size = n()) %>%
      ungroup() %>% 
      tibble::column_to_rownames('row_id')
    
    object$clone_size[is.na(object$clone_id)] <- NA
  }
  
  bigger_thanXcells <- object$clone_id[object$clone_size > cl_size_thresh] %>% unique()
  clone_ids <- clone_ids[clone_ids %in% bigger_thanXcells & !is.na(clone_ids)] %>% unique()
  
  pmhcs <- rownames(scaled_counts)
  
  smoothed_counts <- matrix(0, nrow = nrow(scaled_counts), ncol = ncol(scaled_counts),
                            dimnames = dimnames(scaled_counts))
  
  clones_with_errors <- c()
  
  total_iterations <- length(clone_ids)
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(clone_ids), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "+")   # Character used to create the bar
  
  index <- 0
  object@commands$smooth_pmhc <- list()
  cat(sprintf('\nsmoothing each pMHC within each clone bigger than %d \n', cl_size_thresh))
  for (clone_id in clone_ids) {
    
    clone_cells <- which(clone_assignments == clone_id)
    
    for (pmhc in pmhcs){
      counts <- scaled_counts[pmhc,][clone_cells]
      smoothed_counts[pmhc,][clone_cells] <- counts
      
      error_occurred <- tryCatch({
        # Default parameters
        
        if (!is.null(best_params[[pmhc]])) {
          span_val <- best_params[[pmhc]][1] %>% as.numeric()
          degree_val <- best_params[[pmhc]][2] %>% as.numeric()
          family_val <- best_params[[pmhc]][3]
        }
        
        loess_fit <- loess(counts ~ seq_along(counts), normalize = F, 
                           span = span_val, 
                           degree = degree_val,
                           family = family_val)
        
        smoothed_counts[pmhc,][clone_cells] <- predict(loess_fit)
        object@commands$smooth_pmhc[[clone_id]] <- 'done'
      }, error = function(e) {
        cat(sprintf('\nerror in fitting loess for pmhc %d  in clone %d', pmhc, clone_id))
        object@commands$smooth_pmhc[[clone_id]] <- 'error'
      })
      setTxtProgressBar(pb, index)
      index <- index + 1
    }
  }
  close(pb)
  
  smoothed_counts[is.nan(smoothed_counts)] <- 0
  object@assays[[assay]]@scale.data <- smoothed_counts
  
  return(object)
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
#' @param slot (`character`) within what clot to perform smoothing, we recommend to
#' smooth scaled by ScaleDataNoOutliers counts
#' @param assay (`character`) name of the slot with pMHC counts, default "pMHC"
#' @param entropy_thresh (`numeric`) pMHCs with higher than this threshold will be
#' removed from TCR-pMHC pairing for clones smaller than cl_size_thresh. Meant for
#' extremely noisy pMHCs.
#' @param cl_size_thresh (`numeric`) from what clone size a clone will be smoothed
#' @param delta_threshold what should be a threshold that differentiates 
#' background from specificity signal, recommended from 0.7 to 1
#' 
#' @return (`Seurat`) with the ScaleData slot with scaled pMHC counts smoothed 
#'
#' @examples
#' assign_pmhc(object)
#'
#' @export

assign_pmhc <- function(object, slot='scale.data', assay='pMHC', assign_small_clones=F,
                        cl_size_thresh=3, entropy_thresh=1, delta_threshold = 1){
  
  if (is.null(object@commands$smooth_pmhc)){
    stop('you have to run smooth_pmhc() before running assign_pmhc')
  }
  if (is.null(object@misc$noise_score)){
    stop('you have to run score_pmhc_noise() before running assign_pmhc with assign_small_clones=T')
  }
  if (is.null(object@misc$pmhc)){
    object$pMHC_classification <- 'Negative'
    cat('No pmhc barcode annotation is found in object@misc$pmhc')
    return(object)
  }
  
  clones <- object$clone_id[!is.na(object$clone_id)]
  big_clones <- names(table(clones)[table(clones) > 1])
  
  cat('\naggregating clonotype aggregated pseudobulk pMHC matrix\n')
  clone_bulk <- AverageExpression(object %>% subset(clone_id %in% big_clones), 
                                  assays = assay, group.by = 'clone_id', 
                                  slot = slot, )$pMHC
  
  if (assign_small_clones){
    single_gem_obj <- subset(object, clone_size==1)
    single_gems <- GetAssayData(single_gem_obj, layer='scale.data', assay='pMHC')
    colnames(single_gems) <- single_gem_obj@meta.data$clone_id
    
    rm(single_gem_obj)
    clone_bulk <- cbind(clone_bulk, single_gems)
    
    cat(
      sprintf("\nRemoving features with entropy smaller than %f for clones > than %f gems\n", 
              entropy_thresh, cl_size_thresh))
    high_entropy_pmhc <- rownames(object@misc$noise_score)[object@misc$noise_score$entropy > entropy_thresh]
    small_clones <- colnames(clone_bulk)[colnames(clone_bulk) %in% object$clone_id[object$clone_size < cl_size_thresh]]
    
    clone_bulk[,small_clones][high_entropy_pmhc,] <- -10
  }
  
  sds <- apply(clone_bulk, 2, sd) 
  clone_bulk <- clone_bulk[,!colnames(clone_bulk) %in% names(sds[sds < .2])] 
  
  clone_outliers <- list()
  
  n_iter <- ncol(clone_bulk)
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "+")   # Character used to create the bar
  
  cat('\nassigning pMHC-TCR pairs\n')
  for (i in 1:n_iter) {
    
    cl_size <- object$clone_size[object$clone_id == colnames(clone_bulk)[i]] %>% unique() %>% drop.na()
  
    scores <- sort(clone_bulk[,i], decreasing = F)
    scores[scores<0] <- 0

    diffs <- diff(scores)
    #diffs[1:round(length(diffs)/2)] <- 0
    
    gap_pos <- which(diffs > delta_threshold)[1]
    
    if (!is.na(gap_pos)) {
      outliers <- scores[(gap_pos+1):length(scores)]
      
      entropy <- object@misc$noise_score[names(outliers),]$entropy
      concordance <- calculate_pmhc_concordance(
        clone_obj = subset(object, clone_id == colnames(clone_bulk)[i]),
        assay='pMHC', slot='scale.data')[names(outliers)] 
      
      confidence_scores <- calculate_confidence(entropy = entropy, 
                                                concordances = concordance, 
                                                clone_size = cl_size)
      
      clone_outliers[[colnames(clone_bulk)[i]]] <- list(confidence_scores)
    } else {
      # Handle the case where no gap is found
      clone_outliers[[colnames(clone_bulk)[i]]] <- list(c("Negative"=NA))
    }
    setTxtProgressBar(pb, i)
  }
  
  cat('\nmerging TCR-pMHC information into objects @meta.data\n')
  pmhc_bc <- setNames(object@misc$pmhc$pmhc, object@misc$pmhc$Barcode)
  df <- do.call(rbind, lapply(names(clone_outliers), function(x) {
    data.frame(clone_id = x,
               pMHC_classification = if (clone_outliers[[x]] == "Negative") "Negative" else paste(
                 names(clone_outliers[[x]][[1]]) %>% recode(!!!pmhc_bc), 
                 collapse=":"),
               pMHC_confidence = if (clone_outliers[[x]] == "Negative") NA else paste(
                 clone_outliers[[x]][[1]], 
                 collapse=":"))
  }
  )
  )
  df$clone_id <- gsub('-','_', df$clone_id)
  
  object@meta.data <- object@meta.data %>% 
    tibble::rownames_to_column('row_id') %>% 
    select(-matches('^(pMHC_classification|pMHC_confidence)$')) %>%
    left_join(df, by = 'clone_id') %>% 
    tibble::column_to_rownames('row_id')
  
  return(object)
}








