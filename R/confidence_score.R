#' Calculate Entropy of a Numeric Vector
#'
#' Computes the entropy of a given numeric vector, with considerations for
#' zero proportion and unique values. Entropy is calculated based on the
#' distribution of values across quantile-derived bins.
#'
#' @param vec (`numeric`) A numeric vector for which entropy is to be calculated.
#'
#' @return Returns the calculated entropy of the vector. If the proportion of non-zero
#' values is less than 0.1, or if the vector lacks sufficient unique values, the function
#' returns 0. If an error or warning occurs during calculation, NA is returned with a message.
#'
#' @examples
#' vec <- c(rep(1, 5), rep(2, 5), rep(3, 5))
#' calc_entropy(vec)
#'
#' @export
calc_entropy <- function(vec, quantiles) {
  
  no0propo <- sum(vec!=0)/length(vec)
  if (no0propo < 0.1){
    return(0)
  } else {
    tryCatch({
      if(length(unique(vec[!is.na(vec)])) < 2) {
        return(0) 
      }
      
      quantiles <- quantiles %>% unique()
      
      if(length(quantiles) < 2) { 
        return(0)
      }
      
      bins <- cut(vec, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
      prob_dist <- table(bins) / sum(!is.na(bins))
      entropy <- -sum(prob_dist * log2(prob_dist + 1e-9)) # Adding epsilon to avoid log(0)
      return(entropy)
    }, error = function(e) {
      message("\nError during entropy calculation: ", e$message)
      return(NA) 
    }, warning = function(w) {
      message("\nWarning during entropy calculation: ", w$message)
      return(NA) 
    })
  }
}

#' Calculate Noise Scores for pMHC Data
#'
#' This function calculates noise scores for peptide-MHC (pMHC) data within a Seurat object. 
#' It can operate in two modes: 'per_clone', where calculations are performed for each clone 
#' separately, and 'pseudobulk', where data are aggregated across clones. The function also 
#' supports downsampling to reduce computational load.
#'
#' @param object (`Seurat`) A Seurat object containing pMHC assay data.
#' @param how (`character`) A character string specifying the calculation mode: 'per_clone' or 'pseudobulk'.
#' @param downsample_rate (`numeric`) A numeric value between 0 and 1 indicating the proportion of clones 
#'        to be included in the analysis through downsampling. A value of 1 means no downsampling.
#'
#' @return The input Seurat object with added noise score data in the @misc slot.
#'
#' @details
#' The function computes two metrics for each pMHC feature: the ratio of cells with non-zero 
#' expression (nonzero_ratio) and the entropy of expression levels (entropy). A combined noise 
#' score is also calculated. For 'per_clone' mode, calculations are done for each clone, 
#' and results are averaged. For 'pseudobulk' mode, data from selected clones are aggregated 
#' before calculation. The function updates the input object by storing results in the @misc 
#' slot and returns the modified object.
#'
#' @examples
#' # Assuming 'seurat_obj' is a Seurat object with pMHC assay data
#' seurat_obj <- score_pmhc_noise(seurat_obj, how = 'per_clone', downsample_rate = 0.2)
#'
#' @export
score_pmhc_noise <- function(object, how=c('per_clone', 'pseudobulk'), downsample_rate=.20){
  
  if (downsample_rate < 0 | downsample_rate > 1){
    stop('downsample rate must be between 0 and 1')
  }
  
  pmhc_names <- rownames(object@assays$pMHC)
  object_sub <- subset(object, cells = Cells(object)[!is.na(object$clone_id)])
  
  nonzero <- apply(object@assays$pMHC@counts, 1, function(x) sum(x > 0) / length(Cells(object))) %>% sort()
  object@misc$non0_pmhc <- nonzero
  
  downsampled_clones <- object_sub@meta.data %>% filter(clone_size>1) %>%
    pull(clone_id) %>% unique() %>% sample(., round(length(.) * downsample_rate))
  biggest_clones <- object_sub@meta.data %>% 
    select(clone_id, clone_size) %>%
    unique() %>% slice_max(order_by = clone_size, n = 20) %>% pull(clone_id)
  downsampled_clones <- c(downsampled_clones, biggest_clones) %>% unique()

  entropies <- list()
  
  pmhc_quantiles <- apply(GetAssayData(object, assay = 'pMHC', layer = 'counts'),
                          1, function(vec) quantile(vec, probs = c(.2, .4, .6, .8, 1), 
                                                    na.rm = TRUE, names = FALSE), 
                          simplify = F)
  
  if (how == 'per_clone'){
    # Initializes the progress bar
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = length(downsampled_clones), # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "+")   # Character used to create the bar
    index <- 0
    cat('\ncalculating pMHC entropy on per clone basis\n')
    for (cl in downsampled_clones){
      data_matrix <- GetAssayData(subset(object, clone_id==cl), layer='counts', assay='pMHC')
      rows_list <- split(data_matrix[pmhc_names,], pmhc_names)
      
      entropies <- mapply(calc_entropy, rows_list[pmhc_names], pmhc_quantiles[pmhc_names])
      
      index <- index + 1
      setTxtProgressBar(pb, index)
    }
    close(pb)
    
    entropies_mean <- bind_cols(entropies) %>% rowMeans(na.rm = T)
    names(entropies_mean) <- pmhc_names
  } else if (how == 'pseudobulk') {
    cat('\ncalculating pMHC entropy in pseudobulk mode\n')
    data_matrix <- AverageExpression(subset(object, clone_id %in% downsampled_clones), 
                                     layer='counts', group.by = 'clone_id', assays = 'pMHC')$pMHC
    
    rows_list <- split(data_matrix[pmhc_names,], pmhc_names)
    
    entropies_mean <- mapply(calc_entropy, rows_list[pmhc_names], pmhc_quantiles[pmhc_names])
    names(entropies_mean) <- pmhc_names
  }
  
  object@misc$noise_score <-  data.frame(nonzero_ratio = nonzero[pmhc_names], 
                                         entropy = entropies_mean[pmhc_names],
                                         noise_score = nonzero[pmhc_names] + entropies_mean[pmhc_names])
  
  return(object)
}


#' Calculate pMHC Concordance
#'
#' This function calculates the concordance of pMHC, which for each pMHC is a proportion of gems
#' that have a value for this pMHC highers above all else, the mertic was introtuced in Polvsen et al 2023, elife
#'
#' @param clone_obj (`Seurat`) An Seurat object subset for a single clone. 
#' that includes pMHC assay data.
#' @param preserve_pmhc (`character`) A character vector specifying the names or identifiers of pMHCs to be preserved 
#' for the concordance calculation. If NULL (default), all pMHCs are used.
#' @param slot (`character`) The slot of the clone object from which to retrieve the pMHC data. 
#' Common values include 'counts', 'data', or 'scale.data'. Defaults to 'counts'.
#' @param assay (`character`) The name of the assay to use for extracting the pMHC data. Defaults to 'pMHC'.
#'
#' @return (`Seurat`) A numeric vector representing the concordance for each pMHC. 
#' Concordance is calculated as the proportion of cells (columns) for which a pMHC has 
#' the maximum count value, across all pMHCs considered in the analysis.
#'
#' @examples
#' # Assuming 'clone_obj' is a pre-loaded Seurat object with pMHC data:
#' concordance <- calculate_pmhc_concordance(clone_obj)
#' 
#' # To focus on specific pMHCs:
#' concordance_specific <- calculate_pmhc_concordance(clone_obj, preserve_pmhc = c("pmhc1", "pmhc2"))
#' 
#' @export
calculate_pmhc_concordance <- function(clone_obj, preserve_pmhc=NULL, slot='counts', assay='pMHC'){
  
  pmhc_counts <- GetAssayData(clone_obj, slot = slot, assay = assay)
  
  if (!is.null(preserve_pmhc)){
    pmhc_counts <- pmhc_counts[preserve_pmhc,]
  }
  
  max_vals <- apply(pmhc_counts, 2, max)
  max_counts <- apply(pmhc_counts, 1, function(row) sum(row == max_vals))
  
  concordance <- max_counts/ncol(pmhc_counts)
  
  return(concordance)
}

#' Calculate Confidence Scores for pMHC-TCR Pairs
#'
#' Computes confidence scores for pMHC-TCR pair based on their noise levels,
#' concordance rates, and the size of the clone. The formula used takes into account
#' weighted contributions from each of these factors.
#'
#' @param entropy Numeric vector representing the noise levels for each pMHC-TCR pair.
#' @param concordances Numeric vector representing the concordance rates for each pMHC-TCR pair.
#' @param clone_size A single numeric value representing the size of the clone.
#' @param alpha Weight parameter for the concordance rate contribution to the confidence score.
#' @param beta Weight parameter for the log-transformed clone size contribution to the confidence score.
#' @param gamma Weight parameter for the noise level's contribution to the confidence score.
#'
#' @return A numeric vector of rounded confidence scores for each pMHC-TCR pair, with precision up to two decimal places.
#'
#' @examples
#' entropy <- c(1.2, 0.5, 1.8)
#' concordances <- c(0.9, 0.85, 0.95)
#' clone_size <- 100
#' calculate_confidence(entropy, concordances, clone_size, alpha = 1, beta = 1, gamma = 1)
#'
#' @export
calculate_confidence <- function(entropy, concordances, clone_size, alpha = 1, beta = 1, gamma = 1) {
  if (length(clone_size) > 1) {
    stop("Clone size should be a single digit representing the size of the clone.")
  }
  
  if (length(entropy) != length(concordances)) {
    stop("Length of entropy and concordances must be equal.")
  }
  
  # Calculate confidence score for each pMHC-TCR pair
  confidence_scores <- alpha * concordances + beta * log(clone_size) - gamma * entropy
  
  return(round(confidence_scores, 2))
}

#' Filter pMHC Classifications Based on Confidence Scores
#'
#' This function filters pMHC classifications and confidence scores in a Seurat object's 
#' metadata based on a specified confidence score cutoff. It supports reversible filtering 
#' by creating and using "_unf" backups of the original classifications and confidence scores. 
#' If the function has been run before, it will filter from the "_unf" backups, 
#' allowing for adjustable cutoffs in subsequent analyses.
#'
#' @param object A Seurat object containing pMHC classification and confidence data within its `meta.data`.
#' @param confidence_cutoff A numeric value specifying the confidence score cutoff 
#'        below which pMHC classifications will be considered unreliable and thus filtered out.
#'
#' @return The input Seurat object with `pMHC_classification` and `pMHC_confidence` 
#'         in its `meta.data` filtered based on the specified `confidence_cutoff`. 
#'         The function also updates `object@commands` to indicate that pMHC filtering has been performed.
#'
#' @details
#' The function checks if `pmhc_filter` has been previously executed by inspecting `object@commands`.
#' If not previously executed, original `pMHC_classification` and `pMHC_confidence` data are backed up 
#' with an "_unf" suffix. If `pmhc_filter` was executed before, it starts with these "_unf" backups for filtering.
#' Post-filtering, cells marked as 'Negative' are preserved as such, irrespective of the confidence score.
#'
#' @examples
#' # Assuming 'seurat_obj' is a Seurat object with pMHC classification and confidence information
#' seurat_obj <- filter_pmhc(seurat_obj, confidence_cutoff = 0.5)
#'
#' @export
filter_pmhc <- function(object, confidence_cutoff) {
  
  if (is.null(object@commands$pmhc_filter)) {
    object@meta.data$pMHC_classification_unf <- object@meta.data$pMHC_classification
    object@meta.data$pMHC_confidence_unf <- object@meta.data$pMHC_confidence
  } else {
    object@meta.data$pMHC_classification <- object@meta.data$pMHC_classification_unf
    object@meta.data$pMHC_confidence <- object@meta.data$pMHC_confidence_unf
  }
  
  meta <- object@meta.data 
 
   classifications_list <- strsplit(as.character(meta$pMHC_classification), ":")
  confidence_list <- strsplit(as.character(meta$pMHC_confidence), ":")
  
   for (i in seq_along(classifications_list)) {
  
       confidences <- as.numeric(confidence_list[[i]])
    confidences[is.na(confidences) | confidences <= 0] <- -Inf # Set to -Inf to be filtered out
    
    valid_indices <- which(confidences > confidence_cutoff)
    
    if (length(valid_indices) > 0) {
      meta$pMHC_classification[i] <- paste(classifications_list[[i]][valid_indices], collapse = ":")
      meta$pMHC_confidence[i] <- paste(confidences[valid_indices], collapse = ":")
    } else {
      meta$pMHC_classification[i] <- NA
      meta$pMHC_confidence[i] <- NA
    }
  }
  meta$pMHC_classification[meta$pMHC_classification_unf == 'Negative' &
                               !is.na(meta$pMHC_classification_unf)] <- 'Negative'
  object@meta.data <- meta[Cells(object),]
  
  object@commands$filter_pmhc <- TRUE
  
  return(object)
}

drop.na <- function(vec){
  vec[!is.na(vec)]
}



