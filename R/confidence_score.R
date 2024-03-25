drop.na <- function(vec){
  vec[!is.na(vec)]
}

normalize_to_range <- function(x) {

  valid_x <- x[!is.na(x)]
  
  x_normalized <- (x - min(valid_x)) / (max(valid_x) - min(valid_x))
  
  return(x_normalized)
}


calc_entropy <- function(vec) {
  
  no0propo <- sum(vec!=0)/length(vec)
  if (no0propo < 0.1){
    return(0)
  } else {
    tryCatch({
      if(length(unique(vec[!is.na(vec)])) < 2) {
        warning("Not enough unique non-NA values for entropy calculation.")
        return(NA) 
      }
      
      q <- quantile(vec, probs = 0:4/4, na.rm = TRUE, names = FALSE) %>% 
        unique()
      
      if(length(q) < 2) { 
        return(0)
      }
      
      bins <- cut(vec, breaks = q, include.lowest = TRUE, labels = FALSE)
      prob_dist <- table(bins) / sum(!is.na(bins))
      entropy <- -sum(prob_dist * log2(prob_dist + 1e-9)) # Adding epsilon to avoid log(0)
      return(entropy)
    }, error = function(e) {
      message("Error during entropy calculation: ", e$message)
      return(NA) 
    }, warning = function(w) {
      message("Warning during entropy calculation: ", w$message)
      return(NA) 
    })
  }
}

score_pmhc_noise <- function(object, how=c('per_clone', 'pseudobulk'), downsample_rate=20){
  
  pmhc_names <- rownames(object@assays$pMHC)
  object_sub <- subset(object, cells = Cells(object)[!is.na(object$clone_id)])
  
  nonzero <- apply(object@assays$pMHC@counts, 1, function(x) sum(x > 0) / length(Cells(object))) %>% sort()
  object@misc$non0_pmhc <- nonzero
  
  downsampled_clones <- object_sub@meta.data %>% filter(clone_size>1) %>%
    pull(clone_id) %>% unique() %>% sample(., length(unique(object_sub$clone_id)) / downsample_rate)
  biggest_clones <- object_sub@meta.data %>% 
    select(clone_id, clone_size) %>%
    unique() %>% slice_max(order_by = clone_size, n = 20) %>% pull(clone_id)
  downsampled_clones <- c(downsampled_clones, biggest_clones) %>% unique()

  entropies <- list()
  
  if (how == 'per_clone'){
    for (cl in downsampled_clones){
      print(cl)
      entropies[[cl]] <- apply(GetAssayData(subset(object, clone_id==cl), 
                                            layer='counts', assay='pMHC'), 
                               1, calc_entropy)
    }
    entropies_mean <- bind_cols(entropies) %>% rowMeans(na.rm = T)
    names(entropies_mean) <- pmhc_names
  } else if (how == 'pseudobulk') {
    entropies_mean <- apply(AverageExpression(subset(object, clone_id %in% downsampled_clones), 
                                              layer='counts', group.by = 'clone_id', assays = 'pMHC')$pMHC,
                            1, calc_entropy)
    
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
#' @param clone_obj An Seurat object subset for a single clone. 
#' that includes pMHC assay data.
#' @param preserve_pmhc A character vector specifying the names or identifiers of pMHCs to be preserved 
#' for the concordance calculation. If NULL (default), all pMHCs are used.
#' @param slot The slot of the clone object from which to retrieve the pMHC data. 
#' Common values include 'counts', 'data', or 'scale.data'. Defaults to 'counts'.
#' @param assay The name of the assay to use for extracting the pMHC data. Defaults to 'pMHC'.
#'
#' @return A numeric vector representing the concordance for each pMHC. 
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

calculate_confidence <- function(noises, concordances, clone_size, alpha = 1, beta = 1, gamma = 1) {
  if (length(clone_size) > 1) {
    stop("Clone size should be a single digit representing the size of the clone.")
  }
  
  if (length(noises) != length(concordances)) {
    stop("Length of noises and concordances must be equal.")
  }
  
  # Calculate confidence score for each pMHC-TCR pair
  confidence_scores <- alpha * concordances + beta * log(clone_size) - gamma * noises
  
  return(round(confidence_scores, 2))
}


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

                   



