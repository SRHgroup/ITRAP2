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
scale2 <- function(mat, sd_fun = sd, replace_zero=F, ...) {
  
  if(is.vector(mat)) {
    mean <- mean(mat, na.rm = T)
    sd <- sd_fun(mat, na.rm = T)
    
    if (replace_zero && sd == 0) {
      sd <- .0001
    }
    return( (mat - mean) / sd )
  }
  else {
    means <- rowMeans(mat, na.rm = T)
    sds <- apply(mat, 1, function(x) sd_fun(x, na.rm = T, ...))
    
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

super_scale <- function(umi_matrix, outlier=3, ...) {
  scaled_umi_matrix <- scale2(umi_matrix, ...)  
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

ScaleDataNoOutliers <- function(object, outlier=3) {
  
  # Make sure the provided object is a Seurat object
  if(!is(object, "Seurat")) {
    stop("The provided object is not a Seurat object.")
  }
  
  umi_matrix <- object@assays$pMHC@counts
  iter0 <- apply(umi_matrix, 1, sd, na.rm = T)
  iter <- 0
  sds_iter <- list()
  sds_stats_iter <- list()
  
  while(any(abs(scale2(umi_matrix, sd_fun = sd)) > outlier, na.rm = TRUE)) {
    iter <- iter + 1
    print(paste0('iteration ', iter))
    umi_matrix <- super_scale(umi_matrix, sd_fun = sd, ...)
    nzeroes <- apply(umi_matrix, 1, function(x) sum(x[x != 0], na.rm = T))
    sds_stats_iter[[paste0('iteration=', iter)]] <- nzeroes
    
    sds_iter[[paste0('iteration=', iter)]] <- apply(umi_matrix, 1, sd, na.rm=T)
  }
  
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
  
  # Modify the Seurat object's scale.data slot
  object@assays$pMHC@scale.data <- as.matrix((object@assays$pMHC@data - rowMeans(umi_matrix, na.rm = T)) / sds)
  
  return(object)
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
smooth_pmhc <- function(object, best_params = NULL, slot='scale.data', assay = 'pMHC', cl_size_thresh = 5) {
  
  # Get the scaled count matrix
  scaled_counts <- GetAssayData(object, assay = assay, slot = slot)
  
  # Get clone assignments for each cell
  clone_assignments <- object$clone_id#[!is.na(object$clone_id)]
  
  # Get unique clone IDs
  clone_ids <- unique(clone_assignments)
  
  if (!"clone_size" %in% names(object@meta.data)) {
    object@meta.data <- object@meta.data %>%
      tibble::rownames_to_column('row_id') %>% 
      group_by(clone_id) %>%
      mutate(clone_size = if_else(is.na(clone_id), NA_integer_, n())) %>%
      ungroup() %>% 
      tibble::column_to_rownames('row_id')
  }
  
  bigger_than5cells <- object$clone_id[object$clone_size > cl_size_thresh] %>% unique()
  clone_ids <- clone_ids[clone_ids %in% bigger_than5cells & !is.na(clone_ids)]
  
  # Get list of unique pmhcs
  pmhcs <- rownames(scaled_counts)
  
  # Initialize the matrix to hold smoothed counts
  smoothed_counts <- matrix(0, nrow = nrow(scaled_counts), ncol = ncol(scaled_counts),
                            dimnames = dimnames(scaled_counts))
  
  # Empty vector to store clones with errors
  clones_with_errors <- c()
  
  # add a list to track if the smoothing went through
  object@commands$smooth_pmhc <- list()
  for (clone_id in clone_ids) {
    print(paste0('smoothing all pmhcs for clone ', clone_id))
    
    clone_cells <- which(clone_assignments == clone_id)
    
    for (pmhc in pmhcs){
      counts <- scaled_counts[pmhc,][clone_cells]
      smoothed_counts[pmhc,][clone_cells] <- counts
      
      error_occurred <- tryCatch({
        # Default parameters
        span_val <- 1
        degree_val <- 1
        family_val <- "gaussian"
        
        if (!is.null(best_params[[pmhc]])) {
          span_val <- best_params[[pmhc]][1] %>% as.numeric()
          degree_val <- best_params[[pmhc]][2] %>% as.numeric()
          family_val <- best_params[[pmhc]][3]
        }
        
        # Fit LOESS model
        loess_fit <- loess(counts ~ seq_along(counts), normalize = F, 
                           span = span_val, 
                           degree = degree_val,
                           family = family_val)
        
        # Predict smoothed counts
        smoothed_counts[pmhc,][clone_cells] <- predict(loess_fit)
        #return(FALSE)
      }, error = function(e) {
        print(paste0('error in fitting loess for pmhc ',pmhc, '|clone ', clone_id))
        #return(TRUE)
      })
      
      #if (error_occurred) {
      #  clones_with_errors <- c(clones_with_errors, clone_id)
      #}
    }
    object@commands$smooth_pmhc[[clone_id]] <- 'done'
  }
  
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
#' @param best_params (`list`) ignore it for now
#' @param slot (`character`) within what clot to perform smoothing, we recommend to
#' smooth scaled by ScaleDataNoOutliers counts
#' @param assay (`character`) name of the slot with pMHC counts, default "pMHC"
#' @param cl_size_thresh (`numeric`) from what clone size a clone will be smoothed
#' @param gap_threshold ?????
#' 
#' @return (`Seurat`) with the ScaleData slot with scaled pMHC counts smoothed 
#'
#' @examples
#' assign_pmhc(object)
#'
#' @export

assign_pmhc <- function(object, cl_size_thresh=5, slot='scale.data', assay='pMHC', gap_threshold = 1){
  
  clones <- object$clone_id[!is.na(object$clone_id)]
  big_clones <- names(table(clones)[table(clones) > cl_size_thresh])
  
  if (length(big_clones) == 0){
    object$pMHC_classification <- 'Negative'
    print('No clones of size greater than 5 in the object, returning Negative in pMHC_classification')
    return(object)
  } else if (is.null(object@misc$pmhc)){
    object$pMHC_classification <- 'Negative'
    print('No pmhc barcode annotation is found in object@misc$pmhc')
    return(object)
  } else {
    clone_bulk <- AverageExpression(object %>% subset(clone_id %in% big_clones), 
                                    assays = assay, group.by = 'clone_id', 
                                    slot = slot, )$pMHC
    
    # Initialize a list to store the results
    clone_outliers <- list()
    
    # For each column (clone) in the matrix
    for (i in 1:ncol(clone_bulk)) {
      # Obtain the sorted scores
      scores <- sort(clone_bulk[,i], decreasing = F)
      scores[scores<0] <- 0
      
      # Calculate the differences between consecutive scores
      diffs <- diff(scores)
      #diffs[1:round(length(diffs)/2)] <- 0
      
      # Find the position of the first gap that is larger than the threshold
      gap_pos <- which(diffs > gap_threshold)[1]
      
      # Check if a gap was found
      if (!is.na(gap_pos)) {
        # The outliers are the scores after the largest gap
        outliers <- scores[(gap_pos+1):length(scores)]
        
        # Add the outlier indices to the list, using the column name as the list name
        clone_outliers[[colnames(clone_bulk)[i]]] <- list(names(outliers))
      } else {
        # Handle the case where no gap is found
        clone_outliers[[colnames(clone_bulk)[i]]] <- list("Negative")
      }
    }
    
    # Print or return the result as needed
    clone_outliers
    
    pmhc_bc <- setNames(object@misc$pmhc$pmhc, object@misc$pmhc$Barcode)
    df <- do.call(rbind, lapply(names(clone_outliers), function(x) {
      data.frame(clone_id = x,
                 pMHC_classification = if (clone_outliers[[x]] == "Negative") "Negative" else paste(
                   clone_outliers[[x]][[1]] %>% recode(!!!pmhc_bc), 
                   collapse=":"))
    }
    )
    )
    df$clone_id <- gsub('-','_', df$clone_id)
    
    object@meta.data <- object@meta.data %>% 
      tibble::rownames_to_column('row_id') %>% 
      {
        if ('pMHC_classification' %in% names(.)) {
          dplyr::select(., -pMHC_classification)
        } else {
          .
        }
      } %>%
      left_join(df, by = 'clone_id') %>% 
      tibble::column_to_rownames('row_id')
    
    return(object)
  }
}
