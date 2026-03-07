#' Calculate Entropy of a Numeric Vector
#'
#' Computes the entropy of a given numeric vector, with considerations for
#' zero proportion and unique values. Entropy is calculated based on the
#' distribution of values across quantile-derived bins.
#'
#' @param vec (`numeric`) A numeric vector for which entropy is to be calculated.
#' @param quantiles (`numeric`) vector with the quantiles, expected to receive 
#' pMHC values quantiles or 25, 50, 75% 
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
calc_entropy <- function(vec, quantiles, with_weights=FALSE) {
  
  vec <- vec[!is.na(vec)]
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
      
      if (with_weights) {
        num_bins <- length(prob_dist)   # Number of quantile intervals
        weights <- seq(1, 2, length.out = num_bins)  # Decreasing weights
        #weights <- weights / sum(weights)  # Normalize weights to sum to 1
        
        # Apply weights to the probability distribution
        weighted_prob_dist <- prob_dist * weights[as.numeric(names(prob_dist))]
      } else {
        weighted_prob_dist <- prob_dist
      }
      
      # Calculate weighted entropy
      entropy <- -sum(weighted_prob_dist * log2(prob_dist + 1e-9))  # Adding epsilon to avoid log(0)
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
score_pmhc_noise <- function(object, how=c('per_clone', 'pseudobulk'), 
                             probs = c(.2, .4, .6, .8, 1), with_weights=FALSE,
                             downsample_rate=.20, verbose=T, slot='counts'){
  
  if (downsample_rate < 0 | downsample_rate > 1){
    stop('downsample rate must be between 0 and 1')
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
  
  pmhc_names <- rownames(object@assays$pMHC)
  object_sub <- subset(object, cells = Cells(object)[!is.na(object$clone_id)])
  
  pmhc_mat <- GetAssayData(object, layer = slot)
  
  nonzero <- apply(pmhc_mat, 1, function(x) sum(x > 0) / length(Cells(object))) %>% sort()
  object@misc$non0_pmhc <- nonzero
  
  downsampled_clones <- object_sub@meta.data %>% filter(clone_size>1) %>%
    pull(clone_id) %>% unique() %>% sample(., round(length(.) * downsample_rate))
  biggest_clones <- object_sub@meta.data %>% 
    dplyr::select(clone_id, clone_size) %>%
    unique() %>% slice_max(order_by = clone_size, n = 20) %>% pull(clone_id)
  downsampled_clones <- c(downsampled_clones, biggest_clones) %>% unique()

  entropies <- list()
  
  pmhc_quantiles <- apply(GetAssayData(object, assay = 'pMHC', layer = slot),
                          1, function(vec) quantile(vec, probs = probs, 
                                                    na.rm = TRUE, names = FALSE), 
                          simplify = F)
  
  if (how == 'per_clone'){
    if (verbose){
      pb <- txtProgressBar(min = 0,      
                           max = length(downsampled_clones), 
                           style = 3,    
                           width = 50,  
                           char = "+")   
      index <- 0
      cat('\ncalculating pMHC entropy on per clone basis\n')
    }
    for (cl in downsampled_clones){
      data_matrix <- GetAssayData(subset(object, clone_id==cl), layer=slot, assay='pMHC')
      rows_list <- split(data_matrix[pmhc_names,], pmhc_names)
      
      entropies <- mapply(calc_entropy, rows_list[pmhc_names], pmhc_quantiles[pmhc_names],
                          MoreArgs = list(with_weights = with_weights))
      
      if (verbose) {
        index <- index + 1
        setTxtProgressBar(pb, index)
        }
    }
      if (verbose) close(pb)
    
    entropies_mean <- bind_cols(entropies) %>% rowMeans(na.rm = T)
    names(entropies_mean) <- pmhc_names
  } else if (how == 'pseudobulk') {
    cat('\ncalculating pMHC entropy in pseudobulk mode\n')
    data_matrix <- AverageExpression(subset(object, clone_id %in% downsampled_clones), 
                                     layer=slot, group.by = 'clone_id', assays = 'pMHC')$pMHC
    
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
calculate_pmhc_concordance <- function(clone_obj, exclude_top_pmhc = TRUE, 
                                       z_score_c=1, fraction_c=.5, preserve_pmhc=NULL,
                                       slot = 'counts', assay = 'pMHC') {
  
  if (!is.logical(exclude_top_pmhc) || length(exclude_top_pmhc) != 1) {
    stop("The parameter 'exclude_top_pmhc' must be a single Boolean value (TRUE or FALSE).")
  }
  
  if (is(clone_obj, "Seurat")) {
    pmhc_counts <- GetAssayData(clone_obj, layer = slot, assay = assay)
  } else if (is.matrix(clone_obj) || inherits(clone_obj, "Matrix")) {
    pmhc_counts <- clone_obj
  } else {
    stop("clone_obj must be a Seurat object or a matrix-like pMHC matrix.")
  }
  
  if (!is.null(preserve_pmhc)){
    pmhc_counts <- pmhc_counts[preserve_pmhc,]
  }
  
  calc_concordance <- function(data) {
    max_vals <- apply(data, 2, max, na.rm=T)
    max_counts <- apply(data, 1, function(row) sum(row == max_vals, na.rm=T))
    concordance <- max_counts / ncol(data)
    return(concordance)
  }
  
  concordance <- calc_concordance(pmhc_counts)
  
  if (exclude_top_pmhc) {
    npos <- apply(pmhc_counts > z_score_c, 1, function(x) sum(x) / length(x)) %>% sort(decreasing = TRUE) 
    reactions <- names(npos[npos > fraction_c])
    top_index <- names(concordance)[which.max(concordance)]
    
    if (length(reactions) > 1) {
      if (!top_index %in% reactions){
        reactions <- c(top_index, reactions)
      }
      if (reactions[1] != top_index){
        reactions <- c(top_index, reactions[-which(reactions == top_index)])
      } 
      for (i in 2:length(reactions)) {  
        exclude_indices <- match(reactions[1:(i - 1)], rownames(pmhc_counts))
        remaining_counts <- pmhc_counts[-c(exclude_indices), , drop = FALSE]
        adjusted_concordance <- calc_concordance(remaining_counts)
        
        current_reaction_index <- match(reactions[i], rownames(pmhc_counts))
        concordance[current_reaction_index] <- adjusted_concordance[match(reactions[i], rownames(remaining_counts))]
        }
      }
    }
  
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
calculate_confidence <- function(entropy, concordances, clone_size, deltas, replace.na=T,
                                 alpha = 1, beta = 1, gamma = 1, delta = 1) {
  if (length(clone_size) > 1) {
    stop("Clone size should be a single digit representing the size of the clone.")
  }
  
  if (length(entropy) != length(concordances)) {
    stop("Length of entropy and concordances must be equal.")
  }
  
  safe <- function(x, min_value = 0) {
    if (length(x) == 0) return(min_value)
    x[is.na(x)] <- min_value
    x[x < min_value] <- min_value
    return(x)
  }
  
  confidence_scores <- alpha * safe(concordances) +
    beta  * log(safe(clone_size, min_value = 1)) +
    delta * safe(deltas) -
    gamma * safe(entropy)
  
  return(round(confidence_scores, 2))
}

#' Filter pMHC classifications and metrics in a Seurat object (including deltas and scaled_umis)
#'
#' Filters per-GEM pMHC assignments stored in `object@meta.data` using flexible
#' logical conditions. Supports **reversible** filtering via `*_unf` backups,
#' handles multi-assignment columns that are colon-separated, restores original
#' `"Negative"` labels where appropriate, and warns if an overly strict filter
#' discards most assignments.
#'
#' @param object A Seurat object with per-GEM metadata columns:
#'   `pMHC_classification`, `pMHC_confidence`, `pMHC_pvalues`,
#'   `pMHC_wclone_pvalues`, `pMHC_deltas`, `pMHC_scaled_umis`,
#'   and clone-level columns `clone_id`, `clone_size`.
#' @param condition `character(1)`. Row-wise logical expression that decides which
#'   pMHC rows to **keep** (e.g., `"confidence > 0.9 & pvalues < 0.05 & deltas > 2"`).
#'   You can use **bare tokens** without the `pMHC_` prefix; these are mapped to the
#'   corresponding columns automatically:
#'   - `confidence`   → `pMHC_confidence`
#'   - `pvalues`      → `pMHC_pvalues`
#'   - `wclone_pvalues` → `pMHC_wclone_pvalues`
#'   - `classification` → `pMHC_classification`
#'   - `deltas`       → `pMHC_deltas`
#'   - `scaled_umis`  → `pMHC_scaled_umis`
#'   Set to `NULL` to keep all rows (default).
#' @param condition_scope `character(1)` or `NULL`. Optional logical expression that
#'   defines the **scope** within which `condition` is enforced. Within scope, a row
#'   is kept only if `condition` is TRUE; **outside** scope, rows are kept unchanged.
#'   Accepts the same bare tokens as `condition`. Default `NULL` (apply `condition`
#'   to all rows).
#' @param custom_column `character()`. Optional vector of extra column names from
#'   `meta.data` to carry through the long-format filtering and back-join. Missing
#'   names are silently ignored. Default `character()`.
#'
#' @return The input Seurat object with the six pMHC columns updated (collapsed per
#'   `clone_id` using `:`) according to the filter(s). The function:
#'   1) Backs up originals to `*_unf` on first run and restores from `*_unf` on
#'      subsequent runs (reversible workflow).
#'   2) Drops `"Negative"` rows from filtering but, **after** filtering, restores
#'      `"Negative"` to all pMHC fields for any GEM whose original (`*_unf`)
#'      classification was `"Negative"` and did not receive new assignments.
#'   3) Emits a **warning** if more than **80%** of pMHC rows are removed overall,
#'      reporting how many clones lose all assignments.
#'   4) Records `object@commands$pmhc_filter` to mark that filtering has been run.
#'
#' @details
#' **Data layout & splitting.** pMHC-related columns are expected to be
#' colon-separated strings (e.g., multiple assignments per GEM). Internally,
#' the function expands these columns with `tidyr::separate_rows()`, evaluates
#' your `condition` (optionally inside `condition_scope`), filters rows, and
#' then collapses results **per clone** by joining with `:` again.
#'
#' **Parsing & evaluation.** `condition` and `condition_scope` are parsed once
#' with `rlang::parse_expr()` and evaluated in a dplyr data mask. Bare tokens are
#' mapped to full column names (see @param condition). Typical numeric NA handling:
#' - `pMHC_pvalues`, `pMHC_wclone_pvalues`: NAs are treated as `Inf` so comparisons
#'   like `< 0.05` fail safely.
#' - `pMHC_confidence`: NAs or values `<= 0` are treated as `-Inf` so comparisons
#'   like `> 0` fail safely.
#' - `pMHC_deltas`, `pMHC_scaled_umis`: left as numeric with NA; use `!is.na(...)`
#'   in your expressions if required.
#'
#' **Reversibility.** On the first call, original columns are saved as
#' `pMHC_*_unf`. On subsequent calls, filtering always starts from those backups,
#' enabling you to adjust thresholds without cumulative degradation.
#'
#' **Negative preservation.** GEMs that were originally `"Negative"` are kept as
#' `"Negative"` in all six pMHC columns after the join-back if they receive no new
#' assignments through filtering (prevents unintended NAs).
#'
#' **Warnings on over-filtering.** If more than 80% of long-format pMHC rows are
#' dropped by your filters, the function warns and reports the number of clones that
#' lost all assignments—useful for diagnosing overly strict cut-offs.
#'
#' @examples
#' # Keep confident, significant, and strong-signal assignments
#' obj <- filter_pmhc(
#'   object,
#'   condition = "confidence > 0.9 & pvalues < 0.05 & deltas > 2 & scaled_umis > 1"
#' )
#'
#' # Apply only to singletons; outside scope keep rows unchanged
#' obj <- filter_pmhc(
#'   object,
#'   condition = "pvalues < 0.01 & wclone_pvalues < 0.05",
#'   condition_scope = "clone_size == 1"
#' )
#'
#' # Carry extra columns through (if needed)
#' obj <- filter_pmhc(
#'   object,
#'   condition = "confidence > 0.95",
#'   custom_column = c("donor", "batch")
#' )
#'
#' @export
filter_pmhc <- function(object, condition = NULL, condition_scope = NULL, custom_column = character()) {
  library(dplyr)
  library(tidyr)
  library(rlang)
  library(stringr)
  
   if (is.null(object@commands$pmhc_filter)) {
    object@meta.data <- object@meta.data %>%
      mutate(
        pMHC_classification_unf = pMHC_classification,
        pMHC_confidence_unf     = pMHC_confidence,
        pMHC_pvalues_unf        = pMHC_pvalues,
        pMHC_wclone_pvalues_unf = pMHC_wclone_pvalues,
        pMHC_deltas_unf         = pMHC_deltas,
        pMHC_scaled_umis_unf    = pMHC_scaled_umis
      )
  } else {
    object@meta.data <- object@meta.data %>%
      mutate(
        pMHC_classification     = pMHC_classification_unf,
        pMHC_confidence         = pMHC_confidence_unf,
        pMHC_pvalues            = pMHC_pvalues_unf,
        pMHC_wclone_pvalues     = pMHC_wclone_pvalues_unf,
        pMHC_deltas             = pMHC_deltas_unf,
        pMHC_scaled_umis        = pMHC_scaled_umis_unf
      )
  }
  
  meta <- object@meta.data
  
  custom_column <- intersect(custom_column, colnames(meta))
  
  meta_long <- meta %>%
    dplyr::select(
      pMHC_classification, pMHC_confidence, 
      pMHC_pvalues, pMHC_wclone_pvalues, 
      pMHC_deltas, pMHC_scaled_umis,
      clone_id, clone_size,
      tidyselect::all_of(custom_column)
    ) %>%
    filter(!is.na(pMHC_classification) & pMHC_classification != "Negative") %>%
    separate_rows(
      pMHC_classification, pMHC_confidence, pMHC_pvalues, 
      pMHC_wclone_pvalues, pMHC_deltas, pMHC_scaled_umis,
      sep = ":"
    ) %>%
    mutate(
      pMHC_confidence     = suppressWarnings(as.numeric(pMHC_confidence)),
      pMHC_pvalues        = suppressWarnings(as.numeric(pMHC_pvalues)),
      pMHC_wclone_pvalues = suppressWarnings(as.numeric(pMHC_wclone_pvalues)),
      pMHC_deltas         = suppressWarnings(as.numeric(pMHC_deltas)),
      pMHC_scaled_umis    = suppressWarnings(as.numeric(pMHC_scaled_umis))
    ) %>%
    mutate(
      pMHC_confidence     = ifelse(is.na(pMHC_confidence) | pMHC_confidence <= 0, -Inf, pMHC_confidence),
      pMHC_pvalues        = ifelse(is.na(pMHC_pvalues),        Inf,  pMHC_pvalues),
      pMHC_wclone_pvalues = ifelse(is.na(pMHC_wclone_pvalues), Inf,  pMHC_wclone_pvalues)
        ) %>%
    distinct()
  
  replacements <- c(
    "\\bpvalues\\b"         = "pMHC_pvalues",
    "\\bwclone_pvalues\\b"  = "pMHC_wclone_pvalues",
    "\\bconfidence\\b"      = "pMHC_confidence",
    "\\bclassification\\b"  = "pMHC_classification",
    "\\bdeltas\\b"          = "pMHC_deltas",
    "\\bscaled_umis\\b"     = "pMHC_scaled_umis"
  )
  condition_parsed <- if (is.null(condition)) "TRUE" else str_replace_all(condition, replacements)
  condition_scope_parsed <- if (is.null(condition_scope)) "TRUE" else str_replace_all(condition_scope, replacements)
  
  cond_expr  <- rlang::parse_expr(condition_parsed)        
  scope_expr <- rlang::parse_expr(condition_scope_parsed)   
  
  filtered_long <- meta_long %>%
    mutate(
      .scope = !!scope_expr,
      .valid = dplyr::if_else(.scope, dplyr::coalesce(!!cond_expr, FALSE), TRUE)
    ) %>%
    filter(.valid)
  
  # warning if we drop >80% of pMHC classifications 
  baseline_n <- nrow(meta_long)
  after_n    <- nrow(filtered_long)
  if (baseline_n > 0) {
    loss_prop <- (baseline_n - after_n) / baseline_n
    if (loss_prop >= 0.80) {
      per_clone <- filtered_long %>%
        count(clone_id, name = "kept") %>%
        right_join(count(meta_long, clone_id, name = "total"), by = "clone_id") %>%
        mutate(kept = tidyr::replace_na(kept, 0L),
               kept_ratio = kept / total)
      n_all_dropped <- sum(per_clone$kept == 0)
      
      cond_txt  <- if (is.null(condition)) "NULL" else condition
      scope_txt <- if (is.null(condition_scope)) "NULL" else condition_scope
      
      warning(sprintf(
        "filter_pmhc(): %.1f%% of pMHC classifications were removed (kept %d / %d). \nClones with all assignments dropped: %d. \nConsider relaxing cutoffs. \ncondition='%s'; scope='%s'",
        100 * loss_prop, after_n, baseline_n, n_all_dropped, cond_txt, scope_txt
      ))
    }
  }
  
  meta_filtered_collapsed <- filtered_long %>%
    group_by(clone_id) %>%
    summarise(
      pMHC_classification     = paste(pMHC_classification,     collapse = ":"),
      pMHC_confidence         = paste(pMHC_confidence,         collapse = ":"),
      pMHC_pvalues            = paste(pMHC_pvalues,            collapse = ":"),
      pMHC_wclone_pvalues     = paste(pMHC_wclone_pvalues,     collapse = ":"),
      pMHC_deltas             = paste(pMHC_deltas,             collapse = ":"),
      pMHC_scaled_umis        = paste(pMHC_scaled_umis,        collapse = ":"),
      .groups = "drop"
    )
  
  object@meta.data <- object@meta.data %>%
    tibble::rownames_to_column("rownames") %>%
    dplyr::select(
      -pMHC_classification, -pMHC_confidence, -pMHC_pvalues, -pMHC_wclone_pvalues,
      -pMHC_deltas, -pMHC_scaled_umis
    ) %>%
    dplyr::left_join(meta_filtered_collapsed, by = "clone_id") %>%
    dplyr::mutate(
      across(
        .cols = c(pMHC_classification, pMHC_confidence, pMHC_pvalues,
                  pMHC_wclone_pvalues, pMHC_deltas, pMHC_scaled_umis),
        .fns  = ~ ifelse(pMHC_classification_unf == "Negative" & is.na(.x), "Negative", .x)
      )
    ) %>%
    tibble::column_to_rownames("rownames") %>%
    dplyr::arrange(match(rownames(.), SeuratObject::Cells(object)))
  
  object@commands$pmhc_filter <- list(time = Sys.time(),
                                      condition = condition,
                                      condition_scope = condition_scope)
  
  return(object)
}

is_single_na <- function(x) {
  if (length(x) == 1 && is.na(x)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is_negative <- function(x) {
  length(x) == 1 && !is.na(x) && x == 'Negative'
}


recalculate_confidence <- function(object, assay='pMHC', slot='scale.data', ...){
  
  paired_clones <- object@meta.data %>%
    dplyr::select(clone_id, pMHC_classification, pMHC_confidence) %>%
    filter(!is.na(pMHC_classification) & pMHC_classification != 'Negative') %>% 
    distinct()
  
  uclones <- paired_clones$clone_id %>% unique()
  
  clones_pbulk <- ClonePseudobulk(object %>% subset(clone_id %in% uclones), 
                                  assay = assay, clone_col = 'clone_id', 
                                  slot = slot)
  
  n_iter <- ncol(clones_pbulk)
  pmhc_to_bc <- setNames(object@misc$pmhc$pmhc, object@misc$pmhc$Barcode)
  bc_to_pmhc <- setNames(object@misc$pmhc$Barcode, object@misc$pmhc$pmhc)
  
  new_confidences <- list()
  for (i in 1:n_iter){
    
    specificities_i <- paired_clones %>% 
      filter(clone_id == colnames(clones_pbulk)[i]) %>% 
      separate_rows(c(pMHC_classification, pMHC_confidence), sep = ':') %>%
      mutate(pMHC_classification = recode(pMHC_classification, !!!pmhc_to_bc))
    
    
    cl_size <- object$clone_size[object$clone_id == colnames(clones_pbulk)[i]] %>% unique() %>% drop.na()
    scores <- sort(clones_pbulk[,i], decreasing = F)
    scores[scores < 0] <- 0
    diffs <- diff(scores)
    
    scores_interpolated <- normalize_vector(scores[names(diffs)])
    
    names(scores) <- recode(names(scores), !!!pmhc_to_bc)
    names(scores_interpolated) <- recode(names(scores_interpolated), !!!pmhc_to_bc)
    
    outliers <- scores[names(scores) %in% specificities_i$pMHC_classification]
    
    background <- scores_interpolated[!names(scores_interpolated) %in% specificities_i$pMHC_classification]
    positives <- scores_interpolated[names(scores_interpolated) %in% specificities_i$pMHC_classification]
    
    deltas <- positives - mean(background)
    
    entropy <- object@misc$noise_score[recode(names(outliers), !!!bc_to_pmhc),]$entropy
    concordance <- calculate_pmhc_concordance(clone_obj = subset(object, clone_id == colnames(clones_pbulk)[i]), 
                                              assay='pMHC', slot='scale.data')
    
    names(concordance) <- recode(names(concordance), !!!pmhc_to_bc)
    concordance <- concordance[names(outliers)]
    
    confidences <- calculate_confidence(entropy=entropy, concordances=concordance, 
                                        clone_size=cl_size, deltas=deltas, #...)
                                        gamma=0, beta=0)
    
    specificities_i$pMHC_confidence <- confidences
    
    new_confidences[[colnames(clones_pbulk)[i]]] <- specificities_i 
  }
  new_confidences <- bind_rows(new_confidences) %>%
    group_by(clone_id) %>%
    summarise(across(everything(), ~ paste(., collapse = ":"), .names = "{.col}"))
  
  object@meta.data <- object@meta.data %>%
    left_join(new_confidences, by = "clone_id", suffix = c("", "_new")) %>%
    mutate(
      pMHC_confidence = ifelse(!is.na(pMHC_confidence_new), pMHC_confidence_new, pMHC_confidence),
      pMHC_classification = ifelse(!is.na(pMHC_classification_new), pMHC_classification_new, pMHC_classification)
    ) %>%
    dplyr::select(-c(pMHC_confidence_new, pMHC_classification_new))
    
  return(object)
}
