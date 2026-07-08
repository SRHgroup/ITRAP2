#' Create a Clonal pMHC Heatmap
#'
#' Generates a heatmap to visualize pMHC values within a single clone. This function uses the `ComplexHeatmap` package to create 
#' a heatmap without clustering rows or columns, specifically tailored for displaying pMHC data. It is designed to 
#' visually represent the uniqueness or abundance of pMHC values.
#'
#' @param data_matrix A numeric matrix representing the data to be visualized in the heatmap. The rows should represent 
#' pMHCs and the columns should represent individual cells. If the input is not 
#' already a matrix, the function will attempt to convert it to a matrix.
#'
#' @return A ggplot object representing the heatmap, which can be further modified or directly displayed using ggplot 
#' functions and options. The heatmap includes a legend on the left side by default, indicating the scale of UMI counts 
#' or other metrics used in the data matrix.
#'
#' @examples
#' # Assuming 'data_matrix' is a matrix of UMI counts with rows as pMHCs and columns as cells:
#' heatmap <- clonal_pmhc_heatmap(data_matrix)
#' print(heatmap)
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid gpar unit
#' @export
clonal_pmhc_heatmap <- function(data_matrix, fontsize=8, 
                                heatmap_legend_side = "left", 
                                hm_legent_direction = 'vertical') {

   if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }
  
  heatmap_plot <- ComplexHeatmap::Heatmap(
    data_matrix, 
    cluster_rows = FALSE,
    cluster_columns = FALSE, 
    show_row_names = TRUE,
    row_names_side = 'left',
    show_column_names = FALSE, 
    row_names_gp = grid::gpar(fontsize = fontsize), 
    column_names_gp = grid::gpar(fontsize = 10), 
    rect_gp = grid::gpar(col = "grey", lwd = 0.2), 
    heatmap_legend_param = list(
      title = "UMI counts",
      legend_width = grid::unit(2, "cm"), 
      legend_height = grid::unit(4, "cm"), 
      color_bar = "continuous", 
      position=heatmap_legend_side,
      direction=hm_legent_direction,
      legend_main_gp = grid::gpar(fontsize = 10),
      legend_gp = grid::gpar(fontsize = 8) 
    )
  )
  
  heatmap_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(heatmap_plot, heatmap_legend_side = heatmap_legend_side))
  heatmap_gg <- patchwork::wrap_elements(panel = heatmap_grob)
  
  return(heatmap_gg)
}

#' Distribution of pMHC per Clone
#'
#' This function plots a distribution of pMCH with some extra matrix within a specific clone, visualizing the distribution using 
#' a heatmap and jitter plots. The function aims to provide insights into the diversity and abundance of pMHCs associated 
#' with a particular clone, considering only pMHCs that surpass a specified aggregation threshold.
#'
#' @param object A Seurat object containing the clonotype and pMHC information.
#' @param clone The identifier for the clone to be analyzed within the dataset.
#' @param aggr_threshold A numeric value specifying the minimum aggregate count threshold for pMHCs to be included in the analysis.
#' Default is 10.
#' @param slot The assay slot to use for obtaining pMHC data. Common values include 'counts', 'data', or 'scale.data'.
#' Default is 'counts'.
#' @param xlimits An optional numeric vector of length 2 specifying the x-axis limits for the jitter plots. If NULL (default),
#' x-axis limits are determined automatically.
#'
#' @return A patchwork object combining a heatmap of pMHC distribution, a jitter plot of UMI counts by feature, and a bar plot
#' of concordance values for the pMHCs analyzed. This comprehensive visualization aids in understanding the distribution and
#' concordance of pMHCs per clone.
#'
#' @examples
#' # Assuming 'data' is a Seurat object with pMHC data and clone identifiers:
#' clone_analysis_plot <- pmhc_dist_per_clone(data, clone = "clone1")
#' print(clone_analysis_plot)
#'
#' @importFrom dplyr %>%
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot geom_jitter geom_hline theme_classic labs geom_bar
#' @export
pmhc_dist_per_clone <- function(object, clone, aggr_threshold=10, slot = 'counts', 
                                xlimits=NULL, return_pmhcs=F, display_pvalues=F, 
                                conc_fontsize=8, dot_fontsize=8, hm_fontsize=8, rel_height1=0.05,
                                extra_title=NULL, heatmap_legend_side='left', hm_legent_direction = "vertical",  
                                ...){
  
  clone_obj <- subset(object, clone_id == clone)
  
  pmhc_counts <- GetAssayData(clone_obj, layer = 'counts', assay = 'pMHC')
  pmhc_matrix <- GetAssayData(clone_obj, layer = slot, assay = 'pMHC')
  
  pmhc_labels <- object@misc$pmhc$pmhc
  
  names(pmhc_labels) <- object@misc$pmhc$Barcode
  
  rownames(pmhc_counts) <- recode(rownames(pmhc_counts), !!!pmhc_labels)
  rownames(pmhc_matrix) <- recode(rownames(pmhc_matrix), !!!pmhc_labels)
  
  aggr <- apply(pmhc_counts,1,sum)
  preserve <- names(aggr)[aggr > aggr_threshold] 
  
  pmhc_matrix <- pmhc_matrix[preserve,]
  
  pmhc_matrix2 <- pmhc_matrix
  rown <- unname(apply(pmhc_matrix2 > 2, 1, sum)) 
  rown <- paste('npos=', rown, sep = '')
  rownames(pmhc_matrix2) <- rown 
  
  heatmap_gg <- clonal_pmhc_heatmap(pmhc_matrix2, 
                                    fontsize = hm_fontsize,
                                    heatmap_legend_side=heatmap_legend_side,
                                    hm_legent_direction = hm_legent_direction)
  
  pmhc_matrix_long <- as.data.frame(pmhc_matrix)  %>% 
    rownames_to_column(var = "Feature") %>% 
    gather(key = "gem", value = "UMI", -Feature) %>%
    mutate(Feature = factor(Feature, levels=rev(preserve)))
  
  y_vals <- seq(0.5, nrow(pmhc_matrix) + 0.5)

  dists <- ggplot(pmhc_matrix_long, aes(x = UMI, y = Feature, color = Feature)) +
    geom_jitter() +
    geom_hline(yintercept = y_vals, linetype = "dashed", colour = "grey", size = 0.5) +
    theme_classic() +
    theme(axis.text.y = element_text(size = dot_fontsize)) +
    labs(x = "UMI counts", y = "Feature", color = "Feature") +
    NoLegend()
  
  if (!is.null(xlimits)){
    dists <- dists + xlim(xlimits)
  }
  
  concordance <- calculate_pmhc_concordance(clone_obj = clone_obj, slot=slot, 
                                            preserve_pmhc = names(pmhc_labels)[pmhc_labels %in% preserve], ...)
  names(concordance) <- recode(names(concordance), !!!pmhc_labels)
  concordance_df <- data.frame(pmhc = factor(names(concordance),levels=rev(preserve)),
                               concordance = concordance)
  
  conc_bp <- ggplot(concordance_df,
                    aes(y = pmhc, x = concordance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_hline(yintercept = y_vals, linetype = "dashed", colour = "grey", size = 0.5) +
    theme_classic() +
    labs(x = "concordance", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = conc_fontsize)) +
    xlim(0, 1)
  
  combined_plot <- heatmap_gg + (dists + ylab('')) + conc_bp +
    patchwork::plot_layout(ncol = 3, widths = c(3, 6, 3.5))
  
  title_text <- paste0('clonotype = ', clone, '; #gems=', ncol(pmhc_matrix))
  if (!is.null(extra_title)){
    title_text <- paste0(title_text, ', ', extra_title)
  }
  
  if (display_pvalues){
    pvals <- permutation_test_specific_pair(object = object, bc_or_pmhc = preserve,  
                                            use_pmhc = T, clone = clone)
    
    pvals_df <- data.frame(pmhc = factor(names(pvals),levels=rev(preserve)),
                           pvalues = -log10(pvals+0.0001))
    
    pvals_bp <- ggplot(pvals_df,
                      aes(y = pmhc, x = pvalues)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_hline(yintercept = y_vals, linetype = "dashed", colour = "grey", size = 0.5) +
      theme_classic() +
      geom_vline(xintercept = -log10(0.05), color='red') +
      labs(x = "-log10(p values)", y = "") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = conc_fontsize)) +
      xlim(0, 4)
    
    combined_plot <- heatmap_gg + (dists + ylab('')) + conc_bp + pvals_bp +
      patchwork::plot_layout(ncol = 4, widths = c(3, 4.5, 2.5, 2.5))
  }
  
  if (return_pmhcs){
    return(list(
      'plot'=combined_plot,
      'pmhc_vec'=concordance_df
      )
    )
  } else {
    return(combined_plot)
  }
}

#' Generate a Heatmap of pMHC Distribution Across Clones
#'
#' This function visualizes calues of peptide-MHC (pMHC) across different clones 
#' within a given dataset, offering the option to highlight specific pMHC-TCR pairs. It allows for 
#' customization of the heatmap appearance and can subset the data based on patient identifiers or specific 
#' pMHCs. The heatmap can be ordered based on various criteria to facilitate the identification of patterns 
#' or relationships within the data.
#'
#' @param object Seurat object, containing single-cell data along with 
#' pMHC and clone annotations.
#' @param clones A vector of clone identifiers to include in the analysis.
#' @param patient Optional; a patient identifier to filter the data by a specific patient.
#' @param slot The assay slot from which to retrieve the data (default is 'counts').
#' @param clones_order Optional; an order for the clones to be displayed in the heatmap.
#' @param hm_breaks Numeric vector specifying the breakpoints for heatmap coloring.
#' @param hm_palette A vector of colors corresponding to `hm_breaks` for heatmap coloring.
#' @param condpalette A color palette for conditioning variables, if applicable.
#' @param highlight_pmhc.tcr Optional; specific pMHC-TCR interactions to highlight in the heatmap.
#' @param pmhc_palette A color palette for pMHCs, if differentiating by pMHC is desired.
#' @param rowm.fonts Font size for row names.
#' @param column_title_fonts Font size for the column title.
#' @param column_title_rot Rotation angle for the column title.
#' @param use_original_order Logical; whether to use the original order of clones.
#' @param clean_mat Logical; whether to remove 'Negative' classifications from the matrix.
#' @param add_tcr_cluster Logical; whether to add TCR clustering information, if available.
#' @param show_row_names Logical; whether to display row names in the heatmap.
#' @param pmhc_subset Optional; a subset of pMHCs to include in the heatmap.
#' @param ... Additional arguments passed to the heatmap function.
#' @param lwd thickness of the highliting frame
#'
#' @return A ComplexHeatmap object visualizing the specified pMHC distribution across clones. The heatmap 
#' can be customized extensively via function arguments and supports highlighting of specific pMHC-TCR interactions.
#'
#' @examples
#' # Assuming 'data' is a dataset object with required annotations:
#' heatmap <- pmhc_heatmap(data, clones = c("clone1", "clone2"), patient = "patient1")
#' print(heatmap)
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @export
#' 
pmhc_heatmap <- function(object, clones, patient=NULL, slot='counts', assay = 'pMHC',
                         clones_order=NULL, column_split_var=NULL, split_rows=NULL,
                         highlight_pmhc.tcr=F, hm_breaks = NULL, hm_palette=c("blue", "white", "red"),
                         clean_na_features=F, na_col = 'grey', pmhc_palette=NULL, pmhc_order=NULL,
                         highlight_column="pMHC_classification", highlight_color='green',
                         stop_large_highlighting=T, show_heatmap_legend=T, rowm.fonts=8, column_title_fonts = 10, 
                         column_title_rot = 45, clean_mat=F, add_tcr_cluster=F, 
                         show_row_names=T, pmhc_subset=NULL, clean_na_cells = FALSE, 
                         custom_annotations=c(), custom_ann_palette=list(),  
                         show_legend_ann=FALSE, bugged_width=.6, max_cols=16000, 
                         left_ann_vars=NULL, left_ann_palette=NULL, save_to_disc_highlight=F,
                         verbose=T, lwd=2, flip=FALSE, skip_bugged_frames=F, ...) {

  if(length(clones) == 0) {
    stop("list of clonotypes in `clone` argument is empty")
  }
  missing_cols <- setdiff(custom_annotations, colnames(object@meta.data))
  if(length(missing_cols) > 0){
    stop(paste("Missing columns in object's metadata:", paste(missing_cols, collapse=", ")))
  }
  
  if (!highlight_column %in% colnames(object@meta.data)) {
    object[[highlight_column]] <- NA
  }
  
  cust_cols    <- intersect(custom_annotations, colnames(object@meta.data))
  ann_columns  <- c("clone_id", highlight_column, cust_cols)
  arrange_cols <- c("clone_id", cust_cols)
  
  ann_subset <- object@meta.data %>%
    dplyr::select(all_of(ann_columns)) %>%
    filter(clone_id %in% clones) %>%
    arrange(!!!syms(arrange_cols))
  
  ann_subset$clone_id <- factor(ann_subset$clone_id, levels=unique(clones[clones %in% ann_subset$clone_id]))
  
  pat_pmhc <- object@misc$pmhc %>%
    filter(!is.na(Sequence)) %>%
    { if(!is.null(patient)) filter(., grepl(patient, Patient)) else . } %>%
    pull(Barcode)
  
  pat_pmhc <- pat_pmhc[pat_pmhc %in% rownames(object[[assay]]@counts)]
  
  pmhc_subset_ <- GetAssayData(object, layer = slot, assay = assay) %>% as.data.frame()
  pmhc_subset_ <- pmhc_subset_[,Cells(object)[object$clone_id %in% clones]][pat_pmhc,]
  
  if (!all(rownames(ann_subset) %in% colnames(pmhc_subset_))){
    stop('cells from given clonotypes are not in pmhc matrix, check if you supply right set of clones')
  }
  
  bc_to_pmhc <- setNames(object@misc$pmhc$pmhc, object@misc$pmhc$Barcode)
  rownames(pmhc_subset_) <- recode(rownames(pmhc_subset_), !!!bc_to_pmhc)
  
  if (!is.null(pmhc_subset)){
    pmhc_subset_ <- pmhc_subset_[rownames(pmhc_subset_) %in% pmhc_subset,] 
  }
  
  ncl <- ann_subset$clone_id %>% unique() %>% length()
  
  if (ncl == 1){
    clpalette <- setNames('red' ,ann_subset$clone_id %>% unique())
  } else {
    clpalette <- get_random_grid_colors(ncl, seed=3)
    names(clpalette) <- ann_subset$clone_id %>% unique()
  }
  
  color_list <- list(clone_id = clpalette)
  
  cust_cols_present <- intersect(custom_annotations, colnames(ann_subset))
  
  levels_for <- function(colname, x) {
    pal <- custom_ann_palette[[colname]]
    if (!is.null(pal)) {
      lev <- trimws(names(pal))
      if (is.null(lev)) stop(sprintf("custom_ann_palette[['%s']] must be a named vector", colname))
      lev[lev %in% unique(trimws(as.character(x)))]
    } else {
      v <- trimws(as.character(x))
      v[!duplicated(v)]
    }
  }
  
  if (is.null(clones_order)) {
    cells_order <- ann_subset %>%
      mutate(across(all_of(cust_cols_present),
                    ~ factor(trimws(as.character(.x)), levels = levels_for(cur_column(), .x)))) %>%
      arrange(clone_id, across(all_of(cust_cols_present))) %>%   
      rownames_to_column('cell_id') %>%
      pull('cell_id')
  } else {
    cells_order <- ann_subset %>%
      mutate(clone_order_factor = factor(clone_id, levels = clones_order)) %>%
      mutate(across(all_of(cust_cols_present),
                    ~ factor(trimws(as.character(.x)), levels = levels_for(cur_column(), .x)))) %>%
      arrange(clone_order_factor, across(all_of(cust_cols_present))) %>%
      dplyr::select(-clone_order_factor) %>%
      rownames_to_column('cell_id') %>%
      pull('cell_id')
  }
  
  
  pmhc_mat <- as.matrix(pmhc_subset_[, rownames(ann_subset)])
  
  if (clean_mat){
    positive_bc <- ann_subset[[highlight_column]][ann_subset[[highlight_column]] != 'Negative']
    positive_bc <- positive_bc %>% unique() %>% 
      strsplit(., ":", fixed = TRUE) %>% unlist() %>% unique() %>% drop.na()
    pmhc_mat = pmhc_mat[rownames(pmhc_mat) %in% positive_bc,]
  }
  
  if (is.null(hm_breaks)){
    if (slot == 'scale.data'){
      hm_breaks <- c(-5, 0, 5)
    } else {
      hm_breaks <- c(0, 4, 10)
    }
  }
  
  col_fun = circlize::colorRamp2(hm_breaks, hm_palette)
  
  if (ncol(pmhc_subset_) > max_cols) {
    set.seed(123) 
    sampled_cols <- sample(colnames(pmhc_subset_), max_cols)
    pmhc_subset_ <- pmhc_subset_[, sampled_cols]
    ann_subset <- ann_subset[sampled_cols,]
    message("Number of columns exceeds threshold. Data has been randomly sampled to ", max_cols, " columns.")
  }
  
  if (is.null(clones_order)){
    cells_order <- ann_subset %>% 
      arrange(.data[[highlight_column]], clone_id) %>%
      rownames_to_column('cell_id') %>%
      pull('cell_id')
  } else {
    cells_order <- ann_subset %>%
      mutate(clone_order_factor = factor(clone_id, levels = clones_order)) %>%
      arrange(clone_order_factor) %>%
      dplyr::select(-clone_order_factor) %>% 
      rownames_to_column('cell_id') %>%
      pull(cell_id)
  }
  
  ann_subset_ordered <- ann_subset[cells_order,]
  pmhc_mat_ordered <- pmhc_mat[,cells_order]
  
    if (clean_na_cells) {
    keep_cells <- colnames(pmhc_mat_ordered)[colSums(is.na(pmhc_mat_ordered)) == 0]

    if (length(keep_cells) == 0) {
      stop("clean_na_cells=TRUE removed all cells (every column had at least one NA).")
    }

    pmhc_mat_ordered    <- pmhc_mat_ordered[, keep_cells, drop = FALSE]
    ann_subset_ordered  <- ann_subset_ordered[keep_cells, , drop = FALSE]

    # keep internal ordering consistent for anything downstream
    cells_order <- keep_cells
  }

  
  cust_cols_present <- intersect(custom_annotations, colnames(ann_subset_ordered))
  
  for (cn in cust_cols_present) {
    ann_subset_ordered[[cn]] <- trimws(as.character(ann_subset_ordered[[cn]]))
  }
  
  level_map <- list()
  for (cn in cust_cols_present) {
    if (!is.null(custom_ann_palette[[cn]])) {
      lev <- trimws(names(custom_ann_palette[[cn]]))
      if (is.null(lev)) stop(sprintf("custom_ann_palette[['%s']] must be a named vector", cn))
      present <- unique(ann_subset_ordered[[cn]])
      lev <- lev[lev %in% present]
    } else {
      v <- ann_subset_ordered[[cn]]
      lev <- v[!duplicated(v)]
    }
    level_map[[cn]] <- lev
    ann_subset_ordered[[cn]] <- factor(ann_subset_ordered[[cn]], levels = lev)
  }
  
  for (cn in cust_cols_present) {
    lev <- level_map[[cn]]  # ordered!
    pal <- custom_ann_palette[[cn]]
    if (!is.null(pal)) {
      names(pal) <- trimws(names(pal))
      missing <- setdiff(lev, names(pal))
      if (length(missing) > 0) {
        pal <- c(pal, setNames(randomcoloR::distinctColorPalette(length(missing)), missing))
        warning(sprintf("custom_ann_palette[['%s']] lacked colors for: %s (auto-assigned).",
                        cn, paste(missing, collapse = ", ")))
      }
      color_list[[cn]] <- pal[lev]
    } else {
      color_list[[cn]] <- setNames(randomcoloR::distinctColorPalette(length(lev)), lev)
    }
  }
  
  ann_legend_param <- lapply(cust_cols_present, function(cn) {
    list(at = level_map[[cn]], labels = level_map[[cn]])
  })
  names(ann_legend_param) <- cust_cols_present
  
  
  if (!is.null(pmhc_order)){
    pmhc_order <- pmhc_order[pmhc_order %in% rownames(pmhc_mat_ordered)] 
    pmhc_mat_ordered <- pmhc_mat_ordered[pmhc_order,]
  }
  
  custom_ann_list <- setNames(
    lapply(cust_cols_present, function(cn) ann_subset_ordered[[cn]]),
    cust_cols_present
  )
  
  ann_list <- c(
    list(clone_id = ann_subset_ordered$clone_id),
    custom_ann_list
  )
  
  hm22ann <- do.call(HeatmapAnnotation, c(
    ann_list,
    list(
      col = color_list,
      which = "col",
      show_legend = rep(show_legend_ann, length.out = length(ann_list)),
      # give each track the same width to avoid mismatches
      annotation_width = unit(rep(4, length(ann_list)), "mm"),
      gap = unit(1, "mm")
    )
  ))
  
  if (!is.null(column_split_var)){
    split_cols <- ann_subset[cells_order,][[column_split_var]]
  } else {
    split_cols <- NULL
  }
  
  if (!is.null(left_ann_vars)){
    left_ann_df <- object@misc$pmhc
    
    pats <- object$Patient[object$clone_id %in% clones] %>% unique()
    if (length(pats) == 1){
      left_ann_df <- left_ann_df %>%
        tidyr::separate_rows(Patient, sep=':') %>% dplyr::filter(Patient %in% pats) %>% dplyr::distinct()
    }
    
    if (!is.null(pmhc_subset)){
      left_ann_df <- left_ann_df %>% dplyr::filter(pmhc %in% pmhc_subset)
    }
    
    if (!is.null(pmhc_order)){
      left_ann_df <- left_ann_df %>%
        dplyr::mutate(pmhc = factor(pmhc, levels = pmhc_order)) %>%
        dplyr::arrange(pmhc)
    }
    
    left_ann_df <- left_ann_df %>%
      dplyr::select(dplyr::all_of(c(left_ann_vars, 'pmhc'))) %>%
      dplyr::filter(!is.na(pmhc)) %>%
      dplyr::mutate(pmhc = make.unique(as.character(pmhc), sep = "_"))
    
    for (ann_col in left_ann_vars) {
      left_ann_df[[ann_col]] <- trimws(as.character(left_ann_df[[ann_col]]))
      
      if (!is.null(left_ann_palette) && !is.null(left_ann_palette[[ann_col]])) {
        pal_names <- names(left_ann_palette[[ann_col]])
        if (is.null(pal_names)) {
          stop(sprintf("left_ann_palette[['%s']] must be a named vector", ann_col))
        }
        
        present_vals <- sort(unique(na.omit(left_ann_df[[ann_col]])))
        
        missing_in_palette <- setdiff(present_vals, pal_names)
        if (length(missing_in_palette) > 0) {
          warning(sprintf(
            "Values in '%s' not found in its palette and will not be colored as intended: %s",
            ann_col, paste(missing_in_palette, collapse = ", ")
          ))
        }
        
        left_ann_df[[ann_col]] <- factor(left_ann_df[[ann_col]], levels = pal_names)
      } else {
        left_ann_df[[ann_col]] <- factor(left_ann_df[[ann_col]])
      }
    }
    
    left_ann_df <- tibble::column_to_rownames(left_ann_df, 'pmhc')
    
    left_ann_df <- left_ann_df[rownames(pmhc_mat_ordered), , drop = FALSE]
    
    left_ann <- ComplexHeatmap::rowAnnotation(
      df  = left_ann_df,
      col = left_ann_palette
    )
  } else {
    left_ann <- NULL
  }
  
  if (clean_na_features){
    nonna <- rownames( pmhc_mat_ordered)[rowSums(is.na( pmhc_mat_ordered)) == 0]
    pmhc_mat_ordered <-  pmhc_mat_ordered[nonna,]
  }
  
  if (flip) {
    pmhc_mat_ordered <- t(pmhc_mat_ordered)
    
    tmp <- split_rows; split_rows <- split_cols; split_cols <- tmp
    
    cust_cols_present <- intersect(custom_annotations, colnames(ann_subset_ordered))
    row_ann_df <- ann_subset_ordered[, c("clone_id", cust_cols_present), drop = FALSE]
    stopifnot(identical(rownames(row_ann_df), rownames(pmhc_mat_ordered)))
    
    left_ann <- ComplexHeatmap::rowAnnotation(
      df  = row_ann_df,
      col = color_list,
      show_legend = rep(show_legend_ann, length.out = ncol(row_ann_df))  # hide/show per track
    )
    
    if (!is.null(left_ann_vars)) {
      top_ann_df <- left_ann_df[colnames(pmhc_mat_ordered), , drop = FALSE]
      hm22ann <- ComplexHeatmap::HeatmapAnnotation(
        df   = top_ann_df,
        col  = left_ann_palette,
        which = "column",
        show_legend = rep(show_legend_ann, length.out = ncol(top_ann_df))
      )
    } else {
      hm22ann <- NULL
    }
  }
  
  show_row_names2    <- if (flip) FALSE else show_row_names
  show_column_names <- if (flip) TRUE else FALSE
  pmhc_subset_hmap <- ComplexHeatmap::Heatmap(
    pmhc_mat_ordered, name = "pmhc_tcr_hmap", 
    show_heatmap_legend = show_heatmap_legend, 
    row_split = split_rows, column_split = split_cols,
    show_row_names = show_row_names2, 
    show_column_names = show_column_names,
    col = col_fun, na_col = na_col,
    cluster_rows = F, cluster_columns = F,
    row_names_gp = grid::gpar(fontsize = rowm.fonts),  
    column_title_gp = grid::gpar(fontsize = column_title_fonts), 
    column_title_rot = column_title_rot,
    top_annotation=hm22ann, left_annotation = left_ann)
  
  ### BOOM!
  if (!highlight_pmhc.tcr){
    return(pmhc_subset_hmap)
  } else {
    if (length(clones) > 300 & stop_large_highlighting){
      stop('highlighting more than 300 clones specificity usually crushes the session, \n
           if you want to continue, set stop_large_highlighting=F')
    }
    if (!is.null(split_cols)){
      stop('Highlighting specificities in the clone split heatmap is not supported')
    }
    ComplexHeatmap::draw(pmhc_subset_hmap)
    tcr_pmhc <- ann_subset_ordered %>%
      dplyr::select(clone_id, !!sym(highlight_column)) %>%
      filter(!is.na(.data[[highlight_column]]), .data[[highlight_column]] != 'Negative') %>%
      unique() %>%
      separate_rows(!!sym(highlight_column), sep=':')
    
    n_iter <- nrow(tcr_pmhc)
    if (verbose){
      pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "+") 
    }  
    
    nc_all <- ncol(pmhc_mat_ordered)
    nr_all <- nrow(pmhc_mat_ordered)
    for (i in seq_len(n_iter)) {
      clone_i   <- tcr_pmhc[i, ]$clone_id
      antigen_i <- tcr_pmhc[i, ][[highlight_column]]
      
      cell_ids_i <- which(ann_subset_ordered$clone_id == clone_i)
      if (length(cell_ids_i) == 0) {
        if (verbose) setTxtProgressBar(pb, i)
        next
      }
      
      if (!flip) {
        ## ----- ORIGINAL ORIENTATION: rows = antigens, cols = cells -----
        ## columns belonging to X clone
        col_indices <- match(
          rownames(ann_subset_ordered)[cell_ids_i],
          colnames(pmhc_mat_ordered)
        )
        col_indices <- col_indices[!is.na(col_indices)]
        if (length(col_indices) == 0) {
          if (verbose) setTxtProgressBar(pb, i)
          next
        }
        
        min_col <- min(col_indices)
        max_col <- max(col_indices)
        
        ## row for Y antgen
        pmhc_row <- match(antigen_i, rownames(pmhc_mat_ordered))
        if (is.na(pmhc_row)) {
          if (verbose) setTxtProgressBar(pb, i)
          next
        }
        
        decorate_heatmap_body("pmhc_tcr_hmap", {
          nc <- ncol(pmhc_mat_ordered)
          nr <- nrow(pmhc_mat_ordered)
          
          x_center <- ((min_col + max_col) / 2 - 0.5) / nc
          width    <- (max_col - min_col + 1) / nc
          
          y_center <- 1 - ((pmhc_row - 0.5) / nr)
          height   <- 1 / nr
          
          if (skip_bugged_frames && height > bugged_width) return(NULL)
          
          grid.rect(
            x = x_center, y = y_center,
            width = width, height = height,
            just = c("center", "center"),
            gp = grid::gpar(col = highlight_color, lwd = lwd, fill = NA)
          )
        })
      } else {
        row_indices <- match(rownames(ann_subset_ordered)[cell_ids_i], rownames(pmhc_mat_ordered))
        row_indices <- row_indices[!is.na(row_indices)]
        if (length(row_indices) == 0) {
          if (verbose) setTxtProgressBar(pb, i)
          next
        }
        
        min_row <- min(row_indices)
        max_row <- max(row_indices)
        pmhc_col <- match(antigen_i, colnames(pmhc_mat_ordered))
        if (is.na(pmhc_col)) {
          if (verbose) setTxtProgressBar(pb, i)
          next
        }
        
        decorate_heatmap_body("pmhc_tcr_hmap", {
          x      <- (pmhc_col - 1) / ncol(pmhc_mat_ordered)
          width  <- 1 / ncol(pmhc_mat_ordered)
          y_center <- 1 - (((min_row + max_row) / 2 - 0.5) / nrow(pmhc_mat_ordered))
          height   <- (max_row - min_row + 1) / nrow(pmhc_mat_ordered)
          
          if (skip_bugged_frames && height > bugged_width) return(NULL)
          
          grid.rect(
            x = x, y = y_center, width = width, height = height,
            just = c("left", "center"),
            gp = grid::gpar(col = highlight_color, lwd = lwd, fill = NA)
          )
        })
      }
      if (verbose) setTxtProgressBar(pb, i)
    }
    
  }
  if (verbose) close(pb)
  
  if (save_to_disc_highlight){
    return(pmhc_subset_hmap)
  }
   
}

#' Generate a Heatmap of pMHC Distribution Across Clones
#'
#' This function visualizes calues of peptide-MHC (pMHC) across different clones 
#' within a given dataset, offering the option to highlight specific pMHC-TCR pairs. It allows for 
#' customization of the heatmap appearance and can subset the data based on patient identifiers or specific 
#' pMHCs. The heatmap can be ordered based on various criteria to facilitate the identification of patterns 
#' or relationships within the data.
#'
#' @param object Seurat object, containing single-cell data along with 
#' pMHC and clone annotations.
#' @param clones A vector of clone identifiers to include in the analysis.
#' @param patient Optional; a patient identifier to filter the data by a specific patient.
#' @param slot The assay slot from which to retrieve the data (default is 'counts').
#' @param clones_order Optional; an order for the clones to be displayed in the heatmap.
#' @param hm_breaks Numeric vector specifying the breakpoints for heatmap coloring.
#' @param hm_palette A vector of colors corresponding to `hm_breaks` for heatmap coloring.
#' @param condpalette A color palette for conditioning variables, if applicable.
#' @param highlight_pmhc.tcr Optional; specific pMHC-TCR interactions to highlight in the heatmap.
#' @param pmhc_palette A color palette for pMHCs, if differentiating by pMHC is desired.
#' @param rowm.fonts Font size for row names.
#' @param column_title_fonts Font size for the column title.
#' @param column_title_rot Rotation angle for the column title.
#' @param use_original_order Logical; whether to use the original order of clones.
#' @param clean_mat Logical; whether to remove 'Negative' classifications from the matrix.
#' @param add_tcr_cluster Logical; whether to add TCR clustering information, if available.
#' @param show_row_names Logical; whether to display row names in the heatmap.
#' @param pmhc_subset Optional; a subset of pMHCs to include in the heatmap.
#' @param ... Additional arguments passed to the heatmap function.
#' @param lwd thickness of the highliting frame
#'
#' @return A ComplexHeatmap object visualizing the specified pMHC distribution across clones. The heatmap 
#' can be customized extensively via function arguments and supports highlighting of specific pMHC-TCR interactions.
#'
#' @examples
#' # Assuming 'data' is a dataset object with required annotations:
#' heatmap <- pmhc_heatmap(data, clones = c("clone1", "clone2"), patient = "patient1")
#' print(heatmap)
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @export
#' 
pmhc_heatmap_old <- function(object, clones, patient=NULL, slot='counts', assay = 'pMHC', clones_order=NULL, column_split_var=NULL, split_rows=NULL,
                         highlight_pmhc.tcr=F, hm_breaks = NULL, hm_palette=c("blue", "white", "red"), clean_na_features=F, na_col = 'grey',
                         pmhc_palette=NULL, pmhc_order=NULL, condpalette=NULL, highlight_column="pMHC_classification", 
                         stop_large_highlighting=T, show_heatmap_legend=T, rowm.fonts=8, column_title_fonts = 10, 
                         column_title_rot = 45, annotation_colors=list(), clean_mat=F, add_tcr_cluster=F, 
                         show_row_names=T, pmhc_subset=NULL, custom_annotations=c(), skip_bugged_frames=F, 
                         show_legend_ann=c(F, T, T), bugged_width=.6, max_cols=16000, left_ann_vars=NULL, 
                         left_ann_palette=NULL, verbose=T, lwd=2, ...) {

  if(length(clones) == 0) {
    stop("list of clonotypes in `clone` argument is empty")
  }
  missing_cols <- setdiff(custom_annotations, colnames(object@meta.data))
  if(length(missing_cols) > 0){
    stop(paste("Missing columns in object's metadata:", paste(missing_cols, collapse=", ")))
  }
  
  if (!highlight_column %in% colnames(object@meta.data)) {
    object[[highlight_column]] <- NA
  }
  
  ann_columns <- c("clone_id", highlight_column, custom_annotations)
  ann_subset <- object@meta.data %>%
    dplyr::select(all_of(ann_columns)) %>%
    filter(clone_id %in% clones) %>%
    arrange(clone_id)
  
  ann_subset$clone_id <- factor(ann_subset$clone_id, levels=unique(clones[clones %in% ann_subset$clone_id]))
  
  pat_pmhc <- object@misc$pmhc %>%
    filter(!is.na(Sequence)) %>%
    { if(!is.null(patient)) filter(., grepl(patient, Patient)) else . } %>%
    pull(Barcode)
  
  pat_pmhc <- pat_pmhc[pat_pmhc %in% rownames(object@assays$pMHC@counts)]
  
  pmhc_subset_ <- GetAssayData(object, layer = slot, assay = assay) %>% as.data.frame()
  pmhc_subset_ <- pmhc_subset_[,Cells(object)[object$clone_id %in% clones]][pat_pmhc,]
  
  if (!all(rownames(ann_subset) %in% colnames(pmhc_subset_))){
    stop('cells from given clonotypes are not in pmhc matrix, check if you supply right set of clones')
  }
  
  bc_to_pmhc <- setNames(object@misc$pmhc$pmhc, object@misc$pmhc$Barcode)
  rownames(pmhc_subset_) <- recode(rownames(pmhc_subset_), !!!bc_to_pmhc)
  
  if (!is.null(pmhc_subset)){
    pmhc_subset_ <- pmhc_subset_[rownames(pmhc_subset_) %in% pmhc_subset,] 
  }
  
  ncl <- ann_subset$clone_id %>% unique() %>% length()
  
  if (ncl == 1){
    clpalette <- setNames('red' ,ann_subset$clone_id %>% unique())
  } else {
    clpalette <- get_random_grid_colors(ncl, seed=3)
    names(clpalette) <- ann_subset$clone_id %>% unique()
  }
  
  if (is.null(pmhc_palette)){
    palette_list <- list(clone_id = clpalette,
                         condition = condpalette)
  } else{
    palette_list <- list(clone_id = clpalette,
                         condition = condpalette,
                         pmhc = pmhc_palette)
  }
  
  pmhc_mat <- as.matrix(pmhc_subset_[, rownames(ann_subset)])
  
  if (clean_mat){
    positive_bc <- ann_subset[[highlight_column]][ann_subset[[highlight_column]] != 'Negative']
    positive_bc <- positive_bc %>% unique() %>% 
      strsplit(., ":", fixed = TRUE) %>% unlist() %>% unique() %>% drop.na()
    pmhc_mat = pmhc_mat[rownames(pmhc_mat) %in% positive_bc,]
  }
  
  if (is.null(hm_breaks)){
    if (slot == 'scale.data'){
      hm_breaks <- c(-5, 0, 5)
    } else {
      hm_breaks <- c(0, 4, 10)
    }
  }
  
  col_fun = circlize::colorRamp2(hm_breaks, hm_palette)
  
  if (ncol(pmhc_subset_) > max_cols) {
    set.seed(123) 
    sampled_cols <- sample(colnames(pmhc_subset_), max_cols)
    pmhc_subset_ <- pmhc_subset_[, sampled_cols]
    ann_subset <- ann_subset[sampled_cols,]
    message("Number of columns exceeds threshold. Data has been randomly sampled to ", max_cols, " columns.")
  }
  
  if (is.null(clones_order)){
    cells_order <- ann_subset %>% 
      arrange(.data[[highlight_column]], clone_id) %>%
      rownames_to_column('cell_id') %>%
      pull('cell_id')
  } else {
    cells_order <- ann_subset %>%
      mutate(clone_order_factor = factor(clone_id, levels = clones_order)) %>%
      arrange(clone_order_factor) %>%
      dplyr::select(-clone_order_factor) %>% 
      rownames_to_column('cell_id') %>%
      pull(cell_id)
  }
  
  ann_subset_ordered <- ann_subset[cells_order,]
  pmhc_mat_ordered <- pmhc_mat[,cells_order]
  
  if (!is.null(pmhc_order)){
    pmhc_order <- pmhc_order[pmhc_order %in% rownames(pmhc_mat_ordered)] 
    pmhc_mat_ordered <- pmhc_mat_ordered[pmhc_order,]
  }
  
  color_list <- list()
  
  if (!is.null(custom_annotations)){
    for (ann_col in custom_annotations) {
      if (!is.null(annotation_colors[[ann_col]])) {
        color_list[[ann_col]] <- annotation_colors[[ann_col]]
      } else {
        unique_values <- unique(ann_subset_ordered[[ann_col]]) %>% na.omit()
        color_list[[ann_col]] <- setNames(randomcoloR::distinctColorPalette(length(unique_values)), unique_values)
      }
    }
  }
  
  # Default color palettes redundant?
  #ncl <- ann_subset_ordered$clone_id %>% unique() %>% length()
  #clpalette <- get_random_grid_colors(ncl, seed=3)
  #names(clpalette) <- ann_subset_ordered$clone_id %>% unique()
  
  #palette_list <- list(clone_id = clpalette, condition = condpalette, pmhc = pmhc_palette)
  palette_list <- c(palette_list, color_list)  # Merge custom colors
  
  custom_ann_list <- setNames(lapply(custom_annotations, function(cn) ann_subset_ordered[[cn]]), custom_annotations)
  
  ann_list <- list(
    clone_id = ann_subset_ordered$clone_id#, 
    #pmhc = gsub(':', '\n', ann_subset_ordered[[highlight_column]])
  )
  ann_list <- c(ann_list, custom_ann_list)
  
  hm22ann <- do.call(HeatmapAnnotation, c(ann_list, list(
    col = palette_list, show_legend = show_legend_ann,
    which = 'col',
    annotation_width = unit(c(1, 4), 'cm'),
    gap = unit(1, 'mm')
  )))
  
  if (!is.null(column_split_var)){
    split_cols <- ann_subset[cells_order,][[column_split_var]]
  } else {
    split_cols <- NULL
  }
  
  if (!is.null(left_ann_vars)){
    left_ann_df <- object@misc$pmhc
    
    pats <- object$Patient[object$clone_id %in% clones] %>% unique()
    if (length(pats) == 1){
      left_ann_df <- left_ann_df %>%
        tidyr::separate_rows(Patient, sep=':') %>% dplyr::filter(Patient %in% pats) %>% dplyr::distinct()
    }
    
    if (!is.null(pmhc_subset)){
      left_ann_df <- left_ann_df %>% dplyr::filter(pmhc %in% pmhc_subset)
    }
    
    if (!is.null(pmhc_order)){
      left_ann_df <- left_ann_df %>%
        dplyr::mutate(pmhc = factor(pmhc, levels = pmhc_order)) %>%
        dplyr::arrange(pmhc)
    }
    
    left_ann_df <- left_ann_df %>%
      dplyr::select(dplyr::all_of(c(left_ann_vars, 'pmhc'))) %>%
      dplyr::filter(!is.na(pmhc)) %>%
      dplyr::mutate(pmhc = make.unique(as.character(pmhc), sep = "_"))
    
    for (ann_col in left_ann_vars) {
      left_ann_df[[ann_col]] <- trimws(as.character(left_ann_df[[ann_col]]))
      
      if (!is.null(left_ann_palette) && !is.null(left_ann_palette[[ann_col]])) {
        pal_names <- names(left_ann_palette[[ann_col]])
        if (is.null(pal_names)) {
          stop(sprintf("left_ann_palette[['%s']] must be a named vector", ann_col))
        }
        
        present_vals <- sort(unique(na.omit(left_ann_df[[ann_col]])))
        
        missing_in_palette <- setdiff(present_vals, pal_names)
        if (length(missing_in_palette) > 0) {
          warning(sprintf(
            "Values in '%s' not found in its palette and will not be colored as intended: %s",
            ann_col, paste(missing_in_palette, collapse = ", ")
          ))
        }
        
        left_ann_df[[ann_col]] <- factor(left_ann_df[[ann_col]], levels = pal_names)
      } else {
        left_ann_df[[ann_col]] <- factor(left_ann_df[[ann_col]])
      }
    }
    
    left_ann_df <- tibble::column_to_rownames(left_ann_df, 'pmhc')
    
    left_ann_df <- left_ann_df[rownames(pmhc_mat_ordered), , drop = FALSE]
    
    left_ann <- ComplexHeatmap::rowAnnotation(
      df  = left_ann_df,
      col = left_ann_palette
    )
  } else {
    left_ann <- NULL
  }
  
  if (clean_na_features){
    nonna <- rownames( pmhc_mat_ordered)[rowSums(is.na( pmhc_mat_ordered)) == 0]
    pmhc_mat_ordered <-  pmhc_mat_ordered[nonna,]
  }
  
  pmhc_subset_hmap <- ComplexHeatmap::Heatmap(
    pmhc_mat_ordered, name = "pmhc_tcr_hmap", 
    show_heatmap_legend = show_heatmap_legend, 
    row_split = split_rows, column_split = split_cols,
    show_row_names = show_row_names, show_column_names = FALSE,
    col = col_fun, na_col = na_col,
    cluster_rows = F, cluster_columns = F,
    row_names_gp = grid::gpar(fontsize = rowm.fonts),  
    column_title_gp = grid::gpar(fontsize = column_title_fonts), 
    column_title_rot = column_title_rot,
    top_annotation=hm22ann, left_annotation = left_ann)
  
  if (!highlight_pmhc.tcr){
    return(pmhc_subset_hmap)
  } else {
    if (length(clones) > 300 & stop_large_highlighting){
      stop('highlighting more than 300 clones specificity usually crushes the session, \n
           if you want to continue, set stop_large_highlighting=F')
    }
    if (!is.null(split_cols)){
      stop('Highlighting specificities in the clone split heatmap is not supported')
    }
    ComplexHeatmap::draw(pmhc_subset_hmap)
    tcr_pmhc <- ann_subset_ordered %>%
      dplyr::select(clone_id, !!sym(highlight_column)) %>%
      filter(!is.na(.data[[highlight_column]]), .data[[highlight_column]] != 'Negative') %>%
      unique() %>%
      separate_rows(!!sym(highlight_column), sep=':')
    
    n_iter <- nrow(tcr_pmhc)
    if (verbose){
      pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "+") 
    }  
    
    for (i in 1:n_iter) {
      clone_i <- tcr_pmhc[i,]$clone_id
      antigen_i <- tcr_pmhc[i,][[highlight_column]]
      
      cell_ids_i <- which(ann_subset_ordered$clone_id == clone_i)
      if (length(cell_ids_i) == 0) next
      
      col_indices <- match(rownames(ann_subset_ordered)[cell_ids_i], colnames(pmhc_mat_ordered))
      min_col <- min(col_indices)
      max_col <- max(col_indices)
      
      pmhc_row <- which(rownames(pmhc_mat_ordered) == antigen_i)
      if (length(pmhc_row) == 0) next  # Skip if antigen not found
      if (length(pmhc_row) > 1) {
        pmhc_row <- pmhc_row[1]
      }
      
      decorate_heatmap_body("pmhc_tcr_hmap", {
        x = (min_col - 1) / ncol(pmhc_mat_ordered)
        width = (max_col - min_col + 1) / ncol(pmhc_mat_ordered)
        y = 1 - (pmhc_row - 0.5) / nrow(pmhc_mat_ordered)
        height = 1 / nrow(pmhc_mat_ordered)
        
        if (skip_bugged_frames & width > bugged_width) {
          return(NULL)
        }
        
        grid.rect(x = x, y = y, width = width, height = height,
                  just = c("left", "center"), 
                  gp = grid::gpar(col = "green", lwd = lwd, fill = NA))
      })
      
      if (verbose) setTxtProgressBar(pb, i)
    }
  }
  if (verbose) close(pb)
}






plot_smoothing_changes <- function(scaled, smoothed, raw_counts = NULL, 
                                   point_labels = NULL, pmhc, clone_id, title=NULL) {
  if (is.null(raw_counts)) {
    states <- c("Scaled", "Smoothed")
    values <- list(scaled, smoothed)
  } else {
    states <- c("Raw", "Scaled", "Smoothed")
    values <- list(raw_counts, scaled, smoothed)
  }
  
  data <- data.frame(
    state = rep(states, each = length(scaled)),
    value = unlist(values),
    ids = rep(1:length(scaled), times = length(states))  
  )
  
  if (!is.null(point_labels)) {
    data$label <- rep(point_labels, times = length(states))
  }
  
  data_agg <- data %>%
    group_by(state, value) %>%
    summarise(
      count = n(),          # Count of overlapping points
      ids = list(ids),       # List of IDs for color consistency
      .groups = "drop"
    )
  
  data_expanded <- data_agg %>%
    unnest_longer(ids) 
  
  if (is.null(title)){
    title <- sprintf('Normalisation of pMHC UMI counts\n for %s\n in clone %s', pmhc, clone_id)
  } 
  
  p <- ggplot(data, aes(x = state, y = value, group = ids, color = as.factor(ids))) +
    geom_line(show.legend = FALSE) +  # Lines between points
    geom_point(data = data_expanded, aes(size = count), show.legend = FALSE) +  # Points sized by count
    labs(x = "State", y = "Value", title = "Changes Across States") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.title.position = "plot") +
    ggtitle(title)
  
  
  if (!is.null(point_labels)) {
    p <- p + geom_text(data = subset(data, state == "Smoothed"), 
                       aes(label = label, hjust = -0.1), size = 3)
  }
  
  return(p)
}

plot_smoothing <- function(object, clone_id, pmhc, ...){
  clone_cells <- Cells(object)[which(object$clone_id == clone_id)]
  
  barcode <- object@misc$pmhc$Barcode[object@misc$pmhc$pmhc == pmhc]
  
  counts <- GetAssayData(object, assay = 'pMHC', layer = 'counts')[barcode,][clone_cells]
  scaled <- GetAssayData(object, assay = 'pMHC', layer = 'data')[barcode,][clone_cells]
  smoothed <- GetAssayData(object, assay = 'pMHC', layer = 'scale.data')[barcode,][clone_cells]
  
  plot_smoothing_changes(raw_counts = counts, scaled = scaled, 
                         smoothed = smoothed, pmhc=pmhc, clone_id=clone_id, ...)
}



#' Generate a Set of Random Colors from a Grid
#'
#' This function produces a set of visually distinct colors by generating points within a 3D grid in RGB space. 
#' The approach ensures that the colors are evenly distributed and visually distinct, which is particularly 
#' useful for plotting where a large number of categories must be differentiated by color. The color generation 
#' is based on subdividing the RGB cube into smaller cubes and selecting random colors within these subdivisions.
#'
#' @param ncolor The number of distinct colors to generate.
#' @param seed An optional seed for the random number generator to ensure reproducibility.
#'
#' @return A character vector of colors in hexadecimal format, suitable for use in plotting functions.
#'
#' @examples
#' # Generate 10 random, visually distinct colors:
#' colors <- get_random_grid_colors(10)
#' plot(1:10, pch=19, col=colors, cex=2)
#'
#' @export
get_random_grid_colors <- function(ncolor,seed = 100) {

  set.seed(seed)
  ngrid <- ceiling(ncolor^(1/3))
  x <- seq(0,1,length=ngrid+1)[1:ngrid]
  dx <- (x[2] - x[1])/2
  x <- x + dx
  origins <- expand.grid(x,x,x)
  nbox <- nrow(origins) 
  RGB <- vector("numeric",nbox)
  for(i in seq_len(nbox)) {
    rgb <- runif_in_cube(n=1,d=3,O=as.numeric(origins[i,]),r=dx)
    RGB[i] <- rgb(rgb[1,1],rgb[1,2],rgb[1,3])
  }
  index <- sample(seq(1,nbox),ncolor)
  RGB[index]
}

runif_in_cube <- function(n, d, O = rep(0, d), r = 1) {
  if (length(n) != 1 || length(d) != 1 || n < 1 || d < 1) {
    stop("'n' and 'd' must be positive scalar values")
  }
  if (length(O) != d) {
    stop("'O' must have length 'd'")
  }
  if (length(r) == 1) {
    r <- rep(r, d)
  }
  if (length(r) != d || any(r < 0)) {
    stop("'r' must be a non-negative scalar or a vector of length 'd'")
  }

  sims <- matrix(stats::runif(n * d, min = -1, max = 1), nrow = n, ncol = d)
  sims <- sweep(sims, 2, r, `*`)
  sweep(sims, 2, O, `+`)
}


NoXaxis <- function() {
  theme(
    axis.title.x = element_blank(),   
    axis.text.x = element_blank(),    
    axis.ticks.x = element_blank()    
  )
}


plot_vdj_segment_frequency <- function(data, chain = "beta", segment = "v_call", 
                                       x_var = "pMHC_classification", threshold = 0.05,
                                       metric = "tcr_fraction") {
  
  segment_col <- paste0(segment, "_", chain)
  
  if (!(segment_col %in% colnames(data))) {
    stop(paste("Column", segment_col, "not found in the dataset! Check 'chain' and 'segment' inputs."))
  }
  if (!(x_var %in% colnames(data))) {
    stop(paste("Column", x_var, "not found in the dataset! Check 'x_var' input."))
  }
  if (!(metric %in% c("tcr_fraction", "cell_fraction"))) {
    stop("Invalid 'metric' value. Use 'tcr_fraction' for TCR fraction or 'cell_fraction' for cell fraction.")
  }
  
  if (metric == "tcr_fraction") {
    plot_data <- data %>%
      group_by(!!sym(x_var), !!sym(segment_col)) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(!!sym(x_var)) %>%
      mutate(fraction = count / sum(count)) %>%
      mutate(segment_label = ifelse(fraction < threshold, paste0("Less than ", threshold * 100, "%"), !!sym(segment_col))) %>%
      group_by(!!sym(x_var), segment_label) %>%
      summarise(fraction = sum(fraction), .groups = "drop")
  } else if (metric == "cell_fraction") {
    plot_data <- data %>%
      group_by(!!sym(x_var), !!sym(segment_col)) %>%
      summarise(cell_sum = sum(clone_size, na.rm = TRUE), .groups = "drop") %>%
      group_by(!!sym(x_var)) %>%
      mutate(fraction = cell_sum / sum(cell_sum)) %>%
      mutate(segment_label = ifelse(fraction < threshold, paste0("Less than ", threshold * 100, "%"), !!sym(segment_col))) %>%
      group_by(!!sym(x_var), segment_label) %>%
      summarise(fraction = sum(fraction), .groups = "drop")
  }
  
  ggplot(plot_data, aes_string(x = x_var, y = "fraction", fill = "segment_label")) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = c("grey", setNames(rainbow(length(unique(plot_data$segment_label)) - 1), 
                                  unique(plot_data$segment_label)[unique(plot_data$segment_label) != paste0("Less than ", threshold * 100, "%")]))
    ) +
    labs(
      title = paste0("Fraction of ", toupper(segment %>% gsub('_call', '', .)), " Segment (", toupper(chain), " Chain) by ", x_var),
      x = x_var,
      y = ifelse(metric == "tcr_fraction", "Fraction of TCRs", "Fraction of Cells"),
      fill = paste0(toupper(segment), " Segment")
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}




#' Plot TCR-pMHC Permutation Tests
#'
#' This function creates histograms for permutation tests and highlights, within simulated values per clone, where the observed clone value and stability lie.
#'
#' @param object A Seurat object containing TCR and pMHC data.
#' @param clone The ID of the clonotype to analyze.
#' @param assay The assay to use for retrieving data. Default is 'pMHC'.
#' @param slot The slot to use for retrieving data. Default is 'data'.
#' @param stab_z The stability threshold for determining non-zero UMI counts. Default is 1.
#' @param n_permutations The number of permutations to perform for the test. Default is 1000.
#'
#' @return A combined ggplot object showing histograms of simulated means and stabilities for each pMHC, with observed values highlighted.
#'
#' @examples
#' # Assuming `seurat_obj` is a Seurat object with relevant data
#' plot_tcr_pmhc_permutation(seurat_obj, clone = 'clone_1')
#'
#' @export
plot_tcr_pmhc_permutation <- function(object, clone, assay = 'pMHC', 
                                      slot='data', stab_z=1, n_permutations=1000){

  clone_coords <- which(object$clone_id == clone)
  clone_size <- object$clone_size[object$clone_id == clone] %>% drop.na() %>% unique() 
  
  meta <- object@meta.data %>%
    filter(clone_id == clone) %>% 
    dplyr::select(pMHC_classification, pMHC_confidence, pMHC_pvalues) %>%
    distinct() %>%
    separate_rows(c(pMHC_classification, pMHC_confidence, pMHC_pvalues), sep = ':')
  
  if (sum(meta$pMHC_classification != 'Negative' | !is.na(meta$pMHC_classification)) == 0){
    stop(paste0('clone ', clone, ' doesnt have an assign pmhc specificity'))
  }
  
  pmhcs <- meta %>% pull(pMHC_classification)
  bcs <- object@misc$pmhc %>% filter(pmhc %in% pmhcs) %>% pull(Barcode)
  names(pmhcs) <- bcs
  p_values <- meta %>% pull(pMHC_pvalues)
  names(p_values) <- bcs
  
  pmhc_mat <- GetAssayData(object, assay = assay, layer=slot)[bcs,]
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
  
  for (bc in rownames(simulated_means)) {
    
    current_means <- data.frame(value = simulated_means[bc, ])
    current_obs_mean <- observed_means[bc]
    current_stabs <- data.frame(value = simulated_stabs[bc, ])
    current_obs_stab <- observed_stabs[bc]
    
    means_plot <- ggplot(current_means, aes(x = value)) + 
      geom_histogram(binwidth = 0.1, color = "black", fill = "grey") +
      geom_vline(xintercept = current_obs_mean, linetype = "dashed", color = "red", size = 2) +
      labs(title = paste0('Simulated means\n', 
                          'pmhc=', pmhcs[bc],
                          '\nclone_id=', clone)) +
      theme_minimal()
    
    stabs_plot <- ggplot(current_stabs, aes(x = value)) + 
      geom_histogram(binwidth = 0.01, color = "black", fill = "blue") +
      geom_vline(xintercept = current_obs_stab, linetype = "dashed", color = "green", size = 2) +
      labs(title =  paste0('Simulated pos GEMs ratio',
                           ',\n#GEMS=', clone_size,
                           ', p_value = ', p_values[bc])) +
      theme_minimal()
    
    combined_plot <- means_plot + stabs_plot + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
    plot_list[[bc]] <- combined_plot
    
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
