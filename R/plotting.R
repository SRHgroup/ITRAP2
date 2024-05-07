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
#' @importFrom cowplot ggdraw draw_plot
#' @export
clonal_pmhc_heatmap <- function(data_matrix) {
  
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(grid)
  
   if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }
  
  heatmap_plot <- Heatmap(
    data_matrix, 
    cluster_rows = FALSE,
    cluster_columns = FALSE, 
    show_row_names = TRUE,
    row_names_side = 'left',
    show_column_names = FALSE, 
    row_names_gp = gpar(fontsize = 10), 
    column_names_gp = gpar(fontsize = 10), 
    rect_gp = gpar(col = "grey", lwd = 0.2), 
    heatmap_legend_param = list(
      title = "UMI counts",
      legend_width = unit(2, "cm"), 
      legend_height = unit(4, "cm"), 
      color_bar = "continuous", position='left',
      legend_main_gp = gpar(fontsize = 10),
      legend_gp = gpar(fontsize = 8) 
    )
  )
  
  heatmap_grob <- grid.grabExpr(draw(heatmap_plot, heatmap_legend_side = "left"))
  heatmap_gg <- cowplot::ggdraw() + cowplot::draw_plot(heatmap_grob)
  
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
#' @return A cowplot object combining a heatmap of pMHC distribution, a jitter plot of UMI counts by feature, and a bar plot
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
#' @importFrom cowplot plot_grid draw_label
#' @export
pmhc_dist_per_clone <- function(object, clone, aggr_threshold=10, slot = 'counts', xlimits=NULL){
  
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
  
  heatmap_gg <- clonal_pmhc_heatmap(pmhc_matrix2)
  
  pmhc_matrix_long <- as.data.frame(pmhc_matrix)  %>% 
    rownames_to_column(var = "Feature") %>% 
    gather(key = "gem", value = "UMI", -Feature) %>%
    mutate(Feature = factor(Feature, levels=rev(preserve)))
  
  y_vals <- seq(0.5, nrow(pmhc_matrix) + 0.5)

  dists <- ggplot(pmhc_matrix_long, aes(x = UMI, y = Feature, color = Feature)) +
    geom_jitter() +
    geom_hline(yintercept = y_vals, linetype = "dashed", colour = "grey", size = 0.5) +
    theme_classic() +
    labs(x = "UMI counts", y = "Feature", color = "Feature") +
    NoLegend()
  
  if (!is.null(xlimits)){
    dists <- dists + xlim(xlimits)
  }
  
  concordance <- calculate_pmhc_concordance(clone_obj, slot=slot, preserve_pmhc = names(pmhc_labels)[pmhc_labels %in% preserve])
  names(concordance) <- recode(names(concordance), !!!pmhc_labels)
  concordance_df <- data.frame(pmhc = factor(names(concordance),levels=rev(preserve)),
                               concordance = concordance)
  
  conc_bp <- ggplot(concordance_df,
                    aes(y = pmhc, x = concordance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_hline(yintercept = y_vals, linetype = "dashed", colour = "grey", size = 0.5) +
    theme_classic() +
    labs(x = "concordance", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  combined_plot <- cowplot::plot_grid(
    heatmap_gg, dists+ylab(''), conc_bp,
    align = "h",
    ncol = 3,
    rel_widths = c(3, 6, 3.5), axis = 'l') 
  
  title <- cowplot::ggdraw() + cowplot::draw_label(paste0('clonotype = ', clone, '; #gems=', ncol(pmhc_matrix)), fontface='bold')
  
  return(cowplot::plot_grid(title, combined_plot, ncol = 1, rel_heights=c(0.05, 1)))
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
pmhc_heatmap <- function(object, clones, patient=NULL, slot='counts', clones_order=NULL, split_cols=FALSE,
                         hm_breaks = c(0, 4, 10), hm_palette=c("blue", "white", "red"), pmhc_palette=NULL,
                         condpalette=NULL, highlight_pmhc.tcr=F, stop_large_highlighting=T, rowm.fonts=8,  
                         column_title_fonts = 10, column_title_rot = 45, use_original_order=T, clean_mat=F, 
                         add_tcr_cluster=F, show_row_names=T, pmhc_subset=NULL, custom_annotations=c(), 
                         skip_bugged_frames=F, bugged_width=.6, max_cols=16000, ...){
  
  library(randomcoloR)
  library(circlize)
  library(ComplexHeatmap)
  
  if(length(clones) == 0) {
    stop("list of clonotypes in `clone` argument is empty")
  }
  # Additional steps to ensure custom_annotations are part of the object's metadata
  missing_cols <- setdiff(custom_annotations, colnames(object@meta.data))
  if(length(missing_cols) > 0){
    stop(paste("Missing columns in object's metadata:", paste(missing_cols, collapse=", ")))
  }
  
  if(is.null(object@misc$pmhc)){
    stop("Missing pmhc metadata in object@misc$pmhc")
  }
  
  cells_subset <- Cells(object)[(object$clone_id %in% clones) & !is.na(object$clone_id)]
  
  if (!'pMHC_classification' %in% colnames(object@meta.data)){
    object$pMHC_classification <- NA
  }
  
  ann_columns <- c("clone_id", "pMHC_classification", custom_annotations)
  ann_subset <- object@meta.data %>%
    dplyr::select(all_of(ann_columns)) %>%
    filter(row.names(.) %in% cells_subset) %>%
    arrange(clone_id)
  
  table(ann_subset$clone_id)
  
  ann_subset$clone_id <- factor(ann_subset$clone_id, levels=clones[clones %in% ann_subset$clone_id] %>% unique())
  
  pat_pmhc <- object@misc$pmhc %>%
    filter(!is.na(Sequence)) %>%
    { if(!is.null(patient)) filter(., Patient == patient) else . } %>%
    pull(Barcode)
  
  pat_pmhc <- pat_pmhc[pat_pmhc %in% rownames(object@assays$pMHC@counts)]
  
  pmhc_subset_ <- GetAssayData(object, layer = slot, assay = 'pMHC')
  pmhc_subset_ <- pmhc_subset_[,cells_subset][pat_pmhc,]
  
  pmhc_subset_ann <- object@misc$pmhc %>%
    filter(Barcode %in% row.names(pmhc_subset_)) %>%
    { if(!is.null(patient)) filter(., Patient == patient) else . } %>%
    as.data.frame() %>%
    column_to_rownames('Barcode')
  
  rownames(pmhc_subset_) <- pmhc_subset_ann[rownames(pmhc_subset_),]$pmhc
  
  if (!is.null(pmhc_subset)){
    pmhc_subset_ <- pmhc_subset_[rownames(pmhc_subset_) %in% pmhc_subset,] 
  }
  
  ncl <- ann_subset$clone_id %>%unique() %>% length()
  
  clpalette <- get_random_grid_colors(ncl, seed=3)
  
  names(clpalette) <- ann_subset$clone_id %>% unique()
  
  if (is.null(pmhc_palette)){
    palette_list <- list(clone_id = clpalette,
                         condition = condpalette)
  } else{
    palette_list <- list(clone_id = clpalette,
                         condition = condpalette,
                         pmhc = pmhc_palette)
  }
  
  if (!use_original_order){
    ann_subset$clone_id <- factor(ann_subset$clone_id, 
                                  levels=as.character(ann_subset$clone_id[order(ann_subset$epitope_type, ann_subset$pMHC_classification)]) %>% unique())
    ann_subset <- ann_subset[order(ann_subset$epitope_type, ann_subset$clone_id, ann_subset$condition),]
  } 
  
  pmhc_mat <- as.matrix(pmhc_subset_[, rownames(ann_subset)])
  
  if (clean_mat){
    positive_bc <- ann_subset$pMHC_classification[ann_subset$pMHC_classification != 'Negative']
    positive_bc <- positive_bc %>% unique() %>% 
      strsplit(., ":", fixed = TRUE) %>% unlist() %>% unique()
    pmhc_mat = pmhc_mat[positive_bc,]
  }
  
  col_fun = colorRamp2(hm_breaks, hm_palette)
  
  if (ncol(pmhc_subset_) > max_cols) {
    set.seed(123) 
    sampled_cols <- sample(colnames(pmhc_subset_), max_cols)
    
    pmhc_subset_ <- pmhc_subset_[, sampled_cols]
    ann_subset <- ann_subset[sampled_cols,]
    
    message("Number of columns exceeds threshold. Data has been randomly sampled to ", max_cols, " columns.")
  }
    
  if (is.null(clones_order)){
    cells_order <- ann_subset %>% 
      arrange(pMHC_classification, clone_id) %>%
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
  
  custom_ann_list <- setNames(lapply(custom_annotations, function(cn) ann_subset[cells_order,][[cn]]), custom_annotations)
  
  ann_list <- list(
    clone_id = ann_subset[cells_order,]$clone_id, 
    pmhc = gsub(':', '\n', ann_subset[cells_order,]$pMHC_classification)
  )
  ann_list <- c(ann_list, custom_ann_list) # Combine lists
  
  hm22ann <- do.call(HeatmapAnnotation, c(ann_list, list(
    col = palette_list, show_legend = c(F, T, T),
    which = 'col',
    annotation_width = unit(c(1, 4), 'cm'),
    gap = unit(1, 'mm')
  )))
  
  if (split_cols){
    split_cols <- ann_subset[cells_order,]$clone_id
  } else {
    split_cols <- NULL
  }
  
  pmhc_subset_hmap <- Heatmap(
    pmhc_mat_ordered, name = "pmhc_subset_hmap",
    show_row_names = show_row_names, show_column_names = FALSE,
    col = col_fun, column_split = split_cols,
    cluster_rows = F, cluster_columns = F,
    row_names_gp = grid::gpar(fontsize = rowm.fonts),  
    column_title_gp = gpar(fontsize = column_title_fonts), 
    column_title_rot = column_title_rot,
    top_annotation=hm22ann)
  
  if (!highlight_pmhc.tcr){
    return(pmhc_subset_hmap)
  } else {
    if (length(clones) > 300 & stop_large_highlighting){
      stop('highlighting more than 300 clones specificity usually crushes the session, \n
           if you want to continue, set stop_large_highlighting=F')
    }
    if (!is.null(split_cols)){
      stop('Highlightinh specificities in the clone splited heatmap is not supported')
    }
    draw(pmhc_subset_hmap)
    tcr_pmhc <- ann_subset_ordered %>%
      dplyr::select(clone_id, pMHC_classification) %>%
      filter(!is.na(pMHC_classification) & pMHC_classification != 'Negative') %>%
      unique() %>%
      separate_rows(pMHC_classification, sep=':')
    
    n_iter <- nrow(tcr_pmhc)
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = n_iter, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "+") 
    
    for (i in 1:n_iter){
      cell_ids_i <- ann_subset_ordered$clone_id == tcr_pmhc[i,]$clone_id
      cell_ids_i[is.na(cell_ids_i)] <- FALSE
      cell_ids_i <- match(rownames(ann_subset_ordered)[cell_ids_i], colnames(pmhc_mat_ordered))
      
      pmhc_ids_i <- which(rownames(pmhc_mat_ordered) == tcr_pmhc[i,]$pMHC_classification)
      pmhc_ids_i <- rep(pmhc_ids_i, length(cell_ids_i))
      
      decorate_heatmap_body("pmhc_subset_hmap", {
        for (row in pmhc_ids_i) {
          cols_in_row <- cell_ids_i
          min_col <- min(cols_in_row)
          max_col <- max(cols_in_row)
          
          x = (min_col - 1) / ncol(pmhc_mat_ordered)
          width = (max_col - min_col + 1) / ncol(pmhc_mat_ordered)
          
          y = 1 - (row - 0.5) / nrow(pmhc_mat_ordered)
          height = 1 / nrow(pmhc_mat_ordered)
          
          if (skip_bugged_frames & width > bugged_width){
            next
          }
          
          grid.rect(x = x, y = y, width = width, height = height, just = c("left", "center"), 
                    gp = gpar(col = "green", lwd = 2, fill = NA))
        }
      })
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
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
#' @importFrom uniformly runif_in_cube
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

