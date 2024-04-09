#' Viral screen, data example for pmhcDenoiseR
#'
#' A dataset containing information about something interesting.
#'
#' @format Seurat object with 100 pMHCs
#' \describe{
#'   \item{viral_screen}{PBMC T cells, screen for chronic viral epitopes}
#' }
#' @source Source of your data (Sine Reker Hadrup lab).
#' @examples
#' data(viral_screen)
#' summary(viral_screen)
"viral_screen"

#' 10x Genomics Dataset
#'
#' A dataset containing information from a 10x Genomics experiment.
#'
#' @format Seurat object with 49 pMHCs
#' \describe{
#'   \item{tenx.data}{downsampled 10c genomics public pMHC multimer screen}
#' }
#' @source Source of your data (10c genomics public pMHC multimer screen).
#' @examples
#' data("10x_dataset")
#' summary(tenx.data)
"tenx.data"
