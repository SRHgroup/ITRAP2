ITRAP2 is a Seurat-integrated R workflow for denoising and assignment of TCR–pMHC specificities from single-cell pMHC multimer data.

The package is intended for T-cell scientists and Seurat users, working on single-cell assignments of TCR specificities. 

This repository is public, and the tool is available for free for academic use.
Contributions via pull requests and issue reports are welcome. 

Install the package
```R
devtools::install_github("SRHgroup/ITRAP2")
```


Dependency notes

Some dependencies are Bioconductor packages and may not be installed automatically, depending on the local R setup.
If you encounter errors such as:
```R
Error in library(shades) : there is no package called ‘shades’
Error in library(ComplexHeatmap) : there is no package called ‘ComplexHeatmap’for this version of R
```
Install them manually
```R
install.packages('BiocManager') # if it's not installed already
BiocManager::install('shades')
BiocManager::install('ComplexHeatmap')

#and then try again

devtools::install_github('SRHgroup/ITRAP2')
```


Minimal vignette: denoising and assignment

For a detailed walkthrough, see:
vignettes/Simple_case_vignette.html

Prerequisites

The input object must satisfy the following:

`object` is a Seurat object containing a pMHC assay (object@assays$pMHC)

Clonotype information stored as object$clone_id

A pMHC barcode annotation table stored in object@misc$pmhc
(columns such as: Barcode, HLA, Sequence, Virus, Protein, pmhc)

≥ 20 distinct pMHC barcodes per experiment are recommended for stable panel-based assignments

1) Estimate per-pMHC noise (for downstream confidence scoring)
```R
object <- score_pmhc_noise(object)
# Adds/updates per-pMHC noise metrics used to compute TCR–pMHC confidence later.
```

The package employs a specific way of mean centering and sd scaling, where means and sds are calculated without outliers. 
Used as an alternative to median and mad, since in sparse datasets they often end up as 0.
The scaled pMHC counts matrix will be stored in object@assays$pMHC@scale.data
![Screenshot 2024-03-06 at 12 43 22](https://github.com/SRHgroup/pmhc_denioseR/assets/45093246/8b7fbe80-e8b3-4af4-b496-9133cd167367)

```R
object <- ScaleDataNoOutliers(object)
```

For the next step, we would like you to store your clonotype information in object@metadata$clone_id. We will perform an extra denoising step on every expanded clone.
Next, scaled counts are subjected to smoothing using local regression. The scaled matrix will be updated in object@assays$pMHC@scale.data.

```R
object <- smooth_pmhc(object)
```

Next, we want to label each expanded clone with pMHC specificity. We do that by doing the outlier search within your pMHC panel. Since we look at the relative to your pMHC panel values, 
we expect that you have at least 20 different pMHC barcodes in your experiment; otherwise, the pMHC assignments with this function might be flawed. For the assignments, we also want you to
have a table with pMHC barcode annotations, which you store in object@misc$pmhc. Here is an example
![Screenshot 2024-03-06 at 12 42 55](https://github.com/SRHgroup/pmhc_denioseR/assets/45093246/400afcd4-f1d9-4f91-8ed9-785fa70af79a)
This is done in case the same barcode corresponds to different pMHC within the same experiment, but in different donors, that were demultiplexed by hashing antibodies. In this simple example, row names in the
matrix are the same as the pMHC.


```R
object <- assign_pmhc(object)
```

By the end of this step, your object@meta.data will have a column called  pMHC_classification, it can have multiple pMHCs recorded there, in case the clone is cross-reactive.

```R
object@meta.data %>%
  select(pMHC_classification, pMHC_confidence) %>% 
  head(10)
```
![Screenshot 2024-03-27 at 10 54 10](https://github.com/SRHgroup/pmhc_denioseR/assets/45093246/15395c3c-0970-4605-a404-5484a1e94109)
