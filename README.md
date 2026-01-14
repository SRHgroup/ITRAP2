Seurat-integrated R workflow for TCR-pMHC pairing.

The workflow is intended for T-cell scientists and Seurat users who want to analyze their pMHC data.
This repository is private until we post on bioRxiv. To install, collaborators in SRHgroup need a GitHub Personal Access Token (PAT).

One-time setup (inside R)

```R
# One-time setup
install.packages(c("usethis", "gitcreds", "devtools"), repos = "https://cloud.r-project.org")

# 1) Open GitHub's token page with sensible defaults (select at least the "repo" scope)
usethis::create_github_token()

# 2) Store the token (securely) with your Git credential manager
gitcreds::gitcreds_set()    # paste the token when prompted

# -- For devtools/remotes installs from private repos,
#    ensure a GITHUB_PAT is available in the R session:
usethis::edit_r_environ()   # then add a line like:
# GITHUB_PAT=ghp_your_token_here
# Save, then restart R so it loads.
```

In the file that opens, add a line like:
```
GITHUB_PAT=ghp_your_token_here
```

Save, restart R, then verify:
```R
Sys.getenv("GITHUB_PAT")  # should not be empty
```

Install the package
```R
devtools::install_github("SRHgroup/ITRAP2")
```

Minimal vignette: denoising and assignment
For a more detailed vignette, please check out vignettes/Simple_case_vignette.html

The prerequisites
object is a Seurat object with a pMHC assay. (object@assays$pMHC)
Clonotypes available (store as object$clone_id).
A pMHC barcode annotation table in object@misc$pmhc (columns like: Barcode, HLA, Sequence, Virus, Protein, pmhc).
Aim for ≥ 20 distinct pMHC barcodes per experiment for stable assignments.

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
Next scaled counts are subjected to smoothing using local regression. The scaled matrix will be updated in object@assays$pMHC@scale.data.

```R
object <- smooth_pmhc(object)
```

Next, we want to label each expanded clone with pMHC specificity, we do that by doing the outlier search within your pMHC panel, since we look at the relative to your pMHC panel values, 
we expect that you have at least 20 different pMHC barcodes in your experiment, otherwise, the pMHC assignments with this functions might be flawed. For the assignments we also want you to
have a table with pMHC barcode annotations, that you store in object@misc$pmhc. here is an example
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
