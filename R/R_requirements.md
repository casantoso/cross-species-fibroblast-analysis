# R dependencies

The exact package versions reported from the working analysis environment are
stored in:

```text
environments/r_package_versions.tsv
```

The recorded R version is:

```text
R 4.5.3
```

## Direct packages used by the public R code

These packages are called directly by the scripts/functions in this repository
or are explicitly required by a selected analysis method:

| Package | Recorded version | Why it is needed |
|---|---:|---|
| Seurat | 5.5.1 | Seurat object handling, QC, SCTransform, PCA, clustering, UMAP, markers, plotting |
| Matrix | 1.7.6 | sparse matrix export for SAMap and pySCENIC |
| glmGamPoi | 1.22.0 | `SCTransform(method = "glmGamPoi")` |
| harmony | 2.0.5 | replicate/sample integration for multi-sample species |
| scDblFinder | 1.24.10 | doublet detection |
| SingleCellExperiment | 1.32.0 | Seurat → SingleCellExperiment conversion for scDblFinder |
| dplyr | 1.2.1 | table manipulation |
| tidyr | 1.3.2 | reshaping tables |
| purrr | 1.2.2 | functional iteration in plotting/helpers |
| readr | 2.2.0 | text export |
| stringr | 1.6.0 | enrichment-term cleanup |
| ggplot2 | 4.0.3 | plotting |
| patchwork | 1.3.2 | combining plots |
| enrichR | 3.4 | Enrichr pathway enrichment |
| orthogene | 1.16.1 | pairwise ortholog conversion |
| jsonlite | 2.0.0 | reading/writing shared configuration |
| testthat | install for testing | R unit tests |

`testthat` is required only to run the R tests; it is not part of the analysis
itself.

## Recorded environment packages that are not direct requirements

The exact environment snapshot also records packages such as:

```text
SeuratObject
SummarizedExperiment
S4Vectors
tibble
gprofiler2
reticulate
```

These are useful to retain in the version record because they were present in
the working environment, but the current public R scripts do not call them
directly.

Some are dependencies of the packages above. For example, Seurat uses
SeuratObject, and Bioconductor single-cell packages depend on infrastructure
such as SummarizedExperiment and S4Vectors.

`gprofiler2` is relevant to the g:Profiler backend used through `orthogene`,
but the public helper calls `orthogene` rather than calling `gprofiler2`
directly.

`reticulate` is not currently used by the public R workflow.

## Reproducibility note

This repository records the package versions that were actually present in the
working R environment rather than generating a speculative `renv.lock`.
Package managers will install additional dependencies automatically as needed.
