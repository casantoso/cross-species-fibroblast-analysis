# ============================================================
# Export a Seurat fibroblast object for pySCENIC
#
# Writes:
#   counts.mtx   - raw RNA counts (genes x cells)
#   genes.txt    - gene names
#   cells.txt    - cell IDs
#   metadata.csv - Seurat metadata aligned to the count matrix
#   umap.tsv     - original Seurat UMAP coordinates, when available
# ============================================================

library(Seurat)
library(Matrix)
library(readr)


export_for_pyscenic <- function(
    object,
    output_dir,
    reduction = "umap"
) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  if (!"RNA" %in% Assays(object)) {
    stop("The object does not contain an RNA assay.")
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  DefaultAssay(object) <- "RNA"

  # Join Seurat v5 layers before extracting counts.
  object <- JoinLayers(
    object,
    assay = "RNA"
  )

  counts <- GetAssayData(
    object,
    assay = "RNA",
    layer = "counts"
  )

  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("The RNA count matrix is empty.")
  }

  metadata <- object[[]]
  metadata$cell_id <- rownames(metadata)

  # Force metadata into exactly the same cell order as the matrix.
  metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE
  ]

  Matrix::writeMM(
    counts,
    file.path(output_dir, "counts.mtx")
  )

  readr::write_lines(
    rownames(counts),
    file.path(output_dir, "genes.txt")
  )

  readr::write_lines(
    colnames(counts),
    file.path(output_dir, "cells.txt")
  )

  write.csv(
    metadata,
    file.path(output_dir, "metadata.csv"),
    row.names = FALSE
  )

  if (reduction %in% Reductions(object)) {
    umap_df <- as.data.frame(
      Embeddings(
        object,
        reduction = reduction
      )
    )

    umap_df$cell_id <- rownames(umap_df)

    write.table(
      umap_df,
      file.path(output_dir, "umap.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  } else {
    warning(
      "Reduction '", reduction,
      "' was not found; UMAP coordinates were not exported."
    )
  }

  message(
    "Exported ",
    nrow(counts),
    " genes x ",
    ncol(counts),
    " cells to ",
    normalizePath(output_dir, mustWork = FALSE)
  )

  invisible(object)
}


# Example:
#
# human_fib <- readRDS("outputs/human_fibroblasts.rds")
#
# export_for_pyscenic(
#   human_fib,
#   "outputs/preprocessing/pyscenic_exports/human"
# )
