# ============================================================
# Shared fibroblast subsetting and reclustering functions
# Cross-species intestinal fibroblast analysis
# ============================================================

library(Seurat)
library(dplyr)


# ------------------------------------------------------------
# Validate that broad cell-type annotations are available
# ------------------------------------------------------------

check_cell_type_annotations <- function(
    object,
    metadata_col = "cell_type"
) {
  if (!metadata_col %in% colnames(object[[]])) {
    stop(
      "The object does not contain a `",
      metadata_col,
      "` metadata column. ",
      "Apply the broad cell-type annotations before fibroblast subsetting."
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------
# Subset fibroblasts and perform fibroblast-specific reclustering
# ------------------------------------------------------------

prepare_fibroblast_subset <- function(
    object,
    species_name,
    fibroblast_label = "Fibroblasts",
    metadata_col = "cell_type",
    dims = 1:20,
    resolution = 0.1,
    marker_pct_diff = 0.3,
    seed = 1234
) {
  check_cell_type_annotations(
    object,
    metadata_col = metadata_col
  )

  message(
    "Subsetting fibroblasts for ",
    species_name,
    "..."
  )

  fib <- subset_cell_type(
    object,
    cell_type = fibroblast_label,
    metadata_col = metadata_col
  )

  message(
    species_name,
    ": ",
    ncol(fib),
    " fibroblast cells retained."
  )

  # Re-run the analysis within the fibroblast population.
  fib <- run_sct_clustering(
    fib,
    dims = dims,
    resolution = resolution,
    seed = seed
  )

  markers <- find_cluster_markers(
    fib,
    assay = "SCT",
    min_pct_diff = marker_pct_diff
  )

  list(
    object = fib,
    markers = markers
  )
}
