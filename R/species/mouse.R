# ============================================================
# Mouse intestinal scRNA-seq preprocessing
# Internal dataset
# Input data and unpublished biological results are omitted.
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

source("R/functions/preprocessing_functions.R")
source("R/functions/fibroblast_functions.R")
set.seed(1234)

# Load internal Seurat object
mouse <- readRDS("path/to/private_mouse_dataset.rds")

# Inspect QC metrics already present in the supplied object
VlnPlot(
  mouse,
  features = intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    colnames(mouse[[]])
  ),
  ncol = 3
)

# Final normalization and clustering
mouse <- run_sct_clustering(
  mouse,
  dims = 1:20,
  resolution = 0.1
)

# Marker discovery
mouse_markers <- find_cluster_markers(
  mouse,
  assay = "SCT",
  min_pct_diff = 0.4
)

canonical_markers <- c(
  "Pdgfra", "Dcn", "Lum",
  "Myh11", "Actg2", "Acta2",
  "Vwf", "Pecam1",
  "Rgs5"
)

plot_available_markers(
  mouse,
  canonical_markers
)

# Apply private annotations locally
annotation_file <- "private_annotations/mouse_cluster_annotations.csv"

if (file.exists(annotation_file)) {
  annotations <- read.csv(annotation_file)

  mouse <- apply_cluster_annotations(
    mouse,
    annotations
  )

  mouse_fib <- subset_cell_type(
    mouse,
    cell_type = "Fibroblasts"
  )

  mouse_fib <- run_sct_clustering(
    mouse_fib,
    dims = 1:20,
    resolution = 0.1
  )

  mouse_fib_markers <- find_cluster_markers(
    mouse_fib,
    assay = "SCT",
    min_pct_diff = 0.3
  )
}

# ============================================================
# Fibroblast subsetting and reclustering
# ============================================================

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# Save the broadly annotated object locally.
saveRDS(
  mouse,
  "outputs/mouse_processed.rds"
)

# Continue only when private broad cell-type annotations were loaded.
if ("cell_type" %in% colnames(mouse[[]])) {

  fib_result <- prepare_fibroblast_subset(
    mouse,
    species_name = "Mouse",
    resolution = 0.1,
    marker_pct_diff = 0.3
  )

  mouse_fib <- fib_result$object
  mouse_fib_markers <- fib_result$markers

  saveRDS(
    mouse_fib,
    "outputs/mouse_fibroblasts.rds"
  )

  write.csv(
    mouse_fib_markers,
    "outputs/mouse_fibroblast_cluster_markers.csv",
    row.names = FALSE
  )

  # Fibroblast subtype assignment is performed separately and kept private.
}
