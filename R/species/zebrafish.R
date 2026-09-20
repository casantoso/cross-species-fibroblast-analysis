# ============================================================
# Zebrafish intestinal scRNA-seq preprocessing
# Dataset: GSE271002
# Sample used: GSM8366962
# Description: Control, 5 dpf
#
# This is the non-Mtz control sample used as the baseline
# intestinal condition for the cross-species comparison.
# Public data; unpublished biological results are omitted.
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(scDblFinder)

source("R/functions/preprocessing_functions.R")
source("R/functions/fibroblast_functions.R")
set.seed(1234)

# Load control 5 dpf sample
zf_counts <- Read10X_h5(
  "data/GSE271002/GSM8366962_Control_5dpf_filtered_feature_bc_matrix.h5"
)

zebrafish <- CreateSeuratObject(
  counts = zf_counts,
  project = "zebrafish_control_5dpf",
  min.cells = 3,
  min.features = 200
)

# QC
VlnPlot(
  zebrafish,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 2
)

zebrafish <- filter_cells(
  zebrafish,
  min_features = 200,
  max_features = 5000
)

# Doublet removal
zebrafish <- remove_doublets(zebrafish)

# Normalization and clustering
zebrafish <- run_sct_clustering(
  zebrafish,
  dims = 1:20,
  resolution = 0.8
)

# Marker discovery
zebrafish_markers <- find_cluster_markers(
  zebrafish,
  assay = "SCT",
  min_pct_diff = 0.2
)

canonical_markers <- c(
  "pdgfra", "dcn", "lum",
  "acta2", "myh11a", "tagln",
  "cdh5", "plvapb", "kdrl",
  "pdgfrb", "rgs5",
  "epcam"
)

plot_available_markers(
  zebrafish,
  canonical_markers
)

# Apply private annotations locally
annotation_file <- "private_annotations/zebrafish_cluster_annotations.csv"

if (file.exists(annotation_file)) {
  annotations <- read.csv(annotation_file)

  zebrafish <- apply_cluster_annotations(
    zebrafish,
    annotations
  )

  zebrafish_fib <- subset_cell_type(
    zebrafish,
    cell_type = "Fibroblasts"
  )

  zebrafish_fib <- run_sct_clustering(
    zebrafish_fib,
    dims = 1:20,
    resolution = 0.1
  )

  zebrafish_fib_markers <- find_cluster_markers(
    zebrafish_fib,
    assay = "SCT",
    min_pct_diff = 0.4
  )
}

# ============================================================
# Fibroblast subsetting and reclustering
# ============================================================

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# Save the broadly annotated object locally.
saveRDS(
  zebrafish,
  "outputs/zebrafish_processed.rds"
)

# Continue only when private broad cell-type annotations were loaded.
if ("cell_type" %in% colnames(zebrafish[[]])) {

  fib_result <- prepare_fibroblast_subset(
    zebrafish,
    species_name = "Zebrafish",
    resolution = 0.1,
    marker_pct_diff = 0.3
  )

  zebrafish_fib <- fib_result$object
  zebrafish_fib_markers <- fib_result$markers

  saveRDS(
    zebrafish_fib,
    "outputs/zebrafish_fibroblasts.rds"
  )

  write.csv(
    zebrafish_fib_markers,
    "outputs/zebrafish_fibroblast_cluster_markers.csv",
    row.names = FALSE
  )

  # No fibroblast subtype labels are assigned in this public script.
}
