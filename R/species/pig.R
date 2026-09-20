# ============================================================
# Pig ileum snRNA-seq preprocessing
# Dataset: GSE233285
# Sample: GSM7421580
# Description: Ileum, replicate 1, snRNA-seq
# Public data; unpublished biological results are omitted.
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(scDblFinder)

source("R/functions/preprocessing_functions.R")
source("R/functions/fibroblast_functions.R")
set.seed(1234)

# Load the ileum sample
pig_counts <- Read10X(
  data.dir = "data/GSE233285/GSM7421580_ileum"
)

pig <- CreateSeuratObject(
  counts = pig_counts,
  project = "pig_ileum",
  min.cells = 3,
  min.features = 200
)

# QC
VlnPlot(
  pig,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 2
)

pig <- filter_cells(
  pig,
  min_features = 200,
  max_features = 5000
)

# Doublet removal
pig <- remove_doublets(pig)

# Normalization and clustering
pig <- run_sct_clustering(
  pig,
  dims = 1:20,
  resolution = 0.8
)

# Marker discovery
pig_markers <- find_cluster_markers(
  pig,
  assay = "SCT",
  min_pct_diff = 0.2
)

canonical_markers <- c(
  "PDGFRA", "DCN", "LUM", "COL1A1",
  "MYH11", "ACTG2", "ACTA2",
  "VWF", "PECAM1",
  "PDGFRB", "RGS5",
  "CD3D", "CD3E",
  "CD19", "MS4A1",
  "CDH1", "EPCAM",
  "AIF1", "LYZ", "FCER1G"
)

plot_available_markers(
  pig,
  canonical_markers
)

# Apply private annotations locally
annotation_file <- "private_annotations/pig_cluster_annotations.csv"

if (file.exists(annotation_file)) {
  annotations <- read.csv(annotation_file)

  pig <- apply_cluster_annotations(
    pig,
    annotations
  )

  pig_fib <- subset_cell_type(
    pig,
    cell_type = "Fibroblasts"
  )

  pig_fib <- run_sct_clustering(
    pig_fib,
    dims = 1:20,
    resolution = 1.0
  )

  pig_fib_markers <- find_cluster_markers(
    pig_fib,
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
  pig,
  "outputs/pig_processed.rds"
)

# Continue only when private broad cell-type annotations were loaded.
if ("cell_type" %in% colnames(pig[[]])) {

  fib_result <- prepare_fibroblast_subset(
    pig,
    species_name = "Pig",
    resolution = 1.0,
    marker_pct_diff = 0.3
  )

  pig_fib <- fib_result$object
  pig_fib_markers <- fib_result$markers

  saveRDS(
    pig_fib,
    "outputs/pig_fibroblasts.rds"
  )

  write.csv(
    pig_fib_markers,
    "outputs/pig_fibroblast_cluster_markers.csv",
    row.names = FALSE
  )

  # No fibroblast subtype labels are assigned in this public script.
}
