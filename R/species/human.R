# ============================================================
# Human intestinal scRNA-seq preprocessing
# Dataset: GSE178341
# Tissue used: Normal colon
# Public data; unpublished biological results are omitted.
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(scDblFinder)

source("R/functions/preprocessing_functions.R")
source("R/functions/fibroblast_functions.R")
set.seed(1234)

# Load counts and metadata
counts <- Read10X_h5("data/GSE178341/GSE178341.h5")
human <- CreateSeuratObject(
  counts = counts,
  min.cells = 3,
  min.features = 200
)

meta <- read.csv(
  "data/GSE178341/GSE178341_metatables.csv.gz",
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$cellID
meta$cellID <- NULL
meta <- meta[colnames(human), , drop = FALSE]
human <- AddMetaData(human, metadata = meta)

# Keep normal colon cells
human <- subset(
  human,
  subset = HistologicTypeSimple == "Normal colon"
)

# QC
human <- add_mito_percent(
  human,
  pattern = "^MT-"
)

VlnPlot(
  human,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3
)

human <- filter_cells(
  human,
  min_features = 200,
  max_features = 6000,
  max_percent_mt = 20
)

# Doublet removal
human <- remove_doublets(human)

# Normalization and clustering
human <- run_sct_clustering(
  human,
  dims = 1:20,
  resolution = 0.1
)

# Marker discovery
human_markers <- find_cluster_markers(
  human,
  assay = "SCT",
  min_pct_diff = 0.2
)

canonical_markers <- c(
  "PDGFRA", "DCN", "LUM",
  "MYH11", "ACTG2", "ACTA2",
  "VWF", "PECAM1",
  "RGS5",
  "CD3D", "CD3E",
  "CD19", "MS4A1",
  "CDH1", "EPCAM",
  "AIF1", "LYZ", "FCER1G"
)

plot_available_markers(
  human,
  canonical_markers
)

# Apply private annotations locally
annotation_file <- "private_annotations/human_cluster_annotations.csv"

if (file.exists(annotation_file)) {
  annotations <- read.csv(annotation_file)

  human <- apply_cluster_annotations(
    human,
    annotations
  )

  human_fib <- subset_cell_type(
    human,
    cell_type = "Fibroblasts"
  )

  human_fib <- run_sct_clustering(
    human_fib,
    dims = 1:20,
    resolution = 0.1
  )

  human_fib_markers <- find_cluster_markers(
    human_fib,
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
  human,
  "outputs/human_processed.rds"
)

# Continue only when private broad cell-type annotations were loaded.
if ("cell_type" %in% colnames(human[[]])) {

  fib_result <- prepare_fibroblast_subset(
    human,
    species_name = "Human",
    resolution = 0.1,
    marker_pct_diff = 0.3
  )

  human_fib <- fib_result$object
  human_fib_markers <- fib_result$markers

  saveRDS(
    human_fib,
    "outputs/human_fibroblasts.rds"
  )

  write.csv(
    human_fib_markers,
    "outputs/human_fibroblast_cluster_markers.csv",
    row.names = FALSE
  )

  # Fibroblast subtype assignment is performed separately and kept private.
}
