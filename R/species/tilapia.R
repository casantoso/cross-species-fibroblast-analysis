# ============================================================
# Nile tilapia gut scRNA-seq preprocessing
# Dataset: GSE284663
# Samples used:
#   GSM8689993 - probiotic untreated, biological replicate 1
#   GSM8689994 - probiotic untreated, biological replicate 2
#   GSM8689995 - probiotic untreated, biological replicate 3
#
# Only untreated baseline samples are used.
# Public data; unpublished biological results are omitted.
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(harmony)

source("R/functions/preprocessing_functions.R")
source("R/functions/fibroblast_functions.R")
set.seed(1234)

# Load untreated biological replicates
sample_files <- c(
  GSM8689993 = "data/GSE284663/GSM8689993.txt.gz",
  GSM8689994 = "data/GSE284663/GSM8689994.txt.gz",
  GSM8689995 = "data/GSE284663/GSM8689995.txt.gz"
)

tilapia_samples <- lapply(
  names(sample_files),
  function(sample_name) {
    x <- read.table(
      sample_files[[sample_name]],
      header = TRUE,
      row.names = 1,
      sep = "\t",
      check.names = FALSE
    )

    # GEO matrices include a non-count GeneName column.
    x <- x[, -1, drop = FALSE]

    CreateSeuratObject(
      counts = as.matrix(x),
      project = sample_name,
      min.cells = 3,
      min.features = 200
    )
  }
)

names(tilapia_samples) <- names(sample_files)

# QC per biological replicate
for (sample_name in names(tilapia_samples)) {
  print(
    VlnPlot(
      tilapia_samples[[sample_name]],
      features = c("nCount_RNA", "nFeature_RNA"),
      ncol = 2
    ) +
      ggtitle(sample_name)
  )
}

tilapia_samples <- lapply(
  tilapia_samples,
  filter_cells,
  min_features = 200,
  max_features = 4000,
  max_counts = 20000
)

# Doublet removal separately per replicate
tilapia_samples <- lapply(
  tilapia_samples,
  remove_doublets
)

# Merge untreated replicates
tilapia <- merge(
  tilapia_samples[[1]],
  y = tilapia_samples[-1]
)

# Normalize and run PCA
tilapia <- SCTransform(
  tilapia,
  assay = "RNA",
  method = "glmGamPoi",
  verbose = FALSE
)

tilapia <- RunPCA(
  tilapia,
  verbose = FALSE
)

# Inspect sample structure before Harmony
tilapia <- RunUMAP(
  tilapia,
  dims = 1:20,
  seed.use = 1234,
  verbose = FALSE
)

DimPlot(
  tilapia,
  group.by = "orig.ident"
) +
  ggtitle("Before Harmony")

# Harmony across untreated biological replicates
tilapia <- run_harmony_clustering(
  tilapia,
  batch_var = "orig.ident",
  dims = 1:20,
  resolution = 1.0
)

# Marker discovery
tilapia_markers <- find_cluster_markers(
  tilapia,
  assay = "SCT",
  min_pct_diff = 0.2
)

canonical_markers <- c(
  "pdgfra",
  "pdgfrb",
  "acta2", "myh11",
  "epcam",
  "cd79a", "cd79b",
  "spi1", "irf8", "fcer1g", "cd74",
  "cd3e", "cd2", "lck",
  "vwf"
)

plot_available_markers(
  tilapia,
  canonical_markers
)

# Apply private annotations locally
annotation_file <- "private_annotations/tilapia_cluster_annotations.csv"

if (file.exists(annotation_file)) {
  annotations <- read.csv(annotation_file)

  tilapia <- apply_cluster_annotations(
    tilapia,
    annotations
  )

  tilapia_fib <- subset_cell_type(
    tilapia,
    cell_type = "Fibroblasts"
  )

  tilapia_fib <- run_sct_clustering(
    tilapia_fib,
    dims = 1:20,
    resolution = 0.1
  )

  tilapia_fib_markers <- find_cluster_markers(
    tilapia_fib,
    assay = "SCT",
    min_pct_diff = 0.2
  )
}

# ============================================================
# Fibroblast subsetting and reclustering
# ============================================================

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# Save the broadly annotated object locally.
saveRDS(
  tilapia,
  "outputs/tilapia_processed.rds"
)

# Continue only when private broad cell-type annotations were loaded.
if ("cell_type" %in% colnames(tilapia[[]])) {

  fib_result <- prepare_fibroblast_subset(
    tilapia,
    species_name = "Nile tilapia",
    resolution = 0.1,
    marker_pct_diff = 0.3
  )

  tilapia_fib <- fib_result$object
  tilapia_fib_markers <- fib_result$markers

  saveRDS(
    tilapia_fib,
    "outputs/tilapia_fibroblasts.rds"
  )

  write.csv(
    tilapia_fib_markers,
    "outputs/tilapia_fibroblast_cluster_markers.csv",
    row.names = FALSE
  )

  # No fibroblast subtype labels are assigned in this public script.
}
