# ============================================================
# Chicken intestinal scRNA-seq preprocessing
# Dataset: GSE283090
# Samples used:
#   GSM8655812 - D0 Broiler, replicate 1
#   GSM8655813 - D0 Broiler, replicate 2
#   GSM8655814 - D0 Layer, replicate 1
#   GSM8655815 - D0 Layer, replicate 2
#
# D0 samples correspond to dissociated intestinal villi collected
# before 3 days of enteroid culture.
# Tissue was pooled from duodenum, jejunum, and ileum of ED19 embryos.
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

# Load the four D0 control/source-tissue samples
sample_dirs <- c(
  GSM8655812 = "data/GSE283090/GSM8655812",
  GSM8655813 = "data/GSE283090/GSM8655813",
  GSM8655814 = "data/GSE283090/GSM8655814",
  GSM8655815 = "data/GSE283090/GSM8655815"
)

chicken_samples <- lapply(
  names(sample_dirs),
  function(sample_name) {
    counts <- Read10X(
      data.dir = sample_dirs[[sample_name]]
    )

    CreateSeuratObject(
      counts = counts,
      project = sample_name,
      min.cells = 3,
      min.features = 200
    )
  }
)

names(chicken_samples) <- names(sample_dirs)

# Inspect QC per sample
for (sample_name in names(chicken_samples)) {
  print(
    VlnPlot(
      chicken_samples[[sample_name]],
      features = c("nCount_RNA", "nFeature_RNA"),
      ncol = 2
    ) +
      ggtitle(sample_name)
  )
}

# Apply final sample-specific QC thresholds here if needed.
# Only document thresholds that were actually used.

# Remove doublets separately for each biological replicate
chicken_samples <- lapply(
  chicken_samples,
  remove_doublets
)

# Merge raw-count objects
chicken <- merge(
  chicken_samples[[1]],
  y = chicken_samples[-1]
)

# Normalize and run PCA
chicken <- SCTransform(
  chicken,
  assay = "RNA",
  method = "glmGamPoi",
  verbose = FALSE
)

chicken <- RunPCA(
  chicken,
  verbose = FALSE
)

# Inspect sample structure before batch correction
chicken <- RunUMAP(
  chicken,
  dims = 1:20,
  seed.use = 1234,
  verbose = FALSE
)

DimPlot(
  chicken,
  group.by = "orig.ident"
) +
  ggtitle("Before Harmony")

# Harmony across the four D0 biological samples
chicken <- run_harmony_clustering(
  chicken,
  batch_var = "orig.ident",
  dims = 1:20,
  resolution = 0.8
)

# Marker discovery
chicken_markers <- find_cluster_markers(
  chicken,
  assay = "SCT",
  min_pct_diff = 0.2
)

canonical_markers <- c(
  "PDGFRA", "DCN", "LUM", "COL1A1",
  "ACTA2", "ACTG2", "MYH11", "MYOCD",
  "PDGFRB", "RGS5",
  "CDH1", "EPCAM",
  "PECAM1", "VWF",
  "CD3E",
  "CD68", "CD74"
)

plot_available_markers(
  chicken,
  canonical_markers
)

# Apply private annotations locally
annotation_file <- "private_annotations/chicken_cluster_annotations.csv"

if (file.exists(annotation_file)) {
  annotations <- read.csv(annotation_file)

  chicken <- apply_cluster_annotations(
    chicken,
    annotations
  )

  chicken_fib <- subset_cell_type(
    chicken,
    cell_type = "Fibroblasts"
  )

  chicken_fib <- run_sct_clustering(
    chicken_fib,
    dims = 1:20,
    resolution = 0.5
  )

  chicken_fib_markers <- find_cluster_markers(
    chicken_fib,
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
  chicken,
  "outputs/chicken_processed.rds"
)

# Continue only when private broad cell-type annotations were loaded.
if ("cell_type" %in% colnames(chicken[[]])) {

  fib_result <- prepare_fibroblast_subset(
    chicken,
    species_name = "Chicken",
    resolution = 0.5,
    marker_pct_diff = 0.3
  )

  chicken_fib <- fib_result$object
  chicken_fib_markers <- fib_result$markers

  saveRDS(
    chicken_fib,
    "outputs/chicken_fibroblasts.rds"
  )

  write.csv(
    chicken_fib_markers,
    "outputs/chicken_fibroblast_cluster_markers.csv",
    row.names = FALSE
  )

  # No fibroblast subtype labels are assigned in this public script.
}
