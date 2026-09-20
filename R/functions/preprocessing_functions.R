# ============================================================
# Shared preprocessing functions
# Cross-species intestinal fibroblast sc/snRNA-seq analysis
# ============================================================

library(Seurat)
library(dplyr)


# ------------------------------------------------------------
# Add mitochondrial percentage
# ------------------------------------------------------------

add_mito_percent <- function(
    object,
    pattern = NULL,
    features = NULL,
    column_name = "percent.mt"
) {
  if (is.null(pattern) && is.null(features)) {
    stop("Provide either `pattern` or `features`.")
  }

  object[[column_name]] <- PercentageFeatureSet(
    object,
    pattern = pattern,
    features = features
  )

  object
}


# ------------------------------------------------------------
# Quality-control filtering
# ------------------------------------------------------------

filter_cells <- function(
    object,
    min_features = 200,
    max_features = Inf,
    min_counts = 0,
    max_counts = Inf,
    max_percent_mt = Inf,
    mt_column = "percent.mt"
) {
  meta <- object[[]]

  keep <- meta$nFeature_RNA >= min_features &
    meta$nFeature_RNA <= max_features &
    meta$nCount_RNA >= min_counts &
    meta$nCount_RNA <= max_counts

  if (mt_column %in% colnames(meta) && is.finite(max_percent_mt)) {
    keep <- keep & meta[[mt_column]] <= max_percent_mt
  }

  subset(
    object,
    cells = rownames(meta)[keep]
  )
}


# ------------------------------------------------------------
# Doublet detection with scDblFinder
# ------------------------------------------------------------
# Run after basic QC and before final normalization/clustering.
# scDblFinder uses raw counts and performs its own internal
# normalization/dimensional reduction.

remove_doublets <- function(
    object,
    assay = "RNA",
    seed = 1234
) {
  if (!assay %in% Assays(object)) {
    stop("Assay `", assay, "` is not present in the Seurat object.")
  }

  if (!"counts" %in% Layers(object[[assay]])) {
    stop(
      "Raw counts are required for scDblFinder, but the `",
      assay,
      "` assay does not contain a counts layer."
    )
  }

  set.seed(seed)

  sce <- Seurat::as.SingleCellExperiment(
    object,
    assay = assay
  )

  sce <- scDblFinder::scDblFinder(sce)

  object$scDblFinder.class <- sce$scDblFinder.class
  object$scDblFinder.score <- sce$scDblFinder.score

  subset(
    object,
    subset = scDblFinder.class == "singlet"
  )
}


# ------------------------------------------------------------
# SCTransform + PCA + graph clustering + UMAP
# ------------------------------------------------------------

run_sct_clustering <- function(
    object,
    dims = 1:20,
    resolution = 0.8,
    regress_vars = NULL,
    seed = 1234
) {
  if (!"RNA" %in% Assays(object)) {
    stop("An RNA assay is required for SCTransform.")
  }

  set.seed(seed)
  DefaultAssay(object) <- "RNA"

  object <- SCTransform(
    object,
    assay = "RNA",
    method = "glmGamPoi",
    vars.to.regress = regress_vars,
    verbose = FALSE
  )

  object <- RunPCA(
    object,
    verbose = FALSE
  )

  object <- FindNeighbors(
    object,
    dims = dims,
    verbose = FALSE
  )

  object <- FindClusters(
    object,
    resolution = resolution,
    verbose = FALSE
  )

  object <- RunUMAP(
    object,
    dims = dims,
    seed.use = seed,
    verbose = FALSE
  )

  object
}


# ------------------------------------------------------------
# Harmony batch correction + clustering
# ------------------------------------------------------------
# Assumes PCA has already been calculated.
# Use for true multi-sample datasets when batch correction is
# justified (e.g. biological replicates from the same condition).

run_harmony_clustering <- function(
    object,
    batch_var = "orig.ident",
    dims = 1:20,
    resolution = 0.8,
    seed = 1234
) {
  if (!batch_var %in% colnames(object[[]])) {
    stop("Batch variable `", batch_var, "` is not present in metadata.")
  }

  if (!"pca" %in% Reductions(object)) {
    stop("Run PCA before calling run_harmony_clustering().")
  }

  set.seed(seed)

  object <- harmony::RunHarmony(
    object,
    group.by.vars = batch_var,
    reduction = "pca",
    plot_convergence = FALSE
  )

  object <- FindNeighbors(
    object,
    reduction = "harmony",
    dims = dims,
    verbose = FALSE
  )

  object <- FindClusters(
    object,
    resolution = resolution,
    verbose = FALSE
  )

  object <- RunUMAP(
    object,
    reduction = "harmony",
    dims = dims,
    seed.use = seed,
    verbose = FALSE
  )

  object
}


# ------------------------------------------------------------
# Cluster marker discovery
# ------------------------------------------------------------

find_cluster_markers <- function(
    object,
    assay = "SCT",
    only_pos = TRUE,
    min_log2fc = 0,
    min_pct_diff = NULL
) {
  if (!assay %in% Assays(object)) {
    stop("Assay `", assay, "` is not present in the object.")
  }

  markers <- FindAllMarkers(
    object,
    assay = assay,
    only.pos = only_pos
  )

  markers <- markers %>%
    group_by(cluster) %>%
    arrange(cluster, desc(avg_log2FC))

  if (!is.null(min_log2fc)) {
    markers <- markers %>%
      filter(avg_log2FC > min_log2fc)
  }

  if (!is.null(min_pct_diff)) {
    markers <- markers %>%
      filter((pct.1 - pct.2) > min_pct_diff)
  }

  markers
}


# ------------------------------------------------------------
# Apply a private cluster annotation table
# ------------------------------------------------------------

apply_cluster_annotations <- function(
    object,
    annotation_df,
    cluster_col = "seurat_clusters",
    label_col = "cell_type"
) {
  if (!cluster_col %in% colnames(object[[]])) {
    stop("Cluster column `", cluster_col, "` is not present in metadata.")
  }

  required_cols <- c("cluster", label_col)

  if (!all(required_cols %in% colnames(annotation_df))) {
    stop(
      "Annotation table must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  Idents(object) <- cluster_col

  annotation_map <- setNames(
    as.character(annotation_df[[label_col]]),
    as.character(annotation_df$cluster)
  )

  missing_clusters <- setdiff(
    levels(Idents(object)),
    names(annotation_map)
  )

  if (length(missing_clusters) > 0) {
    stop(
      "Missing annotations for cluster(s): ",
      paste(missing_clusters, collapse = ", ")
    )
  }

  object <- RenameIdents(
    object,
    annotation_map
  )

  object[[label_col]] <- as.character(Idents(object))

  object
}


# ------------------------------------------------------------
# Subset a named cell population
# ------------------------------------------------------------

subset_cell_type <- function(
    object,
    cell_type,
    metadata_col = "cell_type"
) {
  if (!metadata_col %in% colnames(object[[]])) {
    stop("Metadata column `", metadata_col, "` is not present.")
  }

  cells_keep <- rownames(object[[]])[
    object[[]][[metadata_col]] == cell_type
  ]

  if (length(cells_keep) == 0) {
    stop(
      "No cells found with ",
      metadata_col,
      " == '",
      cell_type,
      "'."
    )
  }

  subset(
    object,
    cells = cells_keep
  )
}


# ------------------------------------------------------------
# Plot only markers present in a dataset
# ------------------------------------------------------------

plot_available_markers <- function(
    object,
    markers,
    reduction = "umap",
    order = TRUE
) {
  markers_present <- intersect(
    markers,
    rownames(object)
  )

  markers_missing <- setdiff(
    markers,
    markers_present
  )

  if (length(markers_present) == 0) {
    stop("None of the requested markers are present in the object.")
  }

  if (length(markers_missing) > 0) {
    message(
      "Markers not present and skipped: ",
      paste(markers_missing, collapse = ", ")
    )
  }

  FeaturePlot(
    object,
    features = markers_present,
    reduction = reduction,
    order = order
  )
}
