source(testthat::test_path("..", "..", "R", "functions", "preprocessing_functions.R"))
source(testthat::test_path("..", "..", "R", "functions", "fibroblast_functions.R"))


make_test_seurat <- function() {
  counts <- matrix(
    c(
      5, 0, 1, 0,
      2, 3, 0, 1,
      0, 1, 4, 0,
      1, 2, 1, 3
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("GeneA", "GeneB", "GeneC", "MT-Gene"),
      c("Cell1", "Cell2", "Cell3", "Cell4")
    )
  )

  obj <- Seurat::CreateSeuratObject(counts = counts)

  obj$nFeature_RNA <- c(100, 500, 700, 7000)
  obj$nCount_RNA <- c(500, 1000, 1500, 9000)
  obj$percent.mt <- c(1, 5, 30, 2)

  obj
}


test_that("filter_cells applies feature, count, and mitochondrial thresholds", {
  obj <- make_test_seurat()

  filtered <- filter_cells(
    obj,
    min_features = 200,
    max_features = 6000,
    min_counts = 600,
    max_counts = 5000,
    max_percent_mt = 20
  )

  expect_equal(colnames(filtered), "Cell2")
})


test_that("apply_cluster_annotations assigns labels by cluster", {
  obj <- make_test_seurat()
  obj$seurat_clusters <- factor(c("0", "0", "1", "1"))

  annotations <- data.frame(
    cluster = c("0", "1"),
    cell_type = c("Fibroblasts", "Epithelial"),
    stringsAsFactors = FALSE
  )

  annotated <- apply_cluster_annotations(
    obj,
    annotation_df = annotations
  )

  expect_equal(
    unname(annotated$cell_type),
    c("Fibroblasts", "Fibroblasts", "Epithelial", "Epithelial")
  )
})


test_that("apply_cluster_annotations fails when a cluster is missing", {
  obj <- make_test_seurat()
  obj$seurat_clusters <- factor(c("0", "0", "1", "1"))

  incomplete <- data.frame(
    cluster = "0",
    cell_type = "Fibroblasts",
    stringsAsFactors = FALSE
  )

  expect_error(
    apply_cluster_annotations(obj, incomplete),
    "Missing annotations"
  )
})


test_that("subset_cell_type retains only requested cells", {
  obj <- make_test_seurat()
  obj$cell_type <- c(
    "Fibroblasts",
    "Fibroblasts",
    "Epithelial",
    "Endothelial"
  )

  fib <- subset_cell_type(
    obj,
    cell_type = "Fibroblasts"
  )

  expect_equal(colnames(fib), c("Cell1", "Cell2"))
})


test_that("subset_cell_type fails for an absent population", {
  obj <- make_test_seurat()
  obj$cell_type <- rep("Epithelial", ncol(obj))

  expect_error(
    subset_cell_type(obj, "Fibroblasts"),
    "No cells found"
  )
})


test_that("check_cell_type_annotations validates metadata", {
  obj <- make_test_seurat()

  expect_error(
    check_cell_type_annotations(obj),
    "does not contain"
  )

  obj$cell_type <- rep("Fibroblasts", ncol(obj))
  expect_true(check_cell_type_annotations(obj))
})
