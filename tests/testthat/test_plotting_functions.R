source(testthat::test_path("..", "..", "R", "functions", "plotting_functions.R"))

test_that("gene_has_assay_signal returns FALSE for a missing assay", {
  skip_if_not_installed("Seurat")

  obj <- Seurat::CreateSeuratObject(
    counts = matrix(
      c(1, 0, 2, 0),
      nrow = 2,
      dimnames = list(
        c("GeneA", "GeneB"),
        c("Cell1", "Cell2")
      )
    )
  )

  expect_false(
    gene_has_assay_signal(
      obj,
      gene = "GeneA",
      assay = "not_an_assay"
    )
  )
})
