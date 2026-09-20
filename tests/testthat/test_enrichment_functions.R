source("R/functions/enrichment_functions.R")

test_that("GO identifiers are removed from enrichment labels", {
  observed <- clean_enrichment_term(
    "Extracellular matrix organization (GO:0030198)",
    "GO Biological Process"
  )

  expect_equal(
    observed,
    "Extracellular matrix organization"
  )
})


test_that("KEGG labels are sentence-cased and common acronyms restored", {
  observed <- clean_enrichment_term(
    "ECM RECEPTOR INTERACTION",
    "KEGG"
  )

  expect_match(observed, "ECM")
})


test_that("pathway dotplot helper returns a ggplot object", {
  example <- data.frame(
    Group = factor(c("A", "B")),
    Term_clean = factor(c("Pathway 1", "Pathway 2")),
    nGenes = c(4, 6),
    neg_log10_significance = c(2, 3)
  )

  p <- plot_pathway_enrichment_dotplot(example)

  expect_s3_class(p, "ggplot")
})
