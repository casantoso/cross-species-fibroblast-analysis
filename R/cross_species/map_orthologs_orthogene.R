# ============================================================
# Pairwise and multi-species ortholog mapping with orthogene
#
# Purpose:
#   Convert a gene list from one species to one or more other species.
#
# User-facing functions accept readable species names such as:
#   "human"
#   "mouse"
#   "pig"
#   "chicken"
#   "zebrafish"
#   "tilapia"
#
# These are translated internally to backend organism identifiers.
# ============================================================


SPECIES_SCI <- c(
  human     = "hsapiens",
  mouse     = "mmusculus",
  pig       = "sscrofa",
  chicken   = "ggallus",
  zebrafish = "drerio",
  tilapia   = "oniloticus",
  macaque   = "mfascicularis"
)


resolve_orthogene_species <- function(
    species,
    species_map = SPECIES_SCI
) {
  species <- as.character(species)

  if (species %in% names(species_map)) {
    return(unname(species_map[[species]]))
  }

  if (species %in% unname(species_map)) {
    return(species)
  }

  stop(
    "Unknown species: ", species, ". ",
    "Use one of: ",
    paste(names(species_map), collapse = ", "),
    "; or provide a supported backend identifier."
  )
}


check_orthogene_species <- function(
    species,
    method = "gprofiler"
) {
  if (!requireNamespace("orthogene", quietly = TRUE)) {
    stop(
      "Package `orthogene` is required. Install it with ",
      "BiocManager::install('orthogene')."
    )
  }

  species_backend <- resolve_orthogene_species(species)

  orthogene::map_species(
    species = species_backend,
    method = method
  )
}


.standardize_orthogene_output <- function(
    result,
    input_species,
    output_species
) {
  result <- as.data.frame(result)

  if (ncol(result) < 2) {
    stop(
      "orthogene returned fewer than two columns; ",
      "cannot identify input and ortholog columns."
    )
  }

  # Keep the first two gene columns returned by orthogene and rename them
  # to the readable species names used by this project.
  result <- result[, 1:2, drop = FALSE]
  colnames(result) <- c(
    input_species,
    output_species
  )

  result
}


map_gene_list_orthologs <- function(
    genes,
    input_species,
    output_species,
    method = "gprofiler",
    one_to_one_only = TRUE,
    drop_nonorths = TRUE
) {
  if (!requireNamespace("orthogene", quietly = TRUE)) {
    stop(
      "Package `orthogene` is required. Install it with ",
      "BiocManager::install('orthogene')."
    )
  }

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (length(genes) == 0) {
    stop("`genes` contains no valid gene identifiers.")
  }

  input_backend <- resolve_orthogene_species(input_species)
  output_backend <- resolve_orthogene_species(output_species)

  gene_table <- data.frame(
    gene = genes,
    stringsAsFactors = FALSE
  )

  non121_strategy <- if (one_to_one_only) {
    "drop_both_species"
  } else {
    "keep_both_species"
  }

  result <- orthogene::convert_orthologs(
    gene_df = gene_table,
    gene_input = "gene",
    gene_output = "columns",
    input_species = input_backend,
    output_species = output_backend,
    method = method,
    drop_nonorths = drop_nonorths,
    non121_strategy = non121_strategy
  )

  .standardize_orthogene_output(
    result = result,
    input_species = input_species,
    output_species = output_species
  )
}


map_gene_list_multi_species <- function(
    genes,
    input_species,
    output_species,
    method = "gprofiler",
    one_to_one_only = TRUE,
    drop_nonorths = TRUE
) {
  if (length(output_species) == 0) {
    stop("`output_species` must contain at least one species.")
  }

  output_species <- unique(as.character(output_species))

  base <- data.frame(
    gene = unique(as.character(genes)),
    stringsAsFactors = FALSE
  )

  base <- base[
    !is.na(base$gene) & nzchar(base$gene),
    ,
    drop = FALSE
  ]

  colnames(base) <- input_species

  for (species in output_species) {

    mapped <- map_gene_list_orthologs(
      genes = base[[input_species]],
      input_species = input_species,
      output_species = species,
      method = method,
      one_to_one_only = one_to_one_only,
      drop_nonorths = drop_nonorths
    )

    mapped <- mapped[
      !duplicated(mapped[[input_species]]),
      ,
      drop = FALSE
    ]

    base <- merge(
      base,
      mapped,
      by = input_species,
      all.x = TRUE,
      sort = FALSE
    )
  }

  # Restore the input gene order after successive merges.
  base <- base[
    match(
      unique(as.character(genes)),
      base[[input_species]]
    ),
    ,
    drop = FALSE
  ]

  rownames(base) <- NULL

  base
}


# ------------------------------------------------------------
# Examples
# ------------------------------------------------------------

# Human -> mouse, strict 1:1 mappings
#
# mouse_map <- map_gene_list_orthologs(
#   genes = c("PDGFRA", "COL1A1", "DCN"),
#   input_species = "human",
#   output_species = "mouse",
#   one_to_one_only = TRUE
# )
#
# Output columns:
#   human | mouse


# Human -> zebrafish, preserve one-to-many mappings
#
# zebrafish_map <- map_gene_list_orthologs(
#   genes = c("PDGFRA", "COL1A1", "DCN"),
#   input_species = "human",
#   output_species = "zebrafish",
#   one_to_one_only = FALSE
# )
#
# Important zebrafish mappings should be manually checked in NCBI Gene,
# especially when several paralogs are returned or an expected gene is absent.


# Human -> Nile tilapia, strict 1:1 mappings
#
# tilapia_map <- map_gene_list_orthologs(
#   genes = c("PDGFRA", "COL1A1", "DCN"),
#   input_species = "human",
#   output_species = "tilapia",
#   one_to_one_only = TRUE
# )
#
# Important tilapia mappings should be manually checked in NCBI Gene.
# The tilapia single-cell dataset used in this project contains many genes
# labeled with provisional LOC... identifiers rather than informative symbols,
# so LOC... genes may need manual annotation/ortholog lookup before they are
# interpreted as unmapped or absent.


# Human -> several species in one table
#
# multi_species_map <- map_gene_list_multi_species(
#   genes = c("PDGFRA", "COL1A1", "DCN"),
#   input_species = "human",
#   output_species = c(
#     "mouse",
#     "pig",
#     "chicken",
#     "zebrafish",
#     "tilapia"
#   ),
#   one_to_one_only = TRUE
# )
#
# Output columns:
#   human | mouse | pig | chicken | zebrafish | tilapia
