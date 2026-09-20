
# Generalized pathway enrichment for Seurat marker tables.

library(dplyr)
library(ggplot2)
library(stringr)
library(enrichR)

# ------------------------------------------------------------
# Pathway-library presets and discovery helpers
# ------------------------------------------------------------

# Default libraries used by this project. These are kept in one place so they
# can be changed without editing the enrichment function itself.
DEFAULT_ENRICHR_PATHWAY_DATABASES <- c(
  "GO_Biological_Process_2025",
  "KEGG_2026",
  "Reactome_Pathways_2024"
)

# Optional broader preset. Edit this vector if you want to add or remove
# pathway collections for a particular analysis.
EXTENDED_ENRICHR_PATHWAY_DATABASES <- c(
  DEFAULT_ENRICHR_PATHWAY_DATABASES,
  "WikiPathways_2024_Human"
)


list_enrichr_pathway_libraries <- function(
    keywords = c(
      "GO_Biological_Process",
      "KEGG",
      "Reactome",
      "WikiPathway"
    )
) {
  dbs <- enrichR::listEnrichrDbs()

  keep <- Reduce(
    `|`,
    lapply(
      keywords,
      function(keyword) {
        grepl(
          keyword,
          dbs$libraryName,
          ignore.case = TRUE
        )
      }
    )
  )

  dbs[keep, , drop = FALSE]
}


suggest_latest_enrichr_pathway_libraries <- function(
    families = c(
      "GO_Biological_Process",
      "KEGG",
      "Reactome_Pathways",
      "WikiPathways"
    )
) {
  dbs <- enrichR::listEnrichrDbs()

  picks <- lapply(
    families,
    function(family) {
      candidates <- dbs %>%
        filter(
          grepl(
            family,
            libraryName,
            ignore.case = TRUE
          )
        ) %>%
        mutate(
          year = suppressWarnings(
            as.integer(
              stringr::str_extract(
                libraryName,
                "(?<!\\d)(20\\d{2})(?!\\d)"
              )
            )
          )
        )

      if (nrow(candidates) == 0) {
        return(NULL)
      }

      # Prefer the newest year-tagged version. If no year is encoded,
      # return the first current library name reported by Enrichr.
      if (any(!is.na(candidates$year))) {
        candidates <- candidates %>%
          arrange(
            desc(year),
            libraryName
          )
      }

      candidates %>%
        slice_head(n = 1) %>%
        mutate(
          Family = family
        )
    }
  )

  bind_rows(picks) %>%
    select(
      Family,
      libraryName,
      everything()
    )
}


clean_enrichment_term <- function(term, database = NULL) {
  cleaned <- stringr::str_remove(
    term,
    "\\s*\\(GO:\\d+\\)$"
  )

  if (!is.null(database)) {
    is_kegg <- database == "KEGG"
    cleaned[is_kegg] <- stringr::str_to_sentence(
      stringr::str_to_lower(cleaned[is_kegg])
    )
  }

  stringr::str_replace_all(
    cleaned,
    c(
      "\\bEcm\\b" = "ECM",
      "\\bTgf\\b" = "TGF",
      "\\bMapk\\b" = "MAPK",
      "\\bPi3k\\b" = "PI3K",
      "\\bAkt\\b" = "AKT",
      "\\bJak\\b" = "JAK",
      "\\bStat\\b" = "STAT",
      "\\bWnt\\b" = "WNT",
      "\\bBmp\\b" = "BMP",
      "\\bVegf\\b" = "VEGF",
      "\\bHif\\b" = "HIF"
    )
  )
}


plot_pathway_enrichment_dotplot <- function(
    enrichment_table,
    group_order = NULL
) {
  plot_df <- enrichment_table

  if (!is.null(group_order)) {
    plot_df$Group <- factor(
      plot_df$Group,
      levels = group_order
    )
  }

  ggplot(
    plot_df,
    aes(
      x = Group,
      y = Term_clean
    )
  ) +
    geom_point(
      aes(
        size = nGenes,
        color = neg_log10_significance
      )
    ) +
    scale_size_continuous(name = "Number of genes") +
    scale_color_gradient(
      low = "grey80",
      high = "navy",
      name = expression(-log[10](P))
    ) +
    labs(x = "", y = "") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}


run_marker_pathway_enrichment <- function(
    markers,
    group_col = "cluster",
    group_label = "Markers",
    gene_col = "gene",
    logfc_col = "avg_log2FC",
    marker_padj_col = "p_val_adj",
    marker_padj_cutoff = 0.01,
    top_n_genes = 100,
    only_positive = TRUE,
    databases = DEFAULT_ENRICHR_PATHWAY_DATABASES,
    database_labels = c(
      "GO_Biological_Process_2025" = "GO Biological Process",
      "KEGG_2026" = "KEGG",
      "Reactome_Pathways_2024" = "Reactome",
      "WikiPathways_2024_Human" = "WikiPathways"
    ),
    enrichment_significance_col = "Adjusted.P.value",
    enrichment_p_cutoff = 0.05,
    top_n_terms = 5,
    group_order = NULL
) {
  marker_df <- markers

  # FindMarkers() often stores gene names in row names.
  if (!gene_col %in% colnames(marker_df)) {
    marker_df[[gene_col]] <- rownames(marker_df)
  }

  required_cols <- c(
    gene_col,
    logfc_col,
    marker_padj_col
  )

  missing_cols <- setdiff(required_cols, colnames(marker_df))
  if (length(missing_cols) > 0) {
    stop(
      "Marker table is missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # FindMarkers() has no cluster/group column.
  if (!group_col %in% colnames(marker_df)) {
    marker_df[[group_col]] <- group_label
  }

  marker_sig <- marker_df %>%
    filter(
      !is.na(.data[[marker_padj_col]]),
      .data[[marker_padj_col]] < marker_padj_cutoff
    )

  if (only_positive) {
    marker_sig <- marker_sig %>%
      filter(.data[[logfc_col]] > 0)
  }

  if (nrow(marker_sig) == 0) {
    stop("No markers passed the requested marker filters.")
  }

  top_markers <- marker_sig %>%
    group_by(.data[[group_col]]) %>%
    arrange(desc(.data[[logfc_col]]), .by_group = TRUE) %>%
    slice_head(n = top_n_genes) %>%
    ungroup()

  groups <- unique(as.character(top_markers[[group_col]]))
  enrichment_list <- list()

  for (group_name in groups) {
    genes <- top_markers %>%
      filter(.data[[group_col]] == group_name) %>%
      pull(.data[[gene_col]]) %>%
      unique()

    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0) next

    enrichment <- enrichR::enrichr(
      genes,
      databases
    )

    for (database in databases) {
      if (!database %in% names(enrichment)) next

      database_result <- enrichment[[database]]
      if (nrow(database_result) == 0) next

      enrichment_list[[
        paste(group_name, database, sep = "__")
      ]] <- database_result %>%
        mutate(
          Database = ifelse(
            database %in% names(database_labels),
            database_labels[[database]],
            database
          ),
          Group = group_name,
          SourceDatabase = database
        )
    }
  }

  if (length(enrichment_list) == 0) {
    stop("Enrichr returned no enrichment results.")
  }

  enrichment_all <- bind_rows(enrichment_list)

  if (!enrichment_significance_col %in% colnames(enrichment_all)) {
    stop(
      "Enrichment significance column not found: ",
      enrichment_significance_col
    )
  }

  enrichment_sig <- enrichment_all %>%
    filter(
      !is.na(.data[[enrichment_significance_col]]),
      .data[[enrichment_significance_col]] < enrichment_p_cutoff
    ) %>%
    mutate(
      neg_log10_significance = -log10(
        pmax(
          .data[[enrichment_significance_col]],
          .Machine$double.xmin
        )
      ),
      nGenes = as.numeric(sub("/.*", "", Overlap))
    )

  if (nrow(enrichment_sig) == 0) {
    warning("No enrichment terms passed the enrichment cutoff.")
    return(list(
      top_markers = top_markers,
      enrichment_all = enrichment_all,
      enrichment_significant = enrichment_sig,
      top_enrichment = enrichment_sig,
      plot = NULL
    ))
  }

  top_enrichment <- enrichment_sig %>%
    group_by(Group, Database) %>%
    arrange(.data[[enrichment_significance_col]], .by_group = TRUE) %>%
    slice_head(n = top_n_terms) %>%
    ungroup() %>%
    mutate(
      Term_clean = clean_enrichment_term(
        Term,
        Database
      )
    )

  if (is.null(group_order)) {
    group_order <- groups
  }

  top_enrichment$Group <- factor(
    top_enrichment$Group,
    levels = group_order
  )

  term_order <- top_enrichment %>%
    arrange(
      Group,
      .data[[enrichment_significance_col]]
    ) %>%
    pull(Term_clean) %>%
    unique()

  top_enrichment$Term_clean <- factor(
    top_enrichment$Term_clean,
    levels = rev(term_order)
  )

  plot_obj <- plot_pathway_enrichment_dotplot(
    top_enrichment,
    group_order = group_order
  )

  list(
    top_markers = top_markers,
    enrichment_all = enrichment_all,
    enrichment_significant = enrichment_sig,
    top_enrichment = top_enrichment,
    plot = plot_obj
  )
}


# Example: FindAllMarkers()
#
# enrichment <- run_marker_pathway_enrichment(
#   markers = marker_table,
#   group_col = "cluster",
#   top_n_genes = 100
# )
#
# enrichment$plot
# enrichment$top_enrichment
#
# Example: FindMarkers()
#
# enrichment <- run_marker_pathway_enrichment(
#   markers = one_comparison_markers,
#   group_label = "cluster_A_vs_B"
# )


# ------------------------------------------------------------
# Inspect or change the pathway libraries
# ------------------------------------------------------------
#
# See pathway-related libraries currently reported by Enrichr:
#
# available_pathway_dbs <- list_enrichr_pathway_libraries()
# View(available_pathway_dbs)
#
# Ask the helper to pick the newest year-tagged library for each family:
#
# latest_pathway_dbs <- suggest_latest_enrichr_pathway_libraries()
# latest_pathway_dbs$libraryName
#
# Use a custom set for one analysis:
#
# enrichment <- run_marker_pathway_enrichment(
#   markers = marker_table,
#   databases = c(
#     "GO_Biological_Process_2025",
#     "KEGG_2026"
#   )
# )
#
# Or use whatever the live Enrichr catalog reports as newest:
#
# latest <- suggest_latest_enrichr_pathway_libraries()
#
# enrichment <- run_marker_pathway_enrichment(
#   markers = marker_table,
#   databases = latest$libraryName
# )
