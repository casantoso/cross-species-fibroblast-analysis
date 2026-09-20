
# Reusable plotting helpers for single-cell and cross-species analyses.

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(purrr)

plot_feature_pages_with_dimplot <- function(
    object,
    genes,
    page_title_prefix = "Feature plots",
    assay = NULL,
    reduction = "umap",
    genes_per_page = 8,
    ncol = 3,
    pt.size = NULL,
    label_clusters = TRUE,
    order = TRUE
) {
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.")

  genes <- unique(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (!is.null(assay)) {
    if (!assay %in% Assays(object)) stop("Assay not found: ", assay)
    DefaultAssay(object) <- assay
  }

  genes_present <- intersect(genes, rownames(object[[DefaultAssay(object)]]))
  genes_missing <- setdiff(genes, genes_present)

  if (length(genes_present) == 0) {
    warning("None of the requested genes are present in the object.")
    return(invisible(list(
      pages = list(),
      genes_present = character(0),
      genes_missing = genes_missing
    )))
  }

  gene_chunks <- split(
    genes_present,
    ceiling(seq_along(genes_present) / genes_per_page)
  )

  pages <- vector("list", length(gene_chunks))

  for (i in seq_along(gene_chunks)) {
    genes_chunk <- gene_chunks[[i]]

    p_dim <- DimPlot(
      object,
      reduction = reduction,
      label = label_clusters
    ) +
      ggtitle(paste0(page_title_prefix, " — page ", i)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    feature_plots <- lapply(genes_chunk, function(gene) {
      args <- list(
        object = object,
        features = gene,
        reduction = reduction,
        order = order
      )
      if (!is.null(pt.size)) args$pt.size <- pt.size

      do.call(FeaturePlot, args) +
        ggtitle(gene) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    })

    panel_list <- c(list(p_dim), feature_plots)
    n_target <- genes_per_page + 1
    n_blank <- n_target - length(panel_list)

    if (n_blank > 0) {
      panel_list <- c(
        panel_list,
        replicate(n_blank, patchwork::plot_spacer(), simplify = FALSE)
      )
    }

    pages[[i]] <- wrap_plots(panel_list, ncol = ncol)
  }

  invisible(list(
    pages = pages,
    genes_present = genes_present,
    genes_missing = genes_missing
  ))
}


make_featureplot_pdf <- function(
    object,
    genes,
    output_file,
    assay = "SCT",
    reduction = "umap",
    ncol = 4,
    n_per_page = 12,
    pt.size = 0.3,
    order = TRUE,
    width = 16,
    height = 12
) {
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.")
  if (!assay %in% Assays(object)) stop("Assay not found: ", assay)

  DefaultAssay(object) <- assay

  genes <- unique(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  genes_present <- intersect(genes, rownames(object[[assay]]))
  genes_missing <- setdiff(genes, genes_present)

  if (length(genes_present) == 0) {
    stop("None of the requested genes were found in the object.")
  }

  gene_chunks <- split(
    genes_present,
    ceiling(seq_along(genes_present) / n_per_page)
  )

  grDevices::pdf(output_file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (genes_chunk in gene_chunks) {
    plot_list <- FeaturePlot(
      object,
      features = genes_chunk,
      reduction = reduction,
      pt.size = pt.size,
      order = order,
      combine = FALSE
    )

    plot_list <- lapply(plot_list, function(p) {
      p + theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
    })

    n_blank <- n_per_page - length(plot_list)
    if (n_blank > 0) {
      plot_list <- c(
        plot_list,
        replicate(n_blank, ggplot() + theme_void(), simplify = FALSE)
      )
    }

    print(
      wrap_plots(
        plot_list,
        ncol = ncol,
        nrow = ceiling(n_per_page / ncol)
      )
    )
  }

  invisible(list(
    genes_present = genes_present,
    genes_missing = genes_missing
  ))
}


gene_has_assay_signal <- function(
    object,
    gene,
    assay = "SCT",
    layer = "data",
    min_cells = 5
) {
  if (!assay %in% Assays(object)) return(FALSE)
  if (!gene %in% rownames(object[[assay]])) return(FALSE)

  values <- GetAssayData(
    object,
    assay = assay,
    layer = layer
  )[gene, ]

  values <- values[!is.na(values)]
  if (length(values) == 0) return(FALSE)

  sum(values > 0) >= min_cells
}


make_paralog_panels <- function(
    anchor_gene,
    mapping_df,
    object_list,
    species_order = names(object_list),
    anchor_col = "anchor_gene",
    gene_col = "actual_gene",
    species_col = "species",
    assay = "SCT",
    layer = "data",
    reduction = "umap",
    min_cells = 5,
    order = TRUE
) {
  required_cols <- c(anchor_col, gene_col, species_col)
  missing_cols <- setdiff(required_cols, colnames(mapping_df))
  if (length(missing_cols) > 0) {
    stop("Missing mapping columns: ", paste(missing_cols, collapse = ", "))
  }

  gene_map <- mapping_df %>%
    filter(
      .data[[anchor_col]] == anchor_gene,
      !is.na(.data[[gene_col]]),
      nzchar(.data[[gene_col]])
    )

  panels <- list()

  for (species in species_order) {
    if (!species %in% names(object_list)) next

    rows_species <- gene_map %>%
      filter(.data[[species_col]] == species)

    if (nrow(rows_species) == 0) next

    object <- object_list[[species]]

    for (i in seq_len(nrow(rows_species))) {
      gene <- rows_species[[gene_col]][i]

      if (!gene_has_assay_signal(
        object,
        gene,
        assay = assay,
        layer = layer,
        min_cells = min_cells
      )) next

      DefaultAssay(object) <- assay

      panels[[length(panels) + 1]] <- FeaturePlot(
        object,
        features = gene,
        reduction = reduction,
        order = order
      ) +
        ggtitle(paste0(species, "\n", gene)) +
        theme(
          plot.title = element_text(
            hjust = 0.5,
            face = "bold",
            size = 9
          )
        )
    }
  }

  panels
}


make_gene_pages <- function(
    anchor_gene,
    mapping_df,
    object_list,
    species_order = names(object_list),
    panels_per_page = 9,
    ncol = 3,
    ...
) {
  panels <- make_paralog_panels(
    anchor_gene = anchor_gene,
    mapping_df = mapping_df,
    object_list = object_list,
    species_order = species_order,
    ...
  )

  if (length(panels) == 0) {
    panels <- list(
      ggplot() +
        theme_void() +
        ggtitle(paste0(anchor_gene, "\n(no species passed filter)")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    )
  }

  chunks <- split(
    panels,
    ceiling(seq_along(panels) / panels_per_page)
  )

  pages <- vector("list", length(chunks))

  for (i in seq_along(chunks)) {
    chunk <- chunks[[i]]
    n_blank <- panels_per_page - length(chunk)

    if (n_blank > 0) {
      chunk <- c(
        chunk,
        replicate(n_blank, patchwork::plot_spacer(), simplify = FALSE)
      )
    }

    pages[[i]] <- wrap_plots(chunk, ncol = ncol) +
      plot_annotation(
        title = anchor_gene,
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
        )
      )
  }

  pages
}


prepare_cross_species_marker_heatmap <- function(
    marker_list,
    mapping_df,
    species_order = names(marker_list),
    target_cluster = "Fibroblasts",
    cluster_col = "cluster",
    gene_col = "gene",
    logfc_col = "avg_log2FC",
    anchor_col = "anchor_gene",
    mapped_gene_col = "actual_gene",
    species_col = "species"
) {
  if (is.null(names(marker_list))) stop("`marker_list` must be a named list.")

  scaled_marker_list <- purrr::imap(marker_list, function(markers_species, species) {
    required_cols <- c(cluster_col, gene_col, logfc_col)
    missing_cols <- setdiff(required_cols, colnames(markers_species))

    if (length(missing_cols) > 0) {
      stop(
        "Marker table for ", species, " is missing: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    markers_species %>%
      filter(.data[[cluster_col]] == target_cluster) %>%
      mutate(
        scaled_log2fc = {
          values <- .data[[logfc_col]]
          s <- stats::sd(values, na.rm = TRUE)

          if (is.na(s) || s == 0) {
            rep(0, dplyr::n())
          } else {
            as.numeric(scale(values))
          }
        },
        !!species_col := species
      ) %>%
      select(
        !!species_col,
        !!gene_col,
        scaled_log2fc
      )
  })

  scaled_markers <- bind_rows(scaled_marker_list)

  joined <- mapping_df %>%
    left_join(
      scaled_markers,
      by = setNames(
        c(species_col, gene_col),
        c(species_col, mapped_gene_col)
      )
    ) %>%
    group_by(
      .data[[anchor_col]],
      .data[[species_col]]
    ) %>%
    arrange(desc(scaled_log2fc), .by_group = TRUE) %>%
    mutate(paralog_index = row_number()) %>%
    ungroup() %>%
    mutate(
      row_id = paste(
        .data[[anchor_col]],
        paralog_index,
        sep = "__"
      )
    )

  all_rows <- unique(joined$row_id)

  heatmap_df <- tidyr::expand_grid(
    row_id = all_rows,
    species = species_order
  ) %>%
    left_join(
      joined %>%
        transmute(
          row_id = row_id,
          anchor_gene = .data[[anchor_col]],
          species = .data[[species_col]],
          actual_gene = .data[[mapped_gene_col]],
          scaled_log2fc = scaled_log2fc
        ),
      by = c("row_id", "species")
    )

  heatmap_df$row_id <- factor(
    heatmap_df$row_id,
    levels = rev(all_rows)
  )
  heatmap_df$species <- factor(
    heatmap_df$species,
    levels = species_order
  )

  heatmap_df
}


plot_cross_species_marker_heatmap <- function(
    heatmap_df,
    title = "Marker genes across species",
    na_fill = "grey95"
) {
  ggplot(
    heatmap_df,
    aes(
      x = species,
      y = row_id,
      fill = scaled_log2fc
    )
  ) +
    geom_tile(color = "white") +
    geom_text(
      aes(
        label = ifelse(
          is.na(actual_gene),
          "",
          actual_gene
        )
      ),
      size = 3.5,
      na.rm = TRUE
    ) +
    scale_y_discrete(
      labels = function(x) sub("__.*$", "", x)
    ) +
    scale_fill_gradient2(
      midpoint = 0,
      na.value = na_fill,
      name = "Scaled\nlog2FC"
    ) +
    labs(
      title = title,
      x = "Species",
      y = "Anchor gene"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}
