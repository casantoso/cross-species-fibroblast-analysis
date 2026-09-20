# ============================================================
# Prepare broad cell-type Seurat objects for SAMap
#
# Purpose
# -------
# The broad-cell-type SAMap analysis was used to test whether the same
# major intestinal cell classes align across species before focusing on
# fibroblasts.
#
# Retained broad cell types:
#   Fibroblasts
#   SMCs
#   Endothelial
#   Epithelial
#
# Not every species contains every broad cell type. Extreme cell-count
# imbalances were reduced by reproducible random capping of selected
# abundant populations so that one population did not dominate the
# cross-species alignment.
# ============================================================

library(Seurat)
library(dplyr)
library(jsonlite)

source("R/cross_species/export_for_samap.R")


# ------------------------------------------------------------
# Randomly cap selected cell types at species-specific maxima
# ------------------------------------------------------------

cap_cell_types <- function(
    object,
    caps,
    group_col = "cell_type",
    seed = 123
) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  if (!group_col %in% colnames(object[[]])) {
    stop("Metadata column `", group_col, "` was not found.")
  }

  if (length(caps) == 0) {
    return(object)
  }

  if (is.null(names(caps)) || any(names(caps) == "")) {
    stop("`caps` must be a named numeric vector.")
  }

  set.seed(seed)

  metadata <- object[[]]
  groups <- as.character(metadata[[group_col]])

  selected_cells <- unlist(
    lapply(unique(groups), function(group_name) {

      group_cells <- rownames(metadata)[
        groups == group_name
      ]

      if (
        group_name %in% names(caps) &&
        length(group_cells) > caps[[group_name]]
      ) {
        sample(
          group_cells,
          size = caps[[group_name]],
          replace = FALSE
        )
      } else {
        group_cells
      }
    }),
    use.names = FALSE
  )

  subset(
    object,
    cells = selected_cells
  )
}


# ------------------------------------------------------------
# Keep comparable broad intestinal cell classes
# ------------------------------------------------------------

prepare_samap_celltypes <- function(
    object,
    cell_types = c(
      "Fibroblasts",
      "SMCs",
      "Endothelial",
      "Epithelial"
    ),
    caps = numeric(0),
    group_col = "cell_type",
    seed = 123
) {
  if (!group_col %in% colnames(object[[]])) {
    stop("Metadata column `", group_col, "` was not found.")
  }

  available <- unique(
    as.character(object[[]][[group_col]])
  )

  missing <- setdiff(
    cell_types,
    available
  )

  if (length(missing) > 0) {
    message(
      "Cell types absent from this dataset and skipped: ",
      paste(missing, collapse = ", ")
    )
  }

  keep_types <- intersect(
    cell_types,
    available
  )

  keep_cells <- rownames(object[[]])[
    object[[]][[group_col]] %in% keep_types
  ]

  object <- subset(
    object,
    cells = keep_cells
  )

  object <- cap_cell_types(
    object,
    caps = caps,
    group_col = group_col,
    seed = seed
  )

  object
}


# ------------------------------------------------------------
# Shared SAMap configuration
# ------------------------------------------------------------

load_samap_config <- function(
    config_file = "config/samap_config.json"
) {
  if (!file.exists(config_file)) {
    stop("SAMap config not found: ", config_file)
  }

  jsonlite::fromJSON(
    config_file,
    simplifyVector = FALSE
  )
}


config <- load_samap_config()

shared_cell_types <- unlist(
  config$broad_cell_types,
  use.names = FALSE
)

samap_caps <- lapply(
  config$cell_caps,
  function(x) {
    if (length(x) == 0) {
      return(numeric(0))
    }

    values <- unlist(x, use.names = TRUE)
    as.numeric(values) |>
      setNames(names(values))
  }
)


# ------------------------------------------------------------
# Example workflow
# ------------------------------------------------------------
# Load the species-level processed objects created elsewhere in this repo:
#
# human <- readRDS("outputs/human_processed.rds")
# mouse <- readRDS("outputs/mouse_processed.rds")
# chicken <- readRDS("outputs/chicken_processed.rds")
# pig <- readRDS("outputs/pig_processed.rds")
# zebrafish <- readRDS("outputs/zebrafish_processed.rds")
# tilapia <- readRDS("outputs/tilapia_processed.rds")
#
# samap_objects <- list(
#   human = prepare_samap_celltypes(human, cell_types = shared_cell_types, caps = samap_caps$hu),
#   mouse = prepare_samap_celltypes(mouse, cell_types = shared_cell_types, caps = samap_caps$mo),
#   chicken = prepare_samap_celltypes(chicken, cell_types = shared_cell_types, caps = samap_caps$ch),
#   pig = prepare_samap_celltypes(pig, cell_types = shared_cell_types, caps = samap_caps$pi),
#   zebrafish = prepare_samap_celltypes(zebrafish, cell_types = shared_cell_types, caps = samap_caps$ze),
#   tilapia = prepare_samap_celltypes(tilapia, cell_types = shared_cell_types, caps = samap_caps$ti)
# )
#
# lapply(
#   samap_objects,
#   function(x) table(x$cell_type)
# )
#
# export_for_samap(
#   samap_objects$human,
#   "outputs/samap_exports/celltypes/human"
# )
#
# Repeat export_for_samap() for the other species.
