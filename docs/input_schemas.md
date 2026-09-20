# Input schemas

This document describes the minimum inputs expected by the reusable workflows.

## Species-level Seurat workflows

Species scripts under `R/species/` begin from the source-specific count
matrices/metadata described in `data/README.md`. The public scripts produce
Seurat objects with standard QC metadata and `seurat_clusters`.

## SAMap

SAMap inputs are exported with:

```text
R/cross_species/export_for_samap.R
```

The exporter writes:

```text
matrix.mtx     log-normalized RNA expression, genes × cells
genes.txt      feature names in matrix row order
cells.txt      cell IDs in matrix column order
metadata.csv   metadata aligned to cells
```

`python/samap/make_h5ad.py` converts these files to AnnData (`cells × genes`).

For the broad-cell-type workflow, `adata.obs` must contain:

```text
cell_type
```

For fibroblast-only analyses, `seurat_clusters` can be retained when original
Seurat clusters will later be projected onto the SAMap embedding.

`adata.var_names` must use identifiers consistent with the representative
FASTA/BLAST gene identifiers used by SAMap.

## pySCENIC

pySCENIC inputs are exported with:

```text
R/regulatory_analysis/export_for_pyscenic.R
```

The exporter writes:

```text
counts.mtx     raw RNA counts, genes × cells
genes.txt      feature names
cells.txt      cell IDs
metadata.csv   metadata aligned to cells
umap.tsv       optional original Seurat UMAP
```

`python/pyscenic/make_h5ad.py` converts the raw matrix to an AnnData object.

For group-level RSS/mean-activity analyses, `adata.obs` must contain the
metadata column configured as `group_col` in `config/pyscenic_config.json`.
The public default is `fib_subtype`, but this can be changed to another
grouping variable.

The pySCENIC expression matrix must represent raw counts; the SCENIC UMAP and
regulon activity are computed downstream from the AUCell output.

## Pathway enrichment

`run_marker_pathway_enrichment()` accepts either:

- a `FindAllMarkers()` table containing `cluster`, `gene`, `avg_log2FC`, and
  `p_val_adj`; or
- a single `FindMarkers()` table, where gene names may be stored as row names.

Column names can be changed through function arguments.
