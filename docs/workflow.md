# Workflow overview

```text
public scRNA/snRNA-seq counts
        │
        ├── species-specific QC / doublet removal / normalization
        │       R/species/
        │
        ├── broad cell-type analysis
        │
        └── fibroblast subsetting and reclustering
                │
                ├── marker analysis
                │       ├── pathway enrichment
                │       │       R/functions/enrichment_functions.R
                │       └── reusable visualization
                │               R/functions/plotting_functions.R
                │
                ├── OrthoFinder orthology analysis
                │       └── orthogroups / paralogs / gene-family conservation
                │
                ├── orthogene pairwise ortholog lookup
                │       └── direct species A → species B gene-list conversion
                │
                ├── SAMap cross-species alignment
                │       ├── broad-cell-type workflow
                │       └── fibroblast-only workflow
                │
                └── pySCENIC regulatory analysis
                        ├── human
                        ├── mouse
                        └── human–mouse regulon-target comparison
```

The repository exposes the computational methods while keeping unpublished
cluster mappings, biological findings, internal mouse data, and result files
outside version control.
