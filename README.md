# Cross-species intestinal fibroblast analysis

This repository contains cleaned, reusable code for cross-species intestinal
single-cell/single-nucleus RNA-seq preprocessing and fibroblast-focused
analysis. Unpublished biological findings, internal mouse data, processed
objects, figures, and cluster-to-subtype mappings are intentionally excluded.

## Workflow

```text
public scRNA/snRNA-seq
        ↓
species QC / doublet removal / normalization
        ↓
broad cell annotation
        ↓
fibroblast subsetting and reclustering
        ↓
marker analysis
   ├── pathway enrichment
   ├── reusable plotting
   ├── OrthoFinder orthology analysis
   ├── SAMap cross-species alignment
   └── pySCENIC regulatory-network analysis
```

A longer overview is in `docs/workflow.md`.

## Data sources

Public source accessions and selected samples are documented in
`data/README.md`. The represented datasets are:

- Human: GSE178341, normal colon
- Mouse: internal intestinal dataset; not distributed
- Pig: GSE233285 / GSM7421580
- Chicken: GSE283090 / GSM8655812–GSM8655815
- Zebrafish: GSE271002 / GSM8366962
- Nile tilapia: GSE284663 / GSM8689993–GSM8689995

## Repository structure

```text
config/
├── samap_config.json
└── pyscenic_config.json

data/
└── README.md

docs/
├── workflow.md
├── input_schemas.md
├── output_layout.md
└── PROVENANCE.md

environments/
├── README.md
└── r_package_versions.tsv

R/
├── R_requirements.md
├── functions/
│   ├── preprocessing_functions.R
│   ├── fibroblast_functions.R
│   ├── plotting_functions.R
│   └── enrichment_functions.R
├── species/
│   ├── human.R
│   ├── mouse.R
│   ├── pig.R
│   ├── chicken.R
│   ├── zebrafish.R
│   └── tilapia.R
├── cross_species/
│   ├── prepare_samap_celltypes.R
│   ├── export_for_samap.R
│   └── map_orthologs_orthogene.R
└── regulatory_analysis/
    └── export_for_pyscenic.R

orthology/
├── README.md
└── proteome_manifest.tsv

python/
├── orthology/
│   ├── build_id_symbol_map_from_gtf.py
│   ├── map_orthofinder_ids_to_symbols.py
│   ├── clean_orthogroups.py
│   └── summarize_orthogroups.py
├── samap/
│   ├── samap_analysis_utils.py
│   ├── make_h5ad.py
│   ├── longest_protein_gene_symbols.py
│   └── README.md
└── pyscenic/
    ├── pyscenic_utils.py
    ├── make_h5ad.py
    ├── RESOURCES.md
    ├── resource_manifest.tsv
    └── README.md

notebooks/
├── SAMap_celltype_analysis.ipynb
├── SAMap_fibroblast_analysis.ipynb
├── pySCENIC_human.ipynb
├── pySCENIC_mouse.ipynb
└── pySCENIC_cross_species_regulons.ipynb

scripts/
├── run_orthofinder.sh
├── run_samap_blast.sh
├── run_pyscenic.sh
├── download_pyscenic_resources.sh
└── record_pyscenic_resource_checksums.py

tests/
├── run_r_tests.R
├── test_samap_utils.py
├── test_pyscenic_utils.py
├── test_orthology_utils.py
├── test_orthofinder_mapping.py
└── testthat/
    ├── test_preprocessing_functions.R
    ├── test_plotting_functions.R
    └── test_enrichment_functions.R
```

## Species-level preprocessing

Each species script contains the public preprocessing workflow: count loading,
quality control, doublet removal where applicable, SCTransform normalization,
PCA, Harmony for multi-sample datasets where required, clustering, broad
cell-type annotation hooks, fibroblast subsetting/reclustering, and marker
detection.

The public scripts use a consistent fibroblast marker filter based on:

```text
pct.1 - pct.2 > 0.3
```

Human and mouse subtype mappings are not included while those results remain
unpublished. Pig, chicken, zebrafish, and Nile tilapia are not assigned those
human/mouse subtype labels.


## OrthoFinder orthology analysis

OrthoFinder is used as a complementary evolutionary layer to SAMap. SAMap
provides expression-aware cross-species alignment, whereas OrthoFinder
provides orthogroup and paralog relationships that can be joined back to
fibroblast markers.

The public workflow is documented in:

```text
orthology/README.md
```

Create the environment and run:

```bash
# Activate the OrthoFinder environment used on the analysis machine.

bash scripts/run_orthofinder.sh
```

Prepared one-protein-per-gene FASTA files are kept under
`proteomes/orthofinder/` and are ignored by Git. The large OrthoFinder result
directory is also generated outside version control.

Utility scripts for cleaning and summarizing orthogroup tables are in:

```text
python/orthology/
```

The repository does not include unpublished orthogroup-level biological
results.


For convenient **species A → species B gene-list conversion**, the repository
also includes an optional `orthogene` helper:

```text
R/cross_species/map_orthologs_orthogene.R
```

The helper returns columns named by species (for example `human`, `mouse`) and
also includes `map_gene_list_multi_species()` for adding orthologs from several
species to the same table. Zebrafish mappings should be manually checked in
NCBI Gene when paralogs or missing mappings make the correspondence uncertain.

This can restrict mappings to 1:1 orthologs or preserve one-to-many mappings.
The helper accepts readable names such as `human`, `mouse`, `pig`,
`chicken`, `zebrafish`, and `tilapia`, then translates them internally to
backend identifiers (`hsapiens`, `mmusculus`, `sscrofa`, `ggallus`, `drerio`,
and `oniloticus`). Zebrafish analyses can preserve one-to-many mappings when
paralogs are biologically relevant. The Nile tilapia example uses strict 1:1 mapping. The tilapia single-cell
dataset used in this project contains many provisional `LOC...` gene names, so
important `LOC...` genes should be manually checked in NCBI Gene for updated
annotation and ortholog information before being treated as unmapped, absent,
or non-conserved.

## SAMap

Two workflows are included:

```text
broad intestinal cell types
→ SAMap alignment
→ fibroblast GenePairFinder relationships
→ gene-pair networks / species combinations / paralog analysis
```

and:

```text
fibroblast-only objects
→ SAMap alignment
→ global Leiden clustering
→ species composition / species-specific marker ranking
```

Reusable code is in `python/samap/samap_analysis_utils.py`.

### SAMap environment

The public repository does not include a reconstructed SAMap environment file.
Use the verified SAMap environment from the analysis machine when reproducing
this workflow.

## pySCENIC

The cleaned human and mouse workflow is:

```text
raw fibroblast counts
→ GRNBoost2
→ cisTarget motif pruning
→ AUCell
→ regulon-activity UMAP
→ mean activity / RSS
→ optional human–mouse regulon-target comparison
```

The large cisTarget ranking databases, motif annotations, and TF lists are
external Aerts-lab resources and are not committed. Their provenance and
download instructions are in `python/pyscenic/RESOURCES.md`.

The exact human and mouse resource filenames used in the analysis are tracked
in `python/pyscenic/resource_manifest.tsv` and documented in
`python/pyscenic/RESOURCES.md`. The project uses the cisTarget **ranking**
databases (not the score databases) for both the ±10 kb TSS and
500-bp-upstream/100-bp-downstream search spaces. After download,
`scripts/record_pyscenic_resource_checksums.py` records local SHA256 hashes and
verifies the ranking-database files against the official Aerts Lab
`sha256sum.txt` list.

### pySCENIC environment

The original pySCENIC notebooks used the `scenic_env` kernel. This repository
does not include an unverified reconstructed environment file. The exact
resource filenames used by the analysis are fully documented under
`python/pyscenic/`.

## Plotting and pathway enrichment

Reusable visualization helpers are in:

```text
R/functions/plotting_functions.R
```

Marker-based pathway enrichment is generalized in:

```text
R/functions/enrichment_functions.R
```

For a `FindAllMarkers()` table:

```r
source("R/functions/enrichment_functions.R")

enrichment <- run_marker_pathway_enrichment(
  markers = marker_table,
  group_col = "cluster",
  top_n_genes = 100
)

enrichment$plot
enrichment$top_enrichment
```

The Enrichr database vector can be edited directly or selected from the current
Enrichr library catalog with the included helper functions.

## Input and output conventions

- Expected SAMap, pySCENIC, and enrichment inputs:
  `docs/input_schemas.md`
- Standard generated-output layout:
  `docs/output_layout.md`

Generated data and results are excluded from Git.

## Reproducibility

Python environments are separated under `environments/` because SAMap and
pySCENIC have different dependency stacks.

For R, the exact package versions from the working analysis environment are
already recorded in:

```text
environments/r_package_versions.tsv
```

The recorded R version is `4.5.3`, with exact versions documented for Seurat,
glmGamPoi, harmony, scDblFinder, SingleCellExperiment, orthogene, enrichR, and
the other public-workflow dependencies. `R/R_requirements.md` distinguishes
packages called directly by the public code from transitive/environment
dependencies.



SAMap and pySCENIC runs also save package versions and key run parameters next
to their generated results.

## Tests

Lightweight Python tests cover reusable SAMap and pySCENIC helpers:

```bash
pytest tests/test_samap_utils.py tests/test_pyscenic_utils.py tests/test_orthology_utils.py tests/test_orthofinder_mapping.py
```

R helper tests are provided with `testthat`:

```bash
Rscript tests/run_r_tests.R
```

A GitHub Actions workflow performs Python syntax checks and utility tests on
pushes and pull requests.



## Code and resource provenance

External software, third-party resource files, public datasets, and
project-authored wrappers are distinguished in `docs/PROVENANCE.md`.
