# pySCENIC workflow

This repository contains cleaned pySCENIC workflows for **human** and **mouse**
intestinal fibroblasts.

## Workflow

```text
Seurat fibroblast object
        ↓
R/regulatory_analysis/export_for_pyscenic.R
        ↓
raw counts + metadata + original UMAP
        ↓
python/pyscenic/make_h5ad.py
        ↓
raw-count AnnData
        ↓
GRNBoost2
        ↓
cisTarget motif pruning
        ↓
AUCell
        ↓
regulon-activity UMAP / mean activity / RSS / heatmaps
```

Human and mouse are kept as separate notebooks because pySCENIC uses
species-specific TF lists, cisTarget databases, and motif annotations.

Reusable post-processing functions are in:

```text
python/pyscenic/pyscenic_utils.py
```

The exploratory notebooks contained candidate-TF troubleshooting and
project-specific TF/regulon lists. Those result-specific sections are not
included in the polished public notebooks. The generic human–mouse target-set
comparison method is retained in:

```text
notebooks/pySCENIC_cross_species_regulons.ipynb
```

## External resources

The cisTarget ranking databases, motif annotations, and TF lists are public
Aerts-lab resources but are intentionally **not committed to this repository**.
See `python/pyscenic/RESOURCES.md`.

Download them with:

```bash
bash scripts/download_pyscenic_resources.sh
```

Keeping the external resources out of Git avoids hundreds of MB of third-party
data in the repository and makes their source/provenance explicit.


## Command-line run

After the species input `.h5ad` and external Aerts-lab resources are in place:

```bash
conda activate scenic_env
bash scripts/run_pyscenic.sh human
bash scripts/run_pyscenic.sh mouse
```

The script runs GRNBoost2, cisTarget, and AUCell and writes package versions,
resource paths, and key parameters to each species `metadata/run_metadata.json`.

See `docs/input_schemas.md` for the expected AnnData schema and
`the local pySCENIC environment` for the environment definition.
