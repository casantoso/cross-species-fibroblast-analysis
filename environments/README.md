# Environments

This directory records the software environments used for the main analysis workflows.

## R

`r_package_versions.tsv` records the verified R package versions used for the Seurat-based analysis workflow.

## SAMap

- Environment name: `SAMap`
- Python: 3.12
- SAMap package: `sc-samap==3.0.1`

Files:
- `samap_environment.yml`
- `samap_requirements.txt`

These files contain a curated set of packages directly relevant to the SAMap workflow, using versions verified from the workstation environment.

## pySCENIC

- Environment name: `pyscenic_env`
- Python: 3.10
- pySCENIC: `0.12.1`

Files:
- `pyscenic_environment.yml`
- `pyscenic_requirements.txt`

These files contain a curated set of packages directly relevant to the pySCENIC workflow, using versions verified from the workstation environment.

## OrthoFinder

- Environment name: `of3_env`
- Python: 3.12
- OrthoFinder: `3.1.5`

File:
- `orthofinder_environment.yml`

The OrthoFinder environment is based on the verified conda history from the workstation.

## Reproducibility note

The YAML and requirements files are intentionally curated rather than full package dumps. Jupyter-only packages, transient dependencies, and machine-specific `file:///...` build paths are excluded.
