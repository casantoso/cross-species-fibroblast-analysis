# Code and resource provenance

This repository separates project-authored workflow code from external tools,
third-party resources, and public datasets.

## Project-authored / project-specific code

The reusable wrappers, preprocessing helpers, plotting functions, orthology
table processors, SAMap analysis utilities, pySCENIC post-processing helpers,
tests, configuration files, and workflow notebooks in this repository were
cleaned/generalized for this project.

They orchestrate established external tools; they do not reimplement the core
algorithms of those tools.

## OrthoFinder

OrthoFinder is an external orthology-inference program. The repository includes
a project wrapper (`scripts/run_orthofinder.sh`) and project-specific parsers
for converting its outputs into downstream tables.

The exploratory workflow also used an external utility named
`primary_transcript.py` for selecting one representative/longest protein
isoform per gene. That external script is **not redistributed here as
project-authored code**. The public repository instead documents the
representative-protein step and provides project-specific mapping utilities.

Raw proteome FASTA files and complete OrthoFinder result directories are not
committed.

## SAMap

SAMap is an external cross-species single-cell mapping method. The repository
contains project-specific preparation, configuration, downstream analysis, and
visualization code.

Any upstream SAMap helper such as `map_genes.sh` remains external and is not
presented as code authored by this project.

## pySCENIC / cisTarget resources

pySCENIC, GRNBoost2/arboreto, cisTarget/ctxcore, and AUCell are external
software components.

The TF lists, ranking databases, and motif-to-TF annotation tables are external
Aerts-lab cisTarget resources. They are not redistributed in Git. The exact resource filenames used in the project are listed in:

```text
python/pyscenic/resource_manifest.tsv
```

and an optional local checksum-verification report can be generated with:

```bash
python scripts/record_pyscenic_resource_checksums.py
```

For cisTarget ranking databases, the script also verifies those local hashes
against the official Aerts Lab database checksum list at
`https://resources.aertslab.org/cistarget/databases/sha256sum.txt`.

## Enrichr

Pathway enrichment uses the external Enrichr service through the `enrichR` R
package. Library names can change over time, so the repository includes helpers
that query the live Enrichr catalog rather than assuming all year-tagged
libraries remain unchanged.

## Public sequencing datasets

Public GEO accessions and sample selections are documented in
`data/README.md`. Raw public data are not redistributed by this repository.

The internal mouse dataset is not public and is not included.

## NCBI Gene manual verification

NCBI Gene is used as a manual verification resource for orthology/annotation
cases that are not safely resolved by automated mapping, especially:

- duplicated or ambiguous zebrafish orthologs;
- Nile tilapia genes represented by provisional `LOC...` identifiers;
- unexpected missing mappings in important genes.

Manual verification is treated as a review step, not as a substitute for
recorded automated mappings.
