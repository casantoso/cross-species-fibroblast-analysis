# Data sources

Raw data and processed objects are **not** distributed in this repository.
Public datasets are referenced by accession so the analysis provenance is
transparent.

| Species | Public source used | Selection represented in the workflow |
|---|---|---|
| Human | GSE178341 | Normal colon cells used for the human workflow |
| Mouse | Internal dataset | Internal intestinal dataset; not distributed |
| Pig | GSE233285 / GSM7421580 | Ileum replicate 1, snRNA-seq |
| Chicken | GSE283090 / GSM8655812–GSM8655815 | Day 0 intestinal villus/source-tissue samples: Broiler rep 1–2 and Layer rep 1–2 |
| Zebrafish | GSE271002 / GSM8366962 | Control, 5 dpf |
| Nile tilapia | GSE284663 / GSM8689993–GSM8689995 | Untreated biological replicates 1–3 |

The public repository intentionally excludes unpublished cluster-to-subtype
mappings, internal mouse data, processed count matrices, figures, and result
tables.

Species-specific loading and QC details are documented directly in the scripts
under `R/species/`.


### Nile tilapia annotation note

The Nile tilapia single-cell dataset contains many genes with provisional
`LOC...` identifiers rather than informative gene symbols. For important genes,
these identifiers should be checked manually in NCBI Gene to determine the
current gene annotation and corresponding orthologs before concluding that a
gene is unmapped or absent in tilapia.
