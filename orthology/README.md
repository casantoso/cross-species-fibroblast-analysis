# OrthoFinder orthology workflow

This section documents the evolutionary/orthology layer of the project.

## Purpose

SAMap provides expression-aware cross-species alignment, while OrthoFinder
provides gene-family relationships that are useful for questions such as:

- which fibroblast genes belong to the same orthogroup across species;
- which species contain multiple paralogs in an orthogroup;
- whether a marker is represented by one gene or several paralogs;
- how marker/conservation analyses change when genes are compared at the
  orthogroup rather than strict 1:1-ortholog level.

The public repository includes the **workflow and parsers**, but not the large
proteome FASTA files or the complete OrthoFinder result directory.

## Species

The public workflow is configured for:

```text
human
mouse
pig
chicken
zebrafish
tilapia
```


## Exact protein FASTA sources used

All six protein FASTA files used for the OrthoFinder/SAMap protein-preparation
workflow came from **Ensembl release 116**. The project used the `pep.all`
files, not the `pep.abinitio` files.

The exact filenames, assemblies, and source directories are recorded in:

```text
orthology/proteome_manifest.tsv
```

The files used were:

```text
Human:
Homo_sapiens.GRCh38.pep.all.fa.gz

Mouse:
Mus_musculus.GRCm39.pep.all.fa.gz

Pig:
Sus_scrofa.Sscrofa11.1.pep.all.fa.gz

Chicken:
Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa.gz

Zebrafish:
Danio_rerio.GRCz11.pep.all.fa.gz

Nile tilapia:
Oreochromis_niloticus.O_niloticus_UMD_NMBU.pep.all.fa.gz
```


## 1. Prepare one representative protein sequence per gene

OrthoFinder works best when the input proteome does not contain several protein
isoforms for the same gene.

For the analyses in this repository, use one representative protein per gene,
typically the longest protein isoform.

A reusable gene-symbol-based helper already used elsewhere in this repository
is:

```text
python/samap/longest_protein_gene_symbols.py
```

It is suitable when the Ensembl peptide FASTA headers contain a usable
`gene_symbol` field and you want the FASTA identifiers to match scRNA-seq gene
symbols.

The uploaded `primary_transcript.py` used during exploratory work is an
OrthoFinder-style utility that supports Ensembl/NCBI FASTA headers and keeps
the longest isoform per gene. It is **not copied into this repository as
project-authored code**. If that upstream helper is used, document its source
and version separately.

Prepared proteomes should be placed in:

```text
proteomes/orthofinder/
├── human.fa
├── mouse.fa
├── pig.fa
├── chicken.fa
├── zebrafish.fa
└── tilapia.fa
```

`proteomes/` is ignored by Git.

## 2. Run OrthoFinder

From the repository root:

```bash
# Activate the OrthoFinder environment used on the analysis machine.

bash scripts/run_orthofinder.sh
```

By default the wrapper uses:

```text
input:   proteomes/orthofinder/
output:  outputs/orthology/
```

You can override both:

```bash
bash scripts/run_orthofinder.sh \
  /path/to/prepared_proteomes \
  /path/to/output_directory
```

## 3. Key OrthoFinder outputs

The most useful files for this project are usually found under the generated
OrthoFinder results directory, including:

```text
Orthogroups/Orthogroups.tsv
Orthogroups/Orthogroups.GeneCount.tsv
Orthogroups/Orthogroups_SingleCopyOrthologues.txt
Gene_Duplication_Events/Duplications.tsv
Species_Tree/SpeciesTree_rooted.txt
```

Exact filenames can vary slightly by OrthoFinder version.

If the FASTA headers are gene symbols, `Orthogroups.tsv` is already directly
joinable to gene-level marker tables. If the FASTA headers are stable sequence
or gene IDs, map those IDs to symbols before downstream visualization.


## Reproducible OrthoFinder ID → gene-symbol mapping

Raw OrthoFinder output contains the identifiers present in the input FASTA
headers. If those identifiers are Ensembl gene IDs, protein IDs, or transcript
IDs rather than gene symbols, the symbol-annotation step must be reproducible
rather than performed manually.

### A. Build one mapping table per species

For GTF/GFF-like annotation files:

```bash
python python/orthology/build_id_symbol_map_from_gtf.py \
  --annotation annotations/human.gtf \
  --id-attribute gene_id \
  --output mappings/orthofinder/human.tsv
```

The generated file has:

```text
id    symbol
```

Repeat for each species using the identifier type that matches the FASTA
headers used for OrthoFinder. For example, if the FASTA identifiers are
`protein_id` values, use:

```bash
--id-attribute protein_id
```

The mapping filenames should match the OrthoFinder species-column names:

```text
mappings/orthofinder/
├── human.tsv
├── mouse.tsv
├── pig.tsv
├── chicken.tsv
├── zebrafish.tsv
└── tilapia.tsv
```

Annotation files and generated mapping tables can remain local if licensing or
file size makes redistribution undesirable.

### B. Annotate raw Orthogroups.tsv

```bash
python python/orthology/map_orthofinder_ids_to_symbols.py \
  --orthogroups path/to/Orthogroups.tsv \
  --mapping-dir mappings/orthofinder \
  --output outputs/orthology/tables/Orthogroups_with_symbols.tsv \
  --summary outputs/orthology/tables/orthogroup_mapping_summary.tsv
```

Mapped entries are preserved as:

```text
PDGFRA|ENSG...
```

and unresolved IDs are retained as:

```text
UNMAPPED|ENSG...
```

instead of silently dropping them. This makes it possible to audit mapping
coverage before downstream cleanup.

For tilapia, `LOC...` identifiers in the single-cell dataset may require
manual NCBI Gene inspection even after this automated mapping step. Important
zebrafish mappings should also be checked in NCBI Gene when duplication or
missing annotation makes the correspondence ambiguous.


## 4. Clean a symbol-mapped orthogroup table

`python/orthology/clean_orthogroups.py` generalizes the exploratory
`clean_orthogroups.R/py` step. It removes entries prefixed with `UNMAPPED|`
without hard-coded input/output filenames.

Example:

```bash
python python/orthology/clean_orthogroups.py \
  --input Orthogroups_with_symbols.tsv \
  --output outputs/orthology/tables/Orthogroups_mapped_only.tsv
```

Optionally also write an Excel copy:

```bash
python python/orthology/clean_orthogroups.py \
  --input Orthogroups_with_symbols.tsv \
  --output outputs/orthology/tables/Orthogroups_mapped_only.tsv \
  --xlsx outputs/orthology/tables/Orthogroups_mapped_only.xlsx
```

## 5. Downstream use

The cleaned orthogroup table can be joined to fibroblast-marker tables to
compare conservation at the orthogroup/paralog level rather than assuming
strict one-to-one orthology.

A typical downstream logic is:

```text
species fibroblast markers
        +
OrthoFinder orthogroups
        ↓
marker → orthogroup assignments
        ↓
species representation / paralog counts
        ↓
comparison with SAMap gene-pair relationships
```

This repository intentionally does not include unpublished orthogroup-level
results, specific conserved components, or paralog-substitution findings.


## Direct gene-list conversion with orthogene

OrthoFinder is the preferred workflow for **orthogroups, paralogs, and
gene-family relationships across all species**. However, sometimes the
question is simpler:

```text
I have a gene list from species A.
What are the corresponding orthologs in species B?
```

For that use case, the repository also includes:

```text
R/cross_species/map_orthologs_orthogene.R
```

The helper wraps `orthogene::convert_orthologs()` and can either retain only
strict 1:1 mappings or preserve non-1:1 mappings.

Example, human to mouse:

```r
source("R/cross_species/map_orthologs_orthogene.R")

mouse_map <- map_gene_list_orthologs(
  genes = genes,
  input_species = "human",
  output_species = "mouse",
  one_to_one_only = TRUE
)
```

### Important note for zebrafish

For zebrafish, a human/mammalian gene may correspond to more than one
zebrafish paralog. Therefore, strict 1:1 filtering can remove biologically
meaningful mappings.

For exploratory zebrafish conversion, it may be better to preserve the
one-to-many relationships:

```r
zebrafish_map <- map_gene_list_orthologs(
  genes = genes,
  input_species = "human",
  output_species = "zebrafish",
  one_to_one_only = FALSE
)
```

The resulting paralogs should then be interpreted individually rather than
assuming that one fish paralog is automatically equivalent to the mammalian
gene.

### Species names used in the helper

The public helper intentionally uses readable names:

```text
human
mouse
pig
chicken
zebrafish
tilapia
```

Internally, these are translated through:

```r
SPECIES_SCI <- c(
  human     = "hsapiens",
  mouse     = "mmusculus",
  pig       = "sscrofa",
  chicken   = "ggallus",
  zebrafish = "drerio",
  tilapia   = "oniloticus",
  macaque   = "mfascicularis"
)
```

This means the analysis code can stay readable:

```r
mouse_map <- map_gene_list_orthologs(
  genes = genes,
  input_species = "human",
  output_species = "mouse",
  one_to_one_only = TRUE
)
```

and Nile tilapia can be written simply as:

```r
tilapia_map <- map_gene_list_orthologs(
  genes = genes,
  input_species = "human",
  output_species = "tilapia",
  one_to_one_only = TRUE
)
```

while the helper internally passes `oniloticus` to the selected backend.

For zebrafish, one-to-many relationships can be biologically meaningful
because duplicated fish paralogs may be retained, so the zebrafish example
preserves non-1:1 mappings.

For Nile tilapia, the default example in this repository uses **strict 1:1
mapping** for a simpler gene-list conversion. However, the tilapia annotation
used in this project contains many genes without informative gene symbols or
with incomplete labels. Therefore, important tilapia genes that fail to map,
remain unlabeled, or produce an unexpected ortholog should be checked manually
in **NCBI Gene** (and, when useful, Ensembl) before being treated as missing or
non-conserved.


### Output column names

`map_gene_list_orthologs()` renames the two returned gene columns to the
readable species names used in the call.

For example:

```r
mouse_map <- map_gene_list_orthologs(
  genes = genes,
  input_species = "human",
  output_species = "mouse"
)
```

returns a table with:

```text
human    mouse
```

rather than generic names such as `input_gene` and `ortholog_gene`.

### Multi-species ortholog table

To append orthologs from several species to one table, use:

```r
multi_species_map <- map_gene_list_multi_species(
  genes = genes,
  input_species = "human",
  output_species = c(
    "mouse",
    "pig",
    "chicken",
    "zebrafish",
    "tilapia"
  ),
  one_to_one_only = TRUE
)
```

which returns columns such as:

```text
human | mouse | pig | chicken | zebrafish | tilapia
```

For a strict 1:1 table, keep `one_to_one_only = TRUE`.

For zebrafish, it can also be useful to run the pairwise helper separately with
`one_to_one_only = FALSE` to inspect retained paralogs. Because zebrafish
orthology can be complicated by duplicated genes, important mappings should be
checked manually in **NCBI Gene**, especially when multiple paralogs are
returned or an expected ortholog is missing.

For Nile tilapia, the default project example remains strict 1:1. The tilapia
single-cell dataset used in this project contains many genes labeled with
provisional `LOC...` identifiers rather than informative gene symbols.
Therefore, important `LOC...` genes should be manually checked in **NCBI Gene**
to determine their current annotation and possible orthologs before treating
them as unmapped, absent, or non-conserved. Important missing or ambiguous
non-`LOC...` mappings should also be manually verified.

### Why keep both OrthoFinder and orthogene?

```text
OrthoFinder
    → multi-species orthogroups
    → paralog/gene-family structure
    → evolutionary analyses

orthogene
    → convenient species A → species B gene-list conversion
    → optional strict 1:1 filtering
    → optional retention of one-to-many mappings
```

They are complementary rather than replacements for one another.
