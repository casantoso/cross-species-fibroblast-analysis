# pySCENIC external resources

The pySCENIC workflow requires external cisTarget resources produced and
distributed by the **Aerts lab**. These files are **not project-generated data**
and are not committed to this repository.

Official resource portal:

https://resources.aertslab.org/cistarget/

The portal provides:

- transcription-factor lists (`tf_lists/`);
- gene-based cisTarget ranking databases (`databases/`);
- motif-to-TF annotation tables (`motif2tf/`).

This repository expects the current v10/mc_v10_clust resources used in the
original analysis:

## Human

```text
resources/pyscenic/human/
├── allTFs_hg38.txt
├── hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
├── hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
└── motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

## Mouse

```text
resources/pyscenic/mouse/
├── allTFs_mm.txt
├── mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
├── mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
└── motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
```

Download them with:

```bash
bash scripts/download_pyscenic_resources.sh
```

The ranking databases are hundreds of MB each and the motif annotation tables
are also large, so keeping them outside Git avoids bloating the repository and
makes their external provenance clear.

When using these resources, cite the pySCENIC/cisTarget resources and the
relevant pySCENIC/SCENIC publications according to the Aerts lab documentation.


## Exact resource manifest

The exact filenames, genome builds, annotation releases, and source URLs used
by this project are stored in:

```text
python/pyscenic/resource_manifest.tsv
```

Checksum verification is optional and generated at runtime. Running:

```bash
python scripts/record_pyscenic_resource_checksums.py
```

writes:

```text
outputs/pyscenic/resource_checksums.tsv
```

so the tracked manifest contains no placeholder values.

## Exact resource files used in this project

The following filenames are the exact resources used in the human and mouse
pySCENIC workflows.

### Mouse

```text
allTFs_mm.txt
mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
```

Sources:

```text
https://resources.aertslab.org/cistarget/tf_lists/
https://resources.aertslab.org/cistarget/motif2tf/
https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/
```

### Human

```text
allTFs_hg38.txt
hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

Sources:

```text
https://resources.aertslab.org/cistarget/tf_lists/
https://resources.aertslab.org/cistarget/motif2tf/
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/
```

These filenames are also recorded in `resource_manifest.tsv`.

The SHA256 column remains `` until the actual local files
are hashed. The filename/source information therefore documents exactly which
resources were used even before checksums are recorded.


## Ranking databases and official checksum verification

The Aerts Lab resource page distinguishes:

- `*.rankings.feather`: motif-by-gene ranking matrices used by cisTarget;
- `*.scores.feather`: motif-by-gene CRM score matrices intended for DEM.

This project used the **ranking** databases, not the score databases.

The two regulatory search spaces used were:

- `10kbp_up_10kbp_down`: 10 kb upstream and 10 kb downstream of the TSS;
- `500bp_up_100bp_down`: 500 bp upstream and 100 bp downstream of the TSS.

The Aerts Lab database resource page also provides an official SHA256 list:

```text
https://resources.aertslab.org/cistarget/databases/sha256sum.txt
```

Therefore the repository now distinguishes:

```text
official_sha256
    checksum published in the Aerts database checksum list

local_sha256
    checksum calculated from the exact file on the workstation

checksum_status
    VERIFIED / MISMATCH / LOCAL_HASH_RECORDED / etc.
```

For the cisTarget ranking databases, `VERIFIED` means the locally downloaded
file has exactly the same SHA256 as the value published in the Aerts checksum
list.

For the TF lists and motif-annotation tables, the repository records the local
SHA256 but does not claim an upstream checksum unless one has been explicitly
documented for that resource.
