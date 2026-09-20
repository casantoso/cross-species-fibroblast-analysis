#!/usr/bin/env bash
set -euo pipefail

# Download the exact pySCENIC/cisTarget resources used by this project.
#
# Sources supplied from the Aerts Lab cisTarget resource server:
#   https://resources.aertslab.org/cistarget/tf_lists/
#   https://resources.aertslab.org/cistarget/motif2tf/
#   https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/
#   https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/

ROOT="resources/pyscenic"
MOUSE="${ROOT}/mouse"
HUMAN="${ROOT}/human"

mkdir -p "${MOUSE}" "${HUMAN}"

download_if_missing () {
  local url="$1"
  local output="$2"

  if [[ -f "${output}" ]]; then
    echo "Already present: ${output}"
  else
    echo "Downloading: ${output}"
    curl -L "${url}" -o "${output}"
  fi
}

# ---------------------------
# Mouse
# ---------------------------

download_if_missing \
  "https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt" \
  "${MOUSE}/allTFs_mm.txt"

download_if_missing \
  "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  "${MOUSE}/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

download_if_missing \
  "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  "${MOUSE}/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

download_if_missing \
  "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl" \
  "${MOUSE}/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

# ---------------------------
# Human
# ---------------------------

download_if_missing \
  "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt" \
  "${HUMAN}/allTFs_hg38.txt"

download_if_missing \
  "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  "${HUMAN}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

download_if_missing \
  "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  "${HUMAN}/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

download_if_missing \
  "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
  "${HUMAN}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

echo "Downloading official Aerts cisTarget database checksum list..."
curl -L \
  "https://resources.aertslab.org/cistarget/databases/sha256sum.txt" \
  -o "${ROOT}/sha256sum.txt"

echo "Recording local SHA256 hashes and verifying ranking databases..."
python scripts/record_pyscenic_resource_checksums.py


echo "Done."
