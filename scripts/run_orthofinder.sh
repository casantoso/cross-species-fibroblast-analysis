#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="${1:-proteomes/orthofinder}"
OUTPUT_DIR="${2:-outputs/orthology}"

if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "Input proteome directory not found: ${INPUT_DIR}" >&2
  exit 1
fi

N_FASTA=$(find "${INPUT_DIR}" -maxdepth 1 -type f \
  \( -name "*.fa" -o -name "*.faa" -o -name "*.fasta" -o -name "*.pep" \) \
  | wc -l | tr -d ' ')

if [[ "${N_FASTA}" -lt 2 ]]; then
  echo "Expected at least two prepared proteome FASTA files in ${INPUT_DIR}." >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "Running OrthoFinder"
echo "Input:  ${INPUT_DIR}"
echo "Output: ${OUTPUT_DIR}"

orthofinder \
  -f "${INPUT_DIR}" \
  -o "${OUTPUT_DIR}"

echo "Finished OrthoFinder."
