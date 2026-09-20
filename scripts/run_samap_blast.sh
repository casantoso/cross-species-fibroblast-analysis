#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Run SAMap pairwise sequence mappings for all species.
#
# Dependency:
#   map_genes.sh from the official SAMap repository.
#
# By default this script expects:
#   ./external/SAMap/map_genes.sh
#
# Override with:
#   MAP_GENES=/path/to/map_genes.sh bash scripts/run_samap_blast.sh
# ============================================================

MAP_GENES="${MAP_GENES:-./external/SAMap/map_genes.sh}"
THREADS="${THREADS:-16}"
FASTADIR="${FASTADIR:-./proteomes}"

if [[ ! -f "$MAP_GENES" ]]; then
    echo "SAMap map_genes.sh was not found at: $MAP_GENES" >&2
    echo "Set MAP_GENES to the path of the script from the SAMap repository." >&2
    exit 1
fi

species=(
    "hu:${FASTADIR}/human.fa"
    "mo:${FASTADIR}/mouse.fa"
    "pi:${FASTADIR}/pig.fa"
    "ch:${FASTADIR}/chicken.fa"
    "ze:${FASTADIR}/zebrafish.fa"
    "ti:${FASTADIR}/tilapia.fa"
)

for entry in "${species[@]}"; do
    fasta="${entry#*:}"
    if [[ ! -f "$fasta" ]]; then
        echo "Missing proteome FASTA: $fasta" >&2
        exit 1
    fi
done

# Run each unique species pair exactly once.
for ((i=0; i<${#species[@]}; i++)); do
    for ((j=i+1; j<${#species[@]}; j++)); do

        id1="${species[$i]%%:*}"
        fa1="${species[$i]#*:}"

        id2="${species[$j]%%:*}"
        fa2="${species[$j]#*:}"

        echo "=================================================="
        echo "SAMap sequence mapping: ${id1} vs ${id2}"
        echo "=================================================="

        bash "$MAP_GENES" \
            --tr1 "$fa1" --t1 prot --n1 "$id1" \
            --tr2 "$fa2" --t2 prot --n2 "$id2" \
            --threads "$THREADS"
    done
done

echo "All pairwise SAMap sequence mappings completed."
