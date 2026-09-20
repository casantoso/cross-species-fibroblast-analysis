#!/usr/bin/env bash
set -euo pipefail

# Run the three core pySCENIC stages from config/pyscenic_config.json.
#
# Usage:
#   conda activate scenic_env
#   bash scripts/run_pyscenic.sh human
#   bash scripts/run_pyscenic.sh mouse

SPECIES="${1:-}"
CONFIG="${2:-config/pyscenic_config.json}"

if [[ "${SPECIES}" != "human" && "${SPECIES}" != "mouse" ]]; then
  echo "Usage: bash scripts/run_pyscenic.sh <human|mouse> [config.json]"
  exit 1
fi

readarray -t VALUES < <(
python - "${CONFIG}" "${SPECIES}" <<'PY'
import json
from pathlib import Path
import sys

config_path = Path(sys.argv[1])
species = sys.argv[2]

cfg_all = json.loads(config_path.read_text())
cfg = cfg_all[species]

project = Path(".")
resource_dir = project / cfg["resources_dir"]
result_dir = project / cfg["results_dir"]

values = [
    project / cfg["input_h5ad"],
    resource_dir / cfg["tf_list"],
    resource_dir / cfg["db_10kb"],
    resource_dir / cfg["db_500bp"],
    resource_dir / cfg["motif_annotations"],
    result_dir,
    str(cfg_all["threads"]),
    str(cfg_all["seed"]),
]

for value in values:
    print(value)
PY
)

EXPR="${VALUES[0]}"
TF_LIST="${VALUES[1]}"
DB10="${VALUES[2]}"
DB500="${VALUES[3]}"
ANNOT="${VALUES[4]}"
RESULT_DIR="${VALUES[5]}"
THREADS="${VALUES[6]}"
SEED="${VALUES[7]}"

mkdir -p \
  "${RESULT_DIR}/objects" \
  "${RESULT_DIR}/tables" \
  "${RESULT_DIR}/figures" \
  "${RESULT_DIR}/logs" \
  "${RESULT_DIR}/metadata"

ADJ="${RESULT_DIR}/tables/${SPECIES}_adjacencies.tsv"
REG="${RESULT_DIR}/objects/${SPECIES}_regulons.csv"
SCENIC="${RESULT_DIR}/objects/${SPECIES}_scenic_unmasked.h5ad"

echo "[1/3] GRNBoost2"
pyscenic grn \
  "${EXPR}" \
  "${TF_LIST}" \
  --method grnboost2 \
  --output "${ADJ}" \
  --num_workers "${THREADS}" \
  --seed "${SEED}" \
  2>&1 | tee "${RESULT_DIR}/logs/grn.log"

echo "[2/3] cisTarget"
pyscenic ctx \
  "${ADJ}" \
  "${DB10}" \
  "${DB500}" \
  --annotations_fname "${ANNOT}" \
  --expression_mtx_fname "${EXPR}" \
  --output "${REG}" \
  --num_workers "${THREADS}" \
  --mode dask_multiprocessing \
  2>&1 | tee "${RESULT_DIR}/logs/ctx.log"

echo "[3/3] AUCell"
pyscenic aucell \
  "${EXPR}" \
  "${REG}" \
  --output "${SCENIC}" \
  --num_workers "${THREADS}" \
  --seed "${SEED}" \
  2>&1 | tee "${RESULT_DIR}/logs/aucell.log"

python - "${CONFIG}" "${SPECIES}" <<'PY'
from pathlib import Path
import sys

from python.pyscenic.pyscenic_utils import (
    load_pyscenic_config,
    resolve_species_paths,
    save_run_metadata,
)

config_path = Path(sys.argv[1])
species = sys.argv[2]
config = load_pyscenic_config(config_path)
paths = resolve_species_paths(config, species)

save_run_metadata(
    paths["metadata"] / "run_metadata.json",
    parameters={
        "threads": config["threads"],
        "seed": config["seed"],
        "group_col": config["group_col"],
        "auc_threshold": config["auc_threshold"],
        "umap": config["umap"],
    },
    resources={
        "tf_list": str(paths["tf_list"]),
        "db_10kb": str(paths["db_10kb"]),
        "db_500bp": str(paths["db_500bp"]),
        "motif_annotations": str(paths["motif_annotations"]),
    },
    extra={
        "species": species,
        "expression_file": str(paths["expression"]),
    },
)
PY

echo "Finished pySCENIC for ${SPECIES}"
