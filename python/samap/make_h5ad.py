#!/usr/bin/env python3
"""
Convert a SAMap export directory produced by export_for_samap.R to AnnData.

Expected input files:
    matrix.mtx
    genes.txt
    cells.txt
    metadata.csv

R exports the matrix as genes x cells. AnnData stores cells x genes, so the
matrix is transposed during conversion.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def build_h5ad(input_dir: Path, output_file: Path) -> None:
    required = {
        "matrix": input_dir / "matrix.mtx",
        "genes": input_dir / "genes.txt",
        "cells": input_dir / "cells.txt",
        "metadata": input_dir / "metadata.csv",
    }

    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required SAMap export file(s):\n" + "\n".join(missing)
        )

    genes = pd.read_csv(required["genes"], header=None)[0].astype(str).tolist()
    cells = pd.read_csv(required["cells"], header=None)[0].astype(str).tolist()

    # Matrix Market file from R is genes x cells.
    X = sio.mmread(required["matrix"]).tocsr().T.tocsr()

    if X.shape != (len(cells), len(genes)):
        raise ValueError(
            "Matrix dimensions do not match cells/genes: "
            f"matrix={X.shape}, cells={len(cells)}, genes={len(genes)}"
        )

    metadata = pd.read_csv(required["metadata"])

    if "cell_id" not in metadata.columns:
        raise KeyError("metadata.csv must contain a `cell_id` column.")

    metadata["cell_id"] = metadata["cell_id"].astype(str)
    metadata = metadata.set_index("cell_id")

    missing_metadata = sorted(set(cells) - set(metadata.index))
    if missing_metadata:
        raise ValueError(
            f"{len(missing_metadata)} cells are missing from metadata.csv."
        )

    metadata = metadata.loc[cells].copy()

    adata = ad.AnnData(
        X=X,
        obs=metadata,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )

    adata.obs_names = cells
    adata.var_names = genes
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    adata.X = adata.X.tocsr()

    # Useful sanity check: normalized expression should generally contain
    # non-integer values.
    x = adata.X.data
    noninteger_fraction = (
        float(np.mean(x != np.round(x))) if len(x) else float("nan")
    )

    print(adata)
    print("Sparse matrix:", sp.issparse(adata.X))
    print("Non-integer fraction:", noninteger_fraction)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_file)

    print(f"Saved: {output_file}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Directory containing matrix.mtx, genes.txt, cells.txt, metadata.csv.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output .h5ad path.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_h5ad(args.input_dir, args.output)
