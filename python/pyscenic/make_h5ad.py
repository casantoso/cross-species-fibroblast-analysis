#!/usr/bin/env python3
"""Convert the R pySCENIC export to a raw-count AnnData object."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io as sio


def build_h5ad(input_dir: Path, output_file: Path) -> None:
    genes = pd.read_csv(
        input_dir / "genes.txt",
        header=None,
    )[0].astype(str).tolist()

    cells = pd.read_csv(
        input_dir / "cells.txt",
        header=None,
    )[0].astype(str).tolist()

    # R writes genes x cells; AnnData stores cells x genes.
    X = sio.mmread(
        input_dir / "counts.mtx"
    ).tocsr().T.tocsr()

    if X.shape != (len(cells), len(genes)):
        raise ValueError(
            f"Matrix shape {X.shape} does not match "
            f"{len(cells)} cells x {len(genes)} genes."
        )

    metadata = pd.read_csv(
        input_dir / "metadata.csv"
    )

    if "cell_id" not in metadata:
        raise KeyError("metadata.csv must contain a `cell_id` column.")

    metadata["cell_id"] = metadata["cell_id"].astype(str)
    metadata = metadata.set_index("cell_id").loc[cells].copy()

    adata = ad.AnnData(
        X=X,
        obs=metadata,
        var=pd.DataFrame(index=genes),
    )

    adata.obs_names = cells
    adata.var_names = genes
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    output_file.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_file)

    print(adata)
    print(f"Saved: {output_file}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_h5ad(args.input_dir, args.output)
