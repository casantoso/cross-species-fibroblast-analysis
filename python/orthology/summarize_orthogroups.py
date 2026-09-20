#!/usr/bin/env python3
"""Summarize species representation and paralog counts per orthogroup."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def count_genes(cell: str) -> int:
    if cell is None or not str(cell).strip():
        return 0

    return len([
        gene
        for gene in str(cell).split(",")
        if gene.strip()
    ])


def summarize_orthogroups(df: pd.DataFrame) -> pd.DataFrame:
    if "Orthogroup" not in df.columns:
        raise ValueError("Input table must contain an `Orthogroup` column.")

    species_columns = [
        column
        for column in df.columns
        if column != "Orthogroup"
    ]

    summary = pd.DataFrame({
        "Orthogroup": df["Orthogroup"],
    })

    for species in species_columns:
        summary[f"{species}_n_genes"] = df[species].map(count_genes)
        summary[f"{species}_present"] = summary[f"{species}_n_genes"] > 0

    count_columns = [
        f"{species}_n_genes"
        for species in species_columns
    ]

    present_columns = [
        f"{species}_present"
        for species in species_columns
    ]

    summary["n_species"] = summary[present_columns].sum(axis=1)
    summary["n_genes_total"] = summary[count_columns].sum(axis=1)
    summary["has_paralogs"] = summary[count_columns].gt(1).any(axis=1)

    def species_combo(row):
        present = [
            species
            for species in species_columns
            if row[f"{species}_present"]
        ]
        return "+".join(present)

    summary["species_combo"] = summary.apply(
        species_combo,
        axis=1,
    )

    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Cleaned Orthogroups TSV.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output summary TSV.",
    )
    args = parser.parse_args()

    df = pd.read_csv(
        args.input,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    summary = summarize_orthogroups(df)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(
        args.output,
        sep="\t",
        index=False,
    )

    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
