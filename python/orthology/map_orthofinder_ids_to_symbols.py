#!/usr/bin/env python3
"""Annotate raw OrthoFinder Orthogroups.tsv IDs with gene symbols.

Expected mapping directory:
    mappings/
      human.tsv
      mouse.tsv
      pig.tsv
      chicken.tsv
      zebrafish.tsv
      tilapia.tsv

Each mapping file must contain:
    id    symbol

The species filenames must match the species column names in Orthogroups.tsv.
Each mapped entry is written as:

    SYMBOL|ORIGINAL_ID

Unmapped entries are retained as:

    UNMAPPED|ORIGINAL_ID

This preserves the original identifier while making downstream tables readable.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def split_orthogroup_cell(value: str) -> list[str]:
    """Split a comma-separated OrthoFinder species cell."""
    if value is None or not str(value).strip():
        return []
    return [x.strip() for x in str(value).split(",") if x.strip()]


def load_mapping(path: Path) -> dict[str, str]:
    mapping = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    required = {"id", "symbol"}
    missing = required - set(mapping.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {', '.join(sorted(missing))}"
        )

    mapping = mapping.loc[
        (mapping["id"].str.len() > 0) &
        (mapping["symbol"].str.len() > 0),
        ["id", "symbol"],
    ].drop_duplicates()

    conflicts = (
        mapping.groupby("id")["symbol"]
        .nunique()
        .loc[lambda x: x > 1]
    )
    if len(conflicts):
        raise ValueError(
            f"{path} contains IDs mapped to multiple symbols. "
            f"Example: {conflicts.index[0]}"
        )

    return dict(zip(mapping["id"], mapping["symbol"]))


def annotate_cell(value: str, id_to_symbol: dict[str, str]) -> tuple[str, int, int]:
    """Annotate one species cell and return text, mapped count, unmapped count."""
    ids = split_orthogroup_cell(value)
    annotated: list[str] = []
    n_mapped = 0
    n_unmapped = 0

    for identifier in ids:
        symbol = id_to_symbol.get(identifier)
        if symbol:
            annotated.append(f"{symbol}|{identifier}")
            n_mapped += 1
        else:
            annotated.append(f"UNMAPPED|{identifier}")
            n_unmapped += 1

    return ", ".join(annotated), n_mapped, n_unmapped


def annotate_orthogroups(
    orthogroups: pd.DataFrame,
    mapping_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if "Orthogroup" not in orthogroups.columns:
        raise ValueError("Input must contain an `Orthogroup` column.")

    species_columns = [
        c for c in orthogroups.columns
        if c != "Orthogroup"
    ]

    result = orthogroups.copy()
    summary_rows = []

    for species in species_columns:
        mapping_file = mapping_dir / f"{species}.tsv"
        if not mapping_file.exists():
            raise FileNotFoundError(
                f"Missing mapping file for {species}: {mapping_file}"
            )

        id_to_symbol = load_mapping(mapping_file)

        mapped_total = 0
        unmapped_total = 0
        annotated_values = []

        for value in result[species]:
            annotated, n_mapped, n_unmapped = annotate_cell(
                value,
                id_to_symbol,
            )
            annotated_values.append(annotated)
            mapped_total += n_mapped
            unmapped_total += n_unmapped

        result[species] = annotated_values

        total = mapped_total + unmapped_total
        summary_rows.append({
            "species": species,
            "mapped_ids": mapped_total,
            "unmapped_ids": unmapped_total,
            "total_ids": total,
            "fraction_mapped": mapped_total / total if total else 0.0,
        })

    summary = pd.DataFrame(summary_rows)
    return result, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--mapping-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="Optional TSV reporting mapping success for each species.",
    )
    args = parser.parse_args()

    orthogroups = pd.read_csv(
        args.orthogroups,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    annotated, summary = annotate_orthogroups(
        orthogroups,
        args.mapping_dir,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(args.output, sep="\t", index=False)

    if args.summary:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.summary, sep="\t", index=False)

    print(f"Saved symbol-annotated orthogroups to: {args.output}")
    if args.summary:
        print(f"Saved mapping summary to: {args.summary}")


if __name__ == "__main__":
    main()
