#!/usr/bin/env python3
"""Remove UNMAPPED entries from a symbol-mapped OrthoFinder orthogroup table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


def remove_unmapped(cell: str) -> str:
    """
    Remove `UNMAPPED|...` entries from one species cell.

    Example:
        PDGFRA|ENSG..., UNMAPPED|ENSG..., KIT|ENSG...
    becomes:
        PDGFRA|ENSG..., KIT|ENSG...
    """
    if cell is None or not str(cell).strip():
        return ""

    genes = [
        gene.strip()
        for gene in str(cell).split(",")
        if gene.strip()
    ]

    mapped = [
        gene
        for gene in genes
        if not gene.startswith("UNMAPPED|")
    ]

    return ", ".join(mapped)


def clean_orthogroups(input_file: Path) -> tuple[pd.DataFrame, dict]:
    """Clean a symbol-mapped Orthogroups table and return summary statistics."""
    df = pd.read_csv(
        input_file,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    if "Orthogroup" not in df.columns:
        raise ValueError("Input table must contain an `Orthogroup` column.")

    species_columns = [
        column
        for column in df.columns
        if column != "Orthogroup"
    ]

    removed_genes = 0

    cleaned = df.copy()

    for species in species_columns:
        original_counts = (
            cleaned[species]
            .map(lambda x: len([g for g in str(x).split(",") if g.strip()]))
        )

        cleaned[species] = cleaned[species].map(remove_unmapped)

        cleaned_counts = (
            cleaned[species]
            .map(lambda x: len([g for g in str(x).split(",") if g.strip()]))
        )

        removed_genes += int(
            (original_counts - cleaned_counts).sum()
        )

    has_mapped_gene = (
        cleaned[species_columns]
        .apply(
            lambda row: any(str(value).strip() for value in row),
            axis=1,
        )
    )

    removed_rows = int((~has_mapped_gene).sum())
    cleaned = cleaned.loc[has_mapped_gene].reset_index(drop=True)

    summary = {
        "removed_genes": removed_genes,
        "removed_rows": removed_rows,
        "retained_orthogroups": len(cleaned),
    }

    return cleaned, summary


def write_excel(df: pd.DataFrame, output_file: Path) -> None:
    """Write a readable Excel copy of the cleaned orthogroup table."""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(
        output_file,
        engine="openpyxl",
    ) as writer:
        df.to_excel(
            writer,
            sheet_name="Orthogroups",
            index=False,
        )

        worksheet = writer.sheets["Orthogroups"]
        worksheet.freeze_panes = "B2"
        worksheet.auto_filter.ref = worksheet.dimensions
        worksheet.column_dimensions["A"].width = 16

        for column_number in range(2, len(df.columns) + 1):
            column_letter = worksheet.cell(
                row=1,
                column=column_number,
            ).column_letter
            worksheet.column_dimensions[column_letter].width = 45

        for row in worksheet.iter_rows():
            for cell in row:
                cell.alignment = cell.alignment.copy(
                    vertical="top",
                    wrap_text=True,
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Symbol-mapped Orthogroups TSV.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Cleaned TSV output.",
    )
    parser.add_argument(
        "--xlsx",
        type=Path,
        default=None,
        help="Optional Excel copy of the cleaned table.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    cleaned, summary = clean_orthogroups(args.input)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cleaned.to_csv(
        args.output,
        sep="\t",
        index=False,
    )

    if args.xlsx is not None:
        write_excel(
            cleaned,
            args.xlsx,
        )

    print(f"Clean TSV saved as: {args.output}")
    if args.xlsx is not None:
        print(f"Excel file saved as: {args.xlsx}")
    print(f"Unmapped gene entries removed: {summary['removed_genes']:,}")
    print(f"Completely empty orthogroups removed: {summary['removed_rows']:,}")
    print(f"Orthogroups retained: {summary['retained_orthogroups']:,}")


if __name__ == "__main__":
    main()
