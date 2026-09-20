#!/usr/bin/env python3
"""Build an identifier-to-gene-symbol mapping table from a GTF/GFF-like file.

The output has two columns:

    id    symbol

Typical uses:
- gene_id -> gene_name
- protein_id -> gene_name
- transcript_id -> gene_name

The script intentionally does not guess species-specific naming rules. It
extracts the requested attributes exactly as provided by the annotation file.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_attributes(text: str) -> dict[str, str]:
    """Parse common GTF (`key "value"`) and GFF (`key=value`) attributes."""
    attrs: dict[str, str] = {}

    for field in text.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue

        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
            continue

        match = re.match(r'^(\S+)\s+"?(.*?)"?$', field)
        if match:
            attrs[match.group(1)] = match.group(2).strip().strip('"')

    return attrs


def build_mapping(
    annotation_file: Path,
    id_attribute: str = "gene_id",
    symbol_attributes: tuple[str, ...] = (
        "gene_name",
        "gene_symbol",
        "Name",
        "gene",
    ),
) -> pd.DataFrame:
    """Extract unique identifier -> symbol pairs."""
    records: list[tuple[str, str]] = []

    with annotation_file.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            attrs = parse_attributes(fields[8])
            identifier = attrs.get(id_attribute)

            if not identifier:
                continue

            symbol = None
            for key in symbol_attributes:
                value = attrs.get(key)
                if value:
                    symbol = value
                    break

            if not symbol:
                continue

            records.append((identifier, symbol))

    if not records:
        raise ValueError(
            f"No mappings found for id_attribute={id_attribute!r}. "
            "Check the annotation attribute names."
        )

    mapping = pd.DataFrame(records, columns=["id", "symbol"]).drop_duplicates()

    conflicts = (
        mapping.groupby("id")["symbol"]
        .nunique()
        .loc[lambda x: x > 1]
    )

    if len(conflicts):
        examples = ", ".join(conflicts.index[:10])
        raise ValueError(
            "Some identifiers map to more than one symbol. "
            f"Resolve these before continuing. Example IDs: {examples}"
        )

    mapping = mapping.drop_duplicates("id").sort_values("id").reset_index(drop=True)
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--id-attribute",
        default="gene_id",
        help="Annotation attribute used by the OrthoFinder FASTA IDs.",
    )
    parser.add_argument(
        "--symbol-attributes",
        nargs="+",
        default=["gene_name", "gene_symbol", "Name", "gene"],
        help="Symbol attributes to try in order.",
    )
    args = parser.parse_args()

    mapping = build_mapping(
        annotation_file=args.annotation,
        id_attribute=args.id_attribute,
        symbol_attributes=tuple(args.symbol_attributes),
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(args.output, sep="\t", index=False)

    print(f"Saved {len(mapping):,} ID-to-symbol mappings to {args.output}")


if __name__ == "__main__":
    main()
