#!/usr/bin/env python3

"""
Create a reusable protein FASTA containing one longest protein isoform
per gene, with the FASTA header replaced by the gene symbol.

Designed for Ensembl peptide FASTA headers containing fields such as:

    >ENSP00000354587.4 ... gene:ENSG00000134853.12 ...
    gene_symbol:PDGFRA ...

Output:

    >PDGFRA
    MAAA...

Usage:

    python longest_protein_gene_symbols.py \
        --input Homo_sapiens.GRCh38.pep.all.fa \
        --output hu.fa

    python longest_protein_gene_symbols.py \
        --input Mus_musculus.GRCm39.pep.all.fa \
        --output mo.fa
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import TextIO


GENE_SYMBOL_PATTERNS = (
    re.compile(r"(?:^|\s)gene_symbol:([^\s]+)"),
    re.compile(r"(?:^|\s)gene_symbol=([^\s]+)"),
)

GENE_ID_PATTERNS = (
    re.compile(r"(?:^|\s)gene:([^\s]+)"),
    re.compile(r"(?:^|\s)gene=([^\s]+)"),
)


def open_text(path: Path) -> TextIO:
    """Open an uncompressed or gzip-compressed text file."""

    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt")

    return path.open("r")


def extract_field(
    header: str,
    patterns: tuple[re.Pattern[str], ...],
) -> str | None:
    """Extract the first matching value from a FASTA header."""

    for pattern in patterns:
        match = pattern.search(header)

        if match:
            value = match.group(1).strip()

            if value and value not in {"NA", "None", "."}:
                return value

    return None


def read_fasta(path: Path):
    """
    Yield tuples of:

        full_header, sequence

    The leading '>' is removed from the returned header.
    """

    with open_text(path) as handle:
        header = None
        sequence_parts: list[str] = []

        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)

                header = line[1:]
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(
                        "Sequence encountered before the first FASTA header."
                    )

                sequence_parts.append(line)

        if header is not None:
            yield header, "".join(sequence_parts)


def wrap_sequence(sequence: str, width: int = 60) -> str:
    """Wrap a protein sequence across fixed-width FASTA lines."""

    return "\n".join(
        sequence[i : i + width]
        for i in range(0, len(sequence), width)
    )


def create_longest_protein_fasta(
    input_fasta: Path,
    output_fasta: Path,
    fallback_to_gene_id: bool = False,
) -> None:
    """
    Retain the longest protein isoform for every gene symbol.

    Ties are resolved by retaining the first protein encountered.
    """

    # gene symbol -> (sequence, original header)
    longest_by_gene: dict[str, tuple[str, str]] = {}

    total_proteins = 0
    proteins_with_symbol = 0
    fallback_gene_ids = 0
    unidentified = 0
    replaced_by_longer_isoform = 0

    for header, sequence in read_fasta(input_fasta):
        total_proteins += 1

        gene_symbol = extract_field(
            header,
            GENE_SYMBOL_PATTERNS,
        )

        if gene_symbol is not None:
            proteins_with_symbol += 1

        elif fallback_to_gene_id:
            gene_symbol = extract_field(
                header,
                GENE_ID_PATTERNS,
            )

            if gene_symbol is not None:
                # Strip an Ensembl version suffix:
                # ENSG00000134853.12 -> ENSG00000134853
                gene_symbol = re.sub(
                    r"\.\d+$",
                    "",
                    gene_symbol,
                )

                fallback_gene_ids += 1

        if gene_symbol is None:
            unidentified += 1
            continue

        if not sequence:
            continue

        current = longest_by_gene.get(gene_symbol)

        if current is None:
            longest_by_gene[gene_symbol] = (
                sequence,
                header,
            )

        elif len(sequence) > len(current[0]):
            longest_by_gene[gene_symbol] = (
                sequence,
                header,
            )

            replaced_by_longer_isoform += 1

    if not longest_by_gene:
        raise RuntimeError(
            "No genes were recovered. Inspect the FASTA headers with:\n"
            "grep '^>' input.fa | head\n\n"
            "The script expects a field such as gene_symbol:PDGFRA."
        )

    output_fasta.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_fasta.open("w") as output_handle:
        for gene_symbol in sorted(longest_by_gene):
            sequence, _ = longest_by_gene[gene_symbol]

            # SAMap-compatible header: exact gene symbol only.
            output_handle.write(f">{gene_symbol}\n")
            output_handle.write(
                wrap_sequence(sequence)
            )
            output_handle.write("\n")

    print("\nFinished")
    print("--------")
    print(f"Input FASTA: {input_fasta}")
    print(f"Output FASTA: {output_fasta}")
    print(f"Total protein records: {total_proteins}")
    print(
        "Protein records with gene_symbol:",
        proteins_with_symbol,
    )
    print(
        "Records using gene-ID fallback:",
        fallback_gene_ids,
    )
    print(
        "Records without an identifiable gene:",
        unidentified,
    )
    print(
        "Longer isoforms replacing earlier isoforms:",
        replaced_by_longer_isoform,
    )
    print(
        "Unique genes written:",
        len(longest_by_gene),
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Keep the longest protein isoform per gene and replace "
            "each FASTA header with the Ensembl gene symbol."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input Ensembl peptide FASTA, optionally gzip-compressed.",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output FASTA with one longest protein per gene symbol.",
    )

    parser.add_argument(
        "--fallback-to-gene-id",
        action="store_true",
        help=(
            "Use the Ensembl gene ID when gene_symbol is unavailable. "
            "Do not normally use this for gene-symbol scRNA-seq data."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    if not args.input.exists():
        raise FileNotFoundError(
            f"Input FASTA not found: {args.input}"
        )

    create_longest_protein_fasta(
        input_fasta=args.input,
        output_fasta=args.output,
        fallback_to_gene_id=args.fallback_to_gene_id,
    )


if __name__ == "__main__":
    main()