#!/usr/bin/env python3
"""Verify downloaded cisTarget ranking databases against Aerts Lab checksums.

The tracked resource manifest contains only fixed metadata already known for
this project. This script writes a generated checksum report under `outputs/`
so the repository does not contain placeholder fields waiting to be filled.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import urllib.request

import pandas as pd


OFFICIAL_DB_SHA256_URL = (
    "https://resources.aertslab.org/cistarget/databases/sha256sum.txt"
)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def parse_sha256sum_text(text: str) -> dict[str, str]:
    checksums = {}
    for line in text.splitlines():
        parts = line.strip().split(maxsplit=1)
        if len(parts) != 2:
            continue
        sha256, filename = parts
        checksums[Path(filename.lstrip("*").strip()).name] = sha256
    return checksums


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("python/pyscenic/resource_manifest.tsv"),
    )
    parser.add_argument(
        "--resource-root",
        type=Path,
        default=Path("resources/pyscenic"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/pyscenic/resource_checksums.tsv"),
    )
    args = parser.parse_args()

    manifest = pd.read_csv(
        args.manifest,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    with urllib.request.urlopen(OFFICIAL_DB_SHA256_URL) as response:
        official = parse_sha256sum_text(
            response.read().decode("utf-8")
        )

    records = []

    for row in manifest.itertuples(index=False):
        local_path = args.resource_root / row.species / row.filename

        local_hash = (
            sha256_file(local_path)
            if local_path.exists()
            else ""
        )

        official_hash = (
            official.get(row.filename, "")
            if row.resource_type == "cisTarget ranking database"
            else ""
        )

        if not local_path.exists():
            status = "MISSING_LOCAL_FILE"
        elif official_hash:
            status = "VERIFIED" if local_hash == official_hash else "MISMATCH"
        else:
            status = "LOCAL_HASH_RECORDED"

        records.append({
            "species": row.species,
            "resource_type": row.resource_type,
            "filename": row.filename,
            "local_sha256": local_hash,
            "official_sha256": official_hash,
            "checksum_status": status,
        })

    report = pd.DataFrame(records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.output, sep="\t", index=False)

    print(f"Saved checksum report: {args.output}")


if __name__ == "__main__":
    main()
