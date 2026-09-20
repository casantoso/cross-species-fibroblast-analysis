import pandas as pd

from python.orthology.clean_orthogroups import (
    clean_orthogroups,
    remove_unmapped,
)
from python.orthology.summarize_orthogroups import (
    summarize_orthogroups,
)


def test_remove_unmapped_entries():
    value = "PDGFRA|A, UNMAPPED|B, KIT|C"

    assert remove_unmapped(value) == "PDGFRA|A, KIT|C"


def test_summarize_orthogroups_counts_species_and_paralogs():
    df = pd.DataFrame({
        "Orthogroup": ["OG1", "OG2"],
        "human": ["A", "B, C"],
        "mouse": ["D", ""],
        "zebrafish": ["E, F", ""],
    })

    summary = summarize_orthogroups(df)

    row1 = summary.loc[summary["Orthogroup"] == "OG1"].iloc[0]
    row2 = summary.loc[summary["Orthogroup"] == "OG2"].iloc[0]

    assert row1["n_species"] == 3
    assert bool(row1["has_paralogs"]) is True
    assert row1["species_combo"] == "human+mouse+zebrafish"

    assert row2["n_species"] == 1
    assert bool(row2["has_paralogs"]) is True
