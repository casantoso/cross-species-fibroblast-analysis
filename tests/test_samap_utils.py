import numpy as np
import pandas as pd

from python.samap.samap_analysis_utils import (
    add_species_combo,
    build_gene_pair_network,
    combine_paralog_expression,
    find_gene_candidates_in_source,
    summarize_connected_components,
)


class FakeAdata:
    def __init__(self, var_names):
        self.var_names = pd.Index(var_names)


def test_combine_paralog_expression_takes_cellwise_maximum():
    a = np.array([1, 0, 4])
    b = np.array([0, 3, 2])

    observed = combine_paralog_expression([a, b])

    np.testing.assert_array_equal(
        observed,
        np.array([1, 3, 4]),
    )


def test_gene_matching_detects_exact_and_zebrafish_ab_paralogs():
    adata = FakeAdata(
        ["other", "genea", "geneb", "unrelated"]
    )

    observed = find_gene_candidates_in_source(
        adata,
        species_id="ze",
        anchor_gene="gene",
        paralog_species=("ze",),
    )

    assert observed == ["genea", "geneb"]


def test_species_combo_label_is_exact():
    component_df = pd.DataFrame([{
        "has_hu": True,
        "has_mo": True,
        "has_pi": False,
        "has_ch": False,
        "has_ze": True,
        "has_ti": False,
    }])

    observed = add_species_combo(
        component_df,
        species_order=["hu", "mo", "pi", "ch", "ze", "ti"],
    )

    assert observed.loc[0, "species_combo"] == "hu+mo+ze"


def test_network_construction_and_connected_components():
    pairs = pd.DataFrame([
        {
            "species1": "hu",
            "species2": "mo",
            "gene1": "hu_A",
            "gene2": "mo_A",
            "score": 0.8,
        },
        {
            "species1": "mo",
            "species2": "ze",
            "gene1": "mo_A",
            "gene2": "ze_A",
            "score": 0.7,
        },
        {
            "species1": "hu",
            "species2": "pi",
            "gene1": "hu_B",
            "gene2": "pi_B",
            "score": 0.6,
        },
    ])

    graph = build_gene_pair_network(pairs)

    assert graph.number_of_nodes() == 5
    assert graph.number_of_edges() == 3

    summary = summarize_connected_components(
        graph,
        species_order=["hu", "mo", "pi", "ze"],
    )

    assert len(summary) == 2
    assert set(summary["n_species"]) == {2, 3}
