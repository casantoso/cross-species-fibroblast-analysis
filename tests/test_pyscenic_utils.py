import numpy as np
import pandas as pd

from python.pyscenic.pyscenic_utils import (
    compare_regulon_targets,
    extract_auc_matrix,
    mean_regulon_activity,
)


class FakeScenicAdata:
    def __init__(self):
        self.obs = pd.DataFrame(
            {
                "Regulon(TF_A(+))": [0.1, 0.2, 0.8],
                "Regulon(TF_B(+))": [0.4, 0.5, 0.3],
                "group": ["A", "A", "B"],
            },
            index=["c1", "c2", "c3"],
        )


def test_extract_auc_matrix_keeps_only_regulons():
    adata = FakeScenicAdata()
    auc = extract_auc_matrix(adata)

    assert list(auc.columns) == [
        "Regulon(TF_A(+))",
        "Regulon(TF_B(+))",
    ]
    assert list(auc.index) == ["c1", "c2", "c3"]


def test_mean_regulon_activity_groups_cells():
    adata = FakeScenicAdata()
    auc = extract_auc_matrix(adata)
    groups = adata.obs["group"]

    mean_auc = mean_regulon_activity(auc, groups)

    assert np.isclose(mean_auc.loc["A", "Regulon(TF_A(+))"], 0.15)
    assert np.isclose(mean_auc.loc["B", "Regulon(TF_A(+))"], 0.8)


def test_compare_regulon_targets_direct_overlap():
    human = pd.DataFrame({
        "TF": ["TFX", "TFX", "TFX"],
        "target_gene": ["A", "B", "C"],
    })
    mouse = pd.DataFrame({
        "TF": ["Tfx", "Tfx", "Tfx"],
        "target_gene": ["B", "C", "D"],
    })

    result = compare_regulon_targets(
        human,
        mouse,
        human_tf="TFX",
        mouse_tf="Tfx",
    )

    assert result["shared_targets"] == {"B", "C"}
    assert np.isclose(result["jaccard"], 0.5)


def test_compare_regulon_targets_with_orthology_table():
    human = pd.DataFrame({
        "TF": ["TFX", "TFX"],
        "target_gene": ["H_A", "H_B"],
    })
    mouse = pd.DataFrame({
        "TF": ["Tfx", "Tfx"],
        "target_gene": ["M_A", "M_C"],
    })
    mapping = pd.DataFrame({
        "human_gene": ["H_A", "H_B"],
        "mouse_gene": ["M_A", "M_B"],
    })

    result = compare_regulon_targets(
        human,
        mouse,
        human_tf="TFX",
        mouse_tf="Tfx",
        target_map=mapping,
    )

    assert len(result["mapped_shared_pairs"]) == 1
    assert result["mapped_shared_pairs"].iloc[0]["human_gene"] == "H_A"
