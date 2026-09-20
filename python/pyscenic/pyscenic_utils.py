"""Reusable pySCENIC post-processing utilities."""

from __future__ import annotations

import json
import importlib.metadata
import sys
from pathlib import Path
import re

import numpy as np
import pandas as pd
import scipy.sparse as sp


def load_pyscenic_config(path="config/pyscenic_config.json"):
    """Load the shared human/mouse pySCENIC configuration."""
    with Path(path).open() as handle:
        return json.load(handle)


def resolve_species_paths(config, species, project_dir="."):
    """Resolve project-relative pySCENIC input, resource, and result paths."""
    project_dir = Path(project_dir)
    cfg = config[species]

    resource_dir = project_dir / cfg["resources_dir"]
    result_dir = project_dir / cfg["results_dir"]

    paths = {
        "expression": project_dir / cfg["input_h5ad"],
        "results": result_dir,
        "tf_list": resource_dir / cfg["tf_list"],
        "db_10kb": resource_dir / cfg["db_10kb"],
        "db_500bp": resource_dir / cfg["db_500bp"],
        "motif_annotations": resource_dir / cfg["motif_annotations"],
    }

    result_dir.mkdir(parents=True, exist_ok=True)

    for subdir in [
        "tables",
        "figures",
        "logs",
        "objects",
        "metadata",
    ]:
        (result_dir / subdir).mkdir(exist_ok=True)

    paths.update({
        "tables": result_dir / "tables",
        "figures": result_dir / "figures",
        "logs": result_dir / "logs",
        "objects": result_dir / "objects",
        "metadata": result_dir / "metadata",
    })

    return paths


def validate_pyscenic_inputs(paths, group_col=None):
    """Check required input/resource files and optional metadata before a run."""
    import anndata as ad

    required = [
        "expression",
        "tf_list",
        "db_10kb",
        "db_500bp",
        "motif_annotations",
    ]

    missing = [
        f"{key}: {paths[key]}"
        for key in required
        if not Path(paths[key]).exists()
    ]

    if missing:
        raise FileNotFoundError(
            "Missing pySCENIC input/resource files:\n- "
            + "\n- ".join(missing)
        )

    if group_col:
        adata = ad.read_h5ad(paths["expression"], backed="r")
        if group_col not in adata.obs.columns:
            raise KeyError(
                f"Metadata column {group_col!r} was not found in the input AnnData."
            )
        if getattr(adata, "file", None) is not None:
            adata.file.close()

    return True



def package_versions(packages=None):
    """Return installed versions relevant to a pySCENIC run."""
    packages = packages or [
        "pyscenic",
        "arboreto",
        "ctxcore",
        "scanpy",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "loompy",
        "dask",
        "distributed",
        "umap-learn",
        "matplotlib",
        "pyarrow",
    ]

    versions = {"python": sys.version.split()[0]}

    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None

    return versions


def save_run_metadata(
    output_file,
    parameters,
    resources=None,
    extra=None,
):
    """Save pySCENIC parameters, resources, and package versions as JSON."""
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "parameters": parameters,
        "resources": resources or {},
        "package_versions": package_versions(),
    }

    if extra is not None:
        payload["extra"] = extra

    with output_file.open("w") as handle:
        json.dump(payload, handle, indent=2, default=str)

    return payload


def tf_overlap(adata, tf_list_file):
    """Report overlap between expression features and the supplied TF list."""
    with Path(tf_list_file).open() as handle:
        tfs = {
            line.strip()
            for line in handle
            if line.strip()
        }

    genes = set(adata.var_names.astype(str))
    overlap = genes & tfs

    return {
        "n_genes": len(genes),
        "n_tfs": len(tfs),
        "n_overlap": len(overlap),
        "overlapping_tfs": sorted(overlap),
    }


def extract_auc_matrix(scenic_adata):
    """
    Extract cell x regulon AUCell scores.

    Supports pySCENIC outputs that store regulon scores in obs columns named
    `Regulon(...)`.
    """
    regulon_cols = [
        column
        for column in scenic_adata.obs.columns
        if str(column).startswith("Regulon(")
    ]

    if not regulon_cols:
        raise ValueError(
            "No `Regulon(...)` AUCell columns were found in scenic_adata.obs."
        )

    auc = scenic_adata.obs[regulon_cols].copy()
    auc.index.name = "cell_id"
    return auc


def run_regulon_umap(
    auc_matrix,
    n_neighbors=10,
    min_dist=0.4,
    metric="correlation",
    random_state=123,
):
    """Create a UMAP embedding from cell-by-regulon AUCell activity."""
    import umap

    coordinates = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state,
    ).fit_transform(auc_matrix.to_numpy())

    return pd.DataFrame(
        coordinates,
        columns=["SCENIC_UMAP_1", "SCENIC_UMAP_2"],
        index=auc_matrix.index,
    )


def mean_regulon_activity(auc_matrix, groups):
    """Calculate mean AUCell activity for every regulon within each group."""
    groups = pd.Series(groups, index=auc_matrix.index, name="group")
    return (
        auc_matrix
        .join(groups)
        .groupby("group", observed=True)
        .mean()
    )


def calculate_rss(auc_matrix, groups):
    """Calculate regulon specificity scores after removing invariant regulons."""
    from pyscenic.rss import regulon_specificity_scores

    groups = pd.Series(groups, index=auc_matrix.index).astype(str)

    variable = auc_matrix.nunique(dropna=True) > 1
    filtered = auc_matrix.loc[:, variable].copy()

    return regulon_specificity_scores(
        filtered,
        groups,
    )


def top_rss_regulons(rss, top_n=15):
    """Return a long table of the highest-RSS regulons for every group."""
    rows = []

    for group in rss.index:
        scores = (
            rss.loc[group]
            .sort_values(ascending=False)
            .head(top_n)
        )

        for rank, (regulon, score) in enumerate(scores.items(), start=1):
            rows.append({
                "group": group,
                "rank": rank,
                "regulon": regulon,
                "RSS": score,
            })

    return pd.DataFrame(rows)


def export_regulon_targets(regulon_file, output_file=None):
    """Convert a pySCENIC regulon file into one TF-target row per edge."""
    from pyscenic.cli.utils import load_signatures

    regulons = load_signatures(str(regulon_file))
    rows = []

    for regulon in regulons:
        for gene, weight in regulon.gene2weight.items():
            rows.append({
                "regulon": regulon.name,
                "TF": regulon.transcription_factor,
                "target_gene": gene,
                "weight": weight,
            })

    table = pd.DataFrame(rows)

    if output_file:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(output_file, index=False)

    return table


def plot_embedding_by_group(
    embedding,
    groups,
    title="SCENIC regulon-activity UMAP",
    save_path=None,
):
    """Plot the SCENIC regulon-activity UMAP colored by a metadata group."""
    import matplotlib.pyplot as plt

    plot_df = embedding.copy()
    plot_df["group"] = pd.Series(groups, index=embedding.index).astype(str)

    fig, ax = plt.subplots(figsize=(7, 6), dpi=150)

    for group, sub in plot_df.groupby("group", observed=True):
        ax.scatter(
            sub.iloc[:, 0],
            sub.iloc[:, 1],
            s=8,
            alpha=0.7,
            label=group,
        )

    ax.set_xlabel(embedding.columns[0])
    ax.set_ylabel(embedding.columns[1])
    ax.set_title(title)
    ax.legend(
        title="Group",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax


def plot_regulon_on_embedding(
    embedding,
    auc_matrix,
    regulon,
    title=None,
    save_path=None,
):
    """Color a UMAP/embedding by one regulon's AUCell activity."""
    import matplotlib.pyplot as plt

    if regulon not in auc_matrix.columns:
        raise KeyError(f"{regulon!r} is not present in the AUCell matrix.")

    plot_df = embedding.join(auc_matrix[[regulon]], how="inner")

    fig, ax = plt.subplots(figsize=(7, 6), dpi=150)
    scatter = ax.scatter(
        plot_df.iloc[:, 0],
        plot_df.iloc[:, 1],
        c=plot_df[regulon],
        s=8,
    )

    ax.set_xlabel(embedding.columns[0])
    ax.set_ylabel(embedding.columns[1])
    ax.set_title(title or f"{regulon} activity")
    fig.colorbar(scatter, ax=ax, label="AUCell score")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax


def plot_rss_heatmap(rss, top_n=15, save_path=None):
    """Plot RSS values for the union of the top regulons from all groups."""
    import matplotlib.pyplot as plt

    selected = []
    for group in rss.index:
        selected.extend(
            rss.loc[group]
            .sort_values(ascending=False)
            .head(top_n)
            .index
            .tolist()
        )

    selected = list(dict.fromkeys(selected))
    data = rss.loc[:, selected]

    fig, ax = plt.subplots(
        figsize=(max(8, len(selected) * 0.25), max(3, len(rss) * 0.8)),
        dpi=150,
    )
    image = ax.imshow(data.to_numpy(), aspect="auto")

    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels(data.index)
    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels(data.columns, rotation=90, fontsize=7)
    ax.set_xlabel("Regulon")
    ax.set_ylabel("Group")
    ax.set_title("Top group-specific regulons by RSS")
    fig.colorbar(image, ax=ax, label="RSS")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax, data


def plot_mean_auc_heatmap(
    mean_auc,
    selected_regulons,
    save_path=None,
):
    """Plot mean AUCell activity for a selected regulon set."""
    import matplotlib.pyplot as plt

    selected = [
        regulon
        for regulon in selected_regulons
        if regulon in mean_auc.columns
    ]
    data = mean_auc.loc[:, selected]

    fig, ax = plt.subplots(
        figsize=(max(8, len(selected) * 0.25), max(3, len(mean_auc) * 0.8)),
        dpi=150,
    )
    image = ax.imshow(data.to_numpy(), aspect="auto")

    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels(data.index)
    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels(data.columns, rotation=90, fontsize=7)
    ax.set_xlabel("Regulon")
    ax.set_ylabel("Group")
    ax.set_title("Mean AUCell activity")
    fig.colorbar(image, ax=ax, label="Mean AUCell score")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax, data


def load_original_umap(path):
    """Load the UMAP coordinates exported from Seurat."""
    umap = pd.read_csv(path, sep="\t")
    id_col = "cell_id" if "cell_id" in umap.columns else "CellID"
    umap = umap.set_index(id_col)

    coordinate_cols = [
        column
        for column in umap.columns
        if "UMAP" in str(column).upper()
    ][:2]

    if len(coordinate_cols) < 2:
        coordinate_cols = list(umap.columns[:2])

    result = umap[coordinate_cols].copy()
    result.columns = ["UMAP_1", "UMAP_2"]
    return result


def compare_regulon_targets(
    human_targets,
    mouse_targets,
    human_tf,
    mouse_tf=None,
    target_map=None,
):
    """
    Compare target sets for a homologous TF pair.

    If `target_map` is supplied, it should contain columns:
        human_gene, mouse_gene

    This keeps the comparison method public without hard-coding any
    project-specific TF or marker result.
    """
    mouse_tf = mouse_tf or human_tf

    h = set(
        human_targets.loc[
            human_targets["TF"].astype(str).str.upper() == str(human_tf).upper(),
            "target_gene",
        ].dropna().astype(str)
    )

    m = set(
        mouse_targets.loc[
            mouse_targets["TF"].astype(str).str.upper() == str(mouse_tf).upper(),
            "target_gene",
        ].dropna().astype(str)
    )

    if target_map is None:
        shared = h & m
        union = h | m
        return {
            "human_targets": h,
            "mouse_targets": m,
            "shared_targets": shared,
            "jaccard": len(shared) / len(union) if union else np.nan,
        }

    mapping = target_map[
        target_map["human_gene"].isin(h)
        & target_map["mouse_gene"].isin(m)
    ].drop_duplicates()

    n_union = len(h) + len(m) - len(mapping)

    return {
        "human_targets": h,
        "mouse_targets": m,
        "mapped_shared_pairs": mapping,
        "jaccard_mapped": len(mapping) / n_union if n_union else np.nan,
    }


def find_relevant_regulon_targets(
    adata,
    regulon_targets,
    tf,
    group,
    group_col="fib_subtype",
    auc_threshold=0.01,
):
    """
    Rank one regulon's targets by how often each target appears among the
    top-expressed genes in a selected cell group versus all other cells.

    This generalizes the exploratory human/mouse target-ranking code without
    embedding specific TFs or subtype names.
    """
    tf_targets = (
        regulon_targets.loc[
            regulon_targets["TF"].astype(str).str.upper() == str(tf).upper(),
            ["target_gene", "weight"],
        ]
        .dropna(subset=["target_gene"])
        .drop_duplicates("target_gene")
        .copy()
    )

    if tf_targets.empty:
        raise ValueError(f"No final regulon targets were found for {tf!r}.")

    present = [
        gene
        for gene in tf_targets["target_gene"]
        if gene in adata.var_names
    ]
    if not present:
        raise ValueError("None of the regulon targets are present in the AnnData.")

    X = adata.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)

    n_top = max(1, int(np.floor(adata.n_vars * auc_threshold)))
    group_mask = adata.obs[group_col].astype(str).to_numpy() == str(group)

    target_indices = adata.var_names.get_indexer(present)
    target_set = set(target_indices)
    idx_to_gene = dict(zip(target_indices, present))

    in_group = {gene: 0 for gene in present}
    out_group = {gene: 0 for gene in present}

    for cell_idx in range(X.shape[0]):
        start, end = X.indptr[cell_idx], X.indptr[cell_idx + 1]
        indices = X.indices[start:end]
        values = X.data[start:end]

        if len(values) > n_top:
            keep = np.argpartition(values, -n_top)[-n_top:]
            indices = indices[keep]

        for idx in set(indices) & target_set:
            gene = idx_to_gene[idx]
            if group_mask[cell_idx]:
                in_group[gene] += 1
            else:
                out_group[gene] += 1

    n_in = int(group_mask.sum())
    n_out = int((~group_mask).sum())

    result = tf_targets.copy()
    result["fraction_top_in_group"] = result["target_gene"].map(
        lambda x: in_group.get(x, 0) / max(n_in, 1)
    )
    result["fraction_top_out_group"] = result["target_gene"].map(
        lambda x: out_group.get(x, 0) / max(n_out, 1)
    )
    result["difference"] = (
        result["fraction_top_in_group"]
        - result["fraction_top_out_group"]
    )

    return result.sort_values(
        ["difference", "weight"],
        ascending=[False, False],
    )
