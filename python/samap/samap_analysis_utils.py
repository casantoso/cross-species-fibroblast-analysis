"""
Reusable utilities for cross-species SAMap analysis.

The functions in this module were generalized from exploratory notebooks.
They intentionally avoid project-specific marker names and unpublished
biological conclusions.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import json
from functools import reduce
import gc
import importlib.metadata
import re
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp


DEFAULT_SPECIES_ORDER = ["hu", "mo", "pi", "ch", "ze", "ti"]

DEFAULT_SPECIES_LABELS = {
    "hu": "Human",
    "mo": "Mouse",
    "pi": "Pig",
    "ch": "Chicken",
    "ze": "Zebrafish",
    "ti": "Tilapia",
}



# ====================================================================
# Configuration, validation, and reproducibility helpers
# ====================================================================

def load_samap_config(path="config/samap_config.json"):
    """Load the shared SAMap project configuration JSON."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"SAMap config not found: {path}")
    with path.open() as handle:
        return json.load(handle)


def species_files_from_config(
    config,
    workflow,
    project_dir=".",
):
    """
    Build the species -> .h5ad path dictionary from the shared configuration.

    workflow:
        "celltype" or "fibroblast"
    """
    if workflow not in {"celltype", "fibroblast"}:
        raise ValueError("workflow must be 'celltype' or 'fibroblast'.")

    key = f"{workflow}_h5ad_dir"
    base = Path(project_dir) / config["paths"][key]

    return {
        sid: base / f"{info['slug']}.h5ad"
        for sid, info in config["species"].items()
    }


def result_directories_from_config(
    config,
    workflow,
    project_dir=".",
):
    """
    Return standardized output directories for one SAMap workflow.

    Creates:
        objects/
        tables/
        figures/
        metadata/
    """
    if workflow not in {"celltype", "fibroblast"}:
        raise ValueError("workflow must be 'celltype' or 'fibroblast'.")

    key = f"{workflow}_results_dir"
    root = Path(project_dir) / config["paths"][key]

    dirs = {
        "root": root,
        "objects": root / "objects",
        "tables": root / "tables",
        "figures": root / "figures",
        "metadata": root / "metadata",
    }

    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)

    return dirs


def validate_h5ad_inputs(
    files,
    required_obs_columns=None,
):
    """
    Check that all expected .h5ad files exist and optionally contain required
    metadata columns before a long SAMap run starts.
    """
    import anndata as ad

    required_obs_columns = set(required_obs_columns or [])
    problems = []

    for sid, path in files.items():
        path = Path(path)

        if not path.exists():
            problems.append(f"{sid}: missing file {path}")
            continue

        if required_obs_columns:
            obj = ad.read_h5ad(path, backed="r")
            missing = required_obs_columns - set(obj.obs.columns)
            if missing:
                problems.append(
                    f"{sid}: missing metadata columns {sorted(missing)}"
                )
            if getattr(obj, "file", None) is not None:
                obj.file.close()

    if problems:
        raise ValueError(
            "SAMap input validation failed:\n- "
            + "\n- ".join(problems)
        )

    return True


def expected_map_pairs(species_ids):
    """Return all unique pairwise species IDs expected by SAMap."""
    return list(combinations(species_ids, 2))


def validate_blast_maps(
    maps_dir,
    species_ids,
):
    """
    Check that each expected SAMap map directory and both directional BLAST
    files exist before downstream GenePairFinder/BLAST annotation is run.
    """
    maps_dir = Path(maps_dir)
    problems = []

    for sid1, sid2 in expected_map_pairs(species_ids):
        pair_dir = maps_dir / f"{sid1}{sid2}"
        f12 = pair_dir / f"{sid1}_to_{sid2}.txt"
        f21 = pair_dir / f"{sid2}_to_{sid1}.txt"

        if not pair_dir.exists():
            problems.append(f"missing map directory: {pair_dir}")
            continue

        if not f12.exists():
            problems.append(f"missing BLAST file: {f12}")
        if not f21.exists():
            problems.append(f"missing BLAST file: {f21}")

    if problems:
        raise ValueError(
            "SAMap BLAST-map validation failed:\n- "
            + "\n- ".join(problems)
        )

    return True


def validate_integrated_metadata(
    adata,
    required_columns,
):
    """Check that the integrated AnnData contains required metadata columns."""
    missing = set(required_columns) - set(adata.obs.columns)

    if missing:
        raise ValueError(
            "Integrated AnnData is missing metadata columns: "
            + ", ".join(sorted(missing))
        )

    return True


def package_versions(
    packages=None,
):
    """Return installed package versions for reproducibility metadata."""
    packages = packages or [
        "samap",
        "scanpy",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "networkx",
        "upsetplot",
        "threadpoolctl",
        "igraph",
        "leidenalg",
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
    extra=None,
):
    """
    Save analysis parameters and installed package versions next to the results.
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "parameters": parameters,
        "package_versions": package_versions(),
    }

    if extra is not None:
        payload["extra"] = extra

    with output_file.open("w") as handle:
        json.dump(payload, handle, indent=2, default=str)

    return payload


def combine_paralog_expression(expression_arrays):
    """
    Combine multiple paralogs by taking the maximum expression per cell.

    This is a small pure-NumPy helper so the behavior is easy to test.
    """
    arrays = [
        np.asarray(values, dtype=float).ravel()
        for values in expression_arrays
    ]

    if not arrays:
        raise ValueError("At least one expression array is required.")

    lengths = {len(x) for x in arrays}
    if len(lengths) != 1:
        raise ValueError("All paralog expression arrays must have equal length.")

    if len(arrays) == 1:
        return arrays[0]

    return np.max(np.vstack(arrays), axis=0)


# ====================================================================
# Basic AnnData helpers
# ====================================================================

def load_adata(adata_or_path):
    """Return an AnnData object from an object or .h5ad path."""
    if isinstance(adata_or_path, (str, Path)):
        import scanpy as sc
        return sc.read_h5ad(adata_or_path)
    return adata_or_path


def find_gene_candidates_strict(
    adata,
    species_id,
    gene,
    separator="_",
):
    """
    Find exact/case-insensitive gene candidates for one species.

    SAMap commonly prefixes integrated feature names with species IDs,
    for example `hu_GENE`.
    """
    gene = str(gene)
    expected = f"{species_id}{separator}{gene}"

    exact = [
        x for x in (expected, gene)
        if x in adata.var_names
    ]
    if exact:
        return exact

    matches = [
        feature
        for feature in adata.var_names.astype(str)
        if feature.upper() == expected.upper()
    ]
    return matches


def find_gene_candidates(
    adata,
    species_id,
    gene,
    separator="_",
):
    """Return exact and partial candidate features for a requested gene."""
    strict = find_gene_candidates_strict(
        adata,
        species_id,
        gene,
        separator=separator,
    )
    if strict:
        return strict

    prefix = f"{species_id}{separator}"
    query = str(gene).upper()

    return [
        feature
        for feature in adata.var_names.astype(str)
        if feature.startswith(prefix)
        and query in feature.upper()
    ]


def extract_original_expression(
    adata,
    feature,
):
    """Extract a feature's expression vector from AnnData."""
    values = adata[:, feature].X
    if sp.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


# ====================================================================
# Integrated SAMap plotting: expression
# ====================================================================

def plot_samap_multispecies_expression(
    adata_or_path,
    genes,
    species_col="species",
    species_labels=None,
    species_markers=None,
    species_sizes=None,
    threshold=0,
    shared_scale=False,
    alpha=0.8,
    figsize=(10, 8),
    title=None,
    save_path=None,
    dpi=300,
):
    """
    Plot homologous genes from multiple species on one integrated SAMap UMAP.

    `genes` is a dict such as:
        {"hu": "GENE_A", "mo": "GeneA", "ze": "genea"}
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D

    adata = load_adata(adata_or_path)

    if "X_umap" not in adata.obsm:
        raise KeyError("'X_umap' was not found in adata.obsm.")
    if species_col not in adata.obs:
        raise KeyError(f"{species_col!r} was not found in adata.obs.")

    species_ids = list(genes)
    labels = {
        sid: (species_labels or {}).get(sid, DEFAULT_SPECIES_LABELS.get(sid, sid))
        for sid in species_ids
    }

    markers_default = ["o", "^", "s", "D", "P", "X", "v"]
    species_markers = {
        sid: (species_markers or {}).get(
            sid,
            markers_default[i % len(markers_default)],
        )
        for i, sid in enumerate(species_ids)
    }
    species_sizes = {
        sid: (species_sizes or {}).get(sid, 8)
        for sid in species_ids
    }

    species_vector = adata.obs[species_col].astype(str).to_numpy()
    umap = np.asarray(adata.obsm["X_umap"])

    cmap = LinearSegmentedColormap.from_list(
        "expression_gradient",
        ["blue", "green", "yellow", "red"],
    )

    plot_data = {}
    all_positive = []

    for sid, gene in genes.items():
        candidates = find_gene_candidates_strict(adata, sid, gene)
        if len(candidates) != 1:
            raise KeyError(
                f"Expected one feature for {sid}:{gene}; found {candidates}."
            )
        feature = candidates[0]
        expr = extract_original_expression(adata, feature)
        mask = species_vector == sid
        values = expr[mask]
        positive = values > threshold
        if positive.any():
            all_positive.append(values[positive])

        plot_data[sid] = {
            "gene": gene,
            "feature": feature,
            "xy": umap[mask],
            "values": values,
            "positive": positive,
        }

    shared_vmax = None
    if shared_scale and all_positive:
        shared_vmax = max(
            float(np.quantile(np.concatenate(all_positive), 0.99)),
            threshold + 1e-8,
        )

    fig, ax = plt.subplots(figsize=figsize)

    for sid in species_ids:
        dat = plot_data[sid]
        xy = dat["xy"]
        values = dat["values"]
        positive = dat["positive"]

        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            c="lightgrey",
            s=max(2, species_sizes[sid] * 0.6),
            linewidths=0,
            rasterized=True,
        )

        if positive.any():
            vmax = shared_vmax
            if vmax is None:
                vmax = max(
                    float(np.quantile(values[positive], 0.99)),
                    threshold + 1e-8,
                )

            ax.scatter(
                xy[positive, 0],
                xy[positive, 1],
                c=values[positive],
                cmap=cmap,
                vmin=threshold,
                vmax=vmax,
                marker=species_markers[sid],
                s=species_sizes[sid],
                alpha=alpha,
                linewidths=0,
                rasterized=True,
            )

    handles = [
        Line2D(
            [0], [0],
            marker=species_markers[sid],
            linestyle="none",
            markerfacecolor="black",
            markeredgecolor="none",
            label=labels[sid],
        )
        for sid in species_ids
    ]
    ax.legend(handles=handles, title="Species", frameon=False)

    ax.set_title(title or "SAMap homolog expression")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_box_aspect(1)

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")

    return fig, ax, {
        sid: plot_data[sid]["feature"]
        for sid in species_ids
    }


def plot_samap_expression_separate_panels(
    adata_or_path,
    genes,
    species_col="species",
    species_labels=None,
    species_sizes=None,
    threshold=0,
    shared_scale=False,
    alpha=0.8,
    ncols=3,
    panel_width=5,
    panel_height=5,
    title=None,
    save_path=None,
    dpi=300,
):
    """Plot one homolog/paralog expression panel per species."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    adata = load_adata(adata_or_path)
    species_ids = list(genes)
    labels = {
        sid: (species_labels or {}).get(sid, DEFAULT_SPECIES_LABELS.get(sid, sid))
        for sid in species_ids
    }
    sizes = {
        sid: (species_sizes or {}).get(sid, 6)
        for sid in species_ids
    }

    species_vector = adata.obs[species_col].astype(str).to_numpy()
    umap = np.asarray(adata.obsm["X_umap"])
    cmap = LinearSegmentedColormap.from_list(
        "expression_gradient",
        ["blue", "green", "yellow", "red"],
    )

    data = {}
    positives = []
    features_used = {}

    for sid, gene in genes.items():
        candidates = find_gene_candidates_strict(adata, sid, gene)
        if len(candidates) != 1:
            raise KeyError(
                f"Expected one feature for {sid}:{gene}; found {candidates}."
            )
        feature = candidates[0]
        features_used[sid] = feature
        expr = extract_original_expression(adata, feature)
        mask = species_vector == sid
        values = expr[mask]
        positive = values > threshold
        if positive.any():
            positives.append(values[positive])
        data[sid] = (umap[mask], values, positive)

    shared_vmax = None
    if shared_scale and positives:
        shared_vmax = max(
            float(np.quantile(np.concatenate(positives), 0.99)),
            threshold + 1e-8,
        )

    nrows = int(np.ceil(len(species_ids) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(panel_width * ncols, panel_height * nrows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    axes = axes.ravel()

    xmin, xmax = umap[:, 0].min(), umap[:, 0].max()
    ymin, ymax = umap[:, 1].min(), umap[:, 1].max()

    for i, sid in enumerate(species_ids):
        ax = axes[i]
        xy, values, positive = data[sid]

        ax.scatter(
            xy[:, 0], xy[:, 1],
            c="lightgrey",
            s=max(2, sizes[sid] * 0.6),
            linewidths=0,
            rasterized=True,
        )

        if positive.any():
            vmax = shared_vmax
            if vmax is None:
                vmax = max(
                    float(np.quantile(values[positive], 0.99)),
                    threshold + 1e-8,
                )

            points = ax.scatter(
                xy[positive, 0],
                xy[positive, 1],
                c=values[positive],
                cmap=cmap,
                vmin=threshold,
                vmax=vmax,
                s=sizes[sid],
                alpha=alpha,
                linewidths=0,
                rasterized=True,
            )
            fig.colorbar(points, ax=ax, fraction=0.046, pad=0.04)

        ax.set_title(f"{labels[sid]}: {genes[sid]}", fontweight="bold")
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_box_aspect(1)
        for spine in ax.spines.values():
            spine.set_visible(False)

    for j in range(len(species_ids), len(axes)):
        axes[j].axis("off")

    if title:
        fig.suptitle(title, fontweight="bold")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")

    return fig, axes, features_used


def plot_human_gene_samap(
    adata_or_path,
    mapping_df,
    human_gene,
    human_col="Human",
    species_gene_columns=None,
    **plot_kwargs,
):
    """
    Plot one human-anchored homolog group using a cross-species mapping table.

    No gene names are embedded in the function; the requested homolog group is
    supplied at runtime.
    """
    if species_gene_columns is None:
        species_gene_columns = {
            "hu": "Human",
            "mo": "Mouse",
            "pi": "Pig",
            "ch": "Chicken",
            "ze": "Zebrafish",
            "ti": "Tilapia",
        }

    rows = mapping_df[mapping_df[human_col] == human_gene]
    if rows.empty:
        raise KeyError(f"{human_gene!r} was not found in mapping_df.")

    row = rows.iloc[0]
    genes = {
        sid: row[col]
        for sid, col in species_gene_columns.items()
        if col in row.index and pd.notna(row[col]) and str(row[col]) != ""
    }

    return plot_samap_expression_separate_panels(
        adata_or_path,
        genes=genes,
        **plot_kwargs,
    )


def plot_combined_gene_samap(
    adata_or_path,
    mapping_df,
    human_gene,
    human_col="Human",
    species_gene_columns=None,
    **plot_kwargs,
):
    """Combined-UMAP counterpart of plot_human_gene_samap()."""
    if species_gene_columns is None:
        species_gene_columns = {
            "hu": "Human",
            "mo": "Mouse",
            "pi": "Pig",
            "ch": "Chicken",
            "ze": "Zebrafish",
            "ti": "Tilapia",
        }

    rows = mapping_df[mapping_df[human_col] == human_gene]
    if rows.empty:
        raise KeyError(f"{human_gene!r} was not found in mapping_df.")

    row = rows.iloc[0]
    genes = {
        sid: row[col]
        for sid, col in species_gene_columns.items()
        if col in row.index and pd.notna(row[col]) and str(row[col]) != ""
    }

    return plot_samap_multispecies_expression(
        adata_or_path,
        genes=genes,
        **plot_kwargs,
    )


# ====================================================================
# Integrated SAMap plotting: species and cell types
# ====================================================================

def plot_samap_species(
    adata_or_path,
    species_col="species",
    species_order=None,
    species_labels=None,
    species_colors=None,
    species_markers=None,
    species_sizes=None,
    alpha=0.65,
    figsize=(9, 8),
    title="SAMap alignment by species",
    save_path=None,
    dpi=300,
):
    """Plot all species on the shared SAMap UMAP."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    adata = load_adata(adata_or_path)
    umap = np.asarray(adata.obsm["X_umap"])
    species = adata.obs[species_col].astype(str).to_numpy()

    observed = list(pd.unique(species))
    order = species_order or observed
    order = [sid for sid in order if sid in observed]

    cmap = plt.get_cmap("tab10")
    colors = {
        sid: (species_colors or {}).get(sid, cmap(i % cmap.N))
        for i, sid in enumerate(order)
    }
    markers_default = ["o", "^", "s", "D", "P", "X", "v"]
    markers = {
        sid: (species_markers or {}).get(
            sid,
            markers_default[i % len(markers_default)],
        )
        for i, sid in enumerate(order)
    }
    sizes = {
        sid: (species_sizes or {}).get(sid, 5)
        for sid in order
    }
    labels = {
        sid: (species_labels or {}).get(sid, DEFAULT_SPECIES_LABELS.get(sid, sid))
        for sid in order
    }

    fig, ax = plt.subplots(figsize=figsize)

    for sid in order:
        mask = species == sid
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            c=[colors[sid]],
            marker=markers[sid],
            s=sizes[sid],
            alpha=alpha,
            linewidths=0,
            rasterized=True,
        )

    handles = [
        Line2D(
            [0], [0],
            marker=markers[sid],
            linestyle="none",
            markerfacecolor=colors[sid],
            markeredgecolor="none",
            label=labels[sid],
        )
        for sid in order
    ]
    ax.legend(handles=handles, title="Species", frameon=False)
    ax.set_title(title, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_box_aspect(1)

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")

    return fig, ax, colors


def plot_samap_species_separate(
    adata_or_path,
    species_col="species",
    species_order=None,
    species_labels=None,
    species_colors=None,
    point_size=5,
    alpha=0.65,
    ncols=3,
    figsize=(15, 10),
    save_path=None,
    dpi=300,
):
    """Plot each species separately using shared integrated UMAP coordinates."""
    import matplotlib.pyplot as plt

    adata = load_adata(adata_or_path)
    umap = np.asarray(adata.obsm["X_umap"])
    species = adata.obs[species_col].astype(str).to_numpy()

    observed = list(pd.unique(species))
    order = species_order or observed
    order = [sid for sid in order if sid in observed]

    cmap = plt.get_cmap("tab10")
    colors = {
        sid: (species_colors or {}).get(sid, cmap(i % cmap.N))
        for i, sid in enumerate(order)
    }
    labels = {
        sid: (species_labels or {}).get(sid, DEFAULT_SPECIES_LABELS.get(sid, sid))
        for sid in order
    }

    nrows = int(np.ceil(len(order) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.ravel()

    xlim = (umap[:, 0].min(), umap[:, 0].max())
    ylim = (umap[:, 1].min(), umap[:, 1].max())

    for i, sid in enumerate(order):
        ax = axes[i]
        mask = species == sid
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            c=[colors[sid]],
            s=point_size,
            alpha=alpha,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(labels[sid], fontweight="bold")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_box_aspect(1)
        for spine in ax.spines.values():
            spine.set_visible(False)

    for j in range(len(order), len(axes)):
        axes[j].axis("off")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")

    return fig, axes


def combine_species_metadata(
    adata,
    annotation_columns,
    species_col="species",
    output_col="cell_type_combined",
    default="Unknown",
):
    """
    Combine species-specific annotation columns into one integrated column.

    annotation_columns example:
        {
            "hu": "hu_cell_type",
            "mo": "mo_cell_type",
            ...
        }
    """
    result = pd.Series(default, index=adata.obs_names, dtype="object")

    for sid, column in annotation_columns.items():
        if column not in adata.obs.columns:
            raise KeyError(f"{column!r} was not found in adata.obs.")

        mask = adata.obs[species_col].astype(str) == sid
        result.loc[mask] = adata.obs.loc[mask, column].astype(str)

    adata.obs[output_col] = result
    return adata


def plot_samap_celltypes_combined(
    adata_or_path,
    species_markers,
    species_col="species",
    celltype_col="cell_type_combined",
    species_labels=None,
    species_sizes=None,
    celltype_order=None,
    palette="tab20",
    alpha=0.65,
    figsize=(11, 8),
    title="SAMap alignment by cell type",
    save_path=None,
    dpi=300,
):
    """Plot cell type as color and species as marker on the same SAMap UMAP."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    adata = load_adata(adata_or_path)
    umap = np.asarray(adata.obsm["X_umap"])

    plot_df = pd.DataFrame({
        "x": umap[:, 0],
        "y": umap[:, 1],
        "species": adata.obs[species_col].astype(str).to_numpy(),
        "cell_type": adata.obs[celltype_col].astype(str).to_numpy(),
    }, index=adata.obs_names)

    invalid = {"", "nan", "None", "NA", "<NA>", "Unknown"}
    plot_df = plot_df[~plot_df["cell_type"].isin(invalid)].copy()

    observed_species = list(pd.unique(plot_df["species"]))
    observed_types = set(plot_df["cell_type"])

    cell_types = (
        [x for x in (celltype_order or []) if x in observed_types]
        + sorted(observed_types.difference(celltype_order or []))
    )
    cmap = plt.get_cmap(palette)
    celltype_colors = {
        ct: cmap(i % cmap.N)
        for i, ct in enumerate(cell_types)
    }

    labels = {
        sid: (species_labels or {}).get(sid, DEFAULT_SPECIES_LABELS.get(sid, sid))
        for sid in observed_species
    }
    sizes = {
        sid: (species_sizes or {}).get(sid, 6)
        for sid in observed_species
    }

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.07, 0.10, 0.64, 0.80])
    legend_ax = fig.add_axes([0.75, 0.10, 0.23, 0.80])
    legend_ax.axis("off")

    for sid in observed_species:
        if sid not in species_markers:
            raise KeyError(f"No marker configured for species {sid!r}.")
        for ct in cell_types:
            mask = (
                (plot_df["species"] == sid)
                & (plot_df["cell_type"] == ct)
            )
            if not mask.any():
                continue
            ax.scatter(
                plot_df.loc[mask, "x"],
                plot_df.loc[mask, "y"],
                c=[celltype_colors[ct]],
                marker=species_markers[sid],
                s=sizes[sid],
                alpha=alpha,
                linewidths=0,
                rasterized=True,
            )

    cell_handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor=celltype_colors[ct],
            markeredgecolor="none",
            label=ct,
        )
        for ct in cell_types
    ]
    ct_legend = legend_ax.legend(
        handles=cell_handles,
        title="Cell type",
        loc="upper left",
        frameon=False,
    )
    legend_ax.add_artist(ct_legend)

    species_handles = [
        Line2D(
            [0], [0],
            marker=species_markers[sid],
            linestyle="none",
            markerfacecolor="black",
            markeredgecolor="black",
            label=labels[sid],
        )
        for sid in observed_species
    ]
    legend_ax.legend(
        handles=species_handles,
        title="Species",
        loc="lower left",
        frameon=False,
    )

    ax.set_title(title, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_box_aspect(1)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")

    return fig, ax, celltype_colors


def plot_samap_celltypes_separate_species(
    adata_or_path,
    species_col="species",
    celltype_col="cell_type_combined",
    species_order=None,
    species_labels=None,
    celltype_order=None,
    celltype_colors=None,
    point_sizes=None,
    palette="tab20",
    alpha=0.65,
    ncols=3,
    figsize=(18, 10),
    save_path=None,
    dpi=300,
):
    """Plot cell types in separate species panels with a shared color mapping."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    adata = load_adata(adata_or_path)
    umap = np.asarray(adata.obsm["X_umap"])

    df = pd.DataFrame({
        "x": umap[:, 0],
        "y": umap[:, 1],
        "species": adata.obs[species_col].astype(str).to_numpy(),
        "cell_type": adata.obs[celltype_col].astype(str).to_numpy(),
    }, index=adata.obs_names)

    invalid = {"", "nan", "None", "NA", "<NA>", "Unknown"}
    df = df[~df["cell_type"].isin(invalid)].copy()

    observed_species = list(pd.unique(df["species"]))
    order = species_order or observed_species
    order = [sid for sid in order if sid in observed_species]

    observed_types = set(df["cell_type"])
    cell_types = (
        [x for x in (celltype_order or []) if x in observed_types]
        + sorted(observed_types.difference(celltype_order or []))
    )

    if celltype_colors is None:
        cmap = plt.get_cmap(palette)
        celltype_colors = {
            ct: cmap(i % cmap.N)
            for i, ct in enumerate(cell_types)
        }

    labels = {
        sid: (species_labels or {}).get(sid, DEFAULT_SPECIES_LABELS.get(sid, sid))
        for sid in order
    }
    sizes = {
        sid: (point_sizes or {}).get(sid, 5)
        for sid in order
    }

    nrows = int(np.ceil(len(order) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.ravel()

    xlim = (df["x"].min(), df["x"].max())
    ylim = (df["y"].min(), df["y"].max())

    for i, sid in enumerate(order):
        ax = axes[i]
        sub = df[df["species"] == sid]

        for ct in cell_types:
            mask = sub["cell_type"] == ct
            if not mask.any():
                continue
            ax.scatter(
                sub.loc[mask, "x"],
                sub.loc[mask, "y"],
                c=[celltype_colors[ct]],
                s=sizes[sid],
                alpha=alpha,
                linewidths=0,
                rasterized=True,
            )

        ax.set_title(labels[sid], fontweight="bold")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_box_aspect(1)
        for spine in ax.spines.values():
            spine.set_visible(False)

    for j in range(len(order), len(axes)):
        axes[j].axis("off")

    handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor=celltype_colors[ct],
            markeredgecolor="none",
            label=ct,
        )
        for ct in cell_types
    ]
    fig.legend(
        handles=handles,
        title="Cell type",
        loc="center left",
        bbox_to_anchor=(0.88, 0.5),
        frameon=False,
    )

    plt.tight_layout(rect=[0, 0, 0.87, 1])

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")

    return fig, axes, celltype_colors



# ====================================================================
# Original-expression overlays on the shared SAMap UMAP
# ====================================================================

def find_gene_candidates_in_source(
    species_adata,
    species_id,
    anchor_gene,
    paralog_species=("ze",),
):
    """
    Find the exact gene name used in a species-specific normalized AnnData.

    Matching is case-insensitive.

    For species in `paralog_species`, the function also allows simple a/b
    paralogs such as:
        genea
        geneb
        genea.1
        geneb.2

    This mirrors the exploratory workflow used to check whether a gene of
    interest is present in each species and what exact feature name is stored
    in the source AnnData object.
    """
    anchor = str(anchor_gene).lower()
    matches = []

    for feature in species_adata.var_names.astype(str):
        feature_lower = feature.lower()

        if feature_lower == anchor:
            matches.append(feature)
            continue

        if species_id in set(paralog_species):
            pattern = rf"^{re.escape(anchor)}[ab](?:\.\d+)?$"
            if re.match(pattern, feature_lower):
                matches.append(feature)

    return list(dict.fromkeys(matches))


def check_gene_presence_across_species(
    species_adatas,
    anchor_gene,
    species_order=None,
    species_labels=None,
    paralog_species=("ze",),
):
    """
    Report whether a gene is present in each source AnnData and the exact
    feature name(s) that will be used for plotting.
    """
    order = species_order or list(species_adatas)
    labels = species_labels or DEFAULT_SPECIES_LABELS

    rows = []

    for sid in order:
        if sid not in species_adatas:
            rows.append({
                "species_id": sid,
                "species": labels.get(sid, sid),
                "query_gene": anchor_gene,
                "present": False,
                "matched_features": "",
                "n_matches": 0,
                "note": "source AnnData not provided",
            })
            continue

        matches = find_gene_candidates_in_source(
            species_adatas[sid],
            species_id=sid,
            anchor_gene=anchor_gene,
            paralog_species=paralog_species,
        )

        rows.append({
            "species_id": sid,
            "species": labels.get(sid, sid),
            "query_gene": anchor_gene,
            "present": len(matches) > 0,
            "matched_features": ", ".join(matches),
            "n_matches": len(matches),
            "note": "",
        })

    return pd.DataFrame(rows)


def check_gene_map_presence(
    species_adatas,
    gene_map,
    species_order=None,
    species_labels=None,
):
    """
    Check an explicit homolog/paralog mapping.

    `gene_map` may contain one gene or multiple genes per species:
        {
            "hu": "GENE",
            "mo": "Gene",
            "ze": ["genea", "geneb"],
        }

    Matching is case-insensitive and returns the exact source feature names.
    """
    order = species_order or list(gene_map)
    labels = species_labels or DEFAULT_SPECIES_LABELS
    rows = []

    for sid in order:
        requested = gene_map.get(sid, [])
        if isinstance(requested, str):
            requested = [requested]
        requested = [str(x) for x in requested if pd.notna(x) and str(x) != ""]

        if sid not in species_adatas:
            rows.append({
                "species_id": sid,
                "species": labels.get(sid, sid),
                "requested_genes": ", ".join(requested),
                "present": False,
                "matched_features": "",
                "missing_features": ", ".join(requested),
            })
            continue

        source_features = species_adatas[sid].var_names.astype(str)
        lookup = {x.lower(): x for x in source_features}

        found = []
        missing = []

        for gene in requested:
            hit = lookup.get(gene.lower())
            if hit is None:
                missing.append(gene)
            else:
                found.append(hit)

        rows.append({
            "species_id": sid,
            "species": labels.get(sid, sid),
            "requested_genes": ", ".join(requested),
            "present": len(found) > 0,
            "matched_features": ", ".join(found),
            "missing_features": ", ".join(missing),
        })

    return pd.DataFrame(rows)


def _source_expression_in_integrated_order(
    species_adata,
    feature,
    integrated_cells,
):
    """
    Pull normalized expression from the original species AnnData while
    explicitly matching the cell order used by the integrated SAMap object.
    """
    integrated_cells = list(map(str, integrated_cells))
    source_cells = set(map(str, species_adata.obs_names))

    missing = [
        cell
        for cell in integrated_cells
        if cell not in source_cells
    ]

    if missing:
        raise ValueError(
            f"{len(missing)} integrated SAMap cells are absent from the "
            f"source AnnData. First few: {missing[:5]}"
        )

    ordered = species_adata[
        integrated_cells,
        feature,
    ]

    values = ordered.X
    if sp.issparse(values):
        values = values.toarray()

    return np.asarray(values).ravel()


def _resolve_anchor_gene_map(
    species_adatas,
    anchor_gene,
    species_order,
    paralog_species=("ze",),
):
    """
    Resolve one anchor gene to the exact feature names found in each source
    AnnData object.
    """
    resolved = {}

    for sid in species_order:
        if sid not in species_adatas:
            continue

        matches = find_gene_candidates_in_source(
            species_adatas[sid],
            species_id=sid,
            anchor_gene=anchor_gene,
            paralog_species=paralog_species,
        )

        if matches:
            resolved[sid] = matches

    return resolved


def _normalize_gene_map(gene_map):
    normalized = {}

    for sid, genes in gene_map.items():
        if isinstance(genes, str):
            genes = [genes]

        genes = [
            str(g)
            for g in genes
            if pd.notna(g) and str(g) != ""
        ]

        if genes:
            normalized[sid] = genes

    return normalized


def plot_gene_across_species_panels(
    integrated_adata,
    species_adatas,
    anchor_gene=None,
    gene_map=None,
    species_col="species",
    species_order=None,
    species_labels=None,
    paralog_species=("ze",),
    threshold=0,
    percentile=0.99,
    ncols=3,
    panel_width=5,
    panel_height=5,
    background_color="lightgray",
    background_size=2,
    expression_size=5,
    alpha=1.0,
    save_path=None,
    dpi=300,
):
    """
    Plot a gene/homolog group as separate species panels on the shared SAMap UMAP.

    Coordinates:
        integrated SAMap embedding

    Expression:
        original species-specific normalized AnnData objects

    Scaling:
        each individual gene panel is scaled independently to its own requested
        expression percentile (default 99th percentile).

    Non-expressing cells are gray and expressing cells are drawn from low to
    high expression so highly expressing cells remain visible.

    Supply either:
        `anchor_gene="GENE"` to auto-detect exact names (plus a/b paralogs for
        configured species), or
        `gene_map={"hu": ["GENE"], "ze": ["genea", "geneb"], ...}`.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    if "X_umap" not in integrated_adata.obsm:
        raise KeyError("'X_umap' was not found in integrated_adata.obsm.")
    if species_col not in integrated_adata.obs:
        raise KeyError(f"{species_col!r} was not found in integrated_adata.obs.")

    order = species_order or DEFAULT_SPECIES_ORDER
    labels = species_labels or DEFAULT_SPECIES_LABELS

    if gene_map is None:
        if anchor_gene is None:
            raise ValueError("Provide either `anchor_gene` or `gene_map`.")
        gene_map = _resolve_anchor_gene_map(
            species_adatas,
            anchor_gene=anchor_gene,
            species_order=order,
            paralog_species=paralog_species,
        )
    else:
        gene_map = _normalize_gene_map(gene_map)

    panels = []
    umap = np.asarray(integrated_adata.obsm["X_umap"])
    species_vector = integrated_adata.obs[species_col].astype(str).to_numpy()

    for sid in order:
        if sid not in species_adatas:
            continue

        genes = gene_map.get(sid, [])
        if not genes:
            continue

        mask = species_vector == sid
        integrated_cells = integrated_adata.obs_names[mask]
        coords = umap[mask]

        source_features = {
            x.lower(): x
            for x in species_adatas[sid].var_names.astype(str)
        }

        for requested_gene in genes:
            feature = source_features.get(str(requested_gene).lower())
            if feature is None:
                continue

            expression = _source_expression_in_integrated_order(
                species_adatas[sid],
                feature,
                integrated_cells,
            )

            panels.append({
                "species_id": sid,
                "feature": feature,
                "coords": coords,
                "expression": expression,
            })

    if not panels:
        raise ValueError("None of the requested genes were found.")

    cmap = LinearSegmentedColormap.from_list(
        "expression_gradient",
        ["blue", "green", "yellow", "red"],
    )

    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(panel_width * ncols, panel_height * nrows),
        squeeze=False,
    )
    axes = axes.ravel()

    xmin, xmax = umap[:, 0].min(), umap[:, 0].max()
    ymin, ymax = umap[:, 1].min(), umap[:, 1].max()

    for i, panel in enumerate(panels):
        ax = axes[i]
        coords = panel["coords"]
        expression = panel["expression"]
        positive = expression > threshold

        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=background_color,
            s=background_size,
            linewidths=0,
            rasterized=True,
        )

        if np.any(positive):
            vmax = float(np.quantile(expression[positive], percentile))
            if vmax <= threshold:
                vmax = float(np.max(expression[positive]))

            positive_indices = np.where(positive)[0]
            plot_order = positive_indices[
                np.argsort(expression[positive_indices])
            ]

            points = ax.scatter(
                coords[plot_order, 0],
                coords[plot_order, 1],
                c=expression[plot_order],
                cmap=cmap,
                vmin=threshold,
                vmax=vmax,
                s=expression_size,
                alpha=alpha,
                linewidths=0,
                rasterized=True,
            )

            cbar = fig.colorbar(points, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label("Normalized expression")

        sid = panel["species_id"]
        ax.set_title(
            f"{labels.get(sid, sid)}: {panel['feature']}",
            fontweight="bold",
        )
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_box_aspect(1)

        for spine in ax.spines.values():
            spine.set_visible(False)

    for j in range(len(panels), len(axes)):
        axes[j].axis("off")

    if anchor_gene:
        fig.suptitle(
            f"{anchor_gene} across species",
            fontweight="bold",
        )

    plt.tight_layout()

    if save_path:
        plt.savefig(
            save_path,
            dpi=dpi,
            bbox_inches="tight",
            facecolor="white",
        )

    return fig, axes, gene_map


def plot_gene_combined_relative_expression(
    integrated_adata,
    species_adatas,
    anchor_gene=None,
    gene_map=None,
    species_col="species",
    species_order=None,
    species_labels=None,
    paralog_species=("ze",),
    threshold=0,
    percentile=0.99,
    background_color="lightgray",
    background_size=2,
    expression_size=5,
    alpha=1.0,
    figsize=(8, 8),
    save_path=None,
    dpi=300,
):
    """
    Plot a homolog/paralog group on one combined SAMap UMAP.

    Coordinates:
        integrated SAMap embedding

    Expression:
        original species-specific normalized AnnData objects

    Multiple requested paralogs within one species are combined by taking the
    maximum expression value per cell.

    Scaling:
        expression is scaled independently within each species so its requested
        percentile (default 99th) equals 1. Values above that percentile are
        clipped to 1. This produces one common 0-1 relative-expression scale
        for the combined cross-species UMAP.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    if "X_umap" not in integrated_adata.obsm:
        raise KeyError("'X_umap' was not found in integrated_adata.obsm.")
    if species_col not in integrated_adata.obs:
        raise KeyError(f"{species_col!r} was not found in integrated_adata.obs.")

    order = species_order or DEFAULT_SPECIES_ORDER
    labels = species_labels or DEFAULT_SPECIES_LABELS

    if gene_map is None:
        if anchor_gene is None:
            raise ValueError("Provide either `anchor_gene` or `gene_map`.")
        gene_map = _resolve_anchor_gene_map(
            species_adatas,
            anchor_gene=anchor_gene,
            species_order=order,
            paralog_species=paralog_species,
        )
    else:
        gene_map = _normalize_gene_map(gene_map)

    umap = np.asarray(integrated_adata.obsm["X_umap"])
    species_vector = integrated_adata.obs[species_col].astype(str).to_numpy()

    combined_expression = np.zeros(
        integrated_adata.n_obs,
        dtype=float,
    )
    combined_positive = np.zeros(
        integrated_adata.n_obs,
        dtype=bool,
    )

    diagnostics = []

    for sid in order:
        if sid not in species_adatas:
            continue

        requested_genes = gene_map.get(sid, [])
        if not requested_genes:
            continue

        source_lookup = {
            x.lower(): x
            for x in species_adatas[sid].var_names.astype(str)
        }
        features = [
            source_lookup[g.lower()]
            for g in requested_genes
            if g.lower() in source_lookup
        ]

        if not features:
            continue

        mask = species_vector == sid
        integrated_cells = integrated_adata.obs_names[mask]

        expression_list = [
            _source_expression_in_integrated_order(
                species_adatas[sid],
                feature,
                integrated_cells,
            )
            for feature in features
        ]

        expression = combine_paralog_expression(expression_list)

        positive = expression > threshold
        scaled = np.zeros_like(expression, dtype=float)
        vmax = np.nan

        if np.any(positive):
            vmax = float(np.quantile(expression[positive], percentile))
            if vmax > 0:
                scaled[positive] = expression[positive] / vmax
                scaled = np.clip(scaled, 0, 1)

        combined_expression[mask] = scaled
        combined_positive[mask] = positive

        diagnostics.append({
            "species_id": sid,
            "species": labels.get(sid, sid),
            "features_used": ", ".join(features),
            "n_cells": int(mask.sum()),
            "n_expressing": int(np.count_nonzero(positive)),
            "percentile": percentile,
            "percentile_value": vmax,
        })

    cmap = LinearSegmentedColormap.from_list(
        "expression_gradient",
        ["blue", "green", "yellow", "red"],
    )

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(
        umap[:, 0],
        umap[:, 1],
        c=background_color,
        s=background_size,
        linewidths=0,
        rasterized=True,
    )

    positive_indices = np.where(combined_positive)[0]

    if positive_indices.size > 0:
        plot_order = positive_indices[
            np.argsort(combined_expression[positive_indices])
        ]

        points = ax.scatter(
            umap[plot_order, 0],
            umap[plot_order, 1],
            c=combined_expression[plot_order],
            cmap=cmap,
            vmin=0,
            vmax=1,
            s=expression_size,
            alpha=alpha,
            linewidths=0,
            rasterized=True,
        )

        cbar = fig.colorbar(points, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Relative expression within species (0–1)")

    title = anchor_gene or "Homolog/paralog expression"
    ax.set_title(
        f"{title} across species",
        fontweight="bold",
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_box_aspect(1)

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(
            save_path,
            dpi=dpi,
            bbox_inches="tight",
            facecolor="white",
        )

    return fig, ax, pd.DataFrame(diagnostics), gene_map




def _match_integrated_cells_to_sam_species(
    integrated_cells,
    sam_adata,
    species_id,
):
    """
    Match cells in the integrated SAMap object back to one species' SAM object.

    SAMap may prefix integrated cell names with a species identifier. This
    helper first tries exact matching, then tries removing '<species>_' from
    the integrated cell IDs.
    """
    integrated_cells = [str(x) for x in integrated_cells]
    sam_cells = pd.Index(sam_adata.obs_names.astype(str))
    sam_cell_set = set(sam_cells)

    if all(cell in sam_cell_set for cell in integrated_cells):
        return integrated_cells

    prefix = f"{species_id}_"
    stripped = [
        cell[len(prefix):] if cell.startswith(prefix) else cell
        for cell in integrated_cells
    ]

    if all(cell in sam_cell_set for cell in stripped):
        return stripped

    missing = [
        original
        for original, stripped_cell in zip(integrated_cells, stripped)
        if original not in sam_cell_set and stripped_cell not in sam_cell_set
    ]

    raise ValueError(
        f"Could not match {len(missing)} integrated cells to SAM species "
        f"{species_id!r}. First few: {missing[:5]}"
    )


def check_gene_map_in_sam_species(
    sams,
    gene_map,
    species_order=None,
    species_labels=None,
):
    """
    Verify an explicit homolog/paralog mapping directly in the species-specific
    SAM objects used by SAMap.
    """
    order = species_order or list(gene_map)
    labels = species_labels or DEFAULT_SPECIES_LABELS
    rows = []

    for sid in order:
        requested = gene_map.get(sid, [])
        if isinstance(requested, str):
            requested = [requested]
        requested = [str(x) for x in requested if pd.notna(x) and str(x) != ""]

        if sid not in sams:
            rows.append({
                "species_id": sid,
                "species": labels.get(sid, sid),
                "requested_genes": ", ".join(requested),
                "present": False,
                "matched_features": "",
                "missing_features": ", ".join(requested),
            })
            continue

        sam_adata = sams[sid].adata
        lookup = {
            str(feature).lower(): str(feature)
            for feature in sam_adata.var_names
        }

        found = []
        missing = []

        for gene in requested:
            exact = lookup.get(gene.lower())

            if exact is None:
                prefixed = lookup.get(f"{sid}_{gene}".lower())
                exact = prefixed

            if exact is None:
                missing.append(gene)
            else:
                found.append(exact)

        rows.append({
            "species_id": sid,
            "species": labels.get(sid, sid),
            "requested_genes": ", ".join(requested),
            "present": len(found) > 0,
            "matched_features": ", ".join(found),
            "missing_features": ", ".join(missing),
        })

    return pd.DataFrame(rows)


def plot_gene_across_sam_species_panels(
    integrated_adata,
    sams,
    gene_map,
    species_col="species",
    species_order=None,
    species_labels=None,
    threshold=0,
    percentile=0.99,
    ncols=3,
    panel_width=5,
    panel_height=5,
    background_color="lightgray",
    background_size=2,
    expression_size=5,
    alpha=1.0,
    save_path=None,
    dpi=300,
):
    """
    Plot homolog/paralog expression in separate panels using expression from
    the species-specific SAM objects used by SAMap.

    Coordinates come from the integrated SAMap UMAP.

    Expression comes from:
        sams[species_id].adata

    Each individual gene/paralog panel is scaled independently to its own
    requested expression percentile (default 99th percentile).
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    if "X_umap" not in integrated_adata.obsm:
        raise KeyError("'X_umap' was not found in integrated_adata.obsm.")
    if species_col not in integrated_adata.obs:
        raise KeyError(f"{species_col!r} was not found in integrated_adata.obs.")

    order = species_order or DEFAULT_SPECIES_ORDER
    labels = species_labels or DEFAULT_SPECIES_LABELS
    gene_map = _normalize_gene_map(gene_map)

    umap = np.asarray(integrated_adata.obsm["X_umap"])
    species_vector = integrated_adata.obs[species_col].astype(str).to_numpy()

    panels = []

    for sid in order:
        if sid not in sams:
            continue

        requested_genes = gene_map.get(sid, [])
        if not requested_genes:
            continue

        sam_adata = sams[sid].adata
        feature_lookup = {
            str(feature).lower(): str(feature)
            for feature in sam_adata.var_names
        }

        mask = species_vector == sid
        integrated_cells = integrated_adata.obs_names[mask]
        coords = umap[mask]

        sam_cells = _match_integrated_cells_to_sam_species(
            integrated_cells,
            sam_adata,
            sid,
        )

        for requested_gene in requested_genes:
            feature = feature_lookup.get(requested_gene.lower())

            if feature is None:
                feature = feature_lookup.get(
                    f"{sid}_{requested_gene}".lower()
                )

            if feature is None:
                continue

            values = sam_adata[
                sam_cells,
                feature,
            ].X

            if sp.issparse(values):
                values = values.toarray()

            expression = np.asarray(values).ravel()

            panels.append({
                "species_id": sid,
                "feature": feature,
                "coords": coords,
                "expression": expression,
            })

    if not panels:
        raise ValueError("None of the requested genes were found in the SAM objects.")

    cmap = LinearSegmentedColormap.from_list(
        "expression_gradient",
        ["blue", "green", "yellow", "red"],
    )

    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(panel_width * ncols, panel_height * nrows),
        squeeze=False,
    )
    axes = axes.ravel()

    xmin, xmax = umap[:, 0].min(), umap[:, 0].max()
    ymin, ymax = umap[:, 1].min(), umap[:, 1].max()

    for i, panel in enumerate(panels):
        ax = axes[i]
        coords = panel["coords"]
        expression = panel["expression"]
        positive = expression > threshold

        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=background_color,
            s=background_size,
            linewidths=0,
            rasterized=True,
        )

        if np.any(positive):
            vmax = float(np.quantile(expression[positive], percentile))
            if vmax <= threshold:
                vmax = float(np.max(expression[positive]))

            positive_indices = np.where(positive)[0]
            plot_order = positive_indices[
                np.argsort(expression[positive_indices])
            ]

            points = ax.scatter(
                coords[plot_order, 0],
                coords[plot_order, 1],
                c=expression[plot_order],
                cmap=cmap,
                vmin=threshold,
                vmax=vmax,
                s=expression_size,
                alpha=alpha,
                linewidths=0,
                rasterized=True,
            )

            cbar = fig.colorbar(points, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label("SAM expression")

        sid = panel["species_id"]
        display_feature = panel["feature"]
        if display_feature.startswith(f"{sid}_"):
            display_feature = display_feature[len(sid) + 1:]

        ax.set_title(
            f"{labels.get(sid, sid)}: {display_feature}",
            fontweight="bold",
        )
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_box_aspect(1)

        for spine in ax.spines.values():
            spine.set_visible(False)

    for j in range(len(panels), len(axes)):
        axes[j].axis("off")

    plt.tight_layout()

    if save_path:
        plt.savefig(
            save_path,
            dpi=dpi,
            bbox_inches="tight",
            facecolor="white",
        )

    return fig, axes

# ====================================================================
# Fibroblast-only post-processing: global Leiden and Seurat clusters
# ====================================================================

def run_global_leiden(
    adata,
    resolution=0.7,
    adjacency_key="connectivities",
    output_col="global_leiden",
):
    """Run a single Leiden clustering on the integrated SAMap graph."""
    import scanpy as sc

    if adjacency_key not in adata.obsp:
        raise KeyError(
            f"{adjacency_key!r} was not found in adata.obsp."
        )

    sc.tl.leiden(
        adata,
        adjacency=adata.obsp[adjacency_key],
        resolution=resolution,
        key_added=output_col,
    )
    return adata


def plot_global_leiden_by_species(
    adata,
    species_col="species",
    cluster_col="global_leiden",
    species_order=None,
    species_labels=None,
    point_size=6,
    alpha=0.75,
    ncols=3,
    figsize=(18, 10),
    save_path=None,
    dpi=300,
):
    """Plot global integrated Leiden clusters separately for each species."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    order = species_order or DEFAULT_SPECIES_ORDER
    labels = species_labels or DEFAULT_SPECIES_LABELS
    umap = np.asarray(adata.obsm["X_umap"])

    cluster_ids = sorted(
        adata.obs[cluster_col].dropna().astype(str).unique(),
        key=lambda x: (0, float(x)) if _is_number(x) else (1, x),
    )

    cmap = plt.get_cmap("tab20")
    colors = {
        cluster: cmap(i % cmap.N)
        for i, cluster in enumerate(cluster_ids)
    }

    nrows = int(np.ceil(len(order) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.ravel()

    xlim = (umap[:, 0].min(), umap[:, 0].max())
    ylim = (umap[:, 1].min(), umap[:, 1].max())

    for i, sid in enumerate(order):
        ax = axes[i]
        mask = adata.obs[species_col].astype(str).to_numpy() == sid
        sub_obs = adata.obs.loc[mask].copy()
        sub_umap = umap[mask]

        for cluster in cluster_ids:
            cmask = sub_obs[cluster_col].astype(str).to_numpy() == cluster
            if not cmask.any():
                continue

            ax.scatter(
                sub_umap[cmask, 0],
                sub_umap[cmask, 1],
                c=[colors[cluster]],
                s=point_size,
                alpha=alpha,
                linewidths=0,
                rasterized=True,
            )
            ax.text(
                np.median(sub_umap[cmask, 0]),
                np.median(sub_umap[cmask, 1]),
                cluster,
                ha="center",
                va="center",
                fontweight="bold",
            )

        ax.set_title(labels.get(sid, sid), fontweight="bold")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_box_aspect(1)
        for spine in ax.spines.values():
            spine.set_visible(False)

    for j in range(len(order), len(axes)):
        axes[j].axis("off")

    handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor=colors[cluster],
            markeredgecolor="none",
            label=f"Cluster {cluster}",
        )
        for cluster in cluster_ids
    ]
    fig.legend(
        handles=handles,
        title="Global Leiden",
        loc="center left",
        bbox_to_anchor=(0.88, 0.5),
        frameon=False,
    )
    plt.tight_layout(rect=[0, 0, 0.87, 1])

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")

    return fig, axes, colors


def add_seurat_clusters_from_metadata(
    adata,
    metadata_files,
    species_col="species",
    id_col="cell_id",
    source_cluster_col="seurat_clusters",
    output_template="{sid}_seurat_cluster",
):
    """
    Add original Seurat cluster assignments to an integrated SAMap AnnData.
    """
    for sid, metadata_file in metadata_files.items():
        metadata = pd.read_csv(metadata_file).set_index(id_col)
        mask = adata.obs[species_col].astype(str) == sid
        cells = adata.obs_names[mask]
        out_col = output_template.format(sid=sid)

        adata.obs[out_col] = pd.NA
        adata.obs.loc[mask, out_col] = (
            metadata[source_cluster_col]
            .reindex(cells)
            .to_numpy()
        )

    return adata


def plot_seurat_clusters_by_species(
    adata,
    species_col="species",
    species_order=None,
    species_labels=None,
    cluster_template="{sid}_seurat_cluster",
    point_size=6,
    alpha=0.75,
    ncols=3,
    figsize=(18, 10),
    save_path=None,
    dpi=300,
):
    """Overlay each species' original Seurat clusters on the SAMap UMAP."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    order = species_order or DEFAULT_SPECIES_ORDER
    labels = species_labels or DEFAULT_SPECIES_LABELS
    umap = np.asarray(adata.obsm["X_umap"])

    nrows = int(np.ceil(len(order) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.ravel()

    xlim = (umap[:, 0].min(), umap[:, 0].max())
    ylim = (umap[:, 1].min(), umap[:, 1].max())

    for i, sid in enumerate(order):
        ax = axes[i]
        cluster_col = cluster_template.format(sid=sid)

        if cluster_col not in adata.obs:
            ax.axis("off")
            continue

        mask = adata.obs[species_col].astype(str).to_numpy() == sid
        obs = adata.obs.loc[mask].copy()
        coords = umap[mask]

        valid = obs[cluster_col].notna().to_numpy()
        obs = obs.loc[obs[cluster_col].notna()].copy()
        coords = coords[valid]

        cluster_ids = sorted(
            obs[cluster_col].astype(str).unique(),
            key=lambda x: (0, float(x)) if _is_number(x) else (1, x),
        )
        cmap = plt.get_cmap("tab20")
        colors = {
            cluster: cmap(j % cmap.N)
            for j, cluster in enumerate(cluster_ids)
        }

        for cluster in cluster_ids:
            cmask = obs[cluster_col].astype(str).to_numpy() == cluster
            ax.scatter(
                coords[cmask, 0],
                coords[cmask, 1],
                c=[colors[cluster]],
                s=point_size,
                alpha=alpha,
                linewidths=0,
                rasterized=True,
            )
            ax.text(
                np.median(coords[cmask, 0]),
                np.median(coords[cmask, 1]),
                cluster,
                ha="center",
                va="center",
                fontweight="bold",
            )

        handles = [
            Line2D(
                [0], [0],
                marker="o",
                linestyle="none",
                markerfacecolor=colors[cluster],
                markeredgecolor="none",
                label=f"Cluster {cluster}",
            )
            for cluster in cluster_ids
        ]
        ax.legend(
            handles=handles,
            title="Seurat cluster",
            frameon=False,
            fontsize=7,
        )

        ax.set_title(labels.get(sid, sid), fontweight="bold")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_box_aspect(1)
        for spine in ax.spines.values():
            spine.set_visible(False)

    for j in range(len(order), len(axes)):
        axes[j].axis("off")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")

    return fig, axes


def species_composition_by_cluster(
    adata,
    cluster_col="global_leiden",
    species_col="species",
):
    """Return within-cluster species composition and abundance enrichment."""
    overall = adata.obs[species_col].value_counts(normalize=True)

    composition = pd.crosstab(
        adata.obs[cluster_col],
        adata.obs[species_col],
        normalize="index",
    )

    enrichment = composition.div(overall, axis=1)
    counts = pd.crosstab(
        adata.obs[cluster_col],
        adata.obs[species_col],
    )

    return counts, composition, enrichment


def rank_global_leiden_markers_by_species(
    adata,
    species_ids=None,
    species_col="species",
    cluster_col="global_leiden",
    method="wilcoxon",
    key_added="global_leiden_markers",
):
    """Rank markers for integrated clusters independently within each species."""
    import scanpy as sc

    species_ids = species_ids or DEFAULT_SPECIES_ORDER
    marker_tables = {}

    for sid in species_ids:
        sub = adata[
            adata.obs[species_col].astype(str) == sid
        ].copy()

        if sub.n_obs == 0:
            continue

        sc.tl.rank_genes_groups(
            sub,
            groupby=cluster_col,
            method=method,
            key_added=key_added,
        )

        table = sc.get.rank_genes_groups_df(
            sub,
            group=None,
            key=key_added,
        )
        table["species"] = sid
        marker_tables[sid] = table

    return marker_tables


def summarize_homolog_markers(
    marker_tables,
    cross_species_table,
    species_columns=None,
    p_adj_max=0.05,
    logfc_min=0.5,
):
    """
    Convert species-specific global-Leiden markers into human-anchored groups.
    """
    species_columns = species_columns or {
        "hu": "Human",
        "mo": "Mouse",
        "pi": "Pig",
        "ch": "Chicken",
        "ze": "Zebrafish",
        "ti": "Tilapia",
    }

    lookup = {}
    for gene in cross_species_table["Human"].dropna():
        lookup[("hu", str(gene))] = str(gene)

    for sid, column in species_columns.items():
        if sid == "hu" or column not in cross_species_table:
            continue
        temp = cross_species_table[["Human", column]].dropna()
        for _, row in temp.iterrows():
            lookup[(sid, str(row[column]))] = str(row["Human"])

    rows = []

    for sid, df in marker_tables.items():
        filtered = df[
            (df["pvals_adj"] < p_adj_max)
            & (df["logfoldchanges"] > logfc_min)
        ]

        prefix = f"{sid}_"

        for _, row in filtered.iterrows():
            gene = str(row["names"])
            clean_gene = gene.removeprefix(prefix)
            human_anchor = lookup.get((sid, clean_gene))

            if human_anchor is None:
                continue

            rows.append({
                "global_leiden": str(row["group"]),
                "species": sid,
                "species_name": species_columns.get(sid, sid),
                "gene": clean_gene,
                "Human_anchor": human_anchor,
                "logFC": row["logfoldchanges"],
                "score": row["scores"],
                "pval_adj": row["pvals_adj"],
            })

    homolog_markers = pd.DataFrame(rows)

    if homolog_markers.empty:
        return homolog_markers, pd.DataFrame()

    summary = (
        homolog_markers
        .groupby(["global_leiden", "Human_anchor"])
        .agg(
            n_species=("species", "nunique"),
            species=("species_name", lambda x: ", ".join(sorted(set(x)))),
            genes=("gene", lambda x: ", ".join(sorted(set(x)))),
        )
        .reset_index()
    )

    logfc = (
        homolog_markers
        .pivot_table(
            index=["global_leiden", "Human_anchor"],
            columns="species_name",
            values="logFC",
            aggfunc="first",
        )
        .reset_index()
    )

    summary = summary.merge(
        logfc,
        on=["global_leiden", "Human_anchor"],
        how="left",
    )

    return homolog_markers, summary


def _is_number(value):
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


# ====================================================================
# Mapping scores and GenePairFinder
# ====================================================================

def get_celltype_mapping_scores(
    sm,
    keys,
    n_top=0,
):
    """Wrapper around SAMap get_mapping_scores()."""
    from samap.analysis import get_mapping_scores
    return get_mapping_scores(sm, keys, n_top=n_top)



BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def gene_pair_table_with_blast(
    gpf,
    species1,
    species2,
    celltype1,
    celltype2,
    maps_dir,
):
    """
    Run GenePairFinder for one mapped cell-type pair and append BLAST evidence.
    """
    maps_dir = Path(maps_dir)

    result = gpf.find_genes(celltype1, celltype2)
    gene_pairs = result[0]
    pval1 = result[3]
    pval2 = result[4]
    scores = gpf.gene_pair_scores.reindex(gene_pairs)

    table = pd.DataFrame({
        "gene_pair": gene_pairs,
        "score": scores.values,
        "pval1": pval1,
        "pval2": pval2,
    })

    table[["gene1", "gene2"]] = table["gene_pair"].str.split(
        ";",
        expand=True,
    )

    table["blast_gene1"] = table["gene1"].str.replace(
        rf"^{re.escape(species1)}_",
        "",
        regex=True,
    )
    table["blast_gene2"] = table["gene2"].str.replace(
        rf"^{re.escape(species2)}_",
        "",
        regex=True,
    )

    folder = f"{species1}{species2}"
    blast_file = maps_dir / folder / f"{species1}_to_{species2}.txt"

    if not blast_file.exists():
        raise FileNotFoundError(f"BLAST map not found: {blast_file}")

    blast = pd.read_csv(
        blast_file,
        sep="\t",
        names=BLAST_COLUMNS,
        header=None,
    )

    blast = blast[
        ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
    ].copy()

    blast.columns = [
        "blast_gene1",
        "blast_gene2",
        "blast_pident",
        "blast_length",
        "blast_evalue",
        "blast_bitscore",
    ]

    table = table.merge(
        blast,
        on=["blast_gene1", "blast_gene2"],
        how="left",
    )

    table["species1"] = species1
    table["species2"] = species2
    table["celltype1"] = celltype1
    table["celltype2"] = celltype2

    return table[
        [
            "species1",
            "species2",
            "celltype1",
            "celltype2",
            "gene1",
            "gene2",
            "gene_pair",
            "score",
            "pval1",
            "pval2",
            "blast_pident",
            "blast_length",
            "blast_evalue",
            "blast_bitscore",
        ]
    ]


def run_all_pairwise_gene_pairs(
    gpf,
    celltype_labels,
    maps_dir,
    species_order=None,
    output_dir=None,
):
    """
    Create pairwise GenePairFinder + BLAST CSVs for every species pair.

    celltype_labels example:
        {
            "hu": "hu_Fibroblasts",
            "mo": "mo_Fibroblasts",
            ...
        }
    """
    order = [
        sid
        for sid in (species_order or DEFAULT_SPECIES_ORDER)
        if sid in celltype_labels
    ]

    output_dir = Path(output_dir) if output_dir else None
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)

    tables = []

    for species1, species2 in combinations(order, 2):
        table = gene_pair_table_with_blast(
            gpf,
            species1=species1,
            species2=species2,
            celltype1=celltype_labels[species1],
            celltype2=celltype_labels[species2],
            maps_dir=maps_dir,
        )
        tables.append(table)

        if output_dir:
            table.to_csv(
                output_dir
                / f"{species1}_{species2}_gene_pairs_with_BLAST.csv",
                index=False,
            )

    combined = pd.concat(tables, ignore_index=True)

    if output_dir:
        combined.to_csv(
            output_dir / "all_pairwise_gene_pairs_with_BLAST.csv",
            index=False,
        )

    return combined


def build_human_anchored_gene_table(
    pairwise_files,
    species_prefixes=None,
):
    """
    Merge human-vs-species pairwise tables into one human-anchored table.

    For each human gene and non-human species, the candidate with the highest
    BLAST bit score is retained.
    """
    species_prefixes = species_prefixes or {
        "Mouse": "mo_",
        "Pig": "pi_",
        "Chicken": "ch_",
        "Zebrafish": "ze_",
        "Tilapia": "ti_",
    }

    tables = []
    human_pvalues = []

    for species_name, file in pairwise_files.items():
        df = pd.read_csv(file)
        df["Human"] = df["gene1"].str.replace(r"^hu_", "", regex=True)
        df["Species_gene"] = df["gene2"].str.replace(
            rf"^{re.escape(species_prefixes[species_name])}",
            "",
            regex=True,
        )

        best = (
            df.sort_values(
                "blast_bitscore",
                ascending=False,
                na_position="last",
            )
            .drop_duplicates("Human", keep="first")
            .copy()
        )

        human_pvalues.append(
            best[["Human", "pval1"]].rename(
                columns={"pval1": "Human_pval"}
            )
        )

        tables.append(
            best[
                ["Human", "Species_gene", "score", "blast_bitscore"]
            ].rename(columns={
                "Species_gene": species_name,
                "score": f"{species_name}_score",
                "blast_bitscore": f"{species_name}_bitscore",
            })
        )

    merged = reduce(
        lambda left, right: pd.merge(
            left,
            right,
            on="Human",
            how="outer",
        ),
        tables,
    )

    human_pvalues = (
        pd.concat(human_pvalues, ignore_index=True)
        .sort_values("Human_pval")
        .drop_duplicates("Human", keep="first")
    )

    merged = human_pvalues.merge(
        merged,
        on="Human",
        how="right",
    )

    return merged.sort_values(
        "Human_pval",
        ascending=True,
        kind="stable",
    )


# ====================================================================
# Fibroblast gene-network analysis
# ====================================================================

def build_gene_pair_network(
    all_pairs,
    filter_query=None,
):
    """
    Convert pairwise GenePairFinder results into a NetworkX graph.

    Every species-specific gene is a node and every GenePairFinder-supported
    pair is an edge carrying expression- and BLAST-based evidence.
    """
    import networkx as nx

    pairs = all_pairs.copy()

    if filter_query:
        pairs = pairs.query(filter_query).copy()

    if "species1" not in pairs:
        pairs["species1"] = (
            pairs["gene1"].astype(str).str.split("_", n=1).str[0]
        )
    if "species2" not in pairs:
        pairs["species2"] = (
            pairs["gene2"].astype(str).str.split("_", n=1).str[0]
        )

    graph = nx.Graph()

    for _, row in pairs.iterrows():
        gene1 = row["gene1"]
        gene2 = row["gene2"]

        graph.add_node(gene1, species=row["species1"])
        graph.add_node(gene2, species=row["species2"])

        edge_attrs = {
            "GenePairFinder_score": row.get("score", np.nan),
            "pval1": row.get("pval1", np.nan),
            "pval2": row.get("pval2", np.nan),
            "blast_pident": row.get("blast_pident", np.nan),
            "blast_length": row.get("blast_length", np.nan),
            "blast_evalue": row.get("blast_evalue", np.nan),
            "blast_bitscore": row.get("blast_bitscore", np.nan),
        }

        graph.add_edge(gene1, gene2, **edge_attrs)

    return graph


def summarize_connected_components(
    graph,
    species_order=None,
):
    """Summarize each connected cross-species gene network."""
    import networkx as nx

    species_order = species_order or DEFAULT_SPECIES_ORDER
    rows = []

    for group_id, component in enumerate(
        nx.connected_components(graph),
        start=1,
    ):
        genes = sorted(component)
        species_present = sorted({
            graph.nodes[gene]["species"]
            for gene in genes
        })
        subgraph = graph.subgraph(component)

        row = {
            "group_id": group_id,
            "n_species": len(species_present),
            "species_present": ", ".join(species_present),
            "n_genes": len(genes),
            "n_edges": subgraph.number_of_edges(),
        }

        for sid in species_order:
            sid_genes = [
                gene
                for gene in genes
                if graph.nodes[gene]["species"] == sid
            ]
            row[f"has_{sid}"] = sid in species_present
            row[f"{sid}_genes"] = ", ".join(sorted(sid_genes))

        rows.append(row)

    component_df = pd.DataFrame(rows)

    if component_df.empty:
        return component_df

    return component_df.sort_values(
        ["n_species", "n_genes"],
        ascending=[False, False],
    ).reset_index(drop=True)


def add_species_combo(
    component_df,
    species_order=None,
):
    """Add exact species-combination labels such as hu+mo+ze."""
    species_order = species_order or DEFAULT_SPECIES_ORDER
    out = component_df.copy()

    out["species_combo"] = out.apply(
        lambda row: "+".join(
            sid
            for sid in species_order
            if bool(row.get(f"has_{sid}", False))
        ),
        axis=1,
    )

    return out


def save_components_by_species_combo(
    component_df,
    output_dir,
):
    """Save one connected-component table per exact species combination."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for combo in sorted(component_df["species_combo"].dropna().unique()):
        subset = component_df[
            component_df["species_combo"] == combo
        ].copy()

        subset.to_csv(
            output_dir / f"{combo.replace('+', '_')}.csv",
            index=False,
        )


def plot_species_upset(
    component_df,
    species_labels=None,
    species_order=None,
    title="Species combinations of gene networks",
    save_path=None,
    dpi=300,
):
    """Create an UpSet plot of species membership across connected components."""
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_indicators

    species_order = species_order or DEFAULT_SPECIES_ORDER
    species_labels = species_labels or DEFAULT_SPECIES_LABELS

    presence = pd.DataFrame({
        species_labels.get(sid, sid): component_df[f"has_{sid}"].astype(bool)
        for sid in species_order
        if f"has_{sid}" in component_df
    })

    upset_data = from_indicators(
        presence.columns,
        presence,
    )

    upset = UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        sort_categories_by=None,
    )
    axes = upset.plot()
    plt.suptitle(title, fontweight="bold")

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")

    return axes


def plot_gene_network_component(
    graph,
    component_df,
    group_id,
    species_colors=None,
    species_labels=None,
    figsize=(12, 10),
    save_path=None,
    dpi=300,
):
    """Plot one connected cross-species gene-network component."""
    import matplotlib.pyplot as plt
    import networkx as nx
    from matplotlib.lines import Line2D

    species_labels = species_labels or DEFAULT_SPECIES_LABELS

    if species_colors is None:
        cmap = plt.get_cmap("tab10")
        species_colors = {
            sid: cmap(i % cmap.N)
            for i, sid in enumerate(DEFAULT_SPECIES_ORDER)
        }

    row = component_df.loc[
        component_df["group_id"] == group_id
    ]
    if row.empty:
        raise KeyError(f"group_id {group_id} was not found.")

    genes = set()
    for sid in DEFAULT_SPECIES_ORDER:
        col = f"{sid}_genes"
        if col in row:
            value = row.iloc[0][col]
            if pd.notna(value) and str(value).strip():
                genes.update(
                    x.strip()
                    for x in str(value).split(",")
                    if x.strip()
                )

    subgraph = graph.subgraph(genes).copy()
    pos = nx.spring_layout(subgraph, seed=42, k=1.5)

    node_colors = [
        species_colors.get(
            graph.nodes[node].get("species"),
            "grey",
        )
        for node in subgraph.nodes
    ]

    edge_scores = [
        subgraph[u][v].get("GenePairFinder_score", np.nan)
        for u, v in subgraph.edges
    ]
    valid_scores = [
        x for x in edge_scores
        if pd.notna(x)
    ]

    if valid_scores and max(valid_scores) > min(valid_scores):
        lo, hi = min(valid_scores), max(valid_scores)
        edge_widths = [
            1 + 5 * (score - lo) / (hi - lo)
            if pd.notna(score)
            else 1
            for score in edge_scores
        ]
    else:
        edge_widths = [2] * len(edge_scores)

    fig, ax = plt.subplots(figsize=figsize)

    nx.draw_networkx_nodes(
        subgraph,
        pos,
        node_color=node_colors,
        node_size=900,
        edgecolors="black",
        linewidths=0.8,
        ax=ax,
    )

    nx.draw_networkx_edges(
        subgraph,
        pos,
        width=edge_widths,
        alpha=0.6,
        ax=ax,
    )

    labels = {
        node: node.split("_", 1)[1]
        if "_" in node
        else node
        for node in subgraph.nodes
    }

    nx.draw_networkx_labels(
        subgraph,
        pos,
        labels=labels,
        font_size=9,
        ax=ax,
    )

    present_species = sorted({
        graph.nodes[node].get("species")
        for node in subgraph.nodes
    })

    handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor=species_colors.get(sid, "grey"),
            markeredgecolor="black",
            label=species_labels.get(sid, sid),
        )
        for sid in present_species
    ]

    ax.legend(
        handles=handles,
        title="Species",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )
    ax.set_title(
        f"Cross-species gene network – component {group_id}",
        fontweight="bold",
    )
    ax.axis("off")

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")

    return fig, ax


# ====================================================================
def get_gene_pairs_for_species_set(
    gene_pairs,
    species_ids,
):
    """Restrict a two-column gene-pair array to a requested species set."""
    gene_pairs = np.asarray(gene_pairs)

    if gene_pairs.ndim != 2 or gene_pairs.shape[1] != 2:
        raise ValueError("gene_pairs must be a two-column array.")

    species_ids = set(species_ids)

    first = np.array([
        gene.split("_", 1)[0]
        for gene in gene_pairs[:, 0].astype(str)
    ])
    second = np.array([
        gene.split("_", 1)[0]
        for gene in gene_pairs[:, 1].astype(str)
    ])

    keep = (
        np.isin(first, list(species_ids))
        & np.isin(second, list(species_ids))
        & (first != second)
    )

    result = gene_pairs[keep]
    return np.unique(np.sort(result, axis=1), axis=0)


def count_edges_by_species_pair(gene_pairs):
    """Count pairwise edges by species combination."""
    counts = {}

    for gene1, gene2 in gene_pairs:
        pair = tuple(sorted([
            str(gene1).split("_", 1)[0],
            str(gene2).split("_", 1)[0],
        ]))
        counts[pair] = counts.get(pair, 0) + 1

    return counts


def get_species_gene_pairs(
    gene_pairs,
    species1,
    species2,
):
    """Restrict a two-column gene-pair array to one species pair."""
    gene_pairs = np.asarray(gene_pairs)

    if gene_pairs.ndim != 2 or gene_pairs.shape[1] != 2:
        raise ValueError("gene_pairs must be a two-column array.")

    first = np.array([
        gene.split("_", 1)[0]
        for gene in gene_pairs[:, 0].astype(str)
    ])
    second = np.array([
        gene.split("_", 1)[0]
        for gene in gene_pairs[:, 1].astype(str)
    ])

    keep = (
        ((first == species1) & (second == species2))
        | ((first == species2) & (second == species1))
    )

    result = gene_pairs[keep]
    return np.unique(np.sort(result, axis=1), axis=0)


def run_paralog_substitutions(
    sm,
    ortholog_pairs,
    broad_homolog_pairs,
    species_ids,
    output_dir,
    psub_thr=0.3,
    species_labels=None,
):
    """
    Run SAMap ParalogSubstitutions for all requested pairwise comparisons.

    `ortholog_pairs` should contain the narrower ortholog set and
    `broad_homolog_pairs` the broader homolog/paralog candidate set.
    """
    from samap.analysis import ParalogSubstitutions

    species_labels = species_labels or {
        sid: sid
        for sid in species_ids
    }
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    for species1, species2 in combinations(species_ids, 2):
        direct = get_species_gene_pairs(
            ortholog_pairs,
            species1,
            species2,
        )
        broad = get_species_gene_pairs(
            broad_homolog_pairs,
            species1,
            species2,
        )

        if len(direct) == 0 or len(broad) == 0:
            continue

        result = ParalogSubstitutions(
            sm,
            direct,
            paralog_pairs=broad,
            psub_thr=psub_thr,
        )

        comparison = (
            f"{species_labels.get(species1, species1)}_"
            f"{species_labels.get(species2, species2)}"
        )

        result.to_csv(
            output_dir / f"{comparison}_paralog_substitutions.csv",
            index=False,
        )

        results[comparison] = result

        del direct, broad
        gc.collect()

    if not results:
        return pd.DataFrame()

    combined = []
    for comparison, table in results.items():
        table = table.copy()
        table.insert(0, "species_comparison", comparison)
        combined.append(table)

    combined = pd.concat(combined, ignore_index=True)

    if "corr diff" in combined.columns:
        combined = combined.sort_values(
            "corr diff",
            ascending=False,
        ).reset_index(drop=True)

    combined.to_csv(
        output_dir / "all_paralog_substitutions.csv",
        index=False,
    )

    return combined
