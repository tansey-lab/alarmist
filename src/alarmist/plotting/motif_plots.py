"""
BPTF Motif Visualization Functions

Functions for visualizing BPTF analysis results including:
- LRI clustering and communication heatmaps
- Top LRI interactions per motif
- LRI network graphs
- Pathway-aggregated interactions

Adapted from scripts/bptf_visualization_utils.py for notebook-friendly usage.
"""

import io
import json
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch

from alarmist.constants import (
    COLUMN_NAME_CELL_PAIR,
    COLUMN_NAME_CELLTYPE1,
    COLUMN_NAME_CELLTYPE2,
    COLUMN_NAME_FACTOR,
    COLUMN_NAME_IN_TOP20,
    COLUMN_NAME_LIGAND,
    COLUMN_NAME_LRI_NAME,
    COLUMN_NAME_MODE,
    COLUMN_NAME_MOTIF_IDX,
    COLUMN_NAME_PATHWAY,
    COLUMN_NAME_RECEPTOR,
    COLUMN_NAME_SIGNALING_TYPE,
    COLUMN_NAME_SOURCE,
    COLUMN_NAME_TARGET,
    COLUMN_NAME_WEIGHT,
)

logger = logging.getLogger(__name__)

# Optional dependencies
try:
    from graphviz import Digraph
    from PIL import Image

    GRAPHVIZ_AVAILABLE = True
except ImportError:
    GRAPHVIZ_AVAILABLE = False


# ==============================================================================
# Utility Functions
# ==============================================================================


def parse_lri_full(lri_name: str) -> tuple[str, str, str, str, str]:
    """Parse LRI name into components

    Parameters
    ----------
    lri_name : str
        LRI name in format: celltype1|celltype2|ligand|receptor|mode

    Returns
    -------
    tuple
        (celltype1, celltype2, ligand, receptor, signaling_type)
    """
    parts = lri_name.split("|")
    if len(parts) >= 5:
        return parts[0], parts[1], parts[2], parts[3], parts[4]
    elif len(parts) == 4:
        c1, c2, ligand, receptor = parts
        mode = "autocrine" if c1 == c2 else "paracrine"
        return c1, c2, ligand, receptor, mode
    elif len(parts) == 2:
        return "unknown", "unknown", parts[0], parts[1], "unknown"
    else:
        return "unknown", "unknown", lri_name, lri_name, "unknown"


def get_cell_type_colors(
    unique_ct: list[str], palette: str = "tab20"
) -> dict[str, tuple]:
    """Generate color map for cell types

    Parameters
    ----------
    unique_ct : list of str
        List of unique cell type names
    palette : str, default "tab20"
        Matplotlib colormap name

    Returns
    -------
    dict
        Mapping from cell type to RGBA tuple

    Note
    ----
    Consider using al.set_celltype_colors() instead for global color consistency.
    """
    n = len(unique_ct)
    ct_cmap = plt.get_cmap(palette, n)
    return {ct: ct_cmap(i) for i, ct in enumerate(unique_ct)}


# ==============================================================================
# Data Preprocessing Functions
# ==============================================================================


def add_lri_components(
    lri_motifs_df: pd.DataFrame, lri_col: str = COLUMN_NAME_LRI_NAME
) -> pd.DataFrame:
    """Add parsed LRI components to DataFrame

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with LRI names
    lri_col : str, default 'lri_name'
        Column containing LRI names

    Returns
    -------
    pd.DataFrame
        DataFrame with added columns: celltype1, celltype2, ligand, receptor, signaling_type
    """
    df = lri_motifs_df.copy()

    # Filter out GENE entries
    df = df[~df[lri_col].str.startswith("GENE")]

    # Parse LRI names
    df[
        [
            COLUMN_NAME_CELLTYPE1,
            COLUMN_NAME_CELLTYPE2,
            COLUMN_NAME_LIGAND,
            COLUMN_NAME_RECEPTOR,
            COLUMN_NAME_SIGNALING_TYPE,
        ]
    ] = df[lri_col].apply(lambda x: pd.Series(parse_lri_full(x)))

    return df


def annotate_pathways(
    lri_motifs_df: pd.DataFrame,
    cellchatdb_df: pd.DataFrame,
    ligand_col: str = COLUMN_NAME_LIGAND,
    receptor_col: str = COLUMN_NAME_RECEPTOR,
) -> pd.DataFrame:
    """Annotate LRI interactions with pathway information from CellChatDB

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with ligand and receptor columns
    cellchatdb_df : pd.DataFrame
        CellChatDB DataFrame with 'ligand', 'receptor', 'pathway' columns
    ligand_col : str, default 'ligand'
        Column name for ligand
    receptor_col : str, default 'receptor'
        Column name for receptor

    Returns
    -------
    pd.DataFrame
        DataFrame with added 'pathway' column
    """
    df = lri_motifs_df.copy()

    # Build ligand-receptor to pathway mapping
    lr_to_pathway = {}
    for _, row in cellchatdb_df.iterrows():
        ligands = (
            row[COLUMN_NAME_LIGAND].split("_")
            if pd.notna(row[COLUMN_NAME_LIGAND])
            else []
        )
        receptors = (
            row[COLUMN_NAME_RECEPTOR].split("_")
            if pd.notna(row[COLUMN_NAME_RECEPTOR])
            else []
        )
        pathway = row[COLUMN_NAME_PATHWAY]

        for lig in ligands:
            for rec in receptors:
                lr_to_pathway[(lig, rec)] = pathway
        lr_to_pathway[(row[COLUMN_NAME_LIGAND], row[COLUMN_NAME_RECEPTOR])] = pathway

    # Map pathways
    def get_pathway(row):
        pathway = lr_to_pathway.get((row[ligand_col], row[receptor_col]))
        if pathway:
            return pathway

        # Fuzzy matching
        for (lig, rec), path in lr_to_pathway.items():
            if (row[ligand_col] in lig or lig in row[ligand_col]) and (
                row[receptor_col] in rec or rec in row[receptor_col]
            ):
                return path
        return "Unknown"

    df[COLUMN_NAME_PATHWAY] = df.apply(get_pathway, axis=1)

    return df


# normalize_by_column_mean() has been removed
# Normalization is now done automatically in save_bptf_results()
# The saved lri_motifs.csv already contains 'factor_norm' and 'mean' columns


# ==============================================================================
# Visualization Functions
# ==============================================================================


def plot_lri_clustermap(
    lri_motifs_df: pd.DataFrame,
    factor_col: str = COLUMN_NAME_FACTOR,
    exclude_celltypes: list[str] | None = None,
    save_path: str | None = None,
    figsize: tuple[int, int] = (20, 40),
    cmap: str = "Blues",
    log_scale: bool = True,
) -> object:
    """Create clustermap of LRI cell type pairs across motifs

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with columns: motif_idx, celltype1, celltype2, and factor_col.
        Typically output from process_bptf_results() or load_bptf_results().
    factor_col : str, default "factor"
        Column name for factor values (e.g., 'factor', 'factor_lrnorm', 'score')
    exclude_celltypes : list of str, optional
        Cell types to exclude (case-insensitive substring match).
        Example: ['granulocyte'] to exclude granulocyte-containing pairs.
    save_path : str, optional
        Path to save figure. If None, displays plot.
    figsize : tuple, default (20, 40)
        Figure size
    cmap : str, default "Blues"
        Colormap for heatmap
    log_scale : bool, default True
        Whether to use log scale for color normalization

    Returns
    -------
    seaborn.ClusterGrid
        Clustermap object

    Examples
    --------
    >>> g = plot_lri_clustermap(lri_motifs, factor_col='factor')
    >>> g = plot_lri_clustermap(lri_motifs, factor_col='score', exclude_celltypes=['granulocyte'])
    """
    required_cols = {
        COLUMN_NAME_MOTIF_IDX,
        COLUMN_NAME_CELLTYPE1,
        COLUMN_NAME_CELLTYPE2,
        factor_col,
    }
    missing = required_cols - set(lri_motifs_df.columns)
    if missing:
        raise ValueError(f"Missing columns in lri_motifs_df: {sorted(missing)}")

    df = lri_motifs_df.copy()

    # Filter out excluded cell types
    if exclude_celltypes:
        for ct in exclude_celltypes:
            mask = df[COLUMN_NAME_CELLTYPE1].str.contains(
                ct, case=False, na=False
            ) | df[COLUMN_NAME_CELLTYPE2].str.contains(ct, case=False, na=False)
            df = df[~mask]

    # Create cell_pair column
    df[COLUMN_NAME_CELL_PAIR] = (
        df[COLUMN_NAME_CELLTYPE1] + "|" + df[COLUMN_NAME_CELLTYPE2]
    )

    # Aggregate by cell_pair and motif_idx
    agg = (
        df.groupby([COLUMN_NAME_CELL_PAIR, COLUMN_NAME_MOTIF_IDX])[factor_col]
        .sum()
        .reset_index()
    )

    # Pivot to matrix format (cell_pair × motif)
    pivot = agg.pivot_table(
        index=COLUMN_NAME_CELL_PAIR,
        columns=COLUMN_NAME_MOTIF_IDX,
        values=factor_col,
        fill_value=0,
        aggfunc="sum",
    )

    if pivot.empty:
        raise ValueError("No data to plot after filtering")

    # Create clustermap
    if log_scale:
        # Get min positive value for log scale
        min_pos = pivot[pivot > 0].min().min()
        max_val = pivot.max().max()
        if min_pos > 0 and max_val > min_pos:
            norm = LogNorm(vmin=min_pos, vmax=max_val)
        else:
            norm = None
            log_scale = False
    else:
        norm = None

    g = sns.clustermap(
        pivot,
        cmap=cmap,
        norm=norm,
        linewidths=0.5,
        cbar_kws={"label": f"{factor_col} ({'log scale' if log_scale else 'linear'})"},
        figsize=figsize,
        cbar_pos=(1.02, 0.2, 0.03, 0.6),
        xticklabels=True,
        yticklabels=True,
    )

    g.ax_heatmap.set_xlabel("Motif")
    g.ax_heatmap.set_ylabel("Cell Type Pair (Sender|Receiver)")
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    exclude_text = (
        f" (excluding: {', '.join(exclude_celltypes)})" if exclude_celltypes else ""
    )
    plt.suptitle(
        f"Clustermap of Cell Type Pairs Across Motifs\n{factor_col}{exclude_text}",
        y=0.98,
    )

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, bbox_inches="tight", dpi=300)
        logger.debug(f"Saved: {save_path}")

    return g


def plot_celltype_communication_by_motif(
    lri_motifs_df: pd.DataFrame,
    factor_col: str = COLUMN_NAME_FACTOR,
    save_path: str | None = None,
    n_cols: int = 5,
    figsize_per_motif: tuple[float, float] = (5, 4),
    cmap: str = "Blues",
    exclude_celltypes: list[str] | None = None,
) -> plt.Figure:
    """Create cell-cell communication heatmaps for each motif.

    For each motif, aggregates factor values by (celltype1, celltype2) pairs
    and displays as a heatmap showing sender -> receiver communication strength.

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with columns: motif_idx, celltype1, celltype2, and factor_col.
        Typically output from add_lri_components() or load_bptf_results().
    factor_col : str, default "factor"
        Column name for factor values (e.g., 'factor', 'factor_norm', 'score')
    save_path : str, optional
        Path to save figure. If None, displays plot.
    n_cols : int, default 5
        Number of columns in subplot grid
    figsize_per_motif : tuple, default (5, 4)
        Size per motif subplot
    cmap : str, default "Blues"
        Colormap for heatmaps
    exclude_celltypes : list of str, optional
        Cell types to exclude. If None, no exclusion.
        Pass list like ['granulocyte'] to exclude specific types.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object

    Examples
    --------
    >>> fig = plot_celltype_communication_by_motif(
    ...     lri_motifs,
    ...     factor_col='factor_norm',
    ...     n_cols=5
    ... )

    >>> # Use different colormap and exclude specific cell types
    >>> fig = plot_celltype_communication_by_motif(
    ...     lri_motifs,
    ...     factor_col='score',
    ...     cmap='Reds',
    ...     exclude_celltypes=['granulocyte', 'unknown']
    ... )
    """
    required_cols = {
        COLUMN_NAME_MOTIF_IDX,
        COLUMN_NAME_CELLTYPE1,
        COLUMN_NAME_CELLTYPE2,
        factor_col,
    }
    missing = required_cols - set(lri_motifs_df.columns)
    if missing:
        raise ValueError(f"Missing columns in lri_motifs_df: {sorted(missing)}")

    # Filter out excluded cell types
    df = lri_motifs_df.copy()
    if exclude_celltypes:
        for ct in exclude_celltypes:
            mask = df[COLUMN_NAME_CELLTYPE1].str.contains(
                ct, case=False, na=False
            ) | df[COLUMN_NAME_CELLTYPE2].str.contains(ct, case=False, na=False)
            df = df[~mask]

    # Get all motifs
    motifs = sorted(df[COLUMN_NAME_MOTIF_IDX].unique())
    n_motifs = len(motifs)

    if n_motifs == 0:
        raise ValueError("No motifs found in data after filtering")

    # Setup layout
    n_rows = (n_motifs + n_cols - 1) // n_cols

    # Create figure
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(figsize_per_motif[0] * n_cols, figsize_per_motif[1] * n_rows),
    )

    # Handle single row/col cases
    if n_motifs == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    # Plot each motif
    for idx, motif in enumerate(motifs):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        # Filter to this motif
        dfp = df[df[COLUMN_NAME_MOTIF_IDX] == motif]

        # Aggregate by cell type pairs
        agg = (
            dfp.groupby([COLUMN_NAME_CELLTYPE1, COLUMN_NAME_CELLTYPE2])[factor_col]
            .sum()
            .reset_index()
        )

        if agg.empty:
            ax.text(
                0.5, 0.5, "No Data", ha="center", va="center", transform=ax.transAxes
            )
            ax.set_title(f"Motif {motif}", fontsize=10, fontweight="bold")
            continue

        # Pivot to heatmap format
        heatmap_data = agg.pivot_table(
            index=COLUMN_NAME_CELLTYPE1,
            columns=COLUMN_NAME_CELLTYPE2,
            values=factor_col,
            fill_value=0,
            aggfunc="sum",
        )

        if heatmap_data.empty:
            ax.text(
                0.5, 0.5, "No Data", ha="center", va="center", transform=ax.transAxes
            )
        else:
            sns.heatmap(
                heatmap_data,
                cmap=cmap,
                cbar_kws={"label": factor_col},
                linewidths=0.3,
                square=True,
                xticklabels=True,
                yticklabels=True,
                ax=ax,
            )

        ax.set_title(f"Motif {motif}", fontsize=10, fontweight="bold")
        ax.set_xlabel("Receiver", fontsize=8)
        ax.set_ylabel("Sender", fontsize=8)
        ax.tick_params(axis="x", rotation=45, labelsize=5)
        ax.tick_params(axis="y", rotation=0, labelsize=6)

    # Hide extra subplots
    for idx in range(n_motifs, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)

    # Title
    exclude_text = (
        f" (excluding: {', '.join(exclude_celltypes)})" if exclude_celltypes else ""
    )
    plt.suptitle(
        f"Cell-Cell Communication by Motif\n{factor_col}{exclude_text}",
        fontsize=14,
        fontweight="bold",
        y=1.01,
    )
    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.debug(f"Saved: {save_path}")

    return fig


def plot_top_lri_interactions_dot(
    lri_motifs_df: pd.DataFrame,
    factor_col: str = COLUMN_NAME_FACTOR,
    top_n: int = 35,
    save_path: str | None = None,
    n_cols: int = 3,
    figsize_per_motif: tuple[float, float] = (12, 9),
    ct_colors: dict[str, str] | None = None,
) -> plt.Figure:
    """
    Plot top LRI interactions per motif with lollipop lines + colored dots at endpoints.

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with columns: motif_idx, celltype1, celltype2, ligand, receptor,
        signaling_type, and factor_col
    factor_col : str, default "factor"
        Column name for factor values
    top_n : int, default 35
        Number of top interactions to show per motif
    save_path : str, optional
        Path to save figure
    n_cols : int, default 3
        Number of columns in subplot grid
    figsize_per_motif : tuple, default (12, 9)
        Size per motif subplot
    ct_colors : dict, optional
        Cell type color mapping. If None, uses global colors from set_celltype_colors()
        or auto-generates from data.

    Signaling types by marker:
      ■ autocrine, ● paracrine, ▲ juxtacrine
    """
    from alarmist.plotting.colors import _get_colors_for_plotting

    required_cols = {
        COLUMN_NAME_MOTIF_IDX,
        COLUMN_NAME_CELLTYPE1,
        COLUMN_NAME_CELLTYPE2,
        COLUMN_NAME_LIGAND,
        COLUMN_NAME_RECEPTOR,
        COLUMN_NAME_SIGNALING_TYPE,
        factor_col,
    }
    missing = required_cols - set(lri_motifs_df.columns)
    if missing:
        raise ValueError(f"Missing columns in lri_motifs_df: {sorted(missing)}")

    # Get cell type colors
    all_celltypes = list(
        set(lri_motifs_df[COLUMN_NAME_CELLTYPE1])
        | set(lri_motifs_df[COLUMN_NAME_CELLTYPE2])
    )
    ct_color_map = _get_colors_for_plotting(ct_colors, all_celltypes)

    motifs = sorted(lri_motifs_df[COLUMN_NAME_MOTIF_IDX].unique())
    n_prog = len(motifs)

    n_rows = int(np.ceil(n_prog / n_cols))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * figsize_per_motif[0], n_rows * figsize_per_motif[1]),
        constrained_layout=False,
    )
    axes = np.array(axes).reshape(-1)

    # Give some room for text on both sides
    fig.subplots_adjust(
        right=0.85, left=0.10, top=0.94, bottom=0.06, wspace=0.9, hspace=0.4
    )

    # --- label spacing knobs (in "points") ---
    # More negative => ligand further left; more positive => receptor further right
    ligand_dx_pts = -20
    receptor_dx_pts = 12

    for ax, prog in zip(axes, motifs):
        dfp = lri_motifs_df[lri_motifs_df[COLUMN_NAME_MOTIF_IDX] == prog].copy()
        if dfp.empty:
            ax.set_visible(False)
            continue

        # Filter: keep top per (signaling_type, receptor, direction)
        dfp_sorted = dfp.sort_values(factor_col, ascending=False)
        dfp_filtered = dfp_sorted.drop_duplicates(
            subset=[
                COLUMN_NAME_SIGNALING_TYPE,
                COLUMN_NAME_RECEPTOR,
                COLUMN_NAME_CELLTYPE1,
                COLUMN_NAME_CELLTYPE2,
            ],
            keep="first",
        )

        top_df = dfp_filtered.nlargest(top_n, factor_col).reset_index(drop=True)
        y = np.arange(len(top_df))

        # No y tick labels (avoid huge block)
        ax.set_yticks(y)
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)

        for yi, row in top_df.iterrows():
            total = float(row[factor_col])

            sig_type = str(row[COLUMN_NAME_SIGNALING_TYPE]).lower()
            if sig_type in ["auto", "autocrine"]:
                marker_shape = "s"
            elif sig_type in ["juxta", "juxtacrine"]:
                marker_shape = "^"
            else:
                marker_shape = "o"

            col1 = ct_color_map.get(row[COLUMN_NAME_CELLTYPE1], "gray")
            col2 = ct_color_map.get(row[COLUMN_NAME_CELLTYPE2], "gray")

            # lollipop line
            ax.plot([0, total], [yi, yi], color="gray", linewidth=1.5, zorder=1)

            marker_size = 80
            ax.scatter(
                0,
                yi,
                color=col1,
                s=marker_size,
                marker=marker_shape,
                zorder=2,
                edgecolors="black",
                linewidth=0.5,
            )
            ax.scatter(
                total,
                yi,
                color=col2,
                s=marker_size,
                marker=marker_shape,
                zorder=2,
                edgecolors="black",
                linewidth=0.5,
            )

            # Ligand: further LEFT of the left endpoint (x=0), so it sits left of the axis line
            ax.annotate(
                str(row[COLUMN_NAME_LIGAND]),
                xy=(0, yi),
                xytext=(ligand_dx_pts, 0),
                textcoords="offset points",
                ha="right",
                va="center",
                fontsize=9,
                clip_on=False,
            )

            # Receptor: further RIGHT of the right endpoint
            ax.annotate(
                str(row[COLUMN_NAME_RECEPTOR]),
                xy=(total, yi),
                xytext=(receptor_dx_pts, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=9,
                clip_on=False,
            )

        ax.invert_yaxis()
        ax.set_xlabel(f"{factor_col}", fontsize=10)
        ax.set_title(f"Motif {prog}", fontsize=12, fontweight="bold")

        # More padding so text doesn't feel cramped
        if len(top_df) > 0:
            max_val = float(top_df[factor_col].max())
            ax.set_xlim(-max_val * 0.04, max_val * 1.12)

        # Remove top/right spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Hide unused subplots
    for ax in axes[n_prog:]:
        ax.set_visible(False)

    # Legend
    legend_handles = []
    for ct, col in ct_color_map.items():
        legend_handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=col,
                markersize=8,
                markeredgecolor="black",
                markeredgewidth=0.5,
                label=ct,
            )
        )

    legend_handles.append(plt.Line2D([0], [0], color="none", label=""))
    legend_handles.append(plt.Line2D([0], [0], color="none", label="Signaling Types:"))
    legend_handles.append(
        plt.Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor="lightgray",
            markersize=8,
            markeredgecolor="black",
            markeredgewidth=0.5,
            label="■ = Autocrine",
        )
    )
    legend_handles.append(
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="lightgray",
            markersize=8,
            markeredgecolor="black",
            markeredgewidth=0.5,
            label="● = Paracrine",
        )
    )
    legend_handles.append(
        plt.Line2D(
            [0],
            [0],
            marker="^",
            color="w",
            markerfacecolor="lightgray",
            markersize=8,
            markeredgecolor="black",
            markeredgewidth=0.5,
            label="▲ = Juxtacrine",
        )
    )

    fig.legend(
        handles=legend_handles,
        title="Cell Types & Signaling",
        loc="center right",
        bbox_to_anchor=(0.92, 0.5),
        fontsize=8,
        title_fontsize=10,
        frameon=True,
        fancybox=True,
        shadow=False,
    )

    fig.suptitle(
        "Top LR Interactions per Motif\n"
        "(Left text = Ligand, Right text = Receptor; Shape = Signaling type)",
        fontsize=14,
        y=1,
        fontweight="bold",
    )

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, bbox_inches="tight", dpi=300)
        logger.debug(f"Saved: {save_path}")

    return fig


def plot_single_motif_lri_lollipop(
    lri_motifs_df: pd.DataFrame,
    motif_idx: int,
    factor_col: str = COLUMN_NAME_FACTOR,
    top_n: int = 25,
    sender_type: str | None = None,
    receiver_type: str | None = None,
    save_path: str | None = None,
    figsize: tuple[float, float] = (10, 8),
    title: str | None = None,
    show_legend: bool = True,
    ax: plt.Axes | None = None,
    ct_colors: dict[str, str] | None = None,
) -> plt.Figure:
    """
    Plot top LRI interactions for a SINGLE motif as a lollipop chart.

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with columns: motif_idx, celltype1, celltype2, ligand, receptor,
        signaling_type, and factor_col
    motif_idx : int
        The motif index to plot
    factor_col : str, default "factor"
        Column name for the factor values (e.g., 'factor', 'score', 'factor_lrnorm')
    top_n : int, default 25
        Number of top LRI interactions to show
    sender_type : str, optional
        Filter by sender cell type (celltype1). If None, no filtering.
    receiver_type : str, optional
        Filter by receiver cell type (celltype2). If None, no filtering.
    save_path : str, optional
        Path to save figure. If None, displays plot.
    figsize : tuple, default (10, 8)
        Figure size (width, height)
    title : str, optional
        Custom title. If None, auto-generated based on motif and filters.
    show_legend : bool, default True
        Whether to show the legend
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    ct_colors : dict, optional
        Cell type color mapping. If None, uses global colors from set_celltype_colors()
        or auto-generates from data.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object (or None if ax was provided)

    Examples
    --------
    >>> # All interactions for motif 5
    >>> fig = plot_single_motif_lri_lollipop(lri_motifs, motif_idx=5, factor_col='score', top_n=30)

    >>> # Only interactions where Tumor cells are the sender
    >>> fig = plot_single_motif_lri_lollipop(lri_motifs, motif_idx=5, sender_type='Tumor', top_n=20)

    >>> # Only interactions from Macrophage to T cell
    >>> fig = plot_single_motif_lri_lollipop(
    ...     lri_motifs, motif_idx=5,
    ...     sender_type='Macrophage', receiver_type='T cell', top_n=15
    ... )
    """
    from alarmist.plotting.colors import _get_colors_for_plotting

    required_cols = {
        COLUMN_NAME_MOTIF_IDX,
        COLUMN_NAME_CELLTYPE1,
        COLUMN_NAME_CELLTYPE2,
        COLUMN_NAME_LIGAND,
        COLUMN_NAME_RECEPTOR,
        COLUMN_NAME_SIGNALING_TYPE,
        factor_col,
    }
    missing = required_cols - set(lri_motifs_df.columns)
    if missing:
        raise ValueError(f"Missing columns in lri_motifs_df: {sorted(missing)}")

    # Filter to single motif
    dfp = lri_motifs_df[lri_motifs_df[COLUMN_NAME_MOTIF_IDX] == motif_idx].copy()
    if dfp.empty:
        raise ValueError(f"No data found for motif_idx={motif_idx}")

    # Filter by sender_type (celltype1) if specified
    if sender_type is not None:
        dfp = dfp[dfp[COLUMN_NAME_CELLTYPE1] == sender_type]
        if dfp.empty:
            raise ValueError(
                f"No data found for motif_idx={motif_idx} with sender_type='{sender_type}'"
            )

    # Filter by receiver_type (celltype2) if specified
    if receiver_type is not None:
        dfp = dfp[dfp[COLUMN_NAME_CELLTYPE2] == receiver_type]
        if dfp.empty:
            filter_msg = f"receiver_type='{receiver_type}'"
            if sender_type is not None:
                filter_msg = f"sender_type='{sender_type}' and {filter_msg}"
            raise ValueError(
                f"No data found for motif_idx={motif_idx} with {filter_msg}"
            )

    # Get cell type colors
    all_celltypes = list(
        set(lri_motifs_df[COLUMN_NAME_CELLTYPE1])
        | set(lri_motifs_df[COLUMN_NAME_CELLTYPE2])
    )
    ct_color_map = _get_colors_for_plotting(ct_colors, all_celltypes)

    # Create figure if no ax provided
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        created_fig = True
    else:
        fig = ax.get_figure()

    # Filter: keep top per (signaling_type, receptor, direction)
    dfp_sorted = dfp.sort_values(factor_col, ascending=False)
    dfp_filtered = dfp_sorted.drop_duplicates(
        subset=[
            COLUMN_NAME_SIGNALING_TYPE,
            COLUMN_NAME_RECEPTOR,
            COLUMN_NAME_CELLTYPE1,
            COLUMN_NAME_CELLTYPE2,
        ],
        keep="first",
    )

    top_df = dfp_filtered.nlargest(top_n, factor_col).reset_index(drop=True)
    y = np.arange(len(top_df))

    # Label spacing (in points)
    ligand_dx_pts = -15
    receptor_dx_pts = 10

    # No y tick labels
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)

    for yi, row in top_df.iterrows():
        total = float(row[factor_col])

        # Marker shape by signaling type
        sig_type = str(row[COLUMN_NAME_SIGNALING_TYPE]).lower()
        if sig_type in ["auto", "autocrine"]:
            marker_shape = "s"
        elif sig_type in ["juxta", "juxtacrine"]:
            marker_shape = "^"
        else:
            marker_shape = "o"

        col1 = ct_color_map.get(row[COLUMN_NAME_CELLTYPE1], "gray")
        col2 = ct_color_map.get(row[COLUMN_NAME_CELLTYPE2], "gray")

        # Lollipop line
        ax.plot([0, total], [yi, yi], color="gray", linewidth=1.5, zorder=1)

        # Dots at endpoints
        marker_size = 100
        ax.scatter(
            0,
            yi,
            color=col1,
            s=marker_size,
            marker=marker_shape,
            zorder=2,
            edgecolors="black",
            linewidth=0.5,
        )
        ax.scatter(
            total,
            yi,
            color=col2,
            s=marker_size,
            marker=marker_shape,
            zorder=2,
            edgecolors="black",
            linewidth=0.5,
        )

        # Ligand label (left of x=0)
        ax.annotate(
            str(row[COLUMN_NAME_LIGAND]),
            xy=(0, yi),
            xytext=(ligand_dx_pts, 0),
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize=10,
            clip_on=False,
        )

        # Receptor label (right of endpoint)
        ax.annotate(
            str(row[COLUMN_NAME_RECEPTOR]),
            xy=(total, yi),
            xytext=(receptor_dx_pts, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=10,
            clip_on=False,
        )

    ax.invert_yaxis()
    ax.set_xlabel(factor_col, fontsize=11)

    # Title
    if title is None:
        title = f"Motif {motif_idx} - Top {len(top_df)} LRI Interactions"
        # Add filter info to title
        filter_parts = []
        if sender_type is not None:
            filter_parts.append(f"Sender: {sender_type}")
        if receiver_type is not None:
            filter_parts.append(f"Receiver: {receiver_type}")
        if filter_parts:
            title += f"\n({', '.join(filter_parts)})"
    ax.set_title(title, fontsize=13, fontweight="bold")

    # X-axis limits with padding
    if len(top_df) > 0:
        max_val = float(top_df[factor_col].max())
        ax.set_xlim(-max_val * 0.04, max_val * 1.12)

    # Remove top/right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend
    if show_legend and created_fig:
        legend_handles = []

        # Cell type colors
        for ct, col in ct_color_map.items():
            legend_handles.append(
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=col,
                    markersize=8,
                    markeredgecolor="black",
                    markeredgewidth=0.5,
                    label=ct,
                )
            )

        # Separator and signaling types
        legend_handles.append(plt.Line2D([0], [0], color="none", label=""))
        legend_handles.append(
            plt.Line2D([0], [0], color="none", label="Signaling Types:")
        )
        legend_handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="s",
                color="w",
                markerfacecolor="lightgray",
                markersize=8,
                markeredgecolor="black",
                markeredgewidth=0.5,
                label="■ = Autocrine",
            )
        )
        legend_handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="lightgray",
                markersize=8,
                markeredgecolor="black",
                markeredgewidth=0.5,
                label="● = Paracrine",
            )
        )
        legend_handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="^",
                color="w",
                markerfacecolor="lightgray",
                markersize=8,
                markeredgecolor="black",
                markeredgewidth=0.5,
                label="▲ = Juxtacrine",
            )
        )

        ax.legend(
            handles=legend_handles,
            title="Cell Types & Signaling",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            fontsize=8,
            title_fontsize=9,
            frameon=True,
            fancybox=True,
        )

    if created_fig:
        plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, bbox_inches="tight", dpi=300)
        logger.debug(f"Saved: {save_path}")

    return fig if created_fig else None


def plot_single_motif_cellpair_lollipop(
    lri_motifs_df: pd.DataFrame,
    motif_idx: int,
    factor_col: str = COLUMN_NAME_FACTOR,
    top_n: int = 20,
    sender_type: str | None = None,
    receiver_type: str | None = None,
    save_path: str | None = None,
    figsize: tuple[float, float] = (10, 8),
    title: str | None = None,
    show_legend: bool = True,
    ax: plt.Axes | None = None,
    ct_colors: dict[str, str] | None = None,
) -> plt.Figure:
    """
    Plot aggregated LRI interactions by cell type pairs for a SINGLE motif.

    This function aggregates (sums) the factor values by (celltype1, celltype2) pairs,
    allowing visualization of which cell type pairs have stronger overall interaction.

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with columns: motif_idx, celltype1, celltype2, and factor_col
    motif_idx : int
        The motif index to plot
    factor_col : str, default "factor"
        Column name for the factor values to aggregate (e.g., 'factor', 'factor_lrnorm')
    top_n : int, default 20
        Number of top cell type pairs to show
    sender_type : str, optional
        Filter by sender cell type (celltype1). If None, no filtering.
    receiver_type : str, optional
        Filter by receiver cell type (celltype2). If None, no filtering.
    save_path : str, optional
        Path to save figure. If None, displays plot.
    figsize : tuple, default (10, 8)
        Figure size (width, height)
    title : str, optional
        Custom title. If None, auto-generated based on motif and filters.
    show_legend : bool, default True
        Whether to show the legend
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    ct_colors : dict, optional
        Cell type color mapping. If None, uses global colors from set_celltype_colors()
        or auto-generates from data.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object (or None if ax was provided)

    Examples
    --------
    >>> # All cell type pairs for motif 5
    >>> fig = plot_single_motif_cellpair_lollipop(lri_motifs, motif_idx=5, factor_col='factor', top_n=15)

    >>> # Only pairs where Tumor cells are the sender
    >>> fig = plot_single_motif_cellpair_lollipop(
    ...     lri_motifs, motif_idx=5, sender_type='Tumor', factor_col='factor_lrnorm'
    ... )
    """
    from alarmist.plotting.colors import _get_colors_for_plotting

    required_cols = {
        COLUMN_NAME_MOTIF_IDX,
        COLUMN_NAME_CELLTYPE1,
        COLUMN_NAME_CELLTYPE2,
        factor_col,
    }
    missing = required_cols - set(lri_motifs_df.columns)
    if missing:
        raise ValueError(f"Missing columns in lri_motifs_df: {sorted(missing)}")

    # Filter to single motif
    dfp = lri_motifs_df[lri_motifs_df[COLUMN_NAME_MOTIF_IDX] == motif_idx].copy()
    if dfp.empty:
        raise ValueError(f"No data found for motif_idx={motif_idx}")

    # Filter by sender_type (celltype1) if specified
    if sender_type is not None:
        dfp = dfp[dfp[COLUMN_NAME_CELLTYPE1] == sender_type]
        if dfp.empty:
            raise ValueError(
                f"No data found for motif_idx={motif_idx} with sender_type='{sender_type}'"
            )

    # Filter by receiver_type (celltype2) if specified
    if receiver_type is not None:
        dfp = dfp[dfp[COLUMN_NAME_CELLTYPE2] == receiver_type]
        if dfp.empty:
            filter_msg = f"receiver_type='{receiver_type}'"
            if sender_type is not None:
                filter_msg = f"sender_type='{sender_type}' and {filter_msg}"
            raise ValueError(
                f"No data found for motif_idx={motif_idx} with {filter_msg}"
            )

    # Aggregate by cell type pairs
    agg_df = (
        dfp.groupby([COLUMN_NAME_CELLTYPE1, COLUMN_NAME_CELLTYPE2])
        .agg({factor_col: "sum"})
        .reset_index()
    )

    # Sort and get top N
    agg_df = agg_df.sort_values(factor_col, ascending=False)
    top_df = agg_df.head(top_n).reset_index(drop=True)

    if top_df.empty:
        raise ValueError(f"No aggregated data for motif_idx={motif_idx}")

    # Get cell type colors
    all_celltypes = list(
        set(lri_motifs_df[COLUMN_NAME_CELLTYPE1])
        | set(lri_motifs_df[COLUMN_NAME_CELLTYPE2])
    )
    ct_color_map = _get_colors_for_plotting(ct_colors, all_celltypes)

    # Create figure if no ax provided
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        created_fig = True
    else:
        fig = ax.get_figure()

    y = np.arange(len(top_df))

    # No y tick labels (we'll annotate instead)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)

    # Label spacing (in points)
    sender_dx_pts = -15
    receiver_dx_pts = 10

    for yi, row in top_df.iterrows():
        total = float(row[factor_col])
        ct1 = row[COLUMN_NAME_CELLTYPE1]
        ct2 = row[COLUMN_NAME_CELLTYPE2]

        col1 = ct_color_map.get(ct1, "gray")
        col2 = ct_color_map.get(ct2, "gray")

        # Lollipop line
        ax.plot([0, total], [yi, yi], color="gray", linewidth=1.5, zorder=1)

        # Dots at endpoints
        marker_size = 120
        ax.scatter(
            0,
            yi,
            color=col1,
            s=marker_size,
            marker="o",
            zorder=2,
            edgecolors="black",
            linewidth=0.5,
        )
        ax.scatter(
            total,
            yi,
            color=col2,
            s=marker_size,
            marker="o",
            zorder=2,
            edgecolors="black",
            linewidth=0.5,
        )

        # Sender label (left of x=0)
        ax.annotate(
            ct1,
            xy=(0, yi),
            xytext=(sender_dx_pts, 0),
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize=10,
            clip_on=False,
        )

        # Receiver label (right of endpoint)
        ax.annotate(
            ct2,
            xy=(total, yi),
            xytext=(receiver_dx_pts, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=10,
            clip_on=False,
        )

    ax.invert_yaxis()
    ax.set_xlabel(f"Aggregated {factor_col}", fontsize=11)

    # Title
    if title is None:
        title = f"Motif {motif_idx} - Top {len(top_df)} Cell Type Pairs"
        # Add filter info to title
        filter_parts = []
        if sender_type is not None:
            filter_parts.append(f"Sender: {sender_type}")
        if receiver_type is not None:
            filter_parts.append(f"Receiver: {receiver_type}")
        if filter_parts:
            title += f"\n({', '.join(filter_parts)})"
    ax.set_title(title, fontsize=13, fontweight="bold")

    # X-axis limits with padding
    if len(top_df) > 0:
        max_val = float(top_df[factor_col].max())
        ax.set_xlim(-max_val * 0.03, max_val * 1.12)

    # Remove top/right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend
    if show_legend and created_fig:
        legend_handles = []

        # Cell type colors
        for ct, col in ct_color_map.items():
            legend_handles.append(
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=col,
                    markersize=8,
                    markeredgecolor="black",
                    markeredgewidth=0.5,
                    label=ct,
                )
            )

        ax.legend(
            handles=legend_handles,
            title="Cell Types",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            fontsize=8,
            title_fontsize=9,
            frameon=True,
            fancybox=True,
        )

    if created_fig:
        plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, bbox_inches="tight", dpi=300)
        logger.debug(f"Saved: {save_path}")

    return fig if created_fig else None


def plot_top_lri_interactions_by_pathway(
    lri_motifs_df: pd.DataFrame,
    top_n: int = 35,
    save_path: str | None = None,
    n_cols: int = 3,
    factor_col: str = COLUMN_NAME_FACTOR,
    figsize_per_motif: tuple[float, float] = (12, 9),
    ct_colors: dict[str, str] | None = None,
) -> plt.Figure:
    """Plot top LRI interactions per motif aggregated by pathway

    Same sender, receiver, and pathway combinations are merged (factors summed).

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with parsed components and pathway column
    top_n : int, default 35
        Number of top pathways to show per motif
    save_path : str, optional
        Path to save figure. If None, displays plot.
    n_cols : int, default 3
        Number of columns in subplot grid
    factor_col : str, default "factor"
        Column name for factor values
    figsize_per_motif : tuple, default (12, 9)
        Size per motif subplot
    ct_colors : dict, optional
        Cell type color mapping. If None, uses global colors from set_celltype_colors()
        or auto-generates from data.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object
    """
    from alarmist.plotting.colors import _get_colors_for_plotting

    if COLUMN_NAME_PATHWAY not in lri_motifs_df.columns:
        raise ValueError(
            f"DataFrame must have {COLUMN_NAME_PATHWAY!r} column. Use annotate_pathways() first."
        )

    # Setup colors
    all_celltypes = list(
        set(lri_motifs_df[COLUMN_NAME_CELLTYPE1])
        | set(lri_motifs_df[COLUMN_NAME_CELLTYPE2])
    )
    ct_color_map = _get_colors_for_plotting(ct_colors, all_celltypes)

    motifs = sorted(lri_motifs_df[COLUMN_NAME_MOTIF_IDX].unique())
    n_prog = len(motifs)

    # Layout
    n_rows = int(np.ceil(n_prog / n_cols))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * figsize_per_motif[0], n_rows * figsize_per_motif[1]),
        constrained_layout=False,
    )
    axes = axes.flatten()

    fig.subplots_adjust(
        right=0.85, left=0.15, top=0.94, bottom=0.06, wspace=1.35, hspace=0.4
    )

    for ax, prog in zip(axes, motifs):
        dfp = lri_motifs_df[lri_motifs_df[COLUMN_NAME_MOTIF_IDX] == prog].copy()
        if dfp.empty:
            ax.set_visible(False)
            continue

        # Aggregate by sender, receiver, and pathway
        pathway_agg = (
            dfp.groupby(
                [COLUMN_NAME_CELLTYPE1, COLUMN_NAME_CELLTYPE2, COLUMN_NAME_PATHWAY]
            )
            .agg(
                {
                    factor_col: "sum",
                    COLUMN_NAME_LIGAND: lambda x: ", ".join(x.unique()[:2]),
                    COLUMN_NAME_RECEPTOR: lambda x: ", ".join(x.unique()[:2]),
                }
            )
            .reset_index()
        )

        # Sort and get top N
        pathway_agg = pathway_agg.sort_values(factor_col, ascending=False)
        top_df = pathway_agg.head(top_n).reset_index(drop=True)
        y = np.arange(len(top_df))

        # Draw line bars with dots
        for yi, row in top_df.iterrows():
            total = row[factor_col]

            c1 = row[COLUMN_NAME_CELLTYPE1]
            col1 = ct_color_map.get(c1, "gray")
            c2 = row[COLUMN_NAME_CELLTYPE2]
            col2 = ct_color_map.get(c2, "gray")

            ax.plot([0, total], [yi, yi], color="gray", linewidth=1.5, zorder=1)

            marker_size = 80
            ax.scatter(
                0,
                yi,
                color=col1,
                s=marker_size,
                marker="o",
                zorder=2,
                edgecolors="black",
                linewidth=0.5,
            )
            ax.scatter(
                total,
                yi,
                color=col2,
                s=marker_size,
                marker="o",
                zorder=2,
                edgecolors="black",
                linewidth=0.5,
            )

        # Create labels
        labels = []
        for _, row in top_df.iterrows():
            sender = row[COLUMN_NAME_CELLTYPE1].ljust(12)[:12]
            receiver = row[COLUMN_NAME_CELLTYPE2].ljust(12)[:12]
            pathway = row[COLUMN_NAME_PATHWAY].ljust(20)[:20]
            label = f"{sender} → {receiver} | {pathway}"
            labels.append(label)

        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9, fontfamily="monospace")

        ax.invert_yaxis()
        ax.set_xlabel("Aggregated Normalized Factor", fontsize=10)
        ax.set_title(f"Motif {prog}", fontsize=12, fontweight="bold")

        if len(top_df) > 0:
            max_val = top_df[factor_col].max()
            ax.set_xlim(-max_val * 0.02, max_val * 1.05)

    # Hide unused subplots
    for ax in axes[n_prog:]:
        ax.set_visible(False)

    # Create legend
    legend_handles = []
    for ct, col in ct_color_map.items():
        legend_handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=col,
                markersize=8,
                markeredgecolor="black",
                markeredgewidth=0.5,
                label=ct,
            )
        )

    fig.legend(
        handles=legend_handles,
        title="Cell Types",
        loc="center right",
        bbox_to_anchor=(0.92, 0.5),
        fontsize=8,
        title_fontsize=10,
        frameon=True,
        fancybox=True,
        shadow=False,
    )

    fig.suptitle(
        "Top Pathways per Motif (Aggregated)\n"
        + "Format: Sender → Receiver | Pathway\n"
        + "(Start dot = Sender, End dot = Receiver)",
        fontsize=14,
        y=1,
        fontweight="bold",
    )

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, bbox_inches="tight", dpi=300)
        logger.debug(f"Saved: {save_path}")

    return fig


def build_master_edge_gate(
    lri_motifs_df: pd.DataFrame,
    top_n: int = 200,
    factor_col: str = COLUMN_NAME_FACTOR,
    threshold: float = 1500,
) -> dict[int, set[tuple[str, str]]]:
    """Build edge gates for network filtering

    For each motif, aggregate across all signaling types and keep edges
    (cell1, cell2) whose total weight > threshold.

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with parsed LRI components
    top_n : int, default 200
        Number of top interactions to consider
    threshold : float, default 1500
        Threshold for edge weight

    Returns
    -------
    dict
        Mapping from motif_idx to set of (cell1, cell2) tuples
    """
    gates = {}
    motifs = sorted(lri_motifs_df[COLUMN_NAME_MOTIF_IDX].unique())

    for motif in motifs:
        df = lri_motifs_df[lri_motifs_df[COLUMN_NAME_MOTIF_IDX] == motif].copy()
        if df.empty:
            gates[motif] = set()
            continue

        # Take top N to reduce clutter
        df2 = df.nlargest(top_n, factor_col)

        # Aggregate by cell pair
        agg_all = (
            df2.groupby([COLUMN_NAME_CELLTYPE1, COLUMN_NAME_CELLTYPE2])[factor_col]
            .sum()
            .reset_index(name=COLUMN_NAME_WEIGHT)
        )
        keep = agg_all.loc[
            agg_all[COLUMN_NAME_WEIGHT] > threshold,
            [COLUMN_NAME_CELLTYPE1, COLUMN_NAME_CELLTYPE2],
        ]
        gates[motif] = set(map(tuple, keep.values))

    return gates


def _select_outlier_celltype_pairs(
    agg: pd.DataFrame,
    max_celltypes: int,
    mod_z_thresh: float,
) -> pd.DataFrame:
    """Select edges whose aggregate weight is a statistical outlier (MAD-based
    modified z-score), then greedily keep edges in descending weight order
    while the number of distinct celltype endpoints stays <= max_celltypes.
    """
    if agg.empty or max_celltypes < 2:
        return agg.iloc[0:0]

    w = agg[COLUMN_NAME_WEIGHT].to_numpy(dtype=float)
    med = float(np.median(w))
    mad = float(np.median(np.abs(w - med)))
    if mad > 0:
        score = 0.6745 * (w - med) / mad
    else:
        sd = float(w.std())
        score = (w - float(w.mean())) / sd if sd > 0 else np.zeros_like(w)
    outliers = agg.loc[score > mod_z_thresh].sort_values(
        COLUMN_NAME_WEIGHT, ascending=False
    )
    if outliers.empty:
        return outliers

    kept_idx: list = []
    kept_cells: set[str] = set()
    for idx, row in outliers.iterrows():
        candidate = kept_cells | {
            row[COLUMN_NAME_CELLTYPE1],
            row[COLUMN_NAME_CELLTYPE2],
        }
        if len(candidate) <= max_celltypes:
            kept_idx.append(idx)
            kept_cells = candidate
    return outliers.loc[kept_idx]


def plot_lri_networks(
    lri_motifs_df: pd.DataFrame,
    threshold: float = 0.0,
    top_n: int = 200,
    factor_col: str = COLUMN_NAME_FACTOR,
    annotate_edges: bool = False,
    mode_filter: str | None = None,
    edge_gate: dict[int, set[tuple[str, str]]] | None = None,
    save_path: str | None = None,
    n_cols: int = 5,
    figsize_per_motif: tuple[float, float] = (3, 3),
    ct_colors: dict[str, str] | None = None,
    max_celltypes: int = 4,
    outlier_z: float = 3.5,
) -> plt.Figure:
    """Plot LRI networks for each motif using Graphviz

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with parsed LRI components
    threshold : float, default 0.0
        Hard cutoff on aggregate edge weight. Only used if > 0 and
        ``edge_gate`` is None; otherwise statistical outlier selection is used.
    top_n : int, default 200
        Number of top interactions to consider
    factor_col : str, default "factor"
        Column name for factor values
    annotate_edges : bool, default False
        Whether to annotate edges with L-R pairs
    mode_filter : str, optional
        Filter by signaling mode: None (all), 'paracrine', 'juxtacrine', 'autocrine'
    edge_gate : dict, optional
        Pre-computed edge gates from build_master_edge_gate()
    save_path : str, optional
        Path to save figure. If None, displays plot.
    n_cols : int, default 5
        Number of columns in subplot grid
    figsize_per_motif : tuple, default (3, 3)
        Size per motif subplot
    ct_colors : dict, optional
        Cell type color mapping. If None, uses global colors from set_celltype_colors()
        or auto-generates from data.
    max_celltypes : int, default 4
        When using statistical outlier selection (no ``edge_gate`` and
        ``threshold <= 0``), cap the number of distinct celltypes shown per
        motif. Outlier edges are added in descending weight order until adding
        another would exceed this cap.
    outlier_z : float, default 3.5
        Modified z-score threshold (MAD-based) used to flag a celltype-pair's
        aggregate weight as an outlier. Falls back to a std-based z-score when
        MAD is zero.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object

    Raises
    ------
    ImportError
        If Graphviz is not installed
    """
    from alarmist.plotting.colors import _get_colors_for_plotting

    if not GRAPHVIZ_AVAILABLE:
        raise ImportError(
            "Graphviz required for network plots. Install: pip install graphviz pillow"
        )

    # Setup colors
    all_celltypes = list(
        set(lri_motifs_df[COLUMN_NAME_CELLTYPE1])
        | set(lri_motifs_df[COLUMN_NAME_CELLTYPE2])
    )
    ct_color_map = _get_colors_for_plotting(ct_colors, all_celltypes)

    # Convert RGBA tuples to hex strings for Graphviz compatibility
    def _to_hex(color):
        """Convert color (tuple or string) to hex string for Graphviz."""
        if isinstance(color, str):
            return color
        elif isinstance(color, tuple):
            # RGBA tuple from matplotlib colormap
            r, g, b = int(color[0] * 255), int(color[1] * 255), int(color[2] * 255)
            return f"#{r:02x}{g:02x}{b:02x}"
        else:
            return "#CCCCCC"

    ct_color_map_hex = {ct: _to_hex(col) for ct, col in ct_color_map.items()}

    def _norm_mode(x):
        t = str(x).strip().lower()
        if t in {"auto", "autocrine", "a"}:
            return "autocrine"
        if t in {"juxta", "juxtacrine", "j"}:
            return "juxtacrine"
        return "paracrine"

    motifs = sorted(lri_motifs_df[COLUMN_NAME_MOTIF_IDX].unique())
    n_motifs = len(motifs)
    n_rows = (n_motifs + n_cols - 1) // n_cols

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * figsize_per_motif[0], n_rows * figsize_per_motif[1]),
        constrained_layout=False,
    )
    axes = axes.flatten()
    fig.subplots_adjust(right=0.90, bottom=0.05)

    for i, motif in enumerate(motifs):
        ax = axes[i]
        df = lri_motifs_df[lri_motifs_df[COLUMN_NAME_MOTIF_IDX] == motif].copy()
        if df.empty:
            ax.axis("off")
            continue

        # Normalize mode
        df[COLUMN_NAME_MODE] = df[COLUMN_NAME_SIGNALING_TYPE].apply(_norm_mode)

        # Filter by mode if specified
        if mode_filter is not None:
            mf = _norm_mode(mode_filter)
            df = df[df[COLUMN_NAME_MODE] == mf]
        if df.empty:
            ax.axis("off")
            continue

        df2 = df.nlargest(top_n, factor_col)

        # Aggregate by cell pair
        agg = (
            df2.groupby([COLUMN_NAME_CELLTYPE1, COLUMN_NAME_CELLTYPE2])[factor_col]
            .sum()
            .reset_index(name=COLUMN_NAME_WEIGHT)
        )

        # Apply edge gate, hard threshold, or statistical outlier selection
        if edge_gate is not None:
            keep = edge_gate.get(motif, set())
            if keep:
                mask = [
                    (r[COLUMN_NAME_CELLTYPE1], r[COLUMN_NAME_CELLTYPE2]) in keep
                    for _, r in agg.iterrows()
                ]
                agg = agg.loc[mask]
        elif threshold > 0:
            agg = agg[agg[COLUMN_NAME_WEIGHT] > threshold]
        else:
            agg = _select_outlier_celltype_pairs(agg, max_celltypes, outlier_z)

        if agg.empty:
            ax.axis("off")
            continue

        # Edge annotations
        edge_annotations = {}
        if annotate_edges:
            df_sorted = df.sort_values(factor_col, ascending=False)
            df_filtered = df_sorted.drop_duplicates(
                subset=[COLUMN_NAME_RECEPTOR], keep="first"
            )
            top20_df = df_filtered.nlargest(20, factor_col)[
                [
                    COLUMN_NAME_CELLTYPE1,
                    COLUMN_NAME_CELLTYPE2,
                    COLUMN_NAME_LIGAND,
                    COLUMN_NAME_RECEPTOR,
                ]
            ]
            top20_set = set(
                (
                    r[COLUMN_NAME_CELLTYPE1],
                    r[COLUMN_NAME_CELLTYPE2],
                    r[COLUMN_NAME_LIGAND],
                    r[COLUMN_NAME_RECEPTOR],
                )
                for _, r in top20_df.iterrows()
            )

            for _, row in agg.iterrows():
                key = (row[COLUMN_NAME_CELLTYPE1], row[COLUMN_NAME_CELLTYPE2])
                lr_pairs = df2[
                    (df2[COLUMN_NAME_CELLTYPE1] == row[COLUMN_NAME_CELLTYPE1])
                    & (df2[COLUMN_NAME_CELLTYPE2] == row[COLUMN_NAME_CELLTYPE2])
                ].sort_values(factor_col, ascending=False)
                if not lr_pairs.empty:
                    lr_pairs[COLUMN_NAME_IN_TOP20] = lr_pairs.apply(
                        lambda x: (
                            (
                                x[COLUMN_NAME_CELLTYPE1],
                                x[COLUMN_NAME_CELLTYPE2],
                                x[COLUMN_NAME_LIGAND],
                                x[COLUMN_NAME_RECEPTOR],
                            )
                            in top20_set
                        ),
                        axis=1,
                    )
                    top_sel = lr_pairs[lr_pairs[COLUMN_NAME_IN_TOP20]].head(3)
                    if top_sel.empty:
                        top_sel = lr_pairs.head(3)
                    edge_annotations[key] = "\n".join(
                        f"{r[COLUMN_NAME_LIGAND]}-{r[COLUMN_NAME_RECEPTOR]}"
                        for _, r in top_sel.iterrows()
                    )

        # Build Graphviz graph
        dot = Digraph(engine="neato", format="png")
        dot.graph_attr["dpi"] = "300"
        dot.attr(splines="curved")
        dot.attr(
            "node",
            shape="circle",
            style="filled",
            fontsize="8",
            width="0.5",
            height="0.5",
            fixedsize="true",
        )

        cells = sorted(
            set(agg[COLUMN_NAME_CELLTYPE1]) | set(agg[COLUMN_NAME_CELLTYPE2])
        )
        for c in cells:
            dot.node(c, fillcolor=ct_color_map_hex.get(c, "#CCCCCC"))

        max_w = agg[COLUMN_NAME_WEIGHT].max()
        for _, row in agg.iterrows():
            pen = str(max(0.4, (2.5 * row[COLUMN_NAME_WEIGHT] / max_w)))
            if annotate_edges:
                label = edge_annotations.get(
                    (row[COLUMN_NAME_CELLTYPE1], row[COLUMN_NAME_CELLTYPE2]), ""
                )
                dot.edge(
                    row[COLUMN_NAME_CELLTYPE1],
                    row[COLUMN_NAME_CELLTYPE2],
                    penwidth=pen,
                    color="black",
                    label=label,
                    fontsize="6",
                    fontcolor="darkblue",
                )
            else:
                dot.edge(
                    row[COLUMN_NAME_CELLTYPE1],
                    row[COLUMN_NAME_CELLTYPE2],
                    penwidth=pen,
                    color="black",
                )

        # Render
        try:
            png = dot.pipe()
            img = Image.open(io.BytesIO(png))
            ax.imshow(img)
            ax.set_title(f"Motif {motif}", fontsize=8)
            ax.axis("off")
        except Exception:
            ax.text(
                0.5,
                0.5,
                f"Error rendering\nMotif {motif}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")

    # Hide unused subplots
    for ax in axes[len(motifs) :]:
        ax.axis("off")

    # Legend
    node_handles = [Patch(facecolor=col, label=ct) for ct, col in ct_color_map.items()]
    fig.legend(
        handles=node_handles,
        title="Cell Types",
        loc="center right",
        bbox_to_anchor=(0.98, 0.5),
        fontsize=6,
        title_fontsize=8,
        ncol=1,
    )

    mf_text = "" if mode_filter is None else f" | only {_norm_mode(mode_filter)}"
    edge_text = " (with LR annotations)" if annotate_edges else ""
    if edge_gate is not None:
        gate_text = " | gated by ALL>threshold"
    elif threshold > 0:
        gate_text = f" | weight>{threshold:g}"
    else:
        gate_text = f" | MAD-outliers (z>{outlier_z:g}, ≤{max_celltypes} celltypes)"
    fig.suptitle(
        f"LRI Networks Across Motifs{mf_text}{gate_text}{edge_text}",
        fontsize=12,
        y=0.98,
    )

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=600, bbox_inches="tight")
        logger.debug(f"Saved: {save_path}")

    return fig


def plot_lri_networks_html(
    lri_motifs_df: pd.DataFrame,
    save_path: str,
    top_n: int = 200,
    factor_col: str = COLUMN_NAME_FACTOR,
    mode_filter: str | None = None,
    ct_colors: dict[str, str] | None = None,
    layout_seed: int = 0,
) -> str:
    """Interactive HTML version of :func:`plot_lri_networks`.

    Renders a standalone HTML file (no external libraries) with:
      * a motif selector,
      * a weight-cutoff slider — edges below the cutoff are hidden and any
        celltype with no remaining edges is dropped,
      * click-on-edge → table of every LR pair between those two celltypes,
        sorted by motif loading (``factor_col``) descending.

    Node positions are precomputed in Python via :func:`networkx.spring_layout`
    so layout stays stable as the slider moves.

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        Same shape consumed by :func:`plot_lri_networks` (columns include
        ``motif_idx``, ``celltype1``, ``celltype2``, ``ligand``, ``receptor``,
        ``signaling_type``, and ``factor_col``).
    save_path : str
        Destination ``.html`` path.
    top_n, factor_col, mode_filter, ct_colors
        Same semantics as :func:`plot_lri_networks`.
    layout_seed : int, default 0
        Seed for the spring layout so the file is reproducible.

    Returns
    -------
    str
        The written HTML path.
    """
    import networkx as nx

    from alarmist.plotting.colors import _get_colors_for_plotting

    def _norm_mode(x: str) -> str:
        t = str(x).strip().lower()
        if t in {"auto", "autocrine", "a"}:
            return "autocrine"
        if t in {"juxta", "juxtacrine", "j"}:
            return "juxtacrine"
        return "paracrine"

    def _to_hex(color) -> str:
        if isinstance(color, str):
            return color
        if isinstance(color, tuple):
            r, g, b = int(color[0] * 255), int(color[1] * 255), int(color[2] * 255)
            return f"#{r:02x}{g:02x}{b:02x}"
        return "#cccccc"

    all_celltypes = sorted(
        set(lri_motifs_df[COLUMN_NAME_CELLTYPE1])
        | set(lri_motifs_df[COLUMN_NAME_CELLTYPE2])
    )
    ct_color_map = _get_colors_for_plotting(ct_colors, all_celltypes)
    ct_color_hex = {ct: _to_hex(c) for ct, c in ct_color_map.items()}

    motifs_payload: list[dict] = []
    global_min_w = float("inf")
    global_max_w = 0.0

    for motif in sorted(lri_motifs_df[COLUMN_NAME_MOTIF_IDX].unique()):
        df = lri_motifs_df[lri_motifs_df[COLUMN_NAME_MOTIF_IDX] == motif].copy()
        if df.empty:
            continue
        df[COLUMN_NAME_MODE] = df[COLUMN_NAME_SIGNALING_TYPE].apply(_norm_mode)
        if mode_filter is not None:
            df = df[df[COLUMN_NAME_MODE] == _norm_mode(mode_filter)]
        if df.empty:
            continue

        df_top = df.nlargest(top_n, factor_col)
        agg = (
            df_top.groupby([COLUMN_NAME_CELLTYPE1, COLUMN_NAME_CELLTYPE2])[factor_col]
            .sum()
            .reset_index(name=COLUMN_NAME_WEIGHT)
        )
        if agg.empty:
            continue

        edges = []
        for _, row in agg.iterrows():
            c1, c2, w = (
                row[COLUMN_NAME_CELLTYPE1],
                row[COLUMN_NAME_CELLTYPE2],
                float(row[COLUMN_NAME_WEIGHT]),
            )
            pairs = df[
                (df[COLUMN_NAME_CELLTYPE1] == c1) & (df[COLUMN_NAME_CELLTYPE2] == c2)
            ].sort_values(factor_col, ascending=False)
            lr = [
                {
                    COLUMN_NAME_LIGAND: p[COLUMN_NAME_LIGAND],
                    COLUMN_NAME_RECEPTOR: p[COLUMN_NAME_RECEPTOR],
                    COLUMN_NAME_FACTOR: float(p[factor_col]),
                    COLUMN_NAME_MODE: p[COLUMN_NAME_MODE],
                }
                for _, p in pairs.iterrows()
            ]
            edges.append(
                {
                    COLUMN_NAME_SOURCE: c1,
                    COLUMN_NAME_TARGET: c2,
                    COLUMN_NAME_WEIGHT: w,
                    "lr_pairs": lr,
                }
            )
            global_min_w = min(global_min_w, w)
            global_max_w = max(global_max_w, w)

        nodes_set = sorted(
            set(agg[COLUMN_NAME_CELLTYPE1]) | set(agg[COLUMN_NAME_CELLTYPE2])
        )
        G = nx.DiGraph()
        G.add_nodes_from(nodes_set)
        for e in edges:
            G.add_edge(
                e[COLUMN_NAME_SOURCE],
                e[COLUMN_NAME_TARGET],
                weight=e[COLUMN_NAME_WEIGHT],
            )
        pos = nx.spring_layout(G, seed=layout_seed, k=None, iterations=200)

        nodes = [
            {
                "id": n,
                "x": float(pos[n][0]),
                "y": float(pos[n][1]),
                "color": ct_color_hex.get(n, "#cccccc"),
            }
            for n in nodes_set
        ]
        motifs_payload.append(
            {COLUMN_NAME_MOTIF_IDX: int(motif), "nodes": nodes, "edges": edges}
        )

    if not motifs_payload:
        global_min_w, global_max_w = 0.0, 1.0
    if global_min_w == float("inf"):
        global_min_w = 0.0

    payload = {
        "motifs": motifs_payload,
        "min_weight": global_min_w,
        "max_weight": global_max_w,
        "colors": ct_color_hex,
        "mode_filter": mode_filter,
        "top_n": top_n,
    }

    html = _LRI_NETWORK_HTML_TEMPLATE.replace("__PAYLOAD__", json.dumps(payload))

    if os.path.dirname(save_path):
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
    with open(save_path, "w") as f:
        f.write(html)
    logger.debug(f"Saved: {save_path}")
    return save_path


_LRI_NETWORK_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>LRI Networks</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>
  :root {
    color-scheme: light;
    --bg: #f8fafc;
    --surface: #ffffff;
    --border: #e2e8f0;
    --border-strong: #cbd5e1;
    --muted: #64748b;
    --muted-2: #94a3b8;
    --text: #0f172a;
    --text-2: #334155;
    --accent: #6366f1;
    --accent-2: #4f46e5;
    --accent-soft: #eef2ff;
    --danger: #ef4444;
    --radius: 10px;
    --radius-sm: 6px;
    --shadow-sm: 0 1px 2px rgba(15,23,42,.04), 0 1px 1px rgba(15,23,42,.03);
    --shadow-md: 0 4px 12px -2px rgba(15,23,42,.08), 0 2px 6px -2px rgba(15,23,42,.05);
    --label: 11px;
  }
  * { box-sizing: border-box; }
  html, body { height: 100%; }
  body {
    margin: 0;
    font-family: "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    font-size: 13px; line-height: 1.5;
    color: var(--text); background: var(--bg);
    -webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale;
  }
  header {
    padding: 14px 20px;
    background: var(--surface);
    border-bottom: 1px solid var(--border);
    display: flex; align-items: center; gap: 20px; flex-wrap: wrap;
    box-shadow: var(--shadow-sm);
    position: relative; z-index: 2;
  }
  .brand {
    font-weight: 700; font-size: 13px; letter-spacing: -0.01em;
    color: var(--text);
    display: flex; align-items: center; gap: 8px;
  }
  .brand::before {
    content: ""; width: 8px; height: 8px; border-radius: 2px;
    background: linear-gradient(135deg, var(--accent), #8b5cf6);
  }
  .control { display: flex; flex-direction: column; gap: 4px; }
  .control-label {
    font-size: var(--label); font-weight: 500;
    color: var(--muted); text-transform: uppercase; letter-spacing: 0.04em;
  }
  .control-row { display: flex; align-items: center; gap: 10px; }
  select {
    appearance: none; -webkit-appearance: none;
    background: var(--surface); color: var(--text);
    border: 1px solid var(--border); border-radius: var(--radius-sm);
    padding: 6px 28px 6px 10px; font: inherit; font-weight: 500;
    cursor: pointer; outline: none;
    background-image: url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='12' height='12' viewBox='0 0 24 24' fill='none' stroke='%2364748b' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><polyline points='6 9 12 15 18 9'></polyline></svg>");
    background-repeat: no-repeat; background-position: right 8px center;
    transition: border-color .12s, box-shadow .12s;
  }
  select:hover { border-color: var(--border-strong); }
  select:focus { border-color: var(--accent); box-shadow: 0 0 0 3px rgba(99,102,241,.15); }

  input[type=range] {
    -webkit-appearance: none; appearance: none;
    width: 240px; height: 6px;
    background: var(--border); border-radius: 999px; outline: none; cursor: pointer;
    transition: background .12s;
  }
  input[type=range]:hover { background: var(--border-strong); }
  input[type=range]::-webkit-slider-thumb {
    -webkit-appearance: none; appearance: none;
    width: 16px; height: 16px; border-radius: 50%;
    background: var(--surface); border: 2px solid var(--accent);
    box-shadow: var(--shadow-sm); cursor: pointer;
    transition: transform .12s, box-shadow .12s;
  }
  input[type=range]::-webkit-slider-thumb:hover { transform: scale(1.1); box-shadow: 0 0 0 4px rgba(99,102,241,.15); }
  input[type=range]::-moz-range-thumb {
    width: 14px; height: 14px; border-radius: 50%;
    background: var(--surface); border: 2px solid var(--accent);
    box-shadow: var(--shadow-sm); cursor: pointer;
  }
  #cutoff-val {
    font-family: "JetBrains Mono", ui-monospace, SFMono-Regular, monospace;
    font-size: 12px; font-weight: 500; color: var(--text-2);
    min-width: 56px; padding: 2px 8px;
    background: var(--bg); border: 1px solid var(--border);
    border-radius: var(--radius-sm);
  }

  .legend { display: flex; flex-wrap: wrap; gap: 6px; margin-left: auto; max-width: 50%; }
  .legend-chip {
    display: inline-flex; align-items: center; gap: 6px;
    padding: 3px 8px 3px 6px;
    background: var(--bg); border: 1px solid var(--border);
    border-radius: 999px; font-size: 11px; color: var(--text-2); font-weight: 500;
  }
  .legend-chip .swatch { width: 8px; height: 8px; border-radius: 2px; }

  main { display: grid; grid-template-columns: 1fr 400px; height: calc(100vh - 65px); }
  #plot {
    background: var(--surface);
    border-right: 1px solid var(--border);
    position: relative;
    background-image:
      radial-gradient(circle at 1px 1px, rgba(15,23,42,.04) 1px, transparent 0);
    background-size: 24px 24px;
  }
  #plot svg { width: 100%; height: 100%; display: block; }

  #side {
    background: var(--surface); overflow: auto; padding: 20px;
  }
  #side h2 {
    font-size: var(--label); font-weight: 600; margin: 0 0 12px;
    color: var(--muted); text-transform: uppercase; letter-spacing: 0.06em;
  }
  .hint {
    color: var(--muted); font-size: 12px;
    padding: 16px; background: var(--bg);
    border: 1px dashed var(--border); border-radius: var(--radius);
    text-align: center;
  }
  .edge-header {
    background: var(--accent-soft); border: 1px solid #e0e7ff;
    border-radius: var(--radius); padding: 12px 14px; margin-bottom: 14px;
  }
  .edge-header .pair {
    font-size: 14px; font-weight: 600; color: var(--text);
    display: flex; align-items: center; gap: 8px;
  }
  .edge-header .pair .arrow {
    color: var(--accent); font-weight: 700;
  }
  .edge-header .meta {
    display: flex; gap: 14px; margin-top: 6px;
    font-size: 11px; color: var(--muted); font-weight: 500;
  }
  .edge-header .meta b { color: var(--text-2); font-weight: 600; font-variant-numeric: tabular-nums; }

  table {
    width: 100%; border-collapse: separate; border-spacing: 0;
    font-size: 12px;
    border: 1px solid var(--border); border-radius: var(--radius); overflow: hidden;
  }
  th {
    text-align: left; padding: 8px 10px; background: var(--bg);
    color: var(--muted); font-weight: 600; font-size: 10.5px;
    text-transform: uppercase; letter-spacing: 0.04em;
    border-bottom: 1px solid var(--border);
    position: sticky; top: 0;
  }
  td {
    padding: 8px 10px; border-bottom: 1px solid var(--border);
    color: var(--text-2);
  }
  tr:last-child td { border-bottom: none; }
  tr:hover td { background: var(--bg); }
  td.num { text-align: right; font-variant-numeric: tabular-nums; font-weight: 500; color: var(--text); }
  td.mono { font-family: "JetBrains Mono", ui-monospace, SFMono-Regular, monospace; font-size: 11px; }
  .mode-pill {
    display: inline-block; padding: 1px 6px; border-radius: 999px;
    font-size: 10px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.04em;
  }
  .mode-paracrine { background: #dbeafe; color: #1e40af; }
  .mode-juxtacrine { background: #fce7f3; color: #9d174d; }
  .mode-autocrine { background: #d1fae5; color: #065f46; }

  .node-label {
    font-family: "Inter", sans-serif; font-size: 11px; font-weight: 600;
    pointer-events: none; user-select: none; fill: var(--text);
  }
  .node {
    stroke: rgba(15,23,42,.18); stroke-width: 1.5;
    filter: drop-shadow(0 1px 2px rgba(15,23,42,.08));
    transition: stroke-width .12s, transform .12s;
  }
  .node:hover { stroke-width: 2.5; }
  .edge {
    stroke: #94a3b8; stroke-opacity: 0.55; cursor: pointer; fill: none;
    transition: stroke .12s, stroke-opacity .12s;
  }
  .edge:hover { stroke: var(--accent); stroke-opacity: 1; }
  .edge.selected { stroke: var(--accent-2); stroke-opacity: 1; }
  .empty {
    fill: var(--muted-2); font-size: 13px; font-weight: 500;
  }
</style>
</head>
<body>
<header>
  <div class="brand">LRI Networks</div>
  <div class="control">
    <span class="control-label">Motif</span>
    <div class="control-row"><select id="motif"></select></div>
  </div>
  <div class="control">
    <span class="control-label">Min edge weight</span>
    <div class="control-row">
      <input type="range" id="cutoff" />
      <span id="cutoff-val"></span>
    </div>
  </div>
  <div class="legend" id="legend"></div>
</header>
<main>
  <div id="plot"><svg id="svg" viewBox="0 0 1000 700" preserveAspectRatio="xMidYMid meet"></svg></div>
  <aside id="side">
    <h2>Ligand–Receptor Pairs</h2>
    <div id="table-wrap"><div class="hint">Click an edge to list the ligand-receptor pairs between those two celltypes, ranked by motif loading.</div></div>
  </aside>
</main>
<script>
const DATA = __PAYLOAD__;
const W = 1000, H = 700, PAD = 60;

const motifSel = document.getElementById("motif");
const cutoff = document.getElementById("cutoff");
const cutoffVal = document.getElementById("cutoff-val");
const svg = document.getElementById("svg");
const tableWrap = document.getElementById("table-wrap");

DATA.motifs.forEach(m => {
  const opt = document.createElement("option");
  opt.value = m.motif_idx;
  opt.textContent = "Motif " + m.motif_idx;
  motifSel.appendChild(opt);
});

const lo = DATA.min_weight, hi = DATA.max_weight;
const step = (hi - lo) > 0 ? (hi - lo) / 200 : 1;
cutoff.min = lo; cutoff.max = hi; cutoff.step = step; cutoff.value = lo;

const legend = document.getElementById("legend");
Object.entries(DATA.colors).forEach(([ct, col]) => {
  const s = document.createElement("span");
  s.className = "legend-chip";
  s.innerHTML = `<span class="swatch" style="background:${col}"></span>${ct}`;
  legend.appendChild(s);
});

let selectedEdge = null;

function fmt(v) { return Math.abs(v) >= 100 ? v.toFixed(0) : v.toPrecision(3); }

function render() {
  const motifIdx = +motifSel.value;
  const motif = DATA.motifs.find(m => m.motif_idx === motifIdx);
  const cut = +cutoff.value;
  cutoffVal.textContent = fmt(cut);

  svg.innerHTML = "";
  if (!motif) return;

  const edges = motif.edges.filter(e => e.weight >= cut);
  const liveNodes = new Set();
  edges.forEach(e => { liveNodes.add(e.source); liveNodes.add(e.target); });
  const nodes = motif.nodes.filter(n => liveNodes.has(n.id));
  if (nodes.length === 0) {
    const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
    t.setAttribute("x", W / 2); t.setAttribute("y", H / 2);
    t.setAttribute("text-anchor", "middle"); t.setAttribute("class", "empty");
    t.textContent = "No edges above cutoff";
    svg.appendChild(t);
    return;
  }

  const NODE_PAD = 80;
  const xs = nodes.map(n => n.x), ys = nodes.map(n => n.y);
  const xmin = Math.min(...xs), xmax = Math.max(...xs);
  const ymin = Math.min(...ys), ymax = Math.max(...ys);
  const sx = x => NODE_PAD + (xmax === xmin ? (W - 2*NODE_PAD)/2 : (x - xmin) / (xmax - xmin) * (W - 2*NODE_PAD));
  const sy = y => NODE_PAD + (ymax === ymin ? (H - 2*NODE_PAD)/2 : (y - ymin) / (ymax - ymin) * (H - 2*NODE_PAD));
  const pos = Object.fromEntries(nodes.map(n => [n.id, [sx(n.x), sy(n.y)]]));

  const defs = document.createElementNS("http://www.w3.org/2000/svg", "defs");
  defs.innerHTML = `
    <marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5"
        markerWidth="5" markerHeight="5" orient="auto-start-reverse">
        <path d="M0,0 L10,5 L0,10 z" fill="#64748b"/></marker>
    <marker id="arrow-sel" viewBox="0 0 10 10" refX="8" refY="5"
        markerWidth="5" markerHeight="5" orient="auto-start-reverse">
        <path d="M0,0 L10,5 L0,10 z" fill="#4f46e5"/></marker>`;
  svg.appendChild(defs);

  // Layers so nodes always render above edges
  const edgeLayer = document.createElementNS("http://www.w3.org/2000/svg", "g");
  const nodeLayer = document.createElementNS("http://www.w3.org/2000/svg", "g");
  svg.appendChild(edgeLayer);
  svg.appendChild(nodeLayer);

  // Elide over-long cell-type names so circles stay compact; the full name
  // is preserved in a <title> tooltip on each node.
  const LABEL_MAX = 16;
  const elide = s => {
    s = String(s);
    return s.length > LABEL_MAX ? s.slice(0, LABEL_MAX - 1).trimEnd() + "…" : s;
  };

  // Measure each node's (elided) label first to size circles to fit.
  const radii = {};
  nodes.forEach(n => {
    const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
    t.setAttribute("text-anchor", "middle");
    t.setAttribute("dominant-baseline", "central");
    t.setAttribute("class", "node-label");
    t.textContent = elide(n.id);
    nodeLayer.appendChild(t);
    const bb = t.getBBox();
    const r = Math.max(20, Math.ceil(bb.width / 2) + 10);
    radii[n.id] = r;
    n.__textEl = t;
  });

  // Resolve overlaps: iteratively push apart any pair of circles closer than
  // the sum of their radii (plus a margin), clamping back inside the canvas.
  const MARGIN = 16;
  const ids = nodes.map(n => n.id);
  for (let iter = 0; iter < 400; iter++) {
    let moved = false;
    for (let a = 0; a < ids.length; a++) {
      for (let b = a + 1; b < ids.length; b++) {
        const pa = pos[ids[a]], pb = pos[ids[b]];
        const dx = pb[0] - pa[0], dy = pb[1] - pa[1];
        const d = Math.hypot(dx, dy) || 0.01;
        const need = radii[ids[a]] + radii[ids[b]] + MARGIN;
        if (d < need) {
          const push = (need - d) / 2;
          const ux = dx / d, uy = dy / d;
          pa[0] -= ux * push; pa[1] -= uy * push;
          pb[0] += ux * push; pb[1] += uy * push;
          moved = true;
        }
      }
    }
    ids.forEach(id => {
      const r = radii[id], p = pos[id];
      p[0] = Math.min(W - r, Math.max(r, p[0]));
      p[1] = Math.min(H - r, Math.max(r, p[1]));
    });
    if (!moved) break;
  }

  // Place each label at its node's final (post-relaxation) centre.
  nodes.forEach(n => {
    const [x, y] = pos[n.id];
    n.__textEl.setAttribute("x", x);
    n.__textEl.setAttribute("y", y);
  });

  const wMax = Math.max(...edges.map(e => e.weight));

  const attachEdge = (path, e, i) => {
    path.setAttribute("class", "edge");
    const sw = 1.0 + 5.0 * (e.weight / wMax);
    path.setAttribute("stroke-width", sw.toFixed(2));
    path.setAttribute("stroke-linecap", "round");
    path.setAttribute("marker-end", "url(#arrow)");
    path.addEventListener("click", () => {
      selectedEdge = i;
      svg.querySelectorAll(".edge.selected").forEach(el => {
        el.classList.remove("selected");
        el.setAttribute("marker-end", "url(#arrow)");
      });
      path.classList.add("selected");
      path.setAttribute("marker-end", "url(#arrow-sel)");
      showLR(motif, e);
    });
    edgeLayer.appendChild(path);
  };

  edges.forEach((e, i) => {
    const [x1, y1] = pos[e.source], [x2, y2] = pos[e.target];
    const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
    if (e.source === e.target) {
      // Self / autocrine edge: draw a loop leaving and re-entering the top of
      // the node so it is actually visible (a straight self-edge has length 0).
      const r = radii[e.source];
      const a1 = -Math.PI / 2 - 0.5, a2 = -Math.PI / 2 + 0.5;
      const p1x = x1 + r * Math.cos(a1), p1y = y1 + r * Math.sin(a1);
      const p2x = x1 + r * Math.cos(a2), p2y = y1 + r * Math.sin(a2);
      const loop = r * 2.4;
      const c1x = x1 + loop * Math.cos(a1 - 0.35), c1y = y1 + loop * Math.sin(a1 - 0.35);
      const c2x = x1 + loop * Math.cos(a2 + 0.35), c2y = y1 + loop * Math.sin(a2 + 0.35);
      path.setAttribute("d", `M${p1x},${p1y} C${c1x},${c1y} ${c2x},${c2y} ${p2x},${p2y}`);
    } else {
      const r1 = radii[e.source], r2 = radii[e.target];
      const dx = x2 - x1, dy = y2 - y1;
      const len = Math.hypot(dx, dy) || 1;
      const ux = dx / len, uy = dy / len;
      // shrink endpoints so the line/arrow stops at the node border
      const sx1 = x1 + ux * r1, sy1 = y1 + uy * r1;
      const ex2 = x2 - ux * r2, ey2 = y2 - uy * r2;
      const mx = (sx1 + ex2) / 2 + (ey2 - sy1) * 0.12;
      const my = (sy1 + ey2) / 2 - (ex2 - sx1) * 0.12;
      path.setAttribute("d", `M${sx1},${sy1} Q${mx},${my} ${ex2},${ey2}`);
    }
    attachEdge(path, e, i);
  });

  // Now draw the circles behind the (already-appended) labels.
  nodes.forEach(n => {
    const [x, y] = pos[n.id];
    const c = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    c.setAttribute("cx", x); c.setAttribute("cy", y);
    c.setAttribute("r", radii[n.id]);
    c.setAttribute("fill", n.color);
    c.setAttribute("class", "node");
    const title = document.createElementNS("http://www.w3.org/2000/svg", "title");
    title.textContent = n.id;  // full, un-elided name on hover
    c.appendChild(title);
    // insert before the text so text sits on top
    nodeLayer.insertBefore(c, n.__textEl);
    delete n.__textEl;
  });
}

function showLR(motif, edge) {
  const rows = edge.lr_pairs.map(p => {
    const mc = `mode-${p.mode}`;
    return `<tr>
      <td class="mono">${p.ligand}</td>
      <td class="mono">${p.receptor}</td>
      <td><span class="mode-pill ${mc}">${p.mode}</span></td>
      <td class="num">${fmt(p.factor)}</td>
    </tr>`;
  }).join("");
  tableWrap.innerHTML = `
    <div class="edge-header">
      <div class="pair">${edge.source} <span class="arrow">→</span> ${edge.target}</div>
      <div class="meta">
        <span>Aggregate weight <b>${fmt(edge.weight)}</b></span>
        <span><b>${edge.lr_pairs.length}</b> LR pairs</span>
      </div>
    </div>
    <table>
      <thead><tr><th>Ligand</th><th>Receptor</th><th>Mode</th><th class="num">Loading</th></tr></thead>
      <tbody>${rows}</tbody>
    </table>`;
}

function resetSide() {
  tableWrap.innerHTML = '<div class="hint">Click an edge to list the ligand-receptor pairs between those two celltypes, ranked by motif loading.</div>';
}

motifSel.addEventListener("change", () => { selectedEdge = null; resetSide(); render(); });
cutoff.addEventListener("input", render);
render();
</script>
</body>
</html>
"""
