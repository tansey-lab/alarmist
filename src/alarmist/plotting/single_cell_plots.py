"""
Single cell visualization functions for motif analysis
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from typing import Optional, Dict, Tuple
import os


def plot_motif_celltype_composition(
    df_tidy: pd.DataFrame,
    color_map: Optional[Dict] = None,
    figsize: tuple = (12, 6),
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    cell_type_col: str = "cell_type"
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot stacked bar chart showing cell type composition for each motif.

    Parameters
    ----------
    df_tidy : pd.DataFrame
        Tidy dataframe with columns [motif, cell_type, weight],
        typically from weighted_celltypes_by_motif()
    color_map : dict, optional
        Mapping from cell type to color. If None, uses tab20 colormap.
    figsize : tuple, default (12, 6)
        Figure size
    ylabel : str, optional
        Y-axis label. If None, uses "Fraction of cells"
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure. If None, displays without saving.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    ax : matplotlib.axes.Axes
        Axes object

    Examples
    --------
    >>> from alarmist.core import weighted_celltypes_by_motif
    >>> tidy = weighted_celltypes_by_motif(cell_loadings, cell_meta_df)
    >>> fig, ax = plot_motif_celltype_composition(
    ...     tidy,
    ...     title="Cell type composition per motif",
    ...     save_path="motif_celltype.pdf"
    ... )
    """
    # Pivot to wide format
    wide = df_tidy.pivot_table(
        index="motif",
        columns=cell_type_col,
        values="weight",
        aggfunc="sum",
        fill_value=0.0
    )
    wide = wide.loc[sorted(wide.index)]

    # Order cell types by total weight
    celltype_order = wide.sum(axis=0).sort_values(ascending=False).index.tolist()
    wide = wide[celltype_order]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    bottom = np.zeros(len(wide), dtype=float)
    x = np.arange(len(wide))

    # Determine colors
    if color_map is None:
        cmap = plt.get_cmap("tab20")
        colors = [cmap(i % 20) for i in range(len(wide.columns))]
        color_dict = dict(zip(wide.columns, colors))
    else:
        color_dict = color_map

    # Stack bars
    for ct in wide.columns:
        color = color_dict.get(ct, "lightgray")
        ax.bar(x, wide[ct].to_numpy(), bottom=bottom, label=ct, color=color)
        bottom = bottom + wide[ct].to_numpy()

    # Format plot
    ax.set_xticks(x)
    ax.set_xticklabels([str(m) for m in wide.index])
    ax.set_xlabel("Motif")
    ax.set_ylabel(ylabel if ylabel else "Fraction of cells")
    if title:
        ax.set_title(title)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), title="Cell type")
    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig, ax


def plot_motif_state_counts(
    df_counts: pd.DataFrame,
    figsize: tuple = (12, 6),
    colors: list = ['#66c2a5', '#fc8d62'],
    title: Optional[str] = None,
    save_path: Optional[str] = None
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot stacked bar chart of positive vs negative cells per motif.

    Parameters
    ----------
    df_counts : pd.DataFrame
        Counts dataframe with columns [positive, negative],
        typically from compute_motif_state_counts()
    figsize : tuple, default (12, 6)
        Figure size
    colors : list, default ['#66c2a5', '#fc8d62']
        Colors for [positive, negative] bars
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure. If None, displays without saving.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    ax : matplotlib.axes.Axes
        Axes object

    Examples
    --------
    >>> from alarmist.core import compute_motif_state_counts
    >>> counts = compute_motif_state_counts(adata)
    >>> fig, ax = plot_motif_state_counts(
    ...     counts,
    ...     title="Positive vs Negative cells per motif"
    ... )
    """
    fig, ax = plt.subplots(figsize=figsize)

    df_counts.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=colors
    )

    ax.set_xlabel("Motif")
    ax.set_ylabel("Number of cells")
    if title:
        ax.set_title(title)
    else:
        ax.set_title("Positive vs Negative Cells per Motif")
    ax.legend(title="State")
    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig, ax


def plot_positive_motifs_distribution(
    counts: pd.Series,
    figsize: tuple = (10, 4),
    title: Optional[str] = None,
    save_path: Optional[str] = None
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot distribution of number of positive motifs per cell.

    Parameters
    ----------
    counts : pd.Series
        Series with index = number of positive motifs, values = number of cells,
        typically from compute_positive_motifs_per_cell()
    figsize : tuple, default (10, 4)
        Figure size
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure. If None, displays without saving.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    ax : matplotlib.axes.Axes
        Axes object

    Examples
    --------
    >>> from alarmist.core import compute_positive_motifs_per_cell
    >>> dist = compute_positive_motifs_per_cell(adata)
    >>> fig, ax = plot_positive_motifs_distribution(
    ...     dist,
    ...     title="Distribution of positive motifs per cell"
    ... )
    """
    fig, ax = plt.subplots(figsize=figsize)

    counts.plot(kind='bar', ax=ax)

    ax.set_xlabel('Number of positive motifs per cell')
    ax.set_ylabel('Number of cells')
    if title:
        ax.set_title(title)
    else:
        ax.set_title("Distribution of Positive Motifs per Cell")
    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig, ax


def plot_motif_spatial(
    adata,
    motif_id,
    spot_size: float = 20,
    figsize: tuple = (8, 8),
    palette: list = ["#bdbdbd", "#1f77b4"],
    title: Optional[str] = None,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Plot spatial distribution of motif ON/OFF states.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data object with motif_{motif_id}_state in obs
    motif_id : int or list of int
        Motif ID(s) to visualize.
        - If int: shows ON/OFF states for that single motif
        - If list: shows cells that are ON for ALL motifs in the list
    spot_size : float, default 20
        Size of spots in spatial plot
    figsize : tuple, default (8, 8)
        Figure size
    palette : list, default ["#bdbdbd", "#1f77b4"]
        Colors for [negative, positive] states
    title : str, optional
        Plot title. If None, auto-generates based on motif_id
    save_path : str, optional
        Path to save figure. If None, displays without saving.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object

    Examples
    --------
    >>> # Single motif
    >>> fig = plot_motif_spatial(
    ...     adata,
    ...     motif_id=2,
    ...     title="Motif 2 ON/OFF cells",
    ...     save_path="motif_2_spatial.pdf"
    ... )
    >>>
    >>> # Multiple motifs (cells positive for ALL)
    >>> fig = plot_motif_spatial(
    ...     adata,
    ...     motif_id=[2, 5, 8],
    ...     title="Cells ON for motifs 2, 5, and 8",
    ...     save_path="motifs_2_5_8_spatial.pdf"
    ... )
    """
    # Handle single motif or list of motifs
    if isinstance(motif_id, (list, tuple)):
        # Multiple motifs: find cells positive for ALL
        motif_ids = list(motif_id)
        state_cols = [f'motif_{mid}_state' for mid in motif_ids]

        # Check all required columns exist
        missing = [col for col in state_cols if col not in adata.obs.columns]
        if missing:
            raise ValueError(f"Columns {missing} not found in adata.obs. "
                           f"Run gmm_binarize_all_motifs() first.")

        # Create temporary column: positive only if ALL motifs are positive
        temp_col = '_temp_combined_state'
        all_positive = adata.obs[state_cols[0]] == 'positive'
        for col in state_cols[1:]:
            all_positive = all_positive & (adata.obs[col] == 'positive')

        adata.obs[temp_col] = 'negative'
        adata.obs.loc[all_positive, temp_col] = 'positive'

        state_col = temp_col
        default_title = f"Motifs {motif_ids} (ALL positive)"

    else:
        # Single motif
        state_col = f'motif_{motif_id}_state'

        if state_col not in adata.obs.columns:
            raise ValueError(f"Column '{state_col}' not found in adata.obs. "
                           f"Run gmm_binarize_all_motifs() first.")

        default_title = f"Motif {motif_id} Spatial Distribution"

    fig, ax = plt.subplots(figsize=figsize)

    sc.pl.spatial(
        adata,
        color=state_col,
        spot_size=spot_size,
        ax=ax,
        show=False,
        palette=palette  # gray = negative, blue = positive
    )

    if title:
        ax.set_title(title)
    else:
        ax.set_title(default_title)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {save_path}")

    # Clean up temporary column if created
    if isinstance(motif_id, (list, tuple)) and '_temp_combined_state' in adata.obs.columns:
        adata.obs.drop(columns=['_temp_combined_state'], inplace=True)

    return fig
