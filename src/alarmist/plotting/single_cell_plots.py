"""
Single cell visualization functions for motif analysis
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
from typing import Optional, Dict, Tuple, Union
import os

from alarmist.plotting.colors import _get_colors_for_plotting


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
    adata: Union[anndata.AnnData, Dict[str, anndata.AnnData]],
    motif_idx: Union[int, List[int]],
    sample_column: Optional[str] = None,
    cell_type_column: str = 'cell_type',
    n_cols: int = 4,
    figsize_per_panel: tuple = (6, 6),
    point_size: float = 0.5,
    ct_colors: Optional[Dict] = None,
    negative_color: str = '#d3d3d3',
    show_celltype_legend: bool = False,
    legend_top_n: int = 10,
    output_dir: Optional[str] = None,
) -> Union[plt.Figure, List[plt.Figure]]:
    """
    Plot spatial distribution of motif ON/OFF states.

    Positive cells are colored by cell type, negative cells in gray.

    Supports three input modes:
    1. Single AnnData: single panel
    2. Dict of AnnData: grid of panels (one per sample)
    3. Merged AnnData with sample_column: grid of panels (split by sample)

    Parameters
    ----------
    adata : anndata.AnnData or Dict[str, anndata.AnnData]
        Must contain motif_{k}_state columns in obs (from gmm_binarize_all_motifs)
    motif_idx : int or list of int
        Motif index(es) to visualize.
        - If int: plot single motif, return single Figure
        - If list: plot each motif separately, return list of Figures
    sample_column : str, optional
        Column name in adata.obs identifying samples (for merged AnnData).
        If provided with single AnnData, splits into multiple panels.
    cell_type_column : str, default 'cell_type'
        Column name for cell type annotations
    n_cols : int, default 4
        Number of columns in the grid (for multi-sample)
    figsize_per_panel : tuple, default (6, 6)
        Figure size per panel
    point_size : float, default 0.5
        Size of scatter points
    ct_colors : dict, optional
        Cell type color mapping. If None, uses global color registry.
    negative_color : str, default '#d3d3d3'
        Color for negative (OFF) cells
    show_celltype_legend : bool, default False
        Whether to show cell type legend (can be crowded with many types)
    legend_top_n : int, default 10
        If showing legend, only show top N cell types
    output_dir : str, optional
        Directory to save figures. Files named motif_{k}_spatial.png

    Returns
    -------
    plt.Figure or List[plt.Figure]
        Single Figure if motif_idx is int, list of Figures if motif_idx is list

    Examples
    --------
    >>> # Single motif, single sample
    >>> fig = al.plot_motif_spatial(adata, motif_idx=5)
    >>>
    >>> # Single motif, multi-sample dict
    >>> fig = al.plot_motif_spatial(adata_dict, motif_idx=5, output_dir='results/spatial')
    >>>
    >>> # Multiple motifs
    >>> figs = al.plot_motif_spatial(adata_dict, motif_idx=range(20), output_dir='results/spatial')
    >>>
    >>> # Merged mode
    >>> fig = al.plot_motif_spatial(adata_merged, motif_idx=5, sample_column='patient_id')
    """
    import matplotlib.lines as mlines

    # Handle list of motifs - recursive call
    if isinstance(motif_idx, (list, tuple, range)):
        figs = []
        for k in motif_idx:
            fig = plot_motif_spatial(
                adata=adata,
                motif_idx=k,
                sample_column=sample_column,
                cell_type_column=cell_type_column,
                n_cols=n_cols,
                figsize_per_panel=figsize_per_panel,
                point_size=point_size,
                ct_colors=ct_colors,
                negative_color=negative_color,
                show_celltype_legend=show_celltype_legend,
                legend_top_n=legend_top_n,
                output_dir=output_dir,
            )
            figs.append(fig)
            print(f"Plotted motif {k}")
        return figs

    # Get color map
    color_map = _get_colors_for_plotting(ct_colors=ct_colors)
    if not color_map:
        raise ValueError(
            "No cell type colors available. Either:\n"
            "  1. Call al.set_celltype_colors() first, or\n"
            "  2. Pass ct_colors parameter"
        )

    # Convert colors to hex for matplotlib
    def to_hex(c):
        if isinstance(c, str):
            return c
        elif isinstance(c, tuple) and len(c) >= 3:
            return '#{:02x}{:02x}{:02x}'.format(int(c[0]*255), int(c[1]*255), int(c[2]*255))
        return c

    ct_color_hex = {k: to_hex(v) for k, v in color_map.items()}

    # Prepare adata_dict based on input mode
    if isinstance(adata, dict):
        # Mode 2: Dict of AnnData
        adata_dict = adata
        mode_str = "dict"
    elif sample_column is not None:
        # Mode 3: Merged AnnData - split by sample
        adata_dict = {}
        for sample_id in adata.obs[sample_column].unique():
            mask = adata.obs[sample_column] == sample_id
            adata_dict[sample_id] = adata[mask].copy()
        mode_str = "merged"
    else:
        # Mode 1: Single AnnData
        adata_dict = {'sample': adata}
        mode_str = "single"

    n_samples = len(adata_dict)
    state_col = f'motif_{motif_idx}_state'

    # Determine grid layout
    if n_samples == 1:
        n_rows, n_cols_actual = 1, 1
    else:
        n_cols_actual = min(n_cols, n_samples)
        n_rows = (n_samples + n_cols_actual - 1) // n_cols_actual

    fig, axes = plt.subplots(
        n_rows, n_cols_actual,
        figsize=(figsize_per_panel[0] * n_cols_actual, figsize_per_panel[1] * n_rows),
        squeeze=False
    )
    axes = axes.flatten()

    for i, (sample_id, ad) in enumerate(adata_dict.items()):
        ax = axes[i]
        coords = ad.obsm['spatial'][:, :2]

        if state_col not in ad.obs.columns:
            ax.set_title(f'{sample_id}\n(no motif {motif_idx} data)')
            ax.axis('off')
            continue

        states = ad.obs[state_col].astype(str)
        n_pos = (states == 'positive').sum()
        n_neg = (states == 'negative').sum()
        frac_pos = n_pos / len(states) * 100 if len(states) > 0 else 0

        # Plot negatives (background, gray)
        mask_neg = (states == 'negative').values
        ax.scatter(
            coords[mask_neg, 0],
            coords[mask_neg, 1],
            c=negative_color,
            s=point_size,
            alpha=0.95,
            marker='o',
            edgecolors='none',
            linewidths=0,
            rasterized=True,
        )

        # Plot positives (colored by cell type)
        mask_pos = (states == 'positive').values
        if mask_pos.any():
            ct = ad.obs[cell_type_column].astype(str).values
            pos_colors = np.array([ct_color_hex.get(x, "#000000") for x in ct])[mask_pos]

            ax.scatter(
                coords[mask_pos, 0],
                coords[mask_pos, 1],
                c=pos_colors,
                s=point_size * 6,
                alpha=1,
                marker='o',
                edgecolors='none',
                linewidths=0,
                rasterized=True,
                zorder=3
            )

        # Title
        if mode_str == "single":
            ax.set_title(f'Motif {motif_idx}: {frac_pos:.1f}% ON')
        else:
            ax.set_title(f'{sample_id}\nMotif {motif_idx}: {frac_pos:.1f}% ON')

        ax.set_aspect('equal')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

        # Legend
        off_handle = mlines.Line2D([], [], color=negative_color, marker='o',
                                   linestyle='None', markersize=6, label=f'OFF (n={n_neg:,})')
        on_handle = mlines.Line2D([], [], color='black', marker='o',
                                  linestyle='None', markersize=6, label=f'ON (n={n_pos:,}, {frac_pos:.1f}%)')
        handles = [off_handle, on_handle]

        if show_celltype_legend and mask_pos.any():
            ct_on = ad.obs.loc[mask_pos, cell_type_column].astype(str)
            top_cts = ct_on.value_counts().head(legend_top_n).index.tolist()
            ct_handles = [
                mlines.Line2D([], [], color=ct_color_hex.get(ct, "#000000"), marker='o',
                              linestyle='None', markersize=6, label=ct)
                for ct in top_cts
            ]
            handles += ct_handles

        ax.legend(handles=handles, loc='upper right', fontsize=8, frameon=False)

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    if mode_str != "single":
        plt.suptitle(f'Motif {motif_idx} Spatial Distribution (ON colored by cell type)', fontsize=14, y=1.02)

    plt.tight_layout()

    # Save
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f'motif_{motif_idx}_spatial.png')
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig


def analyze_motif_celltype_composition(
    adata: Union[anndata.AnnData, Dict[str, anndata.AnnData]],
    cell_loadings: np.ndarray,
    cell_type_column: str = 'cell_type',
    sample_column: Optional[str] = None,
    normalize: bool = True,
    ct_colors: Optional[Dict] = None,
    figsize: tuple = (12, 6),
    title: Optional[str] = None,
    output_dir: Optional[str] = None
) -> Tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """
    Analyze and visualize cell type composition for each motif.

    Supports three input modes:
    1. Single AnnData: standard single-sample analysis
    2. Dict of AnnData objects: multi-sample with dict
    3. Single merged AnnData with sample_column: multi-sample merged

    Parameters
    ----------
    adata : anndata.AnnData or Dict[str, anndata.AnnData]
        Single AnnData object, or a dictionary mapping sample_id -> AnnData.
    cell_loadings : np.ndarray
        Cell loadings matrix (n_cells x n_motifs) from project_cell_loadings
    cell_type_column : str, default 'cell_type'
        Column name in adata.obs containing cell type annotations
    sample_column : str, optional
        Column name in adata.obs identifying samples (for merged AnnData).
        If provided with single AnnData, treats it as multi-sample.
    normalize : bool, default True
        If True, per-motif weights sum to 1 (for percent stacked bar)
    ct_colors : dict, optional
        Mapping from cell type to color. If None, uses global color registry
        or auto-generates from tab20 colormap.
    figsize : tuple, default (12, 6)
        Figure size
    title : str, optional
        Plot title. If None, auto-generates based on input mode.
    output_dir : str, optional
        Directory to save figure. If None, figure is not saved.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    ax : matplotlib.axes.Axes
        Axes object
    tidy_df : pd.DataFrame
        Tidy dataframe with columns [motif, cell_type, weight]

    Examples
    --------
    >>> # Single sample
    >>> fig, ax, df = al.analyze_motif_celltype_composition(
    ...     adata, cell_loadings,
    ...     output_dir='results/single_cell'
    ... )
    >>>
    >>> # Multi-sample with dict
    >>> fig, ax, df = al.analyze_motif_celltype_composition(
    ...     {'sample_A': adata_a, 'sample_B': adata_b},
    ...     cell_loadings,
    ...     output_dir='results/single_cell'
    ... )
    >>>
    >>> # Multi-sample with merged AnnData
    >>> fig, ax, df = al.analyze_motif_celltype_composition(
    ...     adata_merged, cell_loadings,
    ...     sample_column='patient_id',
    ...     output_dir='results/single_cell'
    ... )
    """
    from alarmist.core.single_cell import weighted_celltypes_by_motif

    # Build cell metadata DataFrame based on input type
    if isinstance(adata, dict):
        # Mode 2: Dict of AnnData
        cell_meta_dfs = []
        for sample_id, ad in adata.items():
            df = ad.obs[[cell_type_column]].copy()
            df.columns = ['cell_type']
            df['sample_id'] = sample_id
            cell_meta_dfs.append(df)
        cell_meta_df = pd.concat(cell_meta_dfs, axis=0).reset_index(drop=True)
        mode_str = f"Multi-sample ({len(adata)} samples)"
    elif sample_column is not None:
        # Mode 3: Merged AnnData with sample_column
        cell_meta_df = adata.obs[[cell_type_column, sample_column]].copy()
        cell_meta_df.columns = ['cell_type', 'sample_id']
        cell_meta_df = cell_meta_df.reset_index(drop=True)
        n_samples = cell_meta_df['sample_id'].nunique()
        mode_str = f"Multi-sample ({n_samples} samples, merged)"
    else:
        # Mode 1: Single AnnData
        cell_meta_df = adata.obs[[cell_type_column]].copy()
        cell_meta_df.columns = ['cell_type']
        cell_meta_df = cell_meta_df.reset_index(drop=True)
        mode_str = "Single sample"

    # Validate dimensions
    if len(cell_meta_df) != cell_loadings.shape[0]:
        raise ValueError(
            f"Number of cells in adata ({len(cell_meta_df)}) does not match "
            f"cell_loadings rows ({cell_loadings.shape[0]})"
        )

    print(f"Analyzing cell type composition per motif...")
    print(f"  Mode: {mode_str}")
    print(f"  Cells: {len(cell_meta_df):,}")
    print(f"  Motifs: {cell_loadings.shape[1]}")

    # Compute weighted cell types for each motif
    tidy_df = weighted_celltypes_by_motif(
        cell_loadings=cell_loadings,
        metadata_df=cell_meta_df,
        normalize=normalize
    )

    # Get colors
    unique_celltypes = tidy_df['cell_type'].unique().tolist()
    color_map = _get_colors_for_plotting(
        ct_colors=ct_colors,
        df_celltypes=unique_celltypes
    )

    # Determine save path
    save_path = None
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, 'motif_celltype_composition.png')

    # Plot
    fig, ax = plot_motif_celltype_composition(
        tidy_df,
        color_map=color_map,
        figsize=figsize,
        title=title if title else f"Cell Type Composition per Motif ({mode_str})",
        save_path=save_path
    )

    return fig, ax, tidy_df


def analyze_motif_state_counts(
    adata: Union[anndata.AnnData, Dict[str, anndata.AnnData]],
    sample_column: Optional[str] = None,
    figsize: tuple = (12, 6),
    colors: list = ['#66c2a5', '#fc8d62'],
    title: Optional[str] = None,
    output_dir: Optional[str] = None
) -> Tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """
    Analyze and visualize ON/OFF state counts for each motif.

    Supports three input modes:
    1. Single AnnData: standard single-sample analysis
    2. Dict of AnnData objects: multi-sample with dict
    3. Single merged AnnData with sample_column: multi-sample merged

    Parameters
    ----------
    adata : anndata.AnnData or Dict[str, anndata.AnnData]
        Must contain motif_{k}_state columns in obs (from gmm_binarize_all_motifs)
    sample_column : str, optional
        Column name in adata.obs identifying samples (for merged AnnData).
        If provided with single AnnData, treats it as multi-sample merged.
    figsize : tuple, default (12, 6)
        Figure size
    colors : list, default ['#66c2a5', '#fc8d62']
        Colors for [positive, negative] bars
    title : str, optional
        Plot title. If None, auto-generates.
    output_dir : str, optional
        Directory to save figure. If None, figure is not saved.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    ax : matplotlib.axes.Axes
        Axes object
    counts_df : pd.DataFrame
        Counts dataframe with columns [positive, negative] and motif names as index

    Examples
    --------
    >>> # Single sample
    >>> fig, ax, counts = al.analyze_motif_state_counts(adata)
    >>>
    >>> # Multi-sample with dict
    >>> fig, ax, counts = al.analyze_motif_state_counts(adata_dict)
    >>>
    >>> # Multi-sample with merged AnnData
    >>> fig, ax, counts = al.analyze_motif_state_counts(
    ...     adata_merged, sample_column='patient_id'
    ... )
    """
    from alarmist.core.single_cell import compute_motif_state_counts

    # Determine mode for logging
    if isinstance(adata, dict):
        mode_str = f"Multi-sample ({len(adata)} samples)"
        n_cells = sum(ad.n_obs for ad in adata.values())
    elif sample_column is not None:
        n_samples = adata.obs[sample_column].nunique()
        mode_str = f"Multi-sample ({n_samples} samples, merged)"
        n_cells = adata.n_obs
    else:
        mode_str = "Single sample"
        n_cells = adata.n_obs

    print(f"Computing ON/OFF statistics...")
    print(f"  Mode: {mode_str}")
    print(f"  Cells: {n_cells:,}")

    # Compute counts
    counts_df = compute_motif_state_counts(adata, sample_column=sample_column)

    print(f"  Motifs: {len(counts_df)}")

    # Determine save path
    save_path = None
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, 'motif_state_counts.png')

    # Plot
    fig, ax = plot_motif_state_counts(
        counts_df,
        figsize=figsize,
        colors=colors,
        title=title if title else f"Positive vs Negative Cells per Motif ({mode_str})",
        save_path=save_path
    )

    return fig, ax, counts_df
