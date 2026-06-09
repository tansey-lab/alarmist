"""
Single cell visualization functions for motif analysis
"""

import logging
import os

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from alarmist.constants import (
    COLUMN_NAME_CELL_TYPE,
    COLUMN_NAME_MOTIF,
    COLUMN_NAME_SAMPLE_ID,
    COLUMN_NAME_WEIGHT,
)
from alarmist.plotting.colors import _get_colors_for_plotting

logger = logging.getLogger(__name__)


def plot_motif_celltype_composition(
    df_tidy: pd.DataFrame,
    color_map: dict | None = None,
    figsize: tuple = (12, 6),
    ylabel: str | None = None,
    title: str | None = None,
    save_path: str | None = None,
    cell_type_col: str = COLUMN_NAME_CELL_TYPE,
) -> tuple[plt.Figure, plt.Axes]:
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
        index=COLUMN_NAME_MOTIF,
        columns=cell_type_col,
        values=COLUMN_NAME_WEIGHT,
        aggfunc="sum",
        fill_value=0.0,
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
        logger.debug(f"Saved: {save_path}")

    return fig, ax


def plot_motif_state_counts(
    df_counts: pd.DataFrame,
    figsize: tuple = (12, 6),
    colors: list = ["#66c2a5", "#fc8d62"],
    title: str | None = None,
    save_path: str | None = None,
) -> tuple[plt.Figure, plt.Axes]:
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

    df_counts.plot(kind="bar", stacked=True, ax=ax, color=colors)

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
        logger.debug(f"Saved: {save_path}")

    return fig, ax


def plot_positive_motifs_distribution(
    counts: pd.Series,
    figsize: tuple = (10, 4),
    title: str | None = None,
    save_path: str | None = None,
) -> tuple[plt.Figure, plt.Axes]:
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

    counts.plot(kind="bar", ax=ax)

    ax.set_xlabel("Number of positive motifs per cell")
    ax.set_ylabel("Number of cells")
    if title:
        ax.set_title(title)
    else:
        ax.set_title("Distribution of Positive Motifs per Cell")
    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.debug(f"Saved: {save_path}")

    return fig, ax


def plot_motif_spatial(
    adata: anndata.AnnData | dict[str, anndata.AnnData],
    motif_idx: int | list[int],
    sample_column: str | None = None,
    cell_type_column: str = COLUMN_NAME_CELL_TYPE,
    n_cols: int = 4,
    figsize_per_panel: tuple = (6, 6),
    point_size: float = 0.2,
    color_by_celltype: bool = True,
    ct_colors: dict | None = None,
    positive_color: str = "#1f77b4",
    negative_color: str = "#d3d3d3",
    legend_top_n: int = 15,
    wspace: float = 0.3,
    hspace: float = 0.4,
    output_dir: str | None = None,
    intersect: bool = False,
) -> plt.Figure | list[plt.Figure]:
    """
    Plot spatial distribution of motif ON/OFF states.

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
        - If list and intersect=False: plot each motif separately, return list of Figures
        - If list and intersect=True: plot cells positive for ALL motifs, return single Figure
    sample_column : str, optional
        Column name in adata.obs identifying samples (for merged AnnData).
        If provided with single AnnData, splits into multiple panels.
    cell_type_column : str, default 'cell_type'
        Column name for cell type annotations
    n_cols : int, default 4
        Number of columns in the grid (for multi-sample)
    figsize_per_panel : tuple, default (6, 6)
        Figure size per panel
    point_size : float, default 0.2
        Size of scatter points
    color_by_celltype : bool, default True
        If True, positive cells colored by cell type.
        If False, positive cells colored with positive_color (blue).
    ct_colors : dict, optional
        Cell type color mapping. If None, uses global color registry.
        Only used when color_by_celltype=True.
    positive_color : str, default '#1f77b4'
        Color for positive (ON) cells when color_by_celltype=False
    negative_color : str, default '#d3d3d3'
        Color for negative (OFF) cells
    legend_top_n : int, default 15
        When color_by_celltype=True, show top N cell types in legend
    wspace : float, default 0.3
        Horizontal spacing between panels (fraction of panel width)
    hspace : float, default 0.4
        Vertical spacing between panels (fraction of panel height)
    output_dir : str, optional
        Directory to save figures. Files named motif_{k}_spatial.png
    intersect : bool, default False
        If True and motif_idx is a list, plot cells that are positive for ALL
        specified motifs (intersection). Returns a single Figure.
        If False (default), plot each motif separately, returning a list of Figures.

    Returns
    -------
    plt.Figure or List[plt.Figure]
        - Single Figure if motif_idx is int
        - Single Figure if motif_idx is list and intersect=True
        - List of Figures if motif_idx is list and intersect=False

    Examples
    --------
    >>> # Color by cell type (default)
    >>> fig = al.plot_motif_spatial(adata_dict, motif_idx=5, output_dir='results/spatial')
    >>>
    >>> # Simple blue/gray coloring
    >>> fig = al.plot_motif_spatial(adata_dict, motif_idx=5, color_by_celltype=False)
    >>>
    >>> # Multiple motifs - separate figures
    >>> figs = al.plot_motif_spatial(adata_dict, motif_idx=range(20), output_dir='results/spatial')
    >>>
    >>> # Multiple motifs - intersection (cells positive for ALL motifs)
    >>> fig = al.plot_motif_spatial(adata_dict, motif_idx=[0, 5, 10], intersect=True)
    """
    import matplotlib.lines as mlines

    # Handle list of motifs
    if isinstance(motif_idx, list | tuple | range):
        motif_list = list(motif_idx)

        if not intersect:
            # Default behavior: plot each motif separately
            figs = []
            for k in motif_list:
                fig = plot_motif_spatial(
                    adata=adata,
                    motif_idx=k,
                    sample_column=sample_column,
                    cell_type_column=cell_type_column,
                    n_cols=n_cols,
                    figsize_per_panel=figsize_per_panel,
                    point_size=point_size,
                    color_by_celltype=color_by_celltype,
                    ct_colors=ct_colors,
                    positive_color=positive_color,
                    negative_color=negative_color,
                    legend_top_n=legend_top_n,
                    wspace=wspace,
                    hspace=hspace,
                    output_dir=output_dir,
                    intersect=False,
                )
                figs.append(fig)
                logger.debug(f"Plotted motif {k}")
            return figs

        # intersect=True: plot cells positive for ALL motifs
        # We'll handle this below by computing intersection mask

    # Get color map if coloring by cell type
    ct_color_hex = {}
    if color_by_celltype:
        color_map = _get_colors_for_plotting(ct_colors=ct_colors)
        if not color_map:
            raise ValueError(
                "No cell type colors available. Either:\n"
                "  1. Call al.set_celltype_colors() first, or\n"
                "  2. Pass ct_colors parameter, or\n"
                "  3. Set color_by_celltype=False"
            )

        # Convert colors to hex for matplotlib
        def to_hex(c):
            if isinstance(c, str):
                return c
            elif isinstance(c, tuple) and len(c) >= 3:
                return (
                    f"#{int(c[0] * 255):02x}{int(c[1] * 255):02x}{int(c[2] * 255):02x}"
                )
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
        adata_dict = {"sample": adata}
        mode_str = "single"

    n_samples = len(adata_dict)

    # Determine if we're in intersect mode (list of motifs with intersect=True)
    intersect_mode = intersect and isinstance(motif_idx, list | tuple | range)
    if intersect_mode:
        motif_list = list(motif_idx)
        state_cols = [f"motif_{k}_state" for k in motif_list]
        motif_label = f"Motifs {','.join(map(str, motif_list))} (intersection)"
        motif_label_short = f"motifs_{'_'.join(map(str, motif_list))}_intersect"
    else:
        state_cols = [f"motif_{motif_idx}_state"]
        motif_label = f"Motif {motif_idx}"
        motif_label_short = f"motif_{motif_idx}"

    # Determine grid layout
    if n_samples == 1:
        n_rows, n_cols_actual = 1, 1
    else:
        n_cols_actual = min(n_cols, n_samples)
        n_rows = (n_samples + n_cols_actual - 1) // n_cols_actual

    # Add extra width for legend
    fig_width = figsize_per_panel[0] * n_cols_actual + 2.5
    fig_height = figsize_per_panel[1] * n_rows

    fig, axes = plt.subplots(
        n_rows, n_cols_actual, figsize=(fig_width, fig_height), squeeze=False
    )
    axes = axes.flatten()

    # Collect cell type counts across all samples (for legend)
    all_ct_counts = {}
    total_pos = 0
    total_neg = 0

    for i, (sample_id, ad) in enumerate(adata_dict.items()):
        ax = axes[i]
        coords = ad.obsm["spatial"][:, :2]

        # Check if all required state columns exist
        missing_cols = [col for col in state_cols if col not in ad.obs.columns]
        if missing_cols:
            if intersect_mode:
                ax.set_title(f"{sample_id}\n(missing motif data)")
            else:
                ax.set_title(f"{sample_id}\n(no motif {motif_idx} data)")
            ax.axis("off")
            continue

        # Compute positive mask (intersection if multiple motifs)
        if intersect_mode:
            # Intersection: positive for ALL motifs
            mask_pos = np.ones(len(ad), dtype=bool)
            for col in state_cols:
                mask_pos &= (ad.obs[col].astype(str) == "positive").values
        else:
            mask_pos = (ad.obs[state_cols[0]].astype(str) == "positive").values

        mask_neg = ~mask_pos
        n_pos = mask_pos.sum()
        n_neg = mask_neg.sum()
        total_pos += n_pos
        total_neg += n_neg
        frac_pos = n_pos / len(ad) * 100 if len(ad) > 0 else 0

        # Plot negatives (background, gray)
        ax.scatter(
            coords[mask_neg, 0],
            coords[mask_neg, 1],
            c=negative_color,
            s=point_size,
            alpha=0.95,
            marker="o",
            edgecolors="none",
            linewidths=0,
            rasterized=True,
        )

        # Plot positives
        if mask_pos.any():
            if color_by_celltype:
                ct = ad.obs[cell_type_column].astype(str).values
                pos_colors = np.array([ct_color_hex.get(x, "#000000") for x in ct])[
                    mask_pos
                ]

                # Count cell types for legend
                for ct_name in ct[mask_pos]:
                    all_ct_counts[ct_name] = all_ct_counts.get(ct_name, 0) + 1
            else:
                pos_colors = positive_color

            ax.scatter(
                coords[mask_pos, 0],
                coords[mask_pos, 1],
                c=pos_colors,
                s=point_size * 6,
                alpha=1,
                marker="o",
                edgecolors="none",
                linewidths=0,
                rasterized=True,
                zorder=3,
            )

        # Title
        if mode_str == "single":
            ax.set_title(f"{motif_label}: {frac_pos:.1f}% ON")
        else:
            ax.set_title(f"{sample_id}\n{frac_pos:.1f}% ON")

        ax.set_aspect("equal")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    # Create figure-level legend on the right
    total_cells = total_pos + total_neg
    total_frac = total_pos / total_cells * 100 if total_cells > 0 else 0

    handles = []

    # ON/OFF handles
    off_handle = mlines.Line2D(
        [],
        [],
        color=negative_color,
        marker="o",
        linestyle="None",
        markersize=8,
        label=f"OFF (n={total_neg:,})",
    )
    handles.append(off_handle)

    if color_by_celltype:
        on_handle = mlines.Line2D(
            [],
            [],
            color="black",
            marker="o",
            linestyle="None",
            markersize=8,
            label=f"ON (n={total_pos:,}, {total_frac:.1f}%)",
        )
        handles.append(on_handle)

        # Add separator
        handles.append(mlines.Line2D([], [], color="none", label=""))

        # Cell type handles (top N)
        sorted_cts = sorted(all_ct_counts.items(), key=lambda x: -x[1])
        for ct_name, ct_count in sorted_cts[:legend_top_n]:
            ct_handle = mlines.Line2D(
                [],
                [],
                color=ct_color_hex.get(ct_name, "#000000"),
                marker="o",
                linestyle="None",
                markersize=8,
                label=ct_name,
            )
            handles.append(ct_handle)
    else:
        on_handle = mlines.Line2D(
            [],
            [],
            color=positive_color,
            marker="o",
            linestyle="None",
            markersize=8,
            label=f"ON (n={total_pos:,}, {total_frac:.1f}%)",
        )
        handles.append(on_handle)

    # Place legend outside on the right
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        fontsize=10,
        frameon=False,
        title="Cell type" if color_by_celltype else "State",
    )

    # Suptitle
    if color_by_celltype:
        suptitle = f"{motif_label} Spatial Distribution (ON colored by cell type)"
    else:
        suptitle = f"{motif_label} Spatial Distribution"

    if mode_str != "single":
        plt.suptitle(suptitle, fontsize=14)

    plt.subplots_adjust(wspace=wspace, hspace=hspace, right=0.85)
    plt.tight_layout()

    # Save
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"{motif_label_short}_spatial.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.debug(f"Saved: {save_path}")

    return fig


def analyze_motif_celltype_composition(
    adata: anndata.AnnData | dict[str, anndata.AnnData],
    cell_loadings: np.ndarray,
    cell_type_column: str = COLUMN_NAME_CELL_TYPE,
    sample_column: str | None = None,
    normalize: bool = True,
    ct_colors: dict | None = None,
    figsize: tuple = (12, 6),
    title: str | None = None,
    output_dir: str | None = None,
    motif_ids: list[int] | None = None,
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
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
    motif_ids : list of int, optional
        List of motif indices to include in the plot. If None, all motifs
        are plotted. Example: motif_ids=[0, 1, 5] to plot only motifs 0, 1, 5.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    ax : matplotlib.axes.Axes
        Axes object
    tidy_df : pd.DataFrame
        Tidy dataframe with columns [motif, cell_type, weight].
        If motif_ids is specified, only contains the selected motifs.

    Examples
    --------
    >>> # Single sample
    >>> fig, ax, df = al.analyze_motif_celltype_composition(
    ...     adata, cell_loadings,
    ...     output_dir='results/single_cell'
    ... )
    >>>
    >>> # Plot only specific motifs
    >>> fig, ax, df = al.analyze_motif_celltype_composition(
    ...     adata, cell_loadings,
    ...     motif_ids=[0, 1, 5],  # only plot motifs 0, 1, and 5
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
            df.columns = [COLUMN_NAME_CELL_TYPE]
            df[COLUMN_NAME_SAMPLE_ID] = sample_id
            cell_meta_dfs.append(df)
        cell_meta_df = pd.concat(cell_meta_dfs, axis=0).reset_index(drop=True)
        mode_str = f"Multi-sample ({len(adata)} samples)"
    elif sample_column is not None:
        # Mode 3: Merged AnnData with sample_column
        cell_meta_df = adata.obs[[cell_type_column, sample_column]].copy()
        cell_meta_df.columns = [COLUMN_NAME_CELL_TYPE, COLUMN_NAME_SAMPLE_ID]
        cell_meta_df = cell_meta_df.reset_index(drop=True)
        n_samples = cell_meta_df[COLUMN_NAME_SAMPLE_ID].nunique()
        mode_str = f"Multi-sample ({n_samples} samples, merged)"
    else:
        # Mode 1: Single AnnData
        cell_meta_df = adata.obs[[cell_type_column]].copy()
        cell_meta_df.columns = [COLUMN_NAME_CELL_TYPE]
        cell_meta_df = cell_meta_df.reset_index(drop=True)
        mode_str = "Single sample"

    # Validate dimensions
    if len(cell_meta_df) != cell_loadings.shape[0]:
        raise ValueError(
            f"Number of cells in adata ({len(cell_meta_df)}) does not match "
            f"cell_loadings rows ({cell_loadings.shape[0]})"
        )

    logger.debug("Analyzing cell type composition per motif...")
    logger.debug(f"  Mode: {mode_str}")
    logger.debug(f"  Cells: {len(cell_meta_df):,}")
    n_motifs_total = cell_loadings.shape[1]
    if motif_ids is not None:
        logger.debug(f"  Motifs: {len(motif_ids)} selected (of {n_motifs_total} total)")
    else:
        logger.debug(f"  Motifs: {n_motifs_total}")

    # Compute weighted cell types for each motif
    tidy_df = weighted_celltypes_by_motif(
        cell_loadings=cell_loadings, metadata_df=cell_meta_df, normalize=normalize
    )

    # Filter to selected motifs if specified
    if motif_ids is not None:
        tidy_df = tidy_df[tidy_df[COLUMN_NAME_MOTIF].isin(motif_ids)].copy()
        if tidy_df.empty:
            raise ValueError(
                f"No motifs found matching motif_ids={motif_ids}. "
                f"Available motifs: 0-{cell_loadings.shape[1] - 1}"
            )
        logger.debug(f"  Selected motifs: {sorted(motif_ids)}")

    # Get colors
    unique_celltypes = tidy_df[COLUMN_NAME_CELL_TYPE].unique().tolist()
    color_map = _get_colors_for_plotting(
        ct_colors=ct_colors, df_celltypes=unique_celltypes
    )

    # Determine save path
    save_path = None
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, "motif_celltype_weighted.png")

    # Plot
    fig, ax = plot_motif_celltype_composition(
        tidy_df,
        color_map=color_map,
        figsize=figsize,
        title=title if title else "Cell Type Composition (Weighted)",
        save_path=save_path,
    )

    return fig, ax, tidy_df


def analyze_motif_celltype_counts(
    adata: anndata.AnnData | dict[str, anndata.AnnData],
    cell_type_column: str = COLUMN_NAME_CELL_TYPE,
    sample_column: str | None = None,
    normalize: bool = True,
    ct_colors: dict | None = None,
    figsize: tuple = (12, 6),
    title: str | None = None,
    output_dir: str | None = None,
    motif_ids: list[int] | None = None,
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """
    Analyze and visualize cell type composition based on positive cell counts.

    Unlike analyze_motif_celltype_composition which uses cell loadings as weights,
    this function counts the number of positive (ON) cells for each cell type
    per motif based on the binarized motif states from gmm_binarize_all_motifs.

    Supports three input modes:
    1. Single AnnData: standard single-sample analysis
    2. Dict of AnnData objects: multi-sample with dict
    3. Single merged AnnData with sample_column: multi-sample merged

    Parameters
    ----------
    adata : anndata.AnnData or Dict[str, anndata.AnnData]
        Must contain motif_{k}_state columns in obs (from gmm_binarize_all_motifs).
    cell_type_column : str, default 'cell_type'
        Column name in adata.obs containing cell type annotations
    sample_column : str, optional
        Column name in adata.obs identifying samples (for merged AnnData).
        If provided with single AnnData, treats it as multi-sample.
    normalize : bool, default True
        If True, per-motif counts are normalized to sum to 1 (fraction).
        If False, shows raw counts.
    ct_colors : dict, optional
        Mapping from cell type to color. If None, uses global color registry
        or auto-generates from tab20 colormap.
    figsize : tuple, default (12, 6)
        Figure size
    title : str, optional
        Plot title. If None, auto-generates based on input mode.
    output_dir : str, optional
        Directory to save figure. If None, figure is not saved.
    motif_ids : list of int, optional
        List of motif indices to include in the plot. If None, all motifs
        with state columns are plotted.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    ax : matplotlib.axes.Axes
        Axes object
    tidy_df : pd.DataFrame
        Tidy dataframe with columns [motif, cell_type, weight].
        weight represents count (or fraction if normalize=True).

    Examples
    --------
    >>> # After running gmm_binarize_all_motifs
    >>> fig, ax, df = al.analyze_motif_celltype_counts(
    ...     adata,
    ...     cell_type_column='cell_type',
    ...     output_dir='results/single_cell'
    ... )
    >>>
    >>> # Plot only specific motifs
    >>> fig, ax, df = al.analyze_motif_celltype_counts(
    ...     adata,
    ...     motif_ids=[0, 1, 5],
    ...     output_dir='results/single_cell'
    ... )
    >>>
    >>> # Show raw counts instead of fractions
    >>> fig, ax, df = al.analyze_motif_celltype_counts(
    ...     adata,
    ...     normalize=False
    ... )
    """
    # Combine adata if dict
    if isinstance(adata, dict):
        obs_list = []
        for sample_id, ad in adata.items():
            obs_df = ad.obs.copy()
            obs_df["_sample_id"] = sample_id
            obs_list.append(obs_df)
        combined_obs = pd.concat(obs_list, axis=0)
        mode_str = f"Multi-sample ({len(adata)} samples)"
        n_cells = len(combined_obs)
    elif sample_column is not None:
        combined_obs = adata.obs.copy()
        n_samples = combined_obs[sample_column].nunique()
        mode_str = f"Multi-sample ({n_samples} samples, merged)"
        n_cells = len(combined_obs)
    else:
        combined_obs = adata.obs.copy()
        mode_str = "Single sample"
        n_cells = len(combined_obs)

    # Find all motif state columns
    state_cols = [
        c
        for c in combined_obs.columns
        if c.startswith("motif_") and c.endswith("_state")
    ]
    if not state_cols:
        raise ValueError(
            "No motif state columns found in adata.obs. "
            "Run gmm_binarize_all_motifs() first."
        )

    # Extract motif indices from column names
    motif_indices = []
    for col in state_cols:
        # Extract number from 'motif_X_state'
        parts = col.replace("motif_", "").replace("_state", "")
        try:
            motif_indices.append(int(parts))
        except ValueError:
            continue

    motif_indices = sorted(motif_indices)
    n_motifs_total = len(motif_indices)

    # Filter to selected motifs if specified
    if motif_ids is not None:
        motif_indices = [m for m in motif_indices if m in motif_ids]
        if not motif_indices:
            raise ValueError(
                f"No motifs found matching motif_ids={motif_ids}. "
                f"Available motifs: {sorted([int(c.replace('motif_', '').replace('_state', '')) for c in state_cols])}"
            )

    logger.debug("Analyzing cell type counts for positive cells...")
    logger.debug(f"  Mode: {mode_str}")
    logger.debug(f"  Cells: {n_cells:,}")
    if motif_ids is not None:
        logger.debug(
            f"  Motifs: {len(motif_indices)} selected (of {n_motifs_total} total)"
        )
    else:
        logger.debug(f"  Motifs: {n_motifs_total}")

    # Build tidy dataframe: count positive cells by cell type for each motif
    records = []
    cell_types = combined_obs[cell_type_column].astype(str).values

    for k in motif_indices:
        state_col = f"motif_{k}_state"
        if state_col not in combined_obs.columns:
            continue

        states = combined_obs[state_col].astype(str).values
        mask_pos = states == "positive"

        # Count cell types among positive cells
        pos_celltypes = cell_types[mask_pos]
        ct_counts = pd.Series(pos_celltypes).value_counts()

        for ct, count in ct_counts.items():
            records.append(
                {
                    COLUMN_NAME_MOTIF: k,
                    COLUMN_NAME_CELL_TYPE: ct,
                    COLUMN_NAME_WEIGHT: count,
                }
            )

    if not records:
        raise ValueError("No positive cells found for any motif.")

    tidy_df = pd.DataFrame(records)

    # Normalize per motif if requested
    if normalize:
        totals = tidy_df.groupby(COLUMN_NAME_MOTIF)[COLUMN_NAME_WEIGHT].transform("sum")
        tidy_df[COLUMN_NAME_WEIGHT] = tidy_df[COLUMN_NAME_WEIGHT] / totals

    if motif_ids is not None:
        logger.debug(f"  Selected motifs: {sorted(motif_indices)}")

    # Get colors
    unique_celltypes = tidy_df[COLUMN_NAME_CELL_TYPE].unique().tolist()
    color_map = _get_colors_for_plotting(
        ct_colors=ct_colors, df_celltypes=unique_celltypes
    )

    # Determine save path
    save_path = None
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, "motif_celltype_counts.png")

    # Determine ylabel based on normalization
    ylabel = "Fraction of positive cells" if normalize else "Number of positive cells"

    # Plot
    fig, ax = plot_motif_celltype_composition(
        tidy_df,
        color_map=color_map,
        figsize=figsize,
        ylabel=ylabel,
        title=title if title else "Cell Type Composition (Positive Counts)",
        save_path=save_path,
    )

    return fig, ax, tidy_df


def analyze_motif_state_counts(
    adata: anndata.AnnData | dict[str, anndata.AnnData],
    sample_column: str | None = None,
    figsize: tuple = (12, 6),
    colors: list = ["#66c2a5", "#fc8d62"],
    title: str | None = None,
    output_dir: str | None = None,
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
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

    logger.debug("Computing ON/OFF statistics...")
    logger.debug(f"  Mode: {mode_str}")
    logger.debug(f"  Cells: {n_cells:,}")

    # Compute counts
    counts_df = compute_motif_state_counts(adata, sample_column=sample_column)

    logger.debug(f"  Motifs: {len(counts_df)}")

    # Determine save path
    save_path = None
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, "motif_state_counts.png")

    # Plot
    fig, ax = plot_motif_state_counts(
        counts_df,
        figsize=figsize,
        colors=colors,
        title=title if title else f"Positive vs Negative Cells per Motif ({mode_str})",
        save_path=save_path,
    )

    return fig, ax, counts_df
