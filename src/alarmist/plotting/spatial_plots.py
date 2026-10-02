"""
Spatial distribution plotting functions
"""

import logging

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from alarmist.constants import (
    COLUMN_NAME_CELL_TYPE,
    COLUMN_NAME_N_CELLS,
    COLUMN_NAME_PATCH_IDX,
)
from alarmist.plotting.colors import _get_colors_for_plotting

logger = logging.getLogger(__name__)


def plot_cells_per_patch(
    adata: anndata.AnnData | dict[str, anndata.AnnData] | pd.DataFrame,
    bins: int = 50,
    figsize: tuple = (8, 5),
    title: str | None = None,
    color: str = "#4C72B0",
    save_path: str | None = None,
    show: bool = True,
    count_col: str = COLUMN_NAME_N_CELLS,
) -> plt.Figure:
    """
    Plot histogram of cell counts per patch.

    Parameters
    ----------
    adata : AnnData or Dict[str, AnnData]
        AnnData object(s) with 'patch_idx' column in obs (from run_patchify).
        If dict, all samples are combined. A patch-metadata DataFrame (one row
        per patch, e.g. ``results['patch_metadata_df']``) is also accepted; the
        per-patch counts are then read from ``count_col``.
    bins : int, default 50
        Number of histogram bins
    figsize : tuple, default (8, 5)
        Figure size
    title : str, optional
        Plot title. If None, uses "Distribution of Cells per Patch"
    color : str, default '#4C72B0'
        Histogram color
    save_path : str, optional
        Path to save figure. If None, figure is not saved.
    show : bool, default True
        Whether to display the figure
    count_col : str, default 'n_cells'
        Column holding the per-patch cell count. Only used for DataFrame input.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object

    Examples
    --------
    >>> results = analyzer.run_patchify(adata, output_dir=results_dir)
    >>> fig = al.plot_cells_per_patch(adata, save_path='cells_per_patch.png')
    """
    if isinstance(adata, pd.DataFrame):
        if count_col not in adata.columns:
            raise ValueError(f"'{count_col}' not found in patch metadata columns.")
        counts = adata[count_col].to_numpy()
        n_patches, n_cells = len(counts), int(counts.sum())
    elif isinstance(adata, dict):
        # Multi-sample: combine all samples
        all_patch_idx = []
        for sample_id, ad in adata.items():
            if COLUMN_NAME_PATCH_IDX not in ad.obs.columns:
                raise ValueError(
                    f"'patch_idx' not found in adata.obs for sample '{sample_id}'. "
                    "Run run_patchify first."
                )
            all_patch_idx.extend(ad.obs[COLUMN_NAME_PATCH_IDX].values)
        patch_idx = np.array(all_patch_idx)
    else:
        if COLUMN_NAME_PATCH_IDX not in adata.obs.columns:
            raise ValueError(
                "'patch_idx' not found in adata.obs. Run run_patchify first."
            )
        patch_idx = adata.obs[COLUMN_NAME_PATCH_IDX].values

    if not isinstance(adata, pd.DataFrame):
        # Filter out invalid patch indices (-1), then count cells per patch
        valid_patch_idx = patch_idx[patch_idx >= 0]
        unique_patches, counts = np.unique(valid_patch_idx, return_counts=True)
        n_patches, n_cells = len(unique_patches), len(valid_patch_idx)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram
    ax.hist(counts, bins=bins, color=color, edgecolor="white", alpha=0.8)

    # Labels and title
    ax.set_xlabel("Number of Cells per Patch", fontsize=12)
    ax.set_ylabel("Number of Patches", fontsize=12)
    ax.set_title(title or "Distribution of Cells per Patch", fontsize=14)

    # Add statistics as text
    stats_text = (
        f"Total patches: {n_patches}\n"
        f"Total cells: {n_cells}\n"
        f"Mean: {counts.mean():.1f}\n"
        f"Median: {np.median(counts):.1f}\n"
        f"Min: {counts.min()}, Max: {counts.max()}"
    )
    ax.text(
        0.97,
        0.97,
        stats_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()

    # Save if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.debug(f"Saved: {save_path}")

    if show:
        plt.show()
    else:
        plt.close()

    return fig


def _as_sample_dict(
    adata: anndata.AnnData | dict[str, anndata.AnnData],
    sample_column: str | None = None,
) -> dict[str, anndata.AnnData]:
    """Normalise single / dict / merged+sample_column input to {sample: AnnData}."""
    if isinstance(adata, dict):
        return adata
    if sample_column is not None:
        return {
            sid: adata[(adata.obs[sample_column] == sid).values]
            for sid in adata.obs[sample_column].unique()
        }
    return {"sample": adata}


def _panel_grid(
    n_panels: int,
    n_cols: int,
    figsize_per_panel: tuple[float, float],
    extra_width: float = 0.0,
) -> tuple[plt.Figure, np.ndarray]:
    """Create a grid with one axes per panel and hide the unused ones."""
    n_cols = max(1, min(n_cols, n_panels))
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(
            figsize_per_panel[0] * n_cols + extra_width,
            figsize_per_panel[1] * n_rows,
        ),
        squeeze=False,
    )
    axes = axes.flatten()
    for ax in axes[n_panels:]:
        ax.axis("off")
    return fig, axes


def _save(fig: plt.Figure, save_path: str | None, dpi: int = 300) -> None:
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        logger.debug(f"Saved: {save_path}")


def plot_motif_score_spatial(
    adata: anndata.AnnData | dict[str, anndata.AnnData],
    motif_idx: int,
    sample_column: str | None = None,
    value_key: str | None = None,
    log: bool = False,
    spatial_key: str = "spatial",
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    clip_percentiles: tuple[float, float] = (1, 99),
    n_cols: int = 2,
    figsize_per_panel: tuple[float, float] = (6, 6),
    point_size: float = 0.5,
    alpha: float = 0.8,
    invert_y: bool = False,
    title: str | None = None,
    save_path: str | None = None,
) -> plt.Figure:
    """
    Plot a continuous per-cell motif value in space, one panel per sample.

    All panels share one colour scale (and one colourbar), so samples are
    directly comparable. Unless ``vmin``/``vmax`` are given, the scale is
    clipped to percentiles of the values pooled across samples.

    Parameters
    ----------
    adata : AnnData or Dict[str, AnnData]
        Single AnnData, dict of AnnData, or merged AnnData with ``sample_column``.
    motif_idx : int
        Motif to plot.
    sample_column : str, optional
        obs column splitting a merged AnnData into panels.
    value_key : str, optional
        obs column holding the value. Defaults to ``motif_{k}_loading``
        (written by gmm_binarize_all_motifs).
    log : bool, default False
        Plot ``log10(value + eps)`` instead of the raw value, with eps the smallest
        positive value.
    spatial_key : str, default 'spatial'
        obsm key with coordinates.
    cmap : str, default 'viridis'
        Colormap.
    vmin, vmax : float, optional
        Colour-scale limits. Default: ``clip_percentiles`` of the pooled values.
    clip_percentiles : tuple, default (1, 99)
        Percentiles used for the default colour-scale limits.
    n_cols : int, default 2
        Panels per row.
    figsize_per_panel : tuple, default (6, 6)
        Size of each panel.
    point_size : float, default 0.5
        Scatter marker size.
    alpha : float, default 0.8
        Marker alpha.
    invert_y : bool, default False
        Flip the y axis.
    title : str, optional
        Suptitle. Default: ``"Motif k <value>"``.
    save_path : str, optional
        File to save to.

    Returns
    -------
    matplotlib.figure.Figure
    """
    col = value_key or f"motif_{motif_idx}_loading"
    adata_dict = _as_sample_dict(adata, sample_column)

    pooled = [
        ad.obs[col].dropna().to_numpy(dtype=float)
        for ad in adata_dict.values()
        if col in ad.obs.columns
    ]
    if not pooled:
        raise ValueError(f"'{col}' not found in adata.obs of any sample")
    pooled = np.concatenate(pooled)

    eps = 0.0
    if log:
        positive = pooled[pooled > 0]
        eps = float(positive.min()) if positive.size else 1e-10

    def _transform(v):
        return np.log10(v + eps) if log else v

    pooled_t = _transform(pooled)
    if vmin is None:
        vmin = float(np.percentile(pooled_t, clip_percentiles[0]))
    if vmax is None:
        vmax = float(np.percentile(pooled_t, clip_percentiles[1]))

    fig, axes = _panel_grid(len(adata_dict), n_cols, figsize_per_panel, extra_width=1.0)
    mappable = None
    for ax, (sample_id, ad) in zip(axes, adata_dict.items(), strict=False):
        if col not in ad.obs.columns:
            ax.set_title(f"{sample_id}\n(no '{col}')")
            ax.axis("off")
            continue
        coords = np.asarray(ad.obsm[spatial_key])[:, :2]
        values = _transform(ad.obs[col].to_numpy(dtype=float))
        mappable = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=values,
            s=point_size,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            alpha=alpha,
            edgecolors="none",
            rasterized=True,
        )
        ax.set_title(str(sample_id) if len(adata_dict) > 1 else "")
        ax.set_aspect("equal")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        if invert_y:
            ax.invert_yaxis()

    label = f"log10({col})" if log else col
    fig.suptitle(title or f"Motif {motif_idx}: {label}", fontsize=14)
    fig.tight_layout(rect=(0, 0, 0.9, 1))
    if mappable is not None:
        cbar_ax = fig.add_axes((0.92, 0.15, 0.02, 0.7))
        fig.colorbar(mappable, cax=cbar_ax).set_label(label, fontsize=10)

    _save(fig, save_path)
    return fig


def plot_two_motif_spatial(
    adata: anndata.AnnData | dict[str, anndata.AnnData],
    motif_a: int,
    motif_b: int,
    cell_type: str | None = None,
    sample_column: str | None = None,
    cell_type_column: str = COLUMN_NAME_CELL_TYPE,
    spatial_key: str = "spatial",
    a_color: str = "#e74c3c",
    b_color: str = "#3498db",
    on_alpha: float = 0.5,
    background_color: str = "#d3d3d3",
    other_alpha: float = 0.03,
    target_off_alpha: float = 0.08,
    n_cols: int = 2,
    figsize_per_panel: tuple[float, float] = (6, 6),
    point_size: float = 0.5,
    on_size_mult: float = 6.0,
    invert_y: bool = False,
    title: str | None = None,
    save_path: str | None = None,
) -> plt.Figure:
    """
    Show where two motifs are ON, optionally within one cell type.

    Motif A ON cells are drawn in ``a_color`` and motif B ON cells on top in
    ``b_color``, both semi-transparent, so cells ON for both blend into a mixed
    colour (red + blue → purple by default). Panel titles report A-only, B-only
    and A&B counts.

    Parameters
    ----------
    adata : AnnData or Dict[str, AnnData]
        Single AnnData, dict of AnnData, or merged AnnData with ``sample_column``.
        Needs ``motif_{a}_state`` and ``motif_{b}_state`` in obs.
    motif_a, motif_b : int
        The two motifs.
    cell_type : str, optional
        Restrict the comparison to this cell type; other cells form a faint
        background. None compares across all cells.
    sample_column : str, optional
        obs column splitting a merged AnnData into panels.
    cell_type_column : str, default 'cell_type'
        obs column with cell types (only used when ``cell_type`` is set).
    spatial_key : str, default 'spatial'
        obsm key with coordinates.
    a_color, b_color : str
        Colours for motif A and motif B ON cells.
    on_alpha : float, default 0.5
        Alpha of ON cells; must be < 1 for overlaps to blend.
    background_color : str, default '#d3d3d3'
        Colour of cells that are not ON.
    other_alpha : float, default 0.03
        Alpha of cells outside ``cell_type``.
    target_off_alpha : float, default 0.08
        Alpha of cells in ``cell_type`` (or all cells) that are ON for neither motif.
    n_cols : int, default 2
        Panels per row.
    figsize_per_panel : tuple, default (6, 6)
        Size of each panel.
    point_size : float, default 0.5
        Marker size of background cells.
    on_size_mult : float, default 6.0
        ON cells are drawn ``on_size_mult`` times larger.
    invert_y : bool, default False
        Flip the y axis.
    title : str, optional
        Suptitle.
    save_path : str, optional
        File to save to.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.lines as mlines

    state_a, state_b = f"motif_{motif_a}_state", f"motif_{motif_b}_state"
    adata_dict = _as_sample_dict(adata, sample_column)
    target_label = cell_type if cell_type is not None else "Cells"

    fig, axes = _panel_grid(len(adata_dict), n_cols, figsize_per_panel)
    for ax, (sample_id, ad) in zip(axes, adata_dict.items(), strict=False):
        missing = [c for c in (state_a, state_b) if c not in ad.obs.columns]
        if missing:
            ax.set_title(f"{sample_id}\n(missing {', '.join(missing)})")
            ax.axis("off")
            continue

        coords = np.asarray(ad.obsm[spatial_key])[:, :2]
        if cell_type is None:
            is_target = np.ones(ad.n_obs, dtype=bool)
        else:
            is_target = (ad.obs[cell_type_column].astype(str) == str(cell_type)).values
        a_on = is_target & (ad.obs[state_a].astype(str) == "positive").values
        b_on = is_target & (ad.obs[state_b].astype(str) == "positive").values
        neither = is_target & ~a_on & ~b_on

        layers = [
            (~is_target, background_color, other_alpha, point_size),
            (neither, background_color, target_off_alpha, point_size * 1.2),
            (a_on, a_color, on_alpha, point_size * on_size_mult),
            (b_on, b_color, on_alpha, point_size * on_size_mult),
        ]
        for z, (mask, color, a, size) in enumerate(layers, start=1):
            if mask.any():
                ax.scatter(
                    coords[mask, 0],
                    coords[mask, 1],
                    c=color,
                    s=size,
                    alpha=a,
                    edgecolors="none",
                    linewidths=0,
                    rasterized=True,
                    zorder=z,
                )

        n_ab = int((a_on & b_on).sum())
        ax.set_title(
            f"{sample_id}\n"
            f"A only={int((a_on & ~b_on).sum()):,}, B only={int((b_on & ~a_on).sum()):,}, "
            f"A&B={n_ab:,} (n={int(is_target.sum()):,})",
            fontsize=10,
        )
        ax.set_aspect("equal")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        if invert_y:
            ax.invert_yaxis()

    handles = [
        mlines.Line2D(
            [], [], color=background_color, marker="o", ls="None", label="Other cells"
        ),
        mlines.Line2D(
            [], [], color=a_color, marker="o", ls="None", label=f"Motif {motif_a} ON"
        ),
        mlines.Line2D(
            [], [], color=b_color, marker="o", ls="None", label=f"Motif {motif_b} ON"
        ),
    ]
    if cell_type is not None:
        handles.insert(
            1,
            mlines.Line2D(
                [],
                [],
                color=background_color,
                marker="o",
                ls="None",
                label=f"{cell_type} (neither ON)",
            ),
        )
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
        fontsize=9,
    )
    fig.suptitle(
        title or f"{target_label}: Motif {motif_a} vs Motif {motif_b}", fontsize=13
    )
    fig.tight_layout()

    _save(fig, save_path)
    return fig


def plot_spatial_categorical(
    adata: anndata.AnnData | dict[str, anndata.AnnData],
    color_col: str = COLUMN_NAME_CELL_TYPE,
    split_by: str | None = None,
    spatial_key: str = "spatial",
    x_col: str | None = None,
    y_col: str | None = None,
    colors: dict | None = None,
    categories: list[str] | None = None,
    other_color: str = "lightgray",
    n_cols: int = 4,
    figsize_per_panel: tuple[float, float] = (6, 6),
    point_size: float = 0.5,
    alpha: float = 1.0,
    invert_y: bool = False,
    show_axes: bool = True,
    shared_extent: bool = True,
    grid_size: float | None = None,
    grid_origin: str | tuple[float, float] = "min",
    grid_color: str = "k",
    grid_linewidth: float = 0.4,
    grid_alpha: float = 0.25,
    grid_max_lines: int = 5000,
    legend_title: str | None = None,
    title: str | None = None,
    save_path: str | None = None,
) -> plt.Figure:
    """
    Plot a categorical obs column (e.g. cell type) in space, optionally split
    into one panel per sample / core / group, with an optional patch-grid overlay.

    Parameters
    ----------
    adata : AnnData or Dict[str, AnnData]
        AnnData (split into panels by ``split_by`` if given) or a dict of AnnData
        (one panel per entry).
    color_col : str, default 'cell_type'
        Categorical obs column used for colours.
    split_by : str, optional
        obs column giving one panel per value (e.g. 'sample_id', 'tma_id').
    spatial_key : str, default 'spatial'
        obsm key with coordinates. Ignored when ``x_col``/``y_col`` are given.
    x_col, y_col : str, optional
        Read coordinates from these obs columns instead of obsm (e.g. CosMx
        'CenterX_global_px' / 'CenterY_global_px').
    colors : dict, optional
        Category → colour. If None: ``adata.uns[f"{color_col}_colors"]`` when
        present (scanpy convention), else the global cell-type registry, else tab20.
    categories : list of str, optional
        Categories to colour and list in the legend; all others are drawn in
        ``other_color``. Default: all categories.
    other_color : str, default 'lightgray'
        Colour for categories not in ``categories``.
    n_cols : int, default 4
        Panels per row.
    figsize_per_panel : tuple, default (6, 6)
        Size of each panel.
    point_size : float, default 0.5
        Scatter marker size.
    alpha : float, default 1.0
        Marker alpha.
    invert_y : bool, default False
        Flip the y axis (image coordinates).
    show_axes : bool, default True
        Draw axis ticks and labels. False gives clean, tick-free panels.
    shared_extent : bool, default True
        Use the same x/y limits (and grid) in every panel. Only applies to
        ``split_by`` on a single AnnData; dict input always uses per-panel limits.
    grid_size : float, optional
        Draw a square grid with this spacing (same units as the coordinates),
        e.g. the ALARMIST patch size, to show patch boundaries.
    grid_origin : {'min', 'zero'} or (x0, y0), default 'min'
        Grid alignment: snapped to the data minimum, to (0, 0), or an explicit origin.
    grid_color, grid_linewidth, grid_alpha : grid line style.
    grid_max_lines : int, default 5000
        Skip drawing the grid if it would need more lines than this.
    legend_title : str, optional
        Legend title. Defaults to ``color_col``.
    title : str, optional
        Suptitle.
    save_path : str, optional
        File to save to.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(adata, dict):
        panels = dict(adata)
        shared_extent = False
    elif split_by is not None:
        panels = _as_sample_dict(adata, split_by)
        if hasattr(adata.obs[split_by], "cat"):
            order = [c for c in adata.obs[split_by].cat.categories if c in panels]
        else:
            order = sorted(panels, key=str)
        panels = {k: panels[k] for k in order}
    else:
        panels = {"": adata}

    def _coords(ad):
        if x_col is not None and y_col is not None:
            return ad.obs[[x_col, y_col]].to_numpy(dtype=float)
        return np.asarray(ad.obsm[spatial_key])[:, :2].astype(float)

    # Category order: categorical order if available, else sorted
    first = next(iter(panels.values()))
    if hasattr(first.obs[color_col], "cat"):
        all_cats = [str(c) for c in first.obs[color_col].cat.categories]
    else:
        all_cats = []
    seen = set(all_cats)
    for ad in panels.values():
        for c in sorted(ad.obs[color_col].astype(str).unique()):
            if c not in seen:
                all_cats.append(c)
                seen.add(c)
    shown = [str(c) for c in categories] if categories is not None else all_cats

    if colors is None:
        uns_key = f"{color_col}_colors"
        if (
            not isinstance(adata, dict)
            and uns_key in adata.uns
            and hasattr(adata.obs[color_col], "cat")
        ):
            colors = dict(
                zip(
                    adata.obs[color_col].cat.categories.astype(str),
                    adata.uns[uns_key],
                    strict=False,
                )
            )
        else:
            colors = _get_colors_for_plotting(None, all_cats)
    color_lookup = {c: colors.get(c, other_color) for c in shown}

    extent = None
    if shared_extent:
        xy = np.concatenate([_coords(ad) for ad in panels.values()])
        extent = (*xy.min(axis=0), *xy.max(axis=0))

    def _draw_grid(ax, xmin, ymin, xmax, ymax):
        if not grid_size or grid_size <= 0:
            return
        if isinstance(grid_origin, tuple | list):
            x0, y0 = float(grid_origin[0]), float(grid_origin[1])
        elif grid_origin == "zero":
            x0 = y0 = 0.0
        else:
            x0 = np.floor(xmin / grid_size) * grid_size
            y0 = np.floor(ymin / grid_size) * grid_size
        # Grid lines from the origin outwards, covering the extent
        x_start = x0 + np.floor((xmin - x0) / grid_size) * grid_size
        y_start = y0 + np.floor((ymin - y0) / grid_size) * grid_size
        xs = np.arange(x_start, xmax + grid_size, grid_size)
        ys = np.arange(y_start, ymax + grid_size, grid_size)
        if len(xs) + len(ys) > grid_max_lines:
            logger.warning(
                f"Grid needs {len(xs) + len(ys)} lines (> {grid_max_lines}); skipped"
            )
            return
        style = dict(
            colors=grid_color, linewidth=grid_linewidth, alpha=grid_alpha, zorder=10
        )
        ax.vlines(xs, ys[0], ys[-1], **style)
        ax.hlines(ys, xs[0], xs[-1], **style)

    fig, axes = _panel_grid(len(panels), n_cols, figsize_per_panel, extra_width=2.5)
    for ax, (name, ad) in zip(axes, panels.items(), strict=False):
        xy = _coords(ad)
        labels = ad.obs[color_col].astype(str).to_numpy()
        point_colors = [color_lookup.get(lbl, other_color) for lbl in labels]
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            c=point_colors,
            s=point_size,
            alpha=alpha,
            edgecolors="none",
            linewidths=0,
            rasterized=True,
            zorder=2,
        )
        if extent is not None:
            xmin, ymin, xmax, ymax = extent
        else:
            (xmin, ymin), (xmax, ymax) = xy.min(axis=0), xy.max(axis=0)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        _draw_grid(ax, xmin, ymin, xmax, ymax)
        ax.set_aspect("equal", adjustable="box")
        if invert_y:
            ax.invert_yaxis()
        if name != "":
            ax.set_title(str(name), fontsize=10)
        if show_axes:
            ax.set_xlabel(x_col or "X")
            ax.set_ylabel(y_col or "Y")
        else:
            ax.set_xticks([])
            ax.set_yticks([])

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=color_lookup[c],
            markersize=8,
            label=c,
        )
        for c in shown
    ]
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
        fontsize=9,
        title=legend_title or color_col,
    )
    if title:
        fig.suptitle(title, fontsize=14)
    fig.tight_layout()

    _save(fig, save_path)
    return fig
