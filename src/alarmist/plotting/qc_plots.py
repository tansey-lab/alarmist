"""
Input-data quality-control plots.
"""

import logging

import anndata
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import issparse

logger = logging.getLogger(__name__)


def plot_total_counts_per_cell(
    adata: anndata.AnnData,
    layer: str | None = None,
    bins: int = 50,
    log10: bool = False,
    threshold: float | None = None,
    xlim: tuple[float, float] | None = None,
    color: str = "#4C72B0",
    title: str | None = None,
    ax: plt.Axes | None = None,
    figsize: tuple[float, float] = (6, 4),
    save_path: str | None = None,
) -> plt.Figure:
    """
    Histogram of total counts per cell.

    Useful to check that the matrix holds raw counts (ALARMIST needs raw counts
    in ``X`` and ``layers['counts']``) and to pick a low-count filter.

    Parameters
    ----------
    adata : AnnData
        Cells × genes.
    layer : str, optional
        Layer to sum. Default: ``adata.X``.
    bins : int, default 50
        Histogram bins.
    log10 : bool, default False
        Plot ``log10(total + 1)``.
    threshold : float, optional
        Draw a vertical line here (on the raw scale; transformed when ``log10``).
    xlim : tuple, optional
        X-axis limits.
    color : str, default '#4C72B0'
        Bar colour.
    title : str, optional
        Axes title.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into.
    figsize : tuple, default (6, 4)
        Figure size when ``ax`` is None.
    save_path : str, optional
        File to save to.

    Returns
    -------
    matplotlib.figure.Figure
    """
    X = adata.X if layer is None else adata.layers[layer]
    totals = np.asarray(X.sum(axis=1)).ravel() if issparse(X) else np.asarray(X).sum(1)

    values = np.log10(totals + 1) if log10 else totals
    xlabel = "log10(total counts per cell + 1)" if log10 else "Total counts per cell"

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    ax.hist(values, bins=bins, color=color, edgecolor="white")
    if threshold is not None:
        x = np.log10(threshold + 1) if log10 else threshold
        ax.axvline(x, color="red", linewidth=1.5, label=f"threshold = {threshold:g}")
        ax.legend(frameon=False, fontsize=8)
    if xlim is not None:
        ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of cells")
    source = "X" if layer is None else f"layers['{layer}']"
    ax.set_title(title or f"Total counts per cell ({source})")
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.debug(f"Saved: {save_path}")
    return fig
