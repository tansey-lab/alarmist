"""
Motif-similarity heatmaps for comparing two ALARMIST runs.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)


def plot_motif_similarity(
    similarity: np.ndarray,
    row_labels: list | None = None,
    col_labels: list | None = None,
    row_name: str = "Run A motifs",
    col_name: str = "Run B motifs",
    matches: str | list[tuple[int, int]] | None = None,
    min_similarity: float | None = None,
    cmap: str = "Reds",
    vmin: float | None = None,
    vmax: float | None = None,
    annotate: bool = True,
    annotate_min: float | None = None,
    fmt: str = "{:.2f}",
    fontsize: float = 7,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    colorbar_label: str = "Similarity",
    ax: plt.Axes | None = None,
    save_path: str | None = None,
) -> plt.Figure:
    """
    Heatmap of motif-to-motif similarity between two runs, with matches boxed.

    Parameters
    ----------
    similarity : np.ndarray
        (n_a, n_b) matrix, e.g. from :func:`alarmist.core.similarity.row_similarity`.
    row_labels, col_labels : list, optional
        Tick labels. Default: motif indices 0..n-1.
    row_name, col_name : str
        Axis labels.
    matches : 'hungarian' or list of (row, col), optional
        Cells to outline. 'hungarian' computes the one-to-one matching that
        maximises total similarity (:func:`alarmist.core.similarity.match_motifs`).
    min_similarity : float, optional
        With ``matches='hungarian'``, only outline pairs at or above this value.
    cmap : str, default 'Reds'
        Colormap.
    vmin, vmax : float, optional
        Colour limits. Default: data min / max.
    annotate : bool, default True
        Write values in cells.
    annotate_min : float, optional
        Only annotate cells with value >= this.
    fmt : str, default '{:.2f}'
        Annotation format string.
    fontsize : float, default 7
        Annotation font size.
    figsize : tuple, optional
        Figure size; default scales with the matrix shape.
    title : str, optional
        Axes title.
    colorbar_label : str, default 'Similarity'
        Colourbar label.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into.
    save_path : str, optional
        File to save to.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from alarmist.core.similarity import row_similarity
    >>> S = row_similarity(lri_factors_xenium.T, lri_factors_cosmx.T, method="cosine")
    >>> fig = plot_motif_similarity(
    ...     S, row_name="Xenium motifs", col_name="CosMx motifs",
    ...     matches="hungarian", min_similarity=0.4,
    ... )
    """
    from matplotlib.patches import Rectangle

    from alarmist.core.similarity import match_motifs

    S = np.asarray(similarity, dtype=float)
    n_rows, n_cols = S.shape
    row_labels = list(range(n_rows)) if row_labels is None else list(row_labels)
    col_labels = list(range(n_cols)) if col_labels is None else list(col_labels)

    vmin = float(np.nanmin(S)) if vmin is None else vmin
    vmax = float(np.nanmax(S)) if vmax is None else vmax
    if vmax <= vmin:
        vmax = vmin + 1e-12

    if ax is None:
        if figsize is None:
            figsize = (0.35 * n_cols + 2.5, 0.35 * n_rows + 1.5)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    im = ax.imshow(S, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(n_cols))
    ax.set_yticks(np.arange(n_rows))
    ax.set_xticklabels(col_labels, rotation=90)
    ax.set_yticklabels(row_labels)
    ax.set_xlabel(col_name)
    ax.set_ylabel(row_name)
    if title:
        ax.set_title(title)

    if annotate:
        for i in range(n_rows):
            for j in range(n_cols):
                val = S[i, j]
                if np.isnan(val) or (annotate_min is not None and val < annotate_min):
                    continue
                color = "white" if im.norm(val) > 0.6 else "black"
                ax.text(
                    j,
                    i,
                    fmt.format(val),
                    ha="center",
                    va="center",
                    color=color,
                    fontsize=fontsize,
                )

    if isinstance(matches, str):
        if matches != "hungarian":
            raise ValueError(f"Unknown matches '{matches}'")
        matches = match_motifs(S, min_similarity=min_similarity)
    for i, j in matches or []:
        ax.add_patch(
            Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="black", lw=1.5)
        )

    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label(colorbar_label)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.debug(f"Saved: {save_path}")
    return fig
