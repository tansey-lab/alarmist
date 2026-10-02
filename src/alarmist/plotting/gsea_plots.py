"""
Gene-set enrichment (GSEA) result plots.
"""

import logging
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd

logger = logging.getLogger(__name__)


def gsea_barplot(
    gsea_result: Any,
    n_top: int = 10,
    fdr_threshold: float = 0.25,
    term_col: str = "Term",
    nes_col: str = "NES",
    fdr_col: str = "FDR q-val",
    pos_label: str = "Positive NES",
    neg_label: str = "Negative NES",
    pos_color: str = "#E64B35",
    neg_color: str = "#4DBBD5",
    ns_color: str = "#CCCCCC",
    max_term_length: int = 50,
    show_qvalues: bool = True,
    title: str = "GSEA",
    ax: plt.Axes | None = None,
    figsize: tuple[float, float] = (10, 8),
    save_path: str | None = None,
) -> plt.Figure:
    """
    Horizontal NES barplot of the top up- and down-regulated gene sets.

    Takes the ``n_top`` lowest-FDR terms with positive NES and the ``n_top`` with
    negative NES, sorted by NES. Bars at FDR >= ``fdr_threshold`` are grey.

    Parameters
    ----------
    gsea_result : pd.DataFrame or gseapy result
        A gseapy prerank/GSEA result object (its ``res2d`` table is used) or a
        DataFrame with term, NES and FDR columns.
    n_top : int, default 10
        Terms shown per direction.
    fdr_threshold : float, default 0.25
        FDR below which a bar is coloured (0.25 is the GSEA convention).
    term_col, nes_col, fdr_col : str
        Column names (gseapy defaults).
    pos_label, neg_label : str
        Legend labels for positive / negative NES, e.g. 'Up in LUAD' / 'Up in AIS'.
    pos_color, neg_color, ns_color : str
        Bar colours.
    max_term_length : int, default 50
        Truncate term names (underscores become spaces) beyond this length.
    show_qvalues : bool, default True
        Print the FDR next to each bar.
    title : str, default 'GSEA'
        Axes title.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into.
    figsize : tuple, default (10, 8)
        Figure size when ``ax`` is None.
    save_path : str, optional
        File to save to.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> import gseapy as gp
    >>> res = gp.prerank(rnk=ranked_genes, gene_sets="MSigDB_Hallmark_2020")
    >>> fig = gsea_barplot(res, pos_label="Up in LUAD", neg_label="Up in AIS")
    """
    from matplotlib.patches import Patch

    df = gsea_result.res2d if hasattr(gsea_result, "res2d") else gsea_result
    missing = {term_col, nes_col, fdr_col} - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in GSEA results: {sorted(missing)}")
    df = df[[term_col, nes_col, fdr_col]].copy()
    df[nes_col] = df[nes_col].astype(float)
    df[fdr_col] = df[fdr_col].astype(float)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    df_plot = (
        pd.concat(
            [
                df[df[nes_col] > 0].nsmallest(n_top, fdr_col),
                df[df[nes_col] < 0].nsmallest(n_top, fdr_col),
            ]
        )
        .drop_duplicates()
        .sort_values(nes_col)
        .reset_index(drop=True)
    )

    ax.set_title(title, fontsize=12, fontweight="bold")
    if df_plot.empty:
        ax.text(
            0.5, 0.5, "No gene sets", ha="center", va="center", transform=ax.transAxes
        )
        return fig

    sig = df_plot[fdr_col] < fdr_threshold
    colors = [
        (pos_color if nes > 0 else neg_color) if s else ns_color
        for nes, s in zip(df_plot[nes_col], sig, strict=True)
    ]

    def _clean(term):
        term = str(term).replace("_", " ")
        return term[:max_term_length] + "..." if len(term) > max_term_length else term

    y = range(len(df_plot))
    ax.barh(y, df_plot[nes_col], color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(list(y))
    ax.set_yticklabels([_clean(t) for t in df_plot[term_col]], fontsize=9)
    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=11)
    ax.axvline(0, color="black", linewidth=0.8)

    if show_qvalues:
        for i, (nes, q) in enumerate(
            zip(df_plot[nes_col], df_plot[fdr_col], strict=True)
        ):
            ax.text(
                nes + (0.05 if nes > 0 else -0.05),
                i,
                f"q={q:.3f}",
                va="center",
                ha="left" if nes > 0 else "right",
                fontsize=7,
                color="gray",
            )
        lo, hi = ax.get_xlim()
        pad = 0.15 * (hi - lo)
        ax.set_xlim(lo - pad if lo < 0 else lo, hi + pad if hi > 0 else hi)

    ax.legend(
        handles=[
            Patch(facecolor=pos_color, edgecolor="black", label=pos_label),
            Patch(facecolor=neg_color, edgecolor="black", label=neg_label),
            Patch(
                facecolor=ns_color, edgecolor="black", label=f"FDR ≥ {fdr_threshold}"
            ),
        ],
        loc="lower right",
        fontsize=8,
    )
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.debug(f"Saved: {save_path}")
    return fig
