"""
Plot overlap between LRI database genes and input data genes.
"""

import logging
from collections.abc import Iterable

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

logger = logging.getLogger(__name__)


def plot_lri_database_overlap(
    data_genes: Iterable[str],
    database_genes: Iterable[str],
    title: str | None = None,
    data_label: str = "Input data genes",
    database_label: str = "LRI database genes",
    figsize: tuple = (8, 6),
) -> plt.Figure:
    """
    Render a two-circle Venn diagram showing overlap between the gene set in
    the input data and the gene set covered by the ligand-receptor database.

    Parameters
    ----------
    data_genes : Iterable[str]
        Gene names present in the input AnnData (``adata.var_names``).
    database_genes : Iterable[str]
        All ligand/receptor gene names referenced by the LRI database
        (unfiltered, including complex subunits).
    title : str, optional
        Figure title.
    """
    data_set = set(data_genes)
    db_set = set(database_genes)
    only_data = len(data_set - db_set)
    only_db = len(db_set - data_set)
    both = len(data_set & db_set)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal")
    ax.axis("off")

    r = 1.0
    offset = 0.85
    left = Circle((-offset, 0), r, alpha=0.45, color="#4C72B0", linewidth=0)
    right = Circle((offset, 0), r, alpha=0.45, color="#DD8452", linewidth=0)
    ax.add_patch(left)
    ax.add_patch(right)

    ax.text(
        -offset - r * 0.55,
        0,
        f"{only_data:,}",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
    )
    ax.text(
        offset + r * 0.55,
        0,
        f"{only_db:,}",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
    )
    ax.text(0, 0, f"{both:,}", ha="center", va="center", fontsize=16, fontweight="bold")

    ax.text(
        -offset - r * 0.55,
        r + 0.15,
        f"{data_label}\n(n={len(data_set):,})",
        ha="center",
        va="bottom",
        fontsize=11,
    )
    ax.text(
        offset + r * 0.55,
        r + 0.15,
        f"{database_label}\n(n={len(db_set):,})",
        ha="center",
        va="bottom",
        fontsize=11,
    )

    if title is None:
        pct = (both / len(db_set) * 100) if db_set else 0.0
        title = f"LRI database coverage in input data ({both:,}/{len(db_set):,} = {pct:.1f}%)"
    ax.set_title(title, fontsize=12)

    ax.set_xlim(-offset - r - 0.6, offset + r + 0.6)
    ax.set_ylim(-r - 0.3, r + 0.9)
    fig.tight_layout()
    return fig
