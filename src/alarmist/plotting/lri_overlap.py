"""
Plot overlap between LRI database genes and input data genes.
"""

import logging
from collections.abc import Iterable

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle

from alarmist.core.lri import _split_gene_complex

logger = logging.getLogger(__name__)


def _draw_two_circle_venn(
    ax,
    only_left: int,
    only_right: int,
    both: int,
    left_label: str,
    right_label: str,
    left_total: int,
    right_total: int,
    left_color: str = "#4C72B0",
    right_color: str = "#DD8452",
) -> None:
    """Render two overlapping circles with counts inside and labels above."""
    ax.set_aspect("equal")
    ax.axis("off")
    r = 1.0
    offset = 0.85
    ax.add_patch(Circle((-offset, 0), r, alpha=0.45, color=left_color, linewidth=0))
    ax.add_patch(Circle((offset, 0), r, alpha=0.45, color=right_color, linewidth=0))

    ax.text(
        -offset - r * 0.55,
        0,
        f"{only_left:,}",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
    )
    ax.text(
        offset + r * 0.55,
        0,
        f"{only_right:,}",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
    )
    ax.text(0, 0, f"{both:,}", ha="center", va="center", fontsize=16, fontweight="bold")

    ax.text(
        -offset - r * 0.55,
        r + 0.15,
        f"{left_label}\n(n={left_total:,})",
        ha="center",
        va="bottom",
        fontsize=11,
    )
    ax.text(
        offset + r * 0.55,
        r + 0.15,
        f"{right_label}\n(n={right_total:,})",
        ha="center",
        va="bottom",
        fontsize=11,
    )

    ax.set_xlim(-offset - r - 0.6, offset + r + 0.6)
    ax.set_ylim(-r - 0.3, r + 0.9)


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
    _draw_two_circle_venn(
        ax,
        only_left=only_data,
        only_right=only_db,
        both=both,
        left_label=data_label,
        right_label=database_label,
        left_total=len(data_set),
        right_total=len(db_set),
    )

    if title is None:
        pct = (both / len(db_set) * 100) if db_set else 0.0
        title = f"LRI database gene coverage in input data ({both:,}/{len(db_set):,} = {pct:.1f}%)"
    ax.set_title(title, fontsize=12)

    fig.tight_layout()
    return fig


def plot_lr_pair_overlap(
    resource: pd.DataFrame,
    data_genes: Iterable[str],
    title: str | None = None,
    figsize: tuple = (8, 6),
) -> plt.Figure:
    """
    Two-circle Venn at the ligand-receptor *pair* level.

    Left set  = LR pairs whose ligand (all complex subunits) is present in the data.
    Right set = LR pairs whose receptor (all complex subunits) is present in the data.
    Intersection = pairs usable by the analysis (both sides fully covered).
    Outside both = pairs unusable in this dataset (drawn as a count beside the diagram).
    """
    data_set = set(map(str, data_genes))

    total = 0
    ligand_ok_idx: set[int] = set()
    receptor_ok_idx: set[int] = set()
    for idx in range(len(resource)):
        ligand = resource.iloc[idx]["ligand"]
        receptor = resource.iloc[idx]["receptor"]
        if pd.isna(ligand) or pd.isna(receptor):
            continue
        total += 1
        if all(g in data_set for g in _split_gene_complex(str(ligand))):
            ligand_ok_idx.add(idx)
        if all(g in data_set for g in _split_gene_complex(str(receptor))):
            receptor_ok_idx.add(idx)

    both = len(ligand_ok_idx & receptor_ok_idx)
    only_left = len(ligand_ok_idx - receptor_ok_idx)
    only_right = len(receptor_ok_idx - ligand_ok_idx)
    neither = total - (both + only_left + only_right)

    fig, ax = plt.subplots(figsize=figsize)
    _draw_two_circle_venn(
        ax,
        only_left=only_left,
        only_right=only_right,
        both=both,
        left_label="Ligand covered",
        right_label="Receptor covered",
        left_total=len(ligand_ok_idx),
        right_total=len(receptor_ok_idx),
        left_color="#55A868",
        right_color="#C44E52",
    )
    ax.text(
        0,
        -1.15,
        f"Neither side covered: {neither:,}",
        ha="center",
        va="top",
        fontsize=10,
        style="italic",
        color="#444",
    )

    if title is None:
        pct = (both / total * 100) if total else 0.0
        title = (
            f"LR pair coverage in input data "
            f"({both:,}/{total:,} pairs usable = {pct:.1f}%)"
        )
    ax.set_title(title, fontsize=12)

    fig.tight_layout()
    return fig
