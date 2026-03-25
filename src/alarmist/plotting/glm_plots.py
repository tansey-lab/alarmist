"""
GLM results visualization functions

Includes volcano plots and forest plots for differential expression results.
"""

import glob
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
from adjustText import adjust_text
from matplotlib.backends.backend_pdf import PdfPages

from alarmist.constants import (
    COLUMN_NAME_CELL_TYPE,
    COLUMN_NAME_GENE,
    COLUMN_NAME_LOGFC,
    COLUMN_NAME_NEG_LOG10_Q,
    COLUMN_NAME_QVAL,
    COLUMN_NAME_SE,
)

logger = logging.getLogger(__name__)


def volcano_plot(
    df,
    x_col,
    y_col,
    label_col=None,
    fdr=0.1,
    x_threshold=1,
    marker="o",
    n_top=10,
    figsize=(10, 10),
    fontsize=8,
    ax=None,
):
    """
    Draw a volcano plot with customized coloring and labeling

    Parameters
    ----------
    df : pd.DataFrame
        Data frame with DE results
    x_col : str
        Column name for log fold changes
    y_col : str
        Column name for -log10(q-value)
    label_col : str, optional
        Column name for gene labels
    fdr : float, default 0.1
        FDR threshold
    x_threshold : float, default 1
        Log fold change threshold
    marker : str, default 'o'
        Marker style
    n_top : int, default 10
        Number of top genes to label per direction
    figsize : tuple, default (10, 10)
        Figure size
    fontsize : int, default 8
        Font size
    ax : matplotlib.axes.Axes, optional
        Axes to plot on

    Returns
    -------
    matplotlib.axes.Axes
        The axes object
    """
    # Nature journal settings
    plt.rcParams["font.family"] = "Arial"
    colors = [
        "#E64B35FF",
        "#3C5488FF",
        "#00A087FF",
        "#4DBBD5FF",
        "#F39B7FFF",
        "#8491B4FF",
        "#91D1C2FF",
        "#DC0000FF",
        "#7E6148FF",
        "#B09C85FF",
    ]
    sns.set_palette(sns.color_palette(colors))
    sns.set(rc={"figure.figsize": figsize, "font.size": fontsize})
    sns.set_style("white")

    threshold_fdr = -np.log10(fdr)
    x = df[x_col].values
    y = df[y_col].values

    sigmask = y >= threshold_fdr
    effmask = np.abs(x) >= x_threshold
    rightmask = x >= x_threshold
    leftmask = x <= -x_threshold

    # Scoring for top candidate selection
    upper_right = (
        x / np.abs(x).max() + y / y.max()
        if np.abs(x).max() > 0 and y.max() > 0
        else np.zeros_like(x)
    )
    upper_left = (
        -x / np.abs(x).max() + y / y.max()
        if np.abs(x).max() > 0 and y.max() > 0
        else np.zeros_like(x)
    )
    top_right = np.argsort(upper_right)[::-1]
    top_left = np.argsort(upper_left)[::-1]

    indices = np.arange(len(y))
    top_right = top_right[sigmask[top_right] & rightmask[top_right]]
    top_left = top_left[sigmask[top_left] & leftmask[top_left]]
    bottom_right = indices[(~sigmask) & rightmask]
    bottom_left = indices[(~sigmask) & leftmask]
    middle = indices[~effmask]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot the values differently depending on where they lie
    ax.scatter(x[middle], y[middle], color="gray", alpha=0.2, marker=marker)
    ax.scatter(
        x[bottom_right], y[bottom_right], color="black", alpha=0.7, marker=marker
    )
    ax.scatter(x[bottom_left], y[bottom_left], color="black", alpha=0.7, marker=marker)
    ax.scatter(x[top_left], y[top_left], color="blue", alpha=0.7, marker=marker)
    ax.scatter(x[top_right], y[top_right], color="red", alpha=0.7, marker=marker)

    ax.set_xlabel(r"$\log_2$ Fold-Change", fontsize=fontsize + 14)
    ax.set_ylabel(r"-$\log_{10}(q)$", fontsize=fontsize + 14)

    # Set axis limits based on significant hits
    if effmask.any():
        ax.set_ylim([0, df.loc[effmask, y_col].max() * 1.1])
    else:
        ax.set_ylim([0, y.max() * 1.1 if y.max() > 0 else 1])

    n_hits = (sigmask & effmask).sum()
    if n_hits > 0:
        max_abs_sig_x = np.abs(x[sigmask & effmask]).max()
        ax.set_xlim([-max_abs_sig_x * 1.1, max_abs_sig_x * 1.1])
    else:
        ax.set_xlim([-np.abs(x).max() * 1.1, np.abs(x).max() * 1.1])

    # Label top genes with white background boxes
    if label_col is not None:
        head_df = df.iloc[np.concatenate([top_right[:n_top], top_left[:n_top]])]
        texts = []
        for xi, yi, txt in zip(head_df[x_col], head_df[y_col], head_df[label_col]):
            texts.append(
                ax.text(
                    xi,
                    yi,
                    txt,
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    bbox=dict(
                        boxstyle="round,pad=0.2",
                        facecolor="white",
                        edgecolor="gray",
                        alpha=0.9,
                    ),
                )
            )

        adjust_text(
            texts,
            force_text=(0.5, 0.5),
            arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
            ax=ax,
            verbose=False,
        )

    # Threshold lines
    ax.axvline(x_threshold, ls="--", color="black", alpha=0.5)
    ax.axvline(-x_threshold, ls="--", color="black", alpha=0.5)
    ax.axhline(threshold_fdr, ls="--", color="black", alpha=0.5)

    return ax


def forest_plot(
    df,
    ax=None,
    n_top=20,
    effect_col=COLUMN_NAME_LOGFC,
    q_col=COLUMN_NAME_QVAL,
    gene_col=COLUMN_NAME_GENE,
):
    """
    Draw a horizontal forest plot of the strongest absolute logfc signals

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with DE results
    ax : matplotlib.axes.Axes, optional
        Axes to plot on
    n_top : int, default 20
        Number of top genes to show
    effect_col : str, default COLUMN_NAME_LOGFC
        Effect size column
    q_col : str, default COLUMN_NAME_QVAL
        Q-value column
    gene_col : str, default COLUMN_NAME_GENE
        Gene name column

    Returns
    -------
    matplotlib.axes.Axes
        The axes object
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 8))

    # Rank and pick top
    df["_score"] = np.abs(df[effect_col]) * -np.log10(df[q_col].clip(1e-300))
    df = df[np.abs(df[effect_col]) <= 10].copy()
    top = df.sort_values("_score", ascending=False).head(n_top).sort_values(effect_col)

    y = np.arange(len(top))
    colors = top[effect_col].map(lambda x: "red" if x > 0 else "blue")

    # 95% Wald CI: effect ± 1.96 * SE
    se_x = top[COLUMN_NAME_SE]
    lo = top[effect_col] - 1.96 * se_x
    hi = top[effect_col] + 1.96 * se_x
    xerr = np.vstack([top[effect_col] - lo, hi - top[effect_col]])

    # Plot
    ax.errorbar(
        x=top[effect_col], y=y, xerr=xerr, fmt="none", elinewidth=1, capsize=3, zorder=2
    )
    ax.scatter(top[effect_col], y, c=colors, s=25, zorder=2)

    # Gene labels & q-values
    xmax = ax.get_xlim()[1]
    for yi, gene, q in zip(y, top[gene_col], top[q_col]):
        ax.text(xmax * 1.02, yi, f"{gene}  (q={q:.1e})", va="center", fontsize=7)

    ax.axvline(0, color="black", lw=1)
    ax.set_yticks([])
    ax.set_xlabel("log$_2$ Fold-Change")
    sns.despine(ax=ax, left=True)
    ax.set_title(f"Top {n_top} DE genes", fontsize=9)
    plt.tight_layout()

    return ax


def generate_volcano_plots(
    results_dir: str,
    output_dir: str | None,
    n_motifs: int,
    adata,
    cell_types,
    all_genes,
    exclusion_mask,
    min_expression_frac: float,
    fdr_threshold: float,
    lfc_threshold: float,
    n_top_genes: int,
):
    """Generate volcano plots

    Parameters
    ----------
    results_dir : str
        Directory containing DE results
    output_dir : str, optional
        Directory to save plots. If None, returns list of figures.
    n_motifs : int
        Number of motifs
    adata : anndata.AnnData
        Annotated data object
    cell_types : array-like
        Cell type names
    all_genes : array-like
        Gene names
    exclusion_mask : np.ndarray
        Marker gene exclusion mask
    min_expression_frac : float
        Minimum expression fraction
    fdr_threshold : float
        FDR threshold
    lfc_threshold : float
        Log fold change threshold
    n_top_genes : int
        Number of top genes to label

    Returns
    -------
    list of matplotlib.figure.Figure or None
        List of figure objects if output_dir is None, otherwise None
    """
    figures = []

    for motif_id in range(n_motifs):
        files = sorted(
            glob.glob(
                os.path.join(results_dir, f"motif_{motif_id}_celltype_*_de_results.csv")
            )
        )
        if not files:
            logger.debug(f"No results found for motif {motif_id}")
            continue

        fig, axes = plt.subplots(4, 4, figsize=(16, 16))
        axes = axes.flatten()

        for ax, fp in zip(axes, files[:16]):
            df = pd.read_csv(fp)
            ct = (
                os.path.basename(fp)
                .split(f"motif_{motif_id}_celltype_")[1]
                .replace("_de_results.csv", "")
            )

            logger.debug(f"  Motif {motif_id}, {ct}: Original genes = {len(df)}")

            try:
                cidx = list(cell_types).index(ct)
            except ValueError:
                logger.debug(f"Warning: Cell type {ct} not found")
                continue

            # Filter by expression and markers
            df = filter_genes_for_plot(
                df, ct, adata, all_genes, exclusion_mask, cidx, min_expression_frac
            )

            # Compute -log10(q) with jitter
            df[COLUMN_NAME_NEG_LOG10_Q] = -np.log10(df[COLUMN_NAME_QVAL].clip(1e-300))
            m = df[COLUMN_NAME_NEG_LOG10_Q] >= 300
            if m.any():
                df.loc[m, COLUMN_NAME_NEG_LOG10_Q] = 300 + np.random.normal(
                    0, 15, m.sum()
                )

            # Plot
            volcano_plot(
                df,
                COLUMN_NAME_LOGFC,
                COLUMN_NAME_NEG_LOG10_Q,
                label_col=COLUMN_NAME_GENE,
                fdr=fdr_threshold,
                x_threshold=0.2,
                marker="o",
                n_top=30,
                fontsize=5,
                ax=ax,
            )
            ax.set_title(ct, fontsize=10, pad=5)

        # Remove unused subplots
        for j in range(len(files[:16]), 16):
            fig.delaxes(axes[j])

        plt.suptitle(
            f"Motif {motif_id} Volcano Plots (marker genes filtered)", fontsize=14
        )
        plt.tight_layout(rect=[0, 0, 1, 0.95])

        figures.append(fig)
        logger.debug(f"Motif {motif_id} volcano plots completed")

    if output_dir:
        out_volcano = os.path.join(output_dir, "volcano_plots_filtered.pdf")
        with PdfPages(out_volcano) as pdf:
            for fig in figures:
                pdf.savefig(fig)
                plt.close(fig)
        logger.debug(f"Volcano plots saved to: {out_volcano}")
        return None
    else:
        return figures


def generate_forest_plots(
    results_dir: str,
    output_dir: str | None,
    n_motifs: int,
    adata,
    cell_types,
    all_genes,
    exclusion_mask,
    min_expression_frac: float,
):
    """Generate forest plots

    Parameters
    ----------
    results_dir : str
        Directory containing DE results
    output_dir : str, optional
        Directory to save plots. If None, returns list of figures.
    n_motifs : int
        Number of motifs
    adata : anndata.AnnData
        Annotated data object
    cell_types : array-like
        Cell type names
    all_genes : array-like
        Gene names
    exclusion_mask : np.ndarray
        Marker gene exclusion mask
    min_expression_frac : float
        Minimum expression fraction

    Returns
    -------
    list of matplotlib.figure.Figure or None
        List of figure objects if output_dir is None, otherwise None
    """
    figures = []

    for motif_id in range(n_motifs):
        csv_paths = sorted(
            glob.glob(
                os.path.join(results_dir, f"motif_{motif_id}_celltype_*_de_results.csv")
            )
        )
        if not csv_paths:
            continue

        fig, axes = plt.subplots(4, 4, figsize=(16, 24))
        axes = axes.flatten()

        for ax, csv_fp in zip(axes, csv_paths[:16]):
            df = pd.read_csv(csv_fp)
            ct = (
                os.path.basename(csv_fp)
                .split(f"motif_{motif_id}_celltype_")[1]
                .replace("_de_results.csv", "")
            )

            try:
                ct_idx = list(cell_types).index(ct)
            except ValueError:
                continue

            # Filter genes
            df = filter_genes_for_plot(
                df, ct, adata, all_genes, exclusion_mask, ct_idx, min_expression_frac
            )

            forest_plot(df, ax=ax, n_top=30)
            ax.set_title(ct, fontsize=9)

        # Remove unused subplots
        for j in range(len(csv_paths[:16]), 16):
            fig.delaxes(axes[j])

        plt.suptitle(
            f"Motif {motif_id} – Forest plots (marker genes filtered)", fontsize=14
        )
        plt.tight_layout(rect=[0, 0, 1, 0.965])

        figures.append(fig)
        logger.debug(f"Motif {motif_id} forest plots completed")

    if output_dir:
        out_forest = os.path.join(output_dir, "forest_plots_filtered.pdf")
        with PdfPages(out_forest) as pdf:
            for fig in figures:
                pdf.savefig(fig)
                plt.close(fig)
        logger.debug(f"Forest plots saved to: {out_forest}")
        return None
    else:
        return figures


def filter_genes_for_plot(
    df, ct, adata, all_genes, exclusion_mask, cidx, min_expression_frac
):
    """Filter genes by expression and marker status"""
    glm_genes_set = set(df[COLUMN_NAME_GENE])

    # Expression filter
    subX = adata[adata.obs[COLUMN_NAME_CELL_TYPE] == ct].X
    if sp.issparse(subX):
        expr_frac = np.asarray((subX > 0).sum(axis=0)).ravel() / subX.shape[0]
    else:
        expr_frac = np.sum(subX > 0, axis=0) / subX.shape[0]

    gene_to_idx = {gene: i for i, gene in enumerate(all_genes)}
    genes_to_keep = []
    for gene in glm_genes_set:
        if gene in gene_to_idx:
            gene_idx = gene_to_idx[gene]
            if expr_frac[gene_idx] >= min_expression_frac:
                genes_to_keep.append(gene)

    df = df[df[COLUMN_NAME_GENE].isin(genes_to_keep)].copy()
    logger.debug(f"    After expression filter: {len(df)} genes")

    # Marker filter
    other = [i for i in range(len(exclusion_mask)) if i != cidx]
    drop_mask = exclusion_mask[other, :].any(axis=0)
    marker_genes_to_drop = set(all_genes[drop_mask]) & glm_genes_set
    df = df[~df[COLUMN_NAME_GENE].isin(marker_genes_to_drop)]
    logger.debug(f"    After marker filter: {len(df)} genes")

    return df
