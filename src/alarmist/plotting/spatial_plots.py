"""
Spatial distribution plotting functions
"""

import logging

import anndata
import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)


def plot_cells_per_patch(
    adata: anndata.AnnData | dict[str, anndata.AnnData],
    bins: int = 50,
    figsize: tuple = (8, 5),
    title: str | None = None,
    color: str = "#4C72B0",
    save_path: str | None = None,
    show: bool = True,
) -> plt.Figure:
    """
    Plot histogram of cell counts per patch.

    Parameters
    ----------
    adata : AnnData or Dict[str, AnnData]
        AnnData object(s) with 'patch_idx' column in obs (from run_patchify).
        If dict, all samples are combined.
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

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object

    Examples
    --------
    >>> results = analyzer.run_patchify(adata, output_dir=results_dir)
    >>> fig = al.plot_cells_per_patch(adata, save_path='cells_per_patch.png')
    """
    # Handle dict input (multi-sample)
    if isinstance(adata, dict):
        all_patch_idx = []
        for sample_id, ad in adata.items():
            if "patch_idx" not in ad.obs.columns:
                raise ValueError(
                    f"'patch_idx' not found in adata.obs for sample '{sample_id}'. "
                    "Run run_patchify first."
                )
            all_patch_idx.extend(ad.obs["patch_idx"].values)
        patch_idx = np.array(all_patch_idx)
    else:
        if "patch_idx" not in adata.obs.columns:
            raise ValueError(
                "'patch_idx' not found in adata.obs. Run run_patchify first."
            )
        patch_idx = adata.obs["patch_idx"].values

    # Filter out invalid patch indices (-1)
    valid_patch_idx = patch_idx[patch_idx >= 0]

    # Count cells per patch
    unique_patches, counts = np.unique(valid_patch_idx, return_counts=True)

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
        f"Total patches: {len(unique_patches)}\n"
        f"Total cells: {len(valid_patch_idx)}\n"
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
