"""
BPTF visualization functions
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional


def plot_cells_per_patch(patch_metadata_df: pd.DataFrame, save_path: Optional[str] = None):
    """Plot distribution of cells per patch"""
    plt.figure(figsize=(8, 6))
    plt.hist(patch_metadata_df['n_cells'], bins=50, edgecolor='black')
    plt.xlabel('Number of cells per patch')
    plt.ylabel('Number of patches')
    plt.title('Distribution of cells per patch')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return plt.gcf()


def plot_factor_distributions(patch_loadings: np.ndarray,
                              lri_factors: np.ndarray,
                              save_path: Optional[str] = None):
    """Plot LRI participation and patch loading distributions"""
    # Compute data
    lri_participation = lri_factors.sum(axis=0)
    patch_loading_totals = patch_loadings.sum(axis=1)

    # Create side-by-side subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left subplot: LRI participation distribution
    axes[0].hist(lri_participation, bins=100, alpha=0.7, edgecolor='black')
    axes[0].set_title('LRI Factor Distribution')
    axes[0].set_xlabel('Total Factor')
    axes[0].set_ylabel('Number of LRIs (log)')
    axes[0].set_yscale('log')
    axes[0].grid(True)

    # Right subplot: Patch loading distribution
    axes[1].hist(patch_loading_totals, bins=100, alpha=0.7, edgecolor='black')
    axes[1].set_title('Patch Loading Distribution')
    axes[1].set_xlabel('Total Loading')
    axes[1].set_ylabel('Number of Patches')
    axes[1].set_yscale('log')
    axes[1].grid(True)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig


def plot_factor_sparsity(patch_loadings: np.ndarray,
                        lri_factors: np.ndarray,
                        save_path: Optional[str] = None,
                        threshold: float = 1e-6):
    """Plot sparsity of factor matrices"""
    fig, ax = plt.subplots(figsize=(8, 6))

    patch_sparsity = (patch_loadings < threshold).mean()
    lri_sparsity = (lri_factors < threshold).mean()

    ax.bar(['Patch Loadings', 'LRI Factors'], [patch_sparsity, lri_sparsity])
    ax.set_title('Factor Matrix Sparsity')
    ax.set_ylabel(f'Fraction of values < {threshold}')
    ax.set_ylim(0, 1)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig


def plot_motif_activities(patch_loadings: np.ndarray, save_path: Optional[str] = None):
    """Plot motif activity distribution"""
    motif_activities = patch_loadings.sum(axis=0)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(range(len(motif_activities)), motif_activities)
    ax.set_title('Motif Activity Distribution')
    ax.set_xlabel('Motif Index')
    ax.set_ylabel('Total Activity')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig


def plot_bptf_diagnostics(patch_loadings: np.ndarray,
                         lri_factors: np.ndarray,
                         output_dir: Optional[str] = None):
    """Create diagnostic plots for BPTF results

    Parameters
    ----------
    patch_loadings : np.ndarray
        Patch loadings matrix (n_patches × K)
    lri_factors : np.ndarray
        LRI factors matrix (K × n_lris)
    output_dir : str, optional
        Directory to save plots. If None, displays plot without saving.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object with diagnostic plots

    Notes
    -----
    Creates 4 subplots:
    - Row 1: Normalized motif metrics
      - (1,1) Motif coverage across patches
      - (1,2) Overall contribution of each motif
    - Row 2: Factor distributions
      - (2,1) Patch loadings distribution (log-scale)
      - (2,2) LRI factors distribution (log-scale)
    """
    W = patch_loadings  # (n_patches, K)
    V = lri_factors     # (K, n_lri)
    n_patches, K = W.shape

    # ========== Normalization ==========
    # M_k = max_i W_{ik}
    W_max = W.max(axis=0)  # (K,)

    # Avoid division by zero for all-zero motifs
    W_max_safe = W_max.copy()
    W_max_safe[W_max_safe == 0] = 1.0

    # Normalize W: each column's max becomes 1
    W_tilde = W / W_max_safe[None, :]       # (n_patches, K)
    # Scale V accordingly
    V_tilde = V * W_max_safe[:, None]       # (K, n_lri)

    # ========== Compute metrics ==========
    # 1. Coverage: sum_i W_tilde_ik
    S_tilde = W_tilde.sum(axis=0)           # (K,)
    S_tilde_frac = S_tilde / (S_tilde.sum() + 1e-10)

    # 2. Total contribution: sum_{i,j} W_tilde_ik * V_tilde_kj
    V_tilde_sum = V_tilde.sum(axis=1)                        # (K,)
    patch_motif_contrib = W_tilde * V_tilde_sum[None, :]     # (n_patches, K)
    C_k = patch_motif_contrib.sum(axis=0)                    # (K,)
    C_k_frac = C_k / (C_k.sum() + 1e-10)

    # ========== Create figure ==========
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    x = np.arange(K)

    # ========== Row 1: Normalized motif metrics ==========
    # (1,1) Motif coverage
    axes[0, 0].bar(x, S_tilde_frac)
    axes[0, 0].set_xlabel('Motif Index')
    axes[0, 0].set_ylabel('Coverage (fraction)')
    axes[0, 0].set_title('Motif Coverage Across Patches\n(Σᵢ W̃ᵢₖ, normalized)')
    axes[0, 0].set_xticks(x)

    # (1,2) Total contribution
    axes[0, 1].bar(x, C_k_frac)
    axes[0, 1].set_xlabel('Motif Index')
    axes[0, 1].set_ylabel('Total Contribution (fraction)')
    axes[0, 1].set_title('Overall Contribution of Each Motif\n(Σᵢⱼ W̃ᵢₖ Ṽₖⱼ, normalized)')
    axes[0, 1].set_xticks(x)

    # ========== Row 2: Factor distributions ==========
    eps = 1e-10

    # (2,1) Patch loadings distribution
    log_W = np.log(W.flatten() + eps)
    axes[1, 0].hist(log_W, bins=50, edgecolor='black', alpha=0.7)
    axes[1, 0].set_xlabel('log(patch_loadings + 1e-10)')
    axes[1, 0].set_ylabel('Count (log-scale)')
    axes[1, 0].set_yscale('log')
    axes[1, 0].set_title('Distribution of Patch Loadings')

    # (2,2) LRI factors distribution
    log_V = np.log(V.flatten() + eps)
    axes[1, 1].hist(log_V, bins=50, edgecolor='black', alpha=0.7)
    axes[1, 1].set_xlabel('log(lri_factors + 1e-10)')
    axes[1, 1].set_ylabel('Count (log-scale)')
    axes[1, 1].set_yscale('log')
    axes[1, 1].set_title('Distribution of LRI Factors')

    plt.tight_layout()

    if output_dir:
        plots_dir = os.path.join(output_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        save_path = os.path.join(plots_dir, 'bptf_diagnostics.pdf')
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"BPTF diagnostic plots saved to: {plots_dir}")

    return fig
