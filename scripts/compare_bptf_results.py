#!/usr/bin/env python
"""
Script to compare two BPTF results by visualizing their similarity matrices.
Generates two plots: one showing cosine similarity and another with Hungarian alignment.
"""

import pandas as pd
import anndata as ad
import argparse
import os
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import rankdata
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


def row_cosine(A, B):
    """
    Compute row-wise cosine similarity between two matrices.
    
    Args:
        A: Matrix of shape (n, d)
        B: Matrix of shape (m, d)
    
    Returns:
        Cosine similarity matrix of shape (n, m)
    """
    # Normalize rows to unit length
    A_norm = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
    B_norm = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-12)
    
    # Compute cosine similarity
    cos_sim = A_norm @ B_norm.T
    
    return cos_sim


def load_bptf_results(results_dir):
    """Load BPTF results from directory"""
    print(f"Loading BPTF results from: {results_dir}")
    
    # Load factor matrices
    patch_loadings = np.load(os.path.join(results_dir, 'patch_loadings.npy'))
    lri_factors = np.load(os.path.join(results_dir, 'lri_factors.npy'))
    
    # Load detailed analysis
    patch_motifs = pd.read_csv(os.path.join(results_dir, 'patch_motifs.csv'))
    lri_motifs = pd.read_csv(os.path.join(results_dir, 'lri_motifs.csv'))
    
    print(f"Loaded results:")
    print(f"  - Patch loadings: {patch_loadings.shape}")
    print(f"  - LRI factors: {lri_factors.shape}")
    print(f"  - Patch motifs: {patch_motifs.shape}")
    print(f"  - LRI motifs: {lri_motifs.shape}")
    
    return {
        'patch_loadings': patch_loadings,
        'lri_factors': lri_factors,
        'patch_motifs': patch_motifs,
        'lri_motifs': lri_motifs
    }


def subset_lri_factors(result1, result2):
    """
    Subset LRI factors from result1 to match those present in result2.
    
    Args:
        result1: First BPTF result dictionary
        result2: Second BPTF result dictionary
    
    Returns:
        Subset of result1's lri_factors matching LRIs in result2
    """
    # Get LRI names from result2
    result2_lri_names = set(result2['lri_motifs']['lri_name'].unique())
    
    # Subset result1's LRI motifs
    result1_subset = result1['lri_motifs'][
        result1['lri_motifs']['lri_name'].isin(result2_lri_names)
    ]
    
    print(f"Subsetting LRI factors: {result1['lri_factors'].shape} -> ", end="")
    
    # Get indices to keep
    kept_idx = result1_subset['lri_idx'].unique()
    
    # Subset lri_factors
    result1_subset_factors = result1['lri_factors'][:, kept_idx]
    print(f"{result1_subset_factors.shape}")
    
    return result1_subset_factors

def align_lri_factors_consensus(result1, result2, key='lri_name'):
    """
    Align two BPTF results to the consensus set of LRIs (by name),
    and return both lri_factors matrices subsetted and column-aligned.

    Args:
        result1: dict with keys ['lri_motifs', 'lri_factors']
        result2: dict with keys ['lri_motifs', 'lri_factors']
        key:     column name in lri_motifs used to match LRIs (default: 'lri_name')

    Returns:
        f1_aligned: result1['lri_factors'] subsetted to consensus, columns aligned
        f2_aligned: result2['lri_factors'] subsetted to consensus, columns aligned
        common_names: list of LRI names in the aligned order
        idx1: numpy array of column indices taken from result1
        idx2: numpy array of column indices taken from result2
    """
    import numpy as np
    import pandas as pd

    # --- build name -> idx maps (deduplicate by name, keep first) ---
    df1 = (result1['lri_motifs'][[key, 'lri_idx']]
           .drop_duplicates(subset=[key], keep='first'))
    df2 = (result2['lri_motifs'][[key, 'lri_idx']]
           .drop_duplicates(subset=[key], keep='first'))

    names1 = set(df1[key].astype(str))
    names2 = set(df2[key].astype(str))
    common_names = sorted(names1 & names2)

    if len(common_names) == 0:
        raise ValueError("No consensus LRIs found between result1 and result2.")

    # index vectors in the SAME order of common_names
    map1 = df1.set_index(key)['lri_idx']
    map2 = df2.set_index(key)['lri_idx']
    idx1 = map1.loc[common_names].to_numpy()
    idx2 = map2.loc[common_names].to_numpy()

    # subset & align columns
    f1 = result1['lri_factors']
    f2 = result2['lri_factors']

    # basic sanity checks
    if f1.ndim != 2 or f2.ndim != 2:
        raise ValueError("Expected 2D arrays for lri_factors in both results.")
    if np.max(idx1) >= f1.shape[1] or np.max(idx2) >= f2.shape[1]:
        raise IndexError("lri_idx exceeds number of columns in lri_factors.")

    print("[Consensus LRI alignment]")
    print(f"  result1 lri_factors: {f1.shape} -> ({f1.shape[0]}, {len(idx1)})")
    print(f"  result2 lri_factors: {f2.shape} -> ({f2.shape[0]}, {len(idx2)})")
    print(f"  #consensus LRIs: {len(common_names)}")

    f1_aligned = f1[:, idx1]
    f2_aligned = f2[:, idx2]

    return f1_aligned, f2_aligned, common_names, idx1, idx2



def plot_cosine_similarity(matrix1, matrix2, name1, name2, save_path=None, annotate=True, fmt="{:.2f}", fs=8):
    """
    Plot cosine similarity heatmap between two matrices.
    
    Args:
        matrix1: First matrix
        matrix2: Second matrix
        name1: Name for first result
        name2: Name for second result
        save_path: Path to save figure (optional)
        annotate: Whether to annotate cells with values
        fmt: Format string for annotations
        fs: Font size for annotations
    """
    # Compute cosine similarity
    C_cos = row_cosine(matrix1, matrix2)
    
    # Plot
    vmin, vmax = np.nanmin(C_cos), np.nanmax(C_cos)
    if vmin == vmax: 
        vmax = vmin + 1e-12
    
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(C_cos, cmap="Reds", vmin=vmin, vmax=vmax)
    
    ax.set_title(f"Cosine similarity: {name1} vs {name2}")
    ax.set_xlabel(f"{name2} motifs")
    ax.set_ylabel(f"{name1} motifs")
    
    ax.set_xticks(np.arange(C_cos.shape[1]))
    ax.set_yticks(np.arange(C_cos.shape[0]))
    ax.set_xticklabels([f"motif {j}" for j in range(C_cos.shape[1])], rotation=45, ha="right")
    ax.set_yticklabels([f"motif {i}" for i in range(C_cos.shape[0])])
    
    if annotate:
        for i in range(C_cos.shape[0]):
            for j in range(C_cos.shape[1]):
                val = C_cos[i, j]
                color = "white" if (val - vmin)/(vmax - vmin + 1e-12) > 0.6 else "black"
                if val >= 0.4:
                    ax.text(j, i, fmt.format(val), ha="center", va="center", fontsize=fs, color=color)
    
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("Correlation")
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure: {save_path}")
    
    plt.show()
    
    return C_cos


def plot_hungarian_alignment(matrix1, matrix2, name1, name2, save_path=None):
    """
    Plot similarity matrix with Hungarian algorithm alignment.
    
    Args:
        matrix1: First matrix
        matrix2: Second matrix
        name1: Name for first result
        name2: Name for second result
        save_path: Path to save figure (optional)
    """
    # Compute similarity
    C_cos = row_cosine(matrix1, matrix2)
    D = 1 - C_cos  # Convert to distance
    
    # Hungarian matching
    row_ind, col_ind = linear_sum_assignment(D, maximize=False)
    
    # Reorder both rows and columns
    C_cos_reordered = C_cos[np.ix_(row_ind, col_ind)]
    
    # Labels for reordered motifs
    yticklabels = [f"motif {i}" for i in row_ind]
    xticklabels = [f"motif {j}" for j in col_ind]
    
    # Plot
    plt.figure(figsize=(7, 6))
    ax = sns.heatmap(C_cos_reordered, cmap="Reds", cbar=True,
                     xticklabels=xticklabels, yticklabels=yticklabels,
                     cbar_kws={'label': 'Cosine Similarity'})
    
    # Draw black boxes on diagonal for matched pairs
    for i in range(min(C_cos_reordered.shape)):
        ax.add_patch(plt.Rectangle((i, i), 1, 1, fill=False, ec="black", lw=2))
    
    ax.set_xlabel(f"{name2} motifs (reordered)")
    ax.set_ylabel(f"{name1} motifs")
    ax.set_title(f"Hungarian alignment: {name1} vs {name2}")
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure: {save_path}")
    
    plt.show()
    
    return C_cos_reordered, row_ind, col_ind


def main():
    parser = argparse.ArgumentParser(
        description='Compare two BPTF results and generate similarity plots.'
    )
    
    parser.add_argument('--dir1', type=str, required=True,
                        help='Directory containing first BPTF results')
    parser.add_argument('--dir2', type=str, required=True,
                        help='Directory containing second BPTF results')
    parser.add_argument('--name1', type=str, required=True,
                        help='Name for first result (for plot labels)')
    parser.add_argument('--name2', type=str, required=True,
                        help='Name for second result (for plot labels)')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Directory to save output figures')
    parser.add_argument('--subset-lri', action='store_true',
                        help='Subset LRIs from first result to match second result')
    parser.add_argument('--annotate', action='store_true', default=True,
                        help='Annotate heatmap cells with values')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load BPTF results
    print(f"\n{'='*60}")
    print(f"Loading BPTF results...")
    print(f"{'='*60}")
    
    result1 = load_bptf_results(args.dir1)
    result2 = load_bptf_results(args.dir2)
    
    # Get LRI factors
    matrix1 = result1['lri_factors']
    matrix2 = result2['lri_factors']
    
    # Optionally subset LRIs
    if args.subset_lri:
        print(f"\n{'='*60}")
        print(f"Subsetting LRIs from {args.name1} to match {args.name2}...")
        print(f"{'='*60}")
        # matrix1 = subset_lri_factors(result1, result2)
        # matrix2 = subset_lri_factors(result2, result1)
        matrix1, matrix2, lri_names, idx1, idx2 = align_lri_factors_consensus(result1, result2)

    
    # Generate plots
    print(f"\n{'='*60}")
    print(f"Generating plots...")
    print(f"{'='*60}\n")
    
    # Plot 1: Cosine similarity heatmap
    cosine_path = os.path.join(args.output_dir, f'cosine_similarity_{args.name1}_vs_{args.name2}.png')
    C_cos = plot_cosine_similarity(
        matrix1, matrix2, 
        args.name1, args.name2, 
        save_path=cosine_path,
        annotate=args.annotate
    )
    
    # Plot 2: Hungarian alignment
    hungarian_path = os.path.join(args.output_dir, f'hungarian_alignment_{args.name1}_vs_{args.name2}.png')
    C_reordered, row_ind, col_ind = plot_hungarian_alignment(
        matrix1, matrix2,
        args.name1, args.name2,
        save_path=hungarian_path
    )
    
    # Save matching information
    matching_path = os.path.join(args.output_dir, f'matching_{args.name1}_vs_{args.name2}.csv')
    matching_df = pd.DataFrame({
        f'{args.name1}_motif': row_ind,
        f'{args.name2}_motif': col_ind,
        'similarity': [C_cos[i, j] for i, j in zip(row_ind, col_ind)]
    })
    matching_df.to_csv(matching_path, index=False)
    print(f"\nSaved matching information: {matching_path}")
    
    print(f"\n{'='*60}")
    print(f"Analysis complete!")
    print(f"{'='*60}")
    print(f"Output saved to: {args.output_dir}")
    print(f"  - Cosine similarity plot: {os.path.basename(cosine_path)}")
    print(f"  - Hungarian alignment plot: {os.path.basename(hungarian_path)}")
    print(f"  - Matching information: {os.path.basename(matching_path)}")


if __name__ == "__main__":
    main()