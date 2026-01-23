#!/usr/bin/env python3
"""
02 - BPTF Matrix Factorization for Spatial LRI Analysis

Simple implementation using bptf package with detailed output files.
"""

import sys
import os
sys.path.append(os.path.dirname(__file__))

import numpy as np
import pandas as pd
import scipy.sparse as sp
import sparse  # for sparse.COO conversion
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import time
import pickle

# Import BPTF - this should work if bptf is installed
try:
    from bptf import BPTF, save_bptf, load_bptf
    BPTF_AVAILABLE = True
except ImportError:
    BPTF_AVAILABLE = False
    warnings.warn("BPTF not available. Install from: https://github.com/aschein/bptf")

from patch_lri_analysis import load_patch_lri_results


def run_bptf(mat, n_components=50, max_iter=100, verbose=True, random_state=0):
    """Use bptf package for Poisson Matrix Factorization"""
    if isinstance(mat, sp.spmatrix):
        data = sparse.COO(mat)      # Convert scipy sparse to sparse.COO
    else:
        data = mat                  # numpy array or already COO
    
    np.random.seed(random_state)

    model = BPTF(data_shape=data.shape, n_components=n_components)
    model.fit(data, max_iter=max_iter, verbose=verbose)
    return model


def get_top_motifs(patch_loadings, top_k=10):
    """Get top motifs by total activity"""
    motif_activities = patch_loadings.sum(axis=0)
    top_indices = np.argsort(motif_activities)[::-1][:top_k]
    
    return {
        'motif_indices': top_indices.tolist(),
        'activities': motif_activities[top_indices].tolist(),
        'activity_fractions': (motif_activities[top_indices] / motif_activities.sum()).tolist()
    }


def get_top_lris_per_motif(lri_factors, motif_idx, top_k=20):
    """Get top LRIs for a specific motif"""
    if motif_idx >= lri_factors.shape[0]:
        raise ValueError(f"motif index {motif_idx} out of range")
    
    factors = lri_factors[motif_idx, :]
    top_indices = np.argsort(factors)[::-1][:top_k]
    return top_indices.tolist()


def save_bptf_results(model, patch_loadings, lri_factors, column_names, 
                     patch_metadata_df, output_dir, elbo_hist=None, delta_hist=None):
    """Save all BPTF results"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Save model
    model_path = Path(output_dir) / 'bptf_model.npz'
    save_bptf(model, model_path)
    print(f"BPTF model saved to: {model_path}")
    
    # Save factor matrices as numpy arrays
    np.save(os.path.join(output_dir, 'patch_loadings.npy'), patch_loadings)
    np.save(os.path.join(output_dir, 'lri_factors.npy'), lri_factors)
    
    # Save convergence history if available
    if elbo_hist is not None and delta_hist is not None:
        history_df = pd.DataFrame({
            'iteration': range(len(elbo_hist)),
            'elbo': elbo_hist,
            'delta': delta_hist
        })
        history_df.to_csv(os.path.join(output_dir, 'iteration_history.csv'), index=False)
    
    # Analyze patch-motif relationships
    patch_motifs = []
    if 'patch_id' in patch_metadata_df.columns:
        patch_ids = sorted(patch_metadata_df['patch_id'].astype(str).unique().tolist())
    else:
        patch_ids = [f"patch_{i}" for i in range(len(patch_loadings))]
    
    for i, patch_id in enumerate(patch_ids):
        total_loading = patch_loadings[i].sum()
        for k in range(patch_loadings.shape[1]):
            loading = patch_loadings[i, k]
            fraction = loading / total_loading if total_loading > 0 else 0
            
            patch_motifs.append({
                'patch_idx': i,
                'patch_id': patch_id,
                'motif': k,
                'loading': loading,
                'loading_normalized': fraction
            })
    
    # Analyze LRI-motif relationships  
    lri_motifs = []
    for j, column_name in enumerate(column_names):
        motif_factors = lri_factors[:, j]  # factors for this LRI across motifs
        total_factor = motif_factors.sum()
        
        for k in range(len(motif_factors)):
            factor = motif_factors[k]
            fraction = factor / total_factor if total_factor > 0 else 0
            
            lri_motifs.append({
                'lri_idx': j,
                'motif_idx': k,
                'lri_name': column_name,
                'factor': factor,
                'factor_normalized': fraction
            })
    
    # Save analysis results as CSV (matching original format)
    patch_motifs_df = pd.DataFrame(patch_motifs)
    patch_motifs_df.to_csv(os.path.join(output_dir, 'patch_motifs.csv'), index=False)
    
    lri_motifs_df = pd.DataFrame(lri_motifs)
    lri_motifs_df.to_csv(os.path.join(output_dir, 'lri_motifs.csv'), index=False)
    
    # Save motif summary (matching original format)
    motif_activities = patch_loadings.sum(axis=0)
    top_motifs = get_top_motifs(patch_loadings, top_k=len(motif_activities))
    
    top_motifs_df = pd.DataFrame({
        'motif_idx': top_motifs['motif_indices'],
        'total_activity': top_motifs['activities'],
        'activity_fraction': top_motifs['activity_fractions']
    })
    top_motifs_df.to_csv(os.path.join(output_dir, 'top_motifs.csv'), index=False)
    
    # Save model parameters (matching original format)
    params_df = pd.DataFrame({
        'parameter': ['n_components', 'n_patches', 'n_lris', 'method'],
        'value': [
            model.n_components,
            patch_loadings.shape[0],
            lri_factors.shape[1],
            'BPTF'
        ]
    })
    params_df.to_csv(os.path.join(output_dir, 'factorization_parameters.csv'), index=False)
    
    print(f"All results saved to: {output_dir}")


def plot_bptf_diagnostics(patch_loadings, lri_factors, output_dir):
    """Create diagnostic plots"""
    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. motif activity distribution
    motif_activities = patch_loadings.sum(axis=0)
    axes[0, 0].bar(range(len(motif_activities)), motif_activities)
    axes[0, 0].set_title('BPTF motif Activity Distribution')
    axes[0, 0].set_xlabel('motif Index')
    axes[0, 0].set_ylabel('Total Activity')
    
    # 2. LRI participation distribution
    lri_participation = lri_factors.sum(axis=0)
    axes[0, 1].hist(lri_participation, bins=50, alpha=0.7, edgecolor='black')
    axes[0, 1].set_title('LRI Participation Distribution')
    axes[0, 1].set_xlabel('Total Participation')
    axes[0, 1].set_ylabel('Number of LRIs (log)')
    axes[0, 1].set_yscale('log')
    
    # 3. Patch loading distribution
    patch_loading_totals = patch_loadings.sum(axis=1)
    axes[1, 0].hist(patch_loading_totals, bins=50, alpha=0.7, edgecolor='black')
    axes[1, 0].set_title('Patch Loading Distribution')
    axes[1, 0].set_xlabel('Total Loading')
    axes[1, 0].set_ylabel('Number of Patches')
    
    # 4. Factor matrix sparsity
    sparsity_threshold = 1e-6
    patch_sparsity = (patch_loadings < sparsity_threshold).mean()
    lri_sparsity = (lri_factors < sparsity_threshold).mean()
    axes[1, 1].bar(['Patch Loadings', 'LRI Factors'], [patch_sparsity, lri_sparsity])
    axes[1, 1].set_title('Factor Matrix Sparsity')
    axes[1, 1].set_ylabel(f'Fraction of < {sparsity_threshold}')
    axes[1, 1].set_ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'bptf_diagnostics.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"BPTF diagnostic plots saved to: {plots_dir}")


def main():
    """Main BPTF matrix factorization pipeline."""
    parser = argparse.ArgumentParser(description='BPTF matrix factorization for spatial LRI analysis')
    parser.add_argument('--input-dir', default='results/GBM_cellphone',
                       help='Input directory with patch LRI results')
    parser.add_argument('--sparse-matrix-name', default='patch_lri_matrix.npz',
                       help='Name of the sparse matrix file in input directory')
    parser.add_argument('--column-df-name', default='patch_lri_columns.csv',
                       help='Name of the column names file in input directory')
    parser.add_argument('--meta-df-name', default='cell_patch_correspondence.csv',
                       help='Name of the metadata file in input directory')
    parser.add_argument('--output-dir', default='results/GBM_cellphone/bptf',
                       help='Output directory for BPTF results')
    parser.add_argument('--n-components', type=int, default=15,
                       help='Number of latent factors/motifs')
    parser.add_argument('--max-iter', type=int, default=10000,
                       help='Maximum number of iterations')
    parser.add_argument('--random-state', type=int, default=1,
                       help='Random seed for reproducibility')
    parser.add_argument('--spliter', type=str,
                       default='|', help='cell-gene or cell|gene ...')
    parser.add_argument('--neighborhood', type=bool,
                       default=False, help='if neighbood-based or patch-based')
    parser.add_argument('--single-cell', type=bool,
                       default=False, help='if cell-based or patch-based')
    
    # parser.add_argument('--top-k-motifs', type=int, default=10,
    #                    help='Number of top motifs to analyze in detail')
    # parser.add_argument('--top-k-lris', type=int, default=20,
    #                    help='Number of top LRIs per motif to save')
    
    args = parser.parse_args()
    
    print("="*60)
    print("02 - BPTF MATRIX FACTORIZATION")
    print("="*60)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Number of components: {args.n_components}")
    print(f"Max iterations: {args.max_iter}")
    print(f"Random state: {args.random_state}")
    print("="*60)
    
    # Check dependencies
    if not BPTF_AVAILABLE:
        print("Error: BPTF not available!")
        print("Please install from: https://github.com/aschein/bptf")
        print("Or try: pip install git+https://github.com/aschein/bptf.git")
        return
    
    # Check input directory
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory not found: {args.input_dir}")
        print("Please run step 01 first: python scripts/01_run_patch_lri_analysis.py")
        return
    
    # Load patch LRI results
    print("Loading patch-LRI results...")
    try:
        results = load_patch_lri_results(args.input_dir, sparse_matrix_name=args.sparse_matrix_name,
                                         column_df_name = args.column_df_name, meta_df_name = args.meta_df_name,
                                         neighborhood=args.neighborhood, single_cell = args.single_cell)
        mat = results['patch_lri_matrix']
        cols = results['column_names']
        cell_meta_df = results['cell_patch_df']
        if not args.single_cell:
            cell_patch_df = results['patch_tma_df']
        
        print(f"Loaded matrix shape: {mat.shape}")
        print(f"Matrix sparsity: {(1 - mat.nnz / np.prod(mat.shape)) * 100:.2f}%")
        print(f"Number of LRIs: {len(cols)}")
        
    except Exception as e:
        print(f"Error loading patch LRI results: {e}")
        return
    
    # Run BPTF
    print(f"\nRunning BPTF with {args.n_components} components...")
    print("Using simple mode...")
    
    start_time = time.time()
    
    try:
        model = run_bptf(
            mat,
            n_components=args.n_components,
            max_iter=args.max_iter,
            verbose=True,
            random_state=args.random_state
        )
        elbo_hist, delta_hist = None, None
        
        fit_time = time.time() - start_time
        print(f"BPTF completed in {fit_time:.1f} seconds!")
        
    except Exception as e:
        print(f"Error during BPTF fitting: {e}")
        return
    
    # Extract factor matrices
    print("\nExtracting factor matrices...")
    G_DK_M = model.G_DK_M
    patch_loadings = G_DK_M[0]  # patches × motifs
    lri_factors = G_DK_M[1].T   # motifs × LRIs (transposed)
    
    print(f"Patch loadings shape: {patch_loadings.shape}")
    print(f"LRI factors shape: {lri_factors.shape}")
    
    # Save results
    print("\nSaving BPTF results...")
    try:
        save_bptf_results(model, patch_loadings, lri_factors, cols, 
                         cell_meta_df, args.output_dir, elbo_hist, delta_hist)
        
        # Create plots
        print("\nCreating visualizations...")
        plot_bptf_diagnostics(patch_loadings, lri_factors, args.output_dir)
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        return
    
    # Print summary
    print("\n" + "="*60)
    print("BPTF SUMMARY")
    print("="*60)
    
    motif_activities = patch_loadings.sum(axis=0)
    top_motifs = get_top_motifs(patch_loadings, top_k=5)
    
    print(f"Total motifs: {len(motif_activities)}")
    if elbo_hist is not None:
        print(f"Final ELBO: {elbo_hist[-1]:.2f}")
        print(f"Converged after {len(elbo_hist)-1} iterations")
    
    print(f"\nTop 5 motifs by activity:")
    for i, (prog_idx, activity, fraction) in enumerate(zip(
        top_motifs['motif_indices'][:5],
        top_motifs['activities'][:5],
        top_motifs['activity_fractions'][:5]
    )):
        print(f"  {i+1}. motif {prog_idx}: {fraction:.1%} of total activity")
    
    print("\n" + "="*60)
    print("BPTF STEP 02 COMPLETED SUCCESSFULLY!")
    print("="*60)
    print(f"Results saved to: {args.output_dir}")
    print("\nFiles created:")
    print(f"  - bptf_model.npz (BPTF model)")
    print(f"  - patch_loadings.npy, lri_factors.npy (factor matrices)")
    print(f"  - patch_motifs.csv, lri_motifs.csv (detailed analysis)")
    print(f"  - bptf_motif_interpretations.csv (motif summaries)")
    print(f"  - bptf_motif_{{X}}_top_lris.csv (per-motif LRI lists)")
    print(f"  - iteration_history.csv (convergence tracking)")
    print(f"  - plots/ (diagnostic visualizations)")
    
    print("\nNext step:")
    print("Run: python scripts/03_poisson_glm.py")


if __name__ == '__main__':
    main()