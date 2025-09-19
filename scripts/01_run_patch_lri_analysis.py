#!/usr/bin/env python3
"""
01 - Patch-based Ligand-Receptor Interaction Analysis Script

This script runs the patch-based LRI analysis for matrix factorization approaches.
It divides tissue into 50μm patches and counts all-to-all interactions within each patch.
"""

import sys
import os
sys.path.append(os.path.dirname(__file__))

import anndata
import argparse
from pathlib import Path
from os.path import join as pathjoin

import scipy.sparse as sp

from patch_lri_analysis import PatchLRIAnalyzer


def binomial_thinning(sparse_matrix, probs):
    """Perform binomial thinning on a sparse matrix."""
    sparse_matrix = sparse_matrix.tocoo(copy=False)
    data_thinned = np.random.binomial(n=sparse_matrix.data, p=probs[sparse_matrix.col])
    thinned_matrix = sp.coo_matrix((data_thinned, (sparse_matrix.row, sparse_matrix.col)),
                                    shape=sparse_matrix.shape)
    return thinned_matrix.tocsr()


def main():
    """Main patch-based LRI analysis pipeline."""
    parser = argparse.ArgumentParser(description='Patch-based ligand-receptor interaction analysis')
    parser.add_argument('--data-file', default='data/processed/preprocessed_xenium_data_subset.h5ad',
                       help='Processed data file')
    parser.add_argument('--output-dir', default='results/bptf_patch_lri_subset',
                       help='Output directory for results')
    parser.add_argument('--patch-size', type=float, default=50.0,
                       help='Size of spatial patches in micrometers')
    parser.add_argument('--resource', default='cellchatdb',
                       help='LRI database resource name')
    parser.add_argument('--use-batch-processing', action='store_true',
                       help='Use batch processing for large datasets')
    parser.add_argument('--batch-size', type=int, default=500,
                       help='Batch size for processing patches')
    parser.add_argument('--data-thinning', type=float, default=1.0,
                       help='Thin data to which quantile (or not, set as 1)')
    parser.add_argument('--random-state', type=int, default=123,
                       help='Random seed for reproducibility')
    parser.add_argument('--spliter', type=str,
                       default='|', help='cell-gene or cell|gene ...')
    
    args = parser.parse_args()
    
    print("="*60)
    print("01 - PATCH-BASED LIGAND-RECEPTOR INTERACTION ANALYSIS")
    print("="*60)
    print(f"Data file: {args.data_file}")
    print(f"Output directory: {args.output_dir}")
    print(f"Patch size: {args.patch_size} μm")
    print(f"LRI database: {args.resource}")
    print("="*60)
    
    # Check if data file exists
    if not os.path.exists(args.data_file):
        print(f"Error: Data file not found: {args.data_file}")
        return
    
    # Load data
    print("Loading spatial transcriptomics data...")
    adata = anndata.read_h5ad(args.data_file)
    print(f"Data shape: {adata.shape}")
    print(f"Cell types: {adata.obs['cell_type'].cat.categories.tolist()}")
    print(f"TMA IDs: {sorted(adata.obs['tma_id'].unique())}")
    
    # Initialize analyzer
    analyzer = PatchLRIAnalyzer(
        patch_size=args.patch_size,
        resource_name=args.resource,
        spliter=args.spliter
    )
    
    # Run analysis
    results = analyzer.run_analysis(
        adata, 
        args.output_dir, 
        use_batch_processing=args.use_batch_processing,
        batch_size=args.batch_size
    )

    np.random.seed(args.random_state)

    if args.data_thinning != 1.0:
        # Step 1: calculate the average LRI count for each column
        column_means = np.array(results['patch_lri_matrix'].mean(axis=0)).flatten()
        # Step 2: calculate the quantile t of all the means as the threshold
        t = int(args.data_thinning * 100)
        m_t = np.percentile(column_means, t)
        # Step 3: calculate the thinning probability of each LRI
        p_j = np.minimum(1, m_t / column_means)
        patch_lri_matrix_thinned = binomial_thinning(results['patch_lri_matrix'], p_j)
        # Save the thinned matrix
        out_path = pathjoin(args.output_dir, f"patch_lri_matrix_{t}.npz")
        sp.save_npz(out_path, patch_lri_matrix_thinned)
    
    # Print summary
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY")
    print("="*60)
    print(f"Total patches: {results['patch_info']['n_patches']}")
    print(f"Total LRI combinations: {len(results['column_names'])}")
    print(f"Matrix shape: {results['patch_lri_matrix'].shape}")
    print(f"Matrix sparsity: {(1 - results['patch_lri_matrix'].nnz / np.prod(results['patch_lri_matrix'].shape)) * 100:.2f}%")
    print(f"Non-zero interactions: {results['patch_lri_matrix'].nnz:,}")

    if args.data_thinning != 1.0:
        print(f"Thinned matrix sparsity: {(1 - results['patch_lri_matrix'].nnz / np.prod(results['patch_lri_matrix'].shape)) * 100:.2f}%")
    
    # Cell type distribution
    print(f"\nCell type distribution:")
    cell_type_counts = results['cell_patch_df']['cell_type'].value_counts()
    for cell_type, count in cell_type_counts.items():
        print(f"  {cell_type}: {count:,} cells")
    
    # TMA distribution
    print(f"\nTMA distribution:")
    tma_counts = results['patch_tma_df']['tma_id'].value_counts()
    for tma_id, count in tma_counts.items():
        print(f"  TMA {tma_id}: {count} patches")
    
    print("\n" + "="*60)
    print("STEP 01 COMPLETED SUCCESSFULLY!")
    print("="*60)
    print(f"Results saved to: {args.output_dir}")
    print("\nFiles created:")
    print(f"  - patch_lri_matrix.npz (sparse matrix)")
    if args.data_thinning != 1.0:
        print(f"  - patch_lri_matrix_{t}.npz (thinned sparse matrix)")
    print(f"  - patch_lri_columns.csv (column names)")
    print(f"  - patch_tma_correspondence.csv (patch-TMA mapping)")
    print(f"  - cell_patch_correspondence.csv (cell-patch mapping)")
    print(f"  - analysis_parameters.csv (analysis parameters)")
    
    print("\nNext step:")
    print("Run: python scripts/02_bptf_matrix_factorization.py")


if __name__ == '__main__':
    import numpy as np
    main()