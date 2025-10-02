#!/usr/bin/env python3
"""
Cell Neighborhood-based Ligand-Receptor Interaction Analysis Script

This script runs neighborhood-based LRI analysis where each cell is the center
of a square neighborhood and interactions are counted within that neighborhood.
"""

import sys
import os
sys.path.append(os.path.dirname(__file__))

import anndata
import argparse
import numpy as np
import scipy.sparse as sp
from pathlib import Path
from os.path import join as pathjoin

from neighborhood_lri_analysis import NeighborhoodLRIAnalyzer


def binomial_thinning(sparse_matrix, probs):
    """Perform binomial thinning on a sparse matrix."""
    sparse_matrix = sparse_matrix.tocoo(copy=False)
    data_thinned = np.random.binomial(n=sparse_matrix.data, p=probs[sparse_matrix.col])
    thinned_matrix = sp.coo_matrix((data_thinned, (sparse_matrix.row, sparse_matrix.col)),
                                    shape=sparse_matrix.shape)
    return thinned_matrix.tocsr()


def main():
    """Main neighborhood-based LRI analysis pipeline."""
    parser = argparse.ArgumentParser(description='Cell neighborhood-based ligand-receptor interaction analysis')
    parser.add_argument('--data-file', default='/Users/jiayifan/Desktop/Lab/TMA_punch_subfiles/xenium_mm_final_cell_id.h5ad',
                       help='Processed data file')
    parser.add_argument('--output-dir', default='results/single_cell_neighborhood_lri',
                       help='Output directory for results')
    parser.add_argument('--neighborhood-size', type=float, default=50.0,
                       help='Size of square neighborhood in micrometers (edge length)')
    parser.add_argument('--resource', default='cellchatdb',
                       help='LRI database resource name')
    parser.add_argument('--data-thinning', type=float, default=1.0,
                       help='Thin data to which quantile (or not, set as 1)')
    parser.add_argument('--random-state', type=int, default=123,
                       help='Random seed for reproducibility')
    parser.add_argument('--spliter', type=str, default='|',
                       help='Separator for column names (e.g., cell|gene)')
    
    args = parser.parse_args()
    
    print("="*70)
    print("CELL NEIGHBORHOOD-BASED LIGAND-RECEPTOR INTERACTION ANALYSIS")
    print("="*70)
    print(f"Data file: {args.data_file}")
    print(f"Output directory: {args.output_dir}")
    print(f"Neighborhood size: {args.neighborhood_size} µm (square)")
    print(f"LRI database: {args.resource}")
    print("="*70)
    
    # Check if data file exists
    if not os.path.exists(args.data_file):
        print(f"Error: Data file not found: {args.data_file}")
        return
    
    # Load data
    print("\nLoading spatial transcriptomics data...")
    adata = anndata.read_h5ad(args.data_file)
    print(f"Data shape: {adata.shape}")
    print(f"Cell types: {adata.obs['cell_type'].cat.categories.tolist()}")
    print(f"TMA IDs: {sorted(adata.obs['tma_id'].unique())}")
    
    # Initialize analyzer
    analyzer = NeighborhoodLRIAnalyzer(
        neighborhood_size=args.neighborhood_size,
        resource_name=args.resource,
        spliter=args.spliter
    )
    
    # Run analysis
    results = analyzer.run_analysis(adata, args.output_dir)
    
    # Set random seed
    np.random.seed(args.random_state)
    
    # Optional data thinning
    if args.data_thinning != 1.0:
        print(f"\nApplying data thinning (quantile={args.data_thinning})...")
        
        # Calculate column means
        column_means = np.array(results['cell_lri_matrix'].mean(axis=0)).flatten()
        
        # Calculate quantile threshold
        t = int(args.data_thinning * 100)
        m_t = np.percentile(column_means, t)
        
        # Calculate thinning probabilities
        p_j = np.minimum(1, m_t / column_means)
        
        # Apply binomial thinning
        cell_lri_matrix_thinned = binomial_thinning(results['cell_lri_matrix'], p_j)
        
        # Save thinned matrix
        out_path = pathjoin(args.output_dir, f"cell_lri_matrix_{t}.npz")
        sp.save_npz(out_path, cell_lri_matrix_thinned)
        print(f"Thinned matrix saved to: {out_path}")
    
    # Print summary
    print("\n" + "="*70)
    print("ANALYSIS SUMMARY")
    print("="*70)
    print(f"Total cells: {adata.n_obs:,}")
    print(f"Total LRI combinations: {len(results['column_names']):,}")
    print(f"Matrix shape: {results['cell_lri_matrix'].shape}")
    print(f"Matrix sparsity: {(1 - results['cell_lri_matrix'].nnz / np.prod(results['cell_lri_matrix'].shape)) * 100:.2f}%")
    print(f"Non-zero interactions: {results['cell_lri_matrix'].nnz:,}")
    print(f"Average neighborhood size: {results['cell_metadata_df']['neighborhood_size'].mean():.1f} cells")
    print(f"Min neighborhood size: {results['cell_metadata_df']['neighborhood_size'].min()}")
    print(f"Max neighborhood size: {results['cell_metadata_df']['neighborhood_size'].max()}")
    
    if args.data_thinning != 1.0:
        print(f"\nThinned matrix sparsity: {(1 - cell_lri_matrix_thinned.nnz / np.prod(cell_lri_matrix_thinned.shape)) * 100:.2f}%")
        print(f"Thinned non-zero interactions: {cell_lri_matrix_thinned.nnz:,}")
    
    # Cell type distribution
    print(f"\nCell type distribution:")
    cell_type_counts = results['cell_metadata_df']['cell_type'].value_counts()
    for cell_type, count in cell_type_counts.items():
        print(f"  {cell_type}: {count:,} cells")
    
    # TMA distribution
    print(f"\nTMA distribution:")
    tma_counts = results['cell_metadata_df']['tma_id'].value_counts()
    for tma_id, count in tma_counts.items():
        print(f"  TMA {tma_id}: {count:,} cells")
    
    # Neighborhood size distribution by cell type
    print(f"\nAverage neighborhood size by cell type:")
    for cell_type in cell_type_counts.index:
        avg_size = results['cell_metadata_df'][
            results['cell_metadata_df']['cell_type'] == cell_type
        ]['neighborhood_size'].mean()
        print(f"  {cell_type}: {avg_size:.1f} cells")
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETED SUCCESSFULLY!")
    print("="*70)
    print(f"Results saved to: {args.output_dir}")
    print("\nFiles created:")
    print(f"  - cell_lri_matrix.npz (sparse matrix)")
    if args.data_thinning != 1.0:
        print(f"  - cell_lri_matrix_{t}.npz (thinned sparse matrix)")
    print(f"  - cell_lri_columns.csv (column names)")
    print(f"  - cell_metadata.csv (cell metadata with coordinates)")
    print(f"  - analysis_parameters.csv (analysis parameters)")
    
    print("\nNext steps:")
    print("  - Use the cell_lri_matrix for downstream analysis")
    print("  - Compare with patch-based results if available")


if __name__ == '__main__':
    main()