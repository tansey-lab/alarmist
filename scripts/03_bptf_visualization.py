#!/usr/bin/env python3
"""
BPTF Visualization Script

Generate comprehensive visualizations of BPTF analysis results including:
- Factor distributions and sparsity
- Motif activity rankings  
- LRI clustermap and communication heatmaps
- Top LRI interactions per motif
- LRI network graphs
- Spatial motif distributions

Usage:
    python scripts/03_bptf_visualization.py --bptf-dir results/bptf --patch-dir results/patch_lri --data-file data.h5ad --output-dir results/plots
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Import visualization utilities
from bptf_visualization_utils import *


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
    
    return {
        'patch_loadings': patch_loadings,
        'lri_factors': lri_factors,
        'patch_motifs': patch_motifs,
        'lri_motifs': lri_motifs
    }


def load_patch_lri_results(patch_dir, sparse_matrix_name='patch_lri_matrix.npz'):
    """
    Load patch-LRI results from directory
    
    Expected files:
    - patch_lri_matrix.npz: Sparse matrix
    - patch_lri_columns.csv: Column names with 'column_name' column
    - patch_tma_correspondence.csv: Patch metadata
    """
    print(f"Loading patch-LRI results from: {patch_dir}")
    
    # Load sparse matrix
    import scipy.sparse as sp
    mat_path = os.path.join(patch_dir, sparse_matrix_name)
    if os.path.exists(mat_path):
        mat = sp.load_npz(mat_path)
    else:
        # Try default name
        mat = sp.load_npz(os.path.join(patch_dir, 'patch_lri_matrix.npz'))
    
    # Load column names and metadata
    columns_file = os.path.join(patch_dir, 'patch_lri_columns.csv')
    columns_df = pd.read_csv(columns_file)
    column_names = columns_df['column_name'].tolist()
    
    patch_metadata_df = pd.read_csv(os.path.join(patch_dir, 'patch_tma_correspondence.csv'))
    
    print(f"Loaded matrix shape: {mat.shape}")
    print(f"Matrix sparsity: {100 * (1 - mat.nnz / (mat.shape[0] * mat.shape[1])):.2f}%")
    
    return {
        'patch_lri_matrix': mat,
        'column_names': column_names,
        'patch_tma_df': patch_metadata_df
    }


def setup_adata_with_motifs(adata, patch_motifs, patch_dir):
    """Setup AnnData with patch and motif information"""
    print("Setting up AnnData with motif information...")
    
    # Load cell-patch correspondence
    cell_patch_path = os.path.join(patch_dir, 'cell_patch_correspondence.csv')
    cell_patch = pd.read_csv(cell_patch_path)
    
    # Verify alignment using (cell_id, tma_id) pairs for robust mapping
    if len(adata) != len(cell_patch):
        raise ValueError(f"Cell count mismatch: adata {len(adata)}, cell_patch {len(cell_patch)}")
    
    # Verify that the files correspond to the same cells using (cell_id, tma_id) mapping
    adata_keys = [f"{cid}_{tid}" for cid, tid in zip(adata.obs['cell_id'].astype(str), adata.obs['tma_id'].astype(str))]
    patch_keys = [f"{cid}_{tid}" for cid, tid in zip(cell_patch['cell_id'].astype(str), cell_patch['tma_id'].astype(str))]
    
    if adata_keys != patch_keys:
        print("Warning: Row order doesn't match exactly, but checking for complete mapping...")
        if set(adata_keys) != set(patch_keys):
            raise ValueError("Cell-patch correspondence doesn't match adata cells!")
        print("All cells present but order differs - this is acceptable for visualization")
    
    # Add patch_id to adata (using .values preserves row order)
    adata.obs['patch_id'] = cell_patch['patch_id'].values
    
    # Map patches to dominant motifs
    best_idx = patch_motifs.groupby('patch_id')['loading'].idxmax()
    patch2motif = dict(zip(
        patch_motifs.loc[best_idx, 'patch_id'],
        patch_motifs.loc[best_idx, 'motif']
    ))
    
    # Add motif information
    adata.obs['motif'] = adata.obs['patch_id'].map(patch2motif)
    
    # Ensure tma_id is numeric
    adata.obs['tma_id'] = pd.to_numeric(adata.obs['tma_id'], errors='coerce')
    
    print(f"Added motif information for {len(patch2motif)} patches")
    return adata


def get_unique_cell_types(adata=None):
    """Get unique cell types from data"""
    if adata is not None:
        return sorted(set(adata.obs['cell_type'].values))
    else:
        # Default cell types from the notebook
        return ['B cell', 'T cell', 'Tumor', 'dendritic cell', 'endothelial cell', 
                'fat cell', 'fibroblast', 'granulocyte', 'macrophage', 'mast cell',
                'monocyte', 'muscle cell', 'natural killer cell', 'neutrophil', 'pericyte']

def filter_same_celltype_by_lri_name(df, lri_col="lri_name"):
    """
    Remove rows where sender cell type == receiver cell type,
    inferred from 'cell1|cell2|ligand|receptor|signalling' in lri_name.
    """
    if lri_col not in df.columns:
        print(f"Warning: '{lri_col}' not in columns; keeping all rows.")
        return df

    # Split into at most 5 parts to be safe
    parts = df[lri_col].astype(str).str.split("|", n=4, expand=True)
    if parts.shape[1] < 2:
        print("Warning: lri_name does not split into >=2 parts; keeping all rows.")
        return df

    sender = parts[0].str.strip().str.lower()
    receiver = parts[1].str.strip().str.lower()

    keep_mask = sender != receiver
    removed = (~keep_mask).sum()
    print(f"Filtered out {removed} same-cell-type interactions for top LRI plot.")
    return df.loc[keep_mask].copy()



def main():
    """Main visualization function"""
    parser = argparse.ArgumentParser(description='BPTF Visualization Analysis')
    parser.add_argument('--bptf-dir', required=True,
                       help='Directory containing BPTF results')
    parser.add_argument('--patch-dir', required=True,
                       help='Directory containing patch-LRI results')
    parser.add_argument('--data-file', required=True,
                       help='Path to AnnData file')
    parser.add_argument('--output-dir', default='results/bptf_plots',
                       help='Output directory for plots')
    parser.add_argument('--suffix', default='',
                       help='Suffix for output filenames')
    parser.add_argument('--random-state', type=int, default=42,
                       help='Random seed for reproducibility')
    parser.add_argument('--sparse-matrix-name', default='patch_lri_matrix.npz',
                       help='Name of sparse matrix file')
    
    args = parser.parse_args()
    
    # Set random seed
    np.random.seed(args.random_state)
    
    print("=" * 60)
    print("BPTF VISUALIZATION ANALYSIS")
    print("=" * 60)
    print(f"BPTF directory: {args.bptf_dir}")
    print(f"Patch-LRI directory: {args.patch_dir}")
    print(f"Data file: {args.data_file}")
    print(f"Output directory: {args.output_dir}")
    print(f"Random seed: {args.random_state}")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    print("\n1. Loading BPTF results...")
    bptf_results = load_bptf_results(args.bptf_dir)
    
    print("\n2. Loading patch-LRI results...")
    patch_results = load_patch_lri_results(args.patch_dir, args.sparse_matrix_name)
    
    print("\n3. Loading AnnData...")
    adata = ad.read_h5ad(args.data_file)
    
    print("\n4. Setting up motif information...")
    adata = setup_adata_with_motifs(adata, bptf_results['patch_motifs'], args.patch_dir)
    
    # Get unique cell types
    unique_ct = get_unique_cell_types(adata)
    print(f"Found {len(unique_ct)} unique cell types")
    
    # Extract data for plotting
    patch_loadings = bptf_results['patch_loadings']
    lri_factors = bptf_results['lri_factors']
    lri_motifs = bptf_results['lri_motifs']
    column_names = patch_results['column_names']
    patch_metadata_df = patch_results['patch_tma_df']
    
    suffix = f"_{args.suffix}" if args.suffix else ""
    
    print("\n" + "=" * 60)
    print("GENERATING VISUALIZATIONS")
    print("=" * 60)
    
    # 1. Cells per patch distribution
    print("\n1. Plotting cells per patch distribution...")
    save_path = os.path.join(args.output_dir, f"cells_per_patch{suffix}.pdf")
    plot_cells_per_patch(patch_metadata_df, save_path)
    
    # 2. Factor distributions
    print("\n2. Plotting factor distributions...")
    save_path = os.path.join(args.output_dir, f"factor_distributions{suffix}.pdf")
    plot_factor_distributions(patch_loadings, lri_factors, save_path)
    
    # 3. Factor sparsity
    print("\n3. Plotting factor sparsity...")
    save_path = os.path.join(args.output_dir, f"factor_sparsity{suffix}.pdf")
    plot_factor_sparsity(patch_loadings, lri_factors, save_path)
    
    # 4. Motif activities
    print("\n4. Plotting motif activities...")
    save_path = os.path.join(args.output_dir, f"motif_activities{suffix}.pdf")
    plot_motif_activities(patch_loadings, save_path)
    
    # 5. LRI clustermap
    print("\n5. Creating LRI clustermap...")
    save_path = os.path.join(args.output_dir, f"lri_clustermap{suffix}.pdf")
    plot_lri_clustermap(lri_factors, column_names, save_path)
    
    # 6. Cell-cell communication heatmap
    print("\n6. Creating cell-cell communication heatmap...")
    save_path = os.path.join(args.output_dir, f"celltype_communication_heatmap{suffix}.pdf")
    plot_celltype_communication_heatmap(lri_factors, column_names, save_path)
    
    # 7. Cell-cell communication by motif
    print("\n7. Creating cell-cell communication by motif...")
    save_path = os.path.join(args.output_dir, f"celltype_communication_by_motif{suffix}.pdf")
    plot_celltype_communication_by_motif(lri_factors, column_names, args.suffix, save_path)
    
    # # 8. Top LRI interactions
    # print("\n8. Plotting top LRI interactions...")
    # save_path = os.path.join(args.output_dir, f"top_lri_interactions{suffix}.pdf")
    # plot_top_lri_interactions(lri_motifs, unique_ct, args.suffix, save_path)
    # 8. Top LRI interactions (exclude same-cell-type edges for this figure only)
    print("\n8. Plotting top LRI interactions (excluding same-cell-type)...")
    save_path = os.path.join(args.output_dir, f"top_lri_interactions{suffix}.pdf")
    lri_motifs_filtered = filter_same_celltype_by_lri_name(lri_motifs)
    plot_top_lri_interactions(lri_motifs_filtered, unique_ct, args.suffix, save_path)
    
    # 9. LRI networks (without annotations)
    print("\n9. Creating LRI networks...")
    save_path = os.path.join(args.output_dir, f"lri_networks{suffix}.pdf")
    try:
        plot_lri_networks(lri_motifs, unique_ct, args.suffix, 
                         threshold=2000, top_n=200, annotate_edges=False, save_path=save_path)
    except Exception as e:
        print(f"Warning: Could not create LRI networks (requires Graphviz): {e}")
    
    # 10. LRI networks with annotations
    print("\n10. Creating LRI networks with annotations...")
    save_path = os.path.join(args.output_dir, f"lri_networks_annotated{suffix}.pdf")
    try:
        plot_lri_networks(lri_motifs, unique_ct, args.suffix,
                         threshold=500, top_n=200, annotate_edges=True, save_path=save_path)
    except Exception as e:
        print(f"Warning: Could not create annotated LRI networks (requires Graphviz): {e}")
    
    # 11. Spatial motif distribution
    print("\n11. Creating spatial motif distribution...")
    save_path = os.path.join(args.output_dir, f"spatial_motif_all_tmas{suffix}.pdf")
    plot_all_punches_by_cell_type(
        adata=adata, 
        cell_type_column='motif', 
        n_cols=8,
        spot_size=1,
        figsize_per_subplot=(6, 6),
        title='LRI Motifs Distribution Across TMAs',
        save_path=save_path
    )
    
    print("\n" + "=" * 60)
    print("VISUALIZATION COMPLETED SUCCESSFULLY!")
    print("=" * 60)
    print(f"All plots saved to: {args.output_dir}")
    
    # Print summary of generated files
    print("\nGenerated files:")
    for i, filename in enumerate([
        f"cells_per_patch{suffix}.pdf",
        f"factor_distributions{suffix}.pdf", 
        f"factor_sparsity{suffix}.pdf",
        f"motif_activities{suffix}.pdf",
        f"lri_clustermap{suffix}.pdf",
        f"celltype_communication_heatmap{suffix}.pdf",
        f"celltype_communication_by_motif{suffix}.pdf",
        f"top_lri_interactions{suffix}.pdf",
        f"lri_networks{suffix}.pdf",
        f"lri_networks_annotated{suffix}.pdf",
        f"spatial_motif_all_tmas{suffix}.pdf"
    ], 1):
        print(f"  {i:2d}. {filename}")


if __name__ == '__main__':
    main()