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

from patch_lri_analysis import load_patch_lri_results


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

def setup_adata_with_motifs_sc(adata, patch_motifs):
    """
    Add per-patch identifiers to adata.obs by aligning rows in order.

    Parameters
    ----------
    adata : anndata.AnnData
        Has n_obs = number of patches (rows already in the same order as patch_idx).
    patch_motifs : pandas.DataFrame
        Columns: ['patch_idx', 'patch_id', 'motif', 'loading'].
        Likely contains multiple rows per patch (one per motif).

    Returns
    -------
    anndata.AnnData
        adata with adata.obs['patch_idx'] and adata.obs['patch_id'] added.
    """
    import numpy as np
    import pandas as pd

    # Basic checks
    required = {"patch_idx", "patch_id"}
    missing = required - set(patch_motifs.columns)
    if missing:
        raise ValueError(f"patch_motifs missing columns: {missing}")

    # One row per patch_idx (motifs make the table long)
    by_patch = (
        patch_motifs.loc[:, ["patch_idx", "patch_id"]]
        .drop_duplicates(subset="patch_idx", keep="first")
        .sort_values("patch_idx")
        .reset_index(drop=True)
    )

    if by_patch.shape[0] != adata.n_obs:
        raise ValueError(
            f"Length mismatch: {by_patch.shape[0]} unique patch_idx vs {adata.n_obs} rows in adata.obs."
        )

    # Align by row order (you said both are in the same order)
    adata = adata.copy()
    adata.obs = adata.obs.copy()  # avoid view-setting warnings
    adata.obs["patch_idx"] = np.arange(adata.n_obs, dtype=int)
    adata.obs["patch_id"] = pd.Categorical(by_patch["patch_id"].to_numpy())

    return adata


def setup_adata_with_motifs(adata, patch_motifs, patch_dir):
    """Setup AnnData with patch and motif information"""
    print("Setting up AnnData with motif information...")
    
    # Load cell-patch correspondence
    cell_patch_path = os.path.join(patch_dir, 'cell_patch_correspondence.csv')
    cell_patch = pd.read_csv(cell_patch_path)

    adata.obs.index = adata.obs.index.astype(str)
    cell_patch['cell_id'] = cell_patch['cell_id'].astype(str)

    mask = adata.obs.index.isin(cell_patch['cell_id'])
    adata = adata[mask].copy()
    
    # Verify alignment using (cell_id, tma_id) pairs for robust mapping
    if len(adata) != len(cell_patch):
        raise ValueError(f"Cell count mismatch: adata {len(adata)}, cell_patch {len(cell_patch)}")
    
    # # Verify that the files correspond to the same cells using (cell_id, tma_id) mapping
    # adata_keys = [f"{cid}_{tid}" for cid, tid in zip(adata.obs.index.astype(str), adata.obs['tma_id'].astype(str))]
    # patch_keys = [f"{cid}_{tid}" for cid, tid in zip(cell_patch['cell_id'].astype(str), cell_patch['tma_id'].astype(str))]
    
    # if adata_keys != patch_keys:
    #     print("Warning: Row order doesn't match exactly, but checking for complete mapping...")
    #     if set(adata_keys) != set(patch_keys):
    #         raise ValueError("Cell-patch correspondence doesn't match adata cells!")
    #     print("All cells present but order differs - this is acceptable for visualization")
    
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
    # adata.obs['tma_id'] = pd.to_numeric(adata.obs['tma_id'], errors='coerce')
    
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


def plot_top_lri_interactions_by_pathway(lri_motifs, unique_ct, suffix="", save_path=None):
    """Plot top LRI interactions per motif aggregated by pathway.
    Same sender, receiver, and pathway combinations are merged (factors summed).
    """
    # Setup colors
    ct_color_map = get_cell_type_colors(unique_ct)
    
    motifs = sorted(lri_motifs['motif_idx'].unique())
    n_prog = len(motifs)
    top_n = 35

    # Layout
    cols = 3
    rows = int(np.ceil(n_prog / cols))
    fig, axes = plt.subplots(rows, cols,
                             figsize=(cols * 12, rows * 9),
                             constrained_layout=False)
    axes = axes.flatten()

    fig.subplots_adjust(right=0.85, left=0.15, top=0.94, bottom=0.06, wspace=1.35, hspace=0.4)

    # Create a mapping from ligand-receptor pairs to pathways

    for ax, prog in zip(axes, motifs):
        dfp = lri_motifs[lri_motifs['motif_idx'] == prog].copy()
        if dfp.empty:
            ax.set_visible(False)
            continue
        
        # KEY CHANGE: Aggregate by sender, receiver, and pathway (sum factors)
        pathway_agg = dfp.groupby(['celltype1', 'celltype2', 'pathway']).agg({
            'factor_norm': 'sum',  # Sum all factors for same sender-receiver-pathway
            'ligand': lambda x: ', '.join(x.unique()[:2]),  # Keep first 2 unique ligands as examples
            'receptor': lambda x: ', '.join(x.unique()[:2])  # Keep first 2 unique receptors as examples
        }).reset_index()
        
        # Sort by aggregated factor and get top N
        pathway_agg = pathway_agg.sort_values('factor_norm', ascending=False)
        top_df = pathway_agg.head(top_n).reset_index(drop=True)
        y = np.arange(len(top_df))

        # Draw line bars with dots at endpoints
        for yi, row in top_df.iterrows():
            total = row['factor_norm']
            
            # Colors for sender and receiver
            c1 = row['celltype1']
            col1 = ct_color_map.get(c1, 'gray')
            c2 = row['celltype2']
            col2 = ct_color_map.get(c2, 'gray')

            # Draw thin horizontal line (the bar)
            ax.plot([0, total], [yi, yi], color='gray', linewidth=1.5, zorder=1)
            
            # Draw circular markers at start and end (no more signaling type distinction)
            marker_size = 80
            ax.scatter(0, yi, color=col1, s=marker_size, marker='o', 
                      zorder=2, edgecolors='black', linewidth=0.5)
            ax.scatter(total, yi, color=col2, s=marker_size, marker='o',
                      zorder=2, edgecolors='black', linewidth=0.5)

        # Create labels
        labels = []
        for _, row in top_df.iterrows():
            sender = row['celltype1'].ljust(12)[:12]
            receiver = row['celltype2'].ljust(12)[:12]
            pathway = row['pathway'].ljust(20)[:20]
            label = f"{sender} → {receiver} | {pathway}"
            labels.append(label)
        
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9, fontfamily='monospace')
        
        ax.invert_yaxis()
        ax.set_xlabel('Aggregated Normalized Factor', fontsize=10)
        ax.set_title(f'Motif {prog}', fontsize=12, fontweight='bold')
        
        # Set x-axis limits with some padding
        if len(top_df) > 0:
            max_val = top_df['factor_norm'].max()
            ax.set_xlim(-max_val*0.02, max_val*1.05)

    # Hide unused subplots
    for ax in axes[n_prog:]:
        ax.set_visible(False)

    # Create legend (simpler without signaling types)
    legend_handles = []
    
    # Cell type legend (colored dots)
    for ct, col in ct_color_map.items():
        legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor=col, markersize=8, 
                                        markeredgecolor='black', markeredgewidth=0.5,
                                        label=ct))

    fig.legend(handles=legend_handles,
               title='Cell Types',
               loc='center right',
               bbox_to_anchor=(0.92, 0.5),
               fontsize=8,
               title_fontsize=10,
               frameon=True,
               fancybox=True,
               shadow=False)

    fig.suptitle(
        'Top Pathways per Motif (Aggregated)\n' + 
        'Format: Sender → Receiver | Pathway\n' +
        '(Start dot = Sender, End dot = Receiver)',
        fontsize=14, y=1, fontweight='bold'
    )

    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return fig


def plot_top_pathways_only(lri_motifs, unique_ct, cellchatdb, suffix="", save_path=None):
    """Plot top pathways per motif - properly aggregated to avoid duplicates"""
    
    motifs = sorted(lri_motifs['motif_idx'].unique())
    n_prog = len(motifs)
    top_n = 30

    # Layout
    cols = 3
    rows = int(np.ceil(n_prog / cols))
    fig, axes = plt.subplots(rows, cols,
                             figsize=(cols * 8, rows * 7),
                             constrained_layout=False)
    axes = axes.flatten()

    fig.subplots_adjust(right=0.95, left=0.15, top=0.94, bottom=0.06, wspace=0.35, hspace=0.4)

    # Create a mapping from ligand-receptor pairs to pathways
    lr_to_pathway = {}
    for _, row in cellchatdb.iterrows():
        ligands = row['ligand'].split('_') if pd.notna(row['ligand']) else []
        receptors = row['receptor'].split('_') if pd.notna(row['receptor']) else []
        pathway = row['pathway']
        
        for lig in ligands:
            for rec in receptors:
                lr_to_pathway[(lig, rec)] = pathway
        lr_to_pathway[(row['ligand'], row['receptor'])] = pathway

    for ax, prog in zip(axes, motifs):
        dfp = lri_motifs[lri_motifs['motif_idx'] == prog].copy()
        if dfp.empty:
            ax.set_visible(False)
            continue
        
        dfp = dfp[dfp['mean'] > 0.05]
        # Parse interaction names
        dfp = dfp[~dfp['lri_name'].str.startswith('GENE')]
        dfp[['celltype1','celltype2','ligand','receptor','signaling_type']] = (
            dfp['lri_name'].apply(lambda x: pd.Series(parse_lri_full(x)))
        )
        
        # Map to pathways
        def get_pathway(row):
            pathway = lr_to_pathway.get((row['ligand'], row['receptor']))
            if pathway:
                return pathway
            
            for (lig, rec), path in lr_to_pathway.items():
                if (row['ligand'] in lig or lig in row['ligand']) and \
                   (row['receptor'] in rec or rec in row['receptor']):
                    return path
            return 'Unknown'
        
        dfp['pathway'] = dfp.apply(get_pathway, axis=1)
        
        # KEY CHANGE: Aggregate by pathway only (sum all factors for same pathway)
        pathway_agg = dfp.groupby('pathway').agg({
            'factor_norm': 'sum',  # Sum all factors for same pathway
            'celltype1': lambda x: len(x.unique()),  # Count unique sender types
            'celltype2': lambda x: len(x.unique()),  # Count unique receiver types
            'ligand': lambda x: len(x.unique()),     # Count unique ligands
            'receptor': lambda x: len(x.unique())    # Count unique receptors
        }).reset_index()
        pathway_agg.columns = ['pathway', 'factor_norm', 'n_senders', 'n_receivers', 'n_ligands', 'n_receptors']
        
        # Sort and get top pathways
        pathway_agg = pathway_agg.sort_values('factor_norm', ascending=False)
        top_df = pathway_agg.head(top_n).reset_index(drop=True)
        
        # Create horizontal bar plot
        y = np.arange(len(top_df))
        ax.barh(y, top_df['factor_norm'], color='steelblue', alpha=0.7)
        
        # Create labels with counts
        labels = []
        for _, row in top_df.iterrows():
            pathway = row['pathway']
            n_send = row['n_senders']
            n_recv = row['n_receivers']
            n_lig = row['n_ligands']
            n_rec = row['n_receptors']
            label = f"{pathway} ({n_send}S, {n_recv}R, {n_lig}L, {n_rec}R)"
            labels.append(label)
        
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        
        ax.invert_yaxis()
        ax.set_xlabel('Aggregated Normalized Factor', fontsize=10)
        ax.set_title(f'Motif {prog}', fontsize=12, fontweight='bold')

    # Hide unused subplots
    for ax in axes[n_prog:]:
        ax.set_visible(False)

    fig.suptitle(
        'Top Pathways per Motif (Fully Aggregated)\n' + 
        '(S = sender cell types, R = receiver cell types, L = ligands, R = receptors)',
        fontsize=14, y=0.98, fontweight='bold'
    )

    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return fig


def main():
    """Main visualization function"""
    parser = argparse.ArgumentParser(description='BPTF Visualization Analysis')
    parser.add_argument('--bptf-dir', default='results/GBM_cellphone/bptf',
                       help='Directory containing BPTF results')
    parser.add_argument('--patch-dir', default='results/GBM_cellphone',
                       help='Directory containing patch-LRI results')
    parser.add_argument('--data-file', default='/Users/jiayifan/Desktop/Lab/TMA_punch_subfiles/xenium_mm_final_cell_id.h5ad',
                       help='Path to AnnData file')
    parser.add_argument('--output-dir', default='results/GBM_cellphone/bptf_plots',
                       help='Output directory for plots')
    parser.add_argument('--suffix', default='',
                       help='Suffix for output filenames')
    parser.add_argument('--random-state', type=int, default=42,
                       help='Random seed for reproducibility')
    parser.add_argument('--sparse-matrix-name', default='patch_lri_matrix.npz',
                       help='Name of sparse matrix file')
    parser.add_argument('--column-df-name', default='patch_lri_columns.csv',
                       help='Name of the column names file in input directory')
    parser.add_argument('--meta-df-name', default='patch_tma_correspondence.csv',
                       help='Name of the metadata file in input directory')
    parser.add_argument('--neighborhood', type=bool,
                       default=False, help='if neighbood-based or patch-based')
    parser.add_argument('--single-cell', type=bool,
                       default=False, help='if cell-based or patch-based')
        
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
    patch_results = load_patch_lri_results(args.patch_dir, args.sparse_matrix_name, 
                                           column_df_name = args.column_df_name, meta_df_name = args.meta_df_name, 
                                           neighborhood=args.neighborhood, single_cell = args.single_cell)
    
    print("\n3. Loading AnnData...")
    adata = ad.read_h5ad(args.data_file)
    
    print("\n4. Setting up motif information...")
    adata = setup_adata_with_motifs(adata, bptf_results['patch_motifs'], args.patch_dir)
    # adata = setup_adata_with_motifs_sc(adata, bptf_results['patch_motifs'])
    
    # Get unique cell types
    unique_ct = get_unique_cell_types(adata)
    print(f"Found {len(unique_ct)} unique cell types")
    
    # Extract data for plotting
    patch_loadings = bptf_results['patch_loadings']
    lri_factors = bptf_results['lri_factors']
    lri_motifs = bptf_results['lri_motifs']
    column_names = patch_results['column_names']
    if not args.neighborhood:
        patch_metadata_df = patch_results['patch_tma_df']
    
    suffix = f"_{args.suffix}" if args.suffix else ""

    # modify lri_motifs

    lr_to_pathway = {}
    cellchatdb = pd.read_csv('/Users/jiayifan/tansey_lab/alarmist/data/LRdatabase/CellChatDBv2.0.human.csv')
    for _, row in cellchatdb.iterrows():
        ligands = row['ligand'].split('_') if pd.notna(row['ligand']) else []
        receptors = row['receptor'].split('_') if pd.notna(row['receptor']) else []
        pathway = row['pathway']
        
        for lig in ligands:
            for rec in receptors:
                lr_to_pathway[(lig, rec)] = pathway
        lr_to_pathway[(row['ligand'], row['receptor'])] = pathway


    # lri_motifs = lri_motifs[lri_motifs['mean'] > 0.05]
    lri_motifs[['celltype1','celltype2','ligand','receptor','signaling_type']] = (
            lri_motifs['lri_name'].apply(lambda x: pd.Series(parse_lri_full(x)))
        )

    def get_pathway(row):
        pathway = lr_to_pathway.get((row['ligand'], row['receptor']))
        if pathway:
            return pathway
        
        for (lig, rec), path in lr_to_pathway.items():
            if (row['ligand'] in lig or lig in row['ligand']) and \
                (row['receptor'] in rec or rec in row['receptor']):
                return path
        return 'Unknown'

    lri_motifs['pathway'] = lri_motifs.apply(get_pathway, axis=1)

    mat = patch_results['patch_lri_matrix']
    cols = patch_results['column_names']

    column_means = np.array(mat.mean(axis=0)).flatten()
    print(f"Shape of column means: {column_means.shape}")
    print(f"First 10 means: {column_means[:10]}")

    lri_to_mean = dict(zip(cols, column_means))

    lri_motifs['factor_norm'] = lri_motifs.apply(
        lambda row: row['factor'] / lri_to_mean.get(row['lri_name'], 1.0) 
        if lri_to_mean.get(row['lri_name'], 0) > 0 else row['factor'], 
        axis=1
    )

    lri_motifs['mean'] = lri_motifs['lri_name'].map(lri_to_mean)

    print(f"Added factor_norm column. Shape: {lri_motifs.shape}")
    lri_motifs = lri_motifs.drop(columns=['factor_normalized'])

    lri_motifs = lri_motifs[lri_motifs['mean']>0.05]
    
    print("\n" + "=" * 60)
    print("GENERATING VISUALIZATIONS")
    print("=" * 60)
    
    # 1. Cells per patch distribution
    print("\n1. Plotting cells per patch distribution...")
    if not args.neighborhood:
        save_path = os.path.join(args.output_dir, f"cells_per_patch{suffix}.pdf")
        plot_cells_per_patch(patch_metadata_df, save_path)
    else:
        print("  Skipping cells per patch plot for neighborhood-based analysis.")
    
    # 5. LRI clustermap
    print("\n5. Creating LRI clustermap...")
    save_path = os.path.join(args.output_dir, f"lri_clustermap{suffix}.pdf")
    plot_lri_clustermap(lri_factors, column_names, save_path)
    
    # 7. Cell-cell communication by motif
    print("\n7. Creating cell-cell communication by motif...")
    save_path = os.path.join(args.output_dir, f"celltype_communication_by_motif{suffix}.pdf")
    plot_celltype_communication_by_motif(lri_factors, column_names, args.suffix, save_path)
    
    # # 8. Top LRI interactions
    print("\n8. Plotting top LRI interactions...")
    # save_path = os.path.join(args.output_dir, f"top_lri_interactions{suffix}.pdf")
    # plot_top_lri_interactions(lri_motifs, unique_ct, args.suffix, save_path)
    save_path = os.path.join(args.output_dir, f"top_lri_interactions_dot{suffix}.pdf")
    plot_top_lri_interactions_dot(lri_motifs, unique_ct, use_normalized=True, suffix="", save_path = save_path)
    # 8. Top LRI interactions (exclude same-cell-type edges for this figure only)
    # print("\n8. Plotting top LRI interactions (excluding same-cell-type)...")
    # save_path = os.path.join(args.output_dir, f"top_lri_interactions{suffix}.pdf")
    # # lri_motifs_filtered = filter_same_celltype_by_lri_name(lri_motifs)
    # plot_top_lri_interactions(lri_motifs, unique_ct, args.suffix, save_path)
    save_path = os.path.join(args.output_dir, f"pathways_by_motif_normalized{suffix}.pdf")
    plot_top_lri_interactions_by_pathway(lri_motifs, unique_ct, suffix="", save_path=save_path)
    
    # 9. LRI networks (without annotations)
    print("\n9. Creating LRI networks...")
    save_path = os.path.join(args.output_dir, f"lri_networks{suffix}.png")
    # try:
        # plot_lri_networks(lri_motifs, unique_ct, suffix="_all",
        #           threshold=1500, top_n=200,
        #           annotate_edges=False,
        #           save_path="networks_all.png",
        #           mode_filter=None)
        # # 2) 仅 paracrine
        # plot_lri_networks(lri_motifs, unique_ct, suffix="_paracrine",
        #                 threshold=1000, top_n=200,
        #                 annotate_edges=False,
        #                 save_path="networks_paracrine.png",
        #                 mode_filter="paracrine")

        # # 3) 仅 juxtacrine
        # plot_lri_networks(lri_motifs, unique_ct, suffix="_juxtacrine",
        #                 threshold=1000, top_n=200,
        #                 annotate_edges=False,
        #                 save_path="networks_juxtacrine.png",
        #                 mode_filter="juxtacrine")
        
    # 先用 ALL 的聚合构建主门控（>1500）
    master_gate = build_master_edge_gate(lri_motifs, top_n=200, threshold=1500)

    # 1) ALL（自己仍然可设阈值=1500）
    save_path = os.path.join(args.output_dir, f"lri_networks_all{suffix}.png")
    plot_lri_networks(lri_motifs, unique_ct,
                    threshold=10000, top_n=200,
                    annotate_edges=False,
                    save_path=save_path,
                    mode_filter=None,
                    edge_gate=None)  # ALL 图不需要 gate

    # 2) JUxtacrine：不再阈值，使用 ALL>1500 的门控集合
    save_path = os.path.join(args.output_dir, f"lri_networks_jux{suffix}.png")
    plot_lri_networks(lri_motifs, unique_ct,
                    threshold=0, top_n=200,
                    annotate_edges=False,
                    save_path=save_path,
                    mode_filter="juxtacrine",
                    edge_gate=master_gate)

    # 3) PARacrine：同理
    save_path = os.path.join(args.output_dir, f"lri_networks_par{suffix}.png")
    plot_lri_networks(lri_motifs, unique_ct,
                    threshold=0, top_n=200,
                    annotate_edges=False,
                    save_path=save_path,
                    mode_filter="paracrine",
                    edge_gate=master_gate)

    # except Exception as e:
        # print(f"Warning: Could not create LRI networks (requires Graphviz): {e}")
    
    # 10. LRI networks with annotations
    print("\n10. Creating LRI networks with annotations...")
    save_path = os.path.join(args.output_dir, f"lri_networks_annotated{suffix}.png")
    try:
        plot_lri_networks(lri_motifs, unique_ct, args.suffix,
                         threshold=1500, top_n=200, annotate_edges=True, save_path=save_path)
    except Exception as e:
        print(f"Warning: Could not create annotated LRI networks (requires Graphviz): {e}")
    
    # 11. Spatial motif distribution
    print("\n11. Creating spatial motif distribution...")
    save_path = os.path.join(args.output_dir, f"spatial_motif_all_tmas{suffix}.pdf")
    # plot_all_punches_by_cell_type(
    #     adata=adata, 
    #     cell_type_column='motif', 
    #     n_cols=8,
    #     spot_size=1,
    #     figsize_per_subplot=(6, 6),
    #     title='LRI Motifs Distribution Across TMAs',
    #     save_path=save_path
    # )
    
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