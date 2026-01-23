import alarmist as al
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import os

from alarmist.plotting import (
    add_lri_components,
    annotate_pathways,
    plot_lri_clustermap,
    plot_celltype_communication_by_motif,
    plot_top_lri_interactions_dot,
    plot_top_lri_interactions_by_pathway,
    build_master_edge_gate,
    plot_lri_networks
)

results_dir = "/Users/jiayifan/tansey_lab/alarmist/results/SAHA_IBD_more"

n_components = 20

adata = sc.read_h5ad('/Users/jiayifan/tansey_lab/saha_ibd/SAHA_IBD_RNA.h5ad')
adata = adata[adata.obs['cell_type_general'].notna()].copy()

adata_dict = {batch: adata[adata.obs['section_ID'] == batch].copy() for batch in adata.obs['section_ID'].unique()}


patch_size_px = 50/0.12028

analyzer = al.PatchLRIAnalyzer(
    patch_size=patch_size_px,
    resource_name='cellchatdb',
    cell_type_column='cell_type_general'
)
analyzer.cellchatdb_path = 'data/LRdatabase/CellChatDBv2.0.human.csv'

results = analyzer.run_patchify(
    adata_dict,
    output_dir=results_dir
)

# results = al.load_patch_lri_results(results_dir)

# Run BPTF
model = al.run_bptf(
    results['patch_lri_matrix'],
    n_components=n_components,
    max_iter=500,
    verbose=True,
    random_state=42
)

# Extract factors
patch_loadings, lri_factors = al.extract_factors(model)

print(f"Patch loadings shape: {patch_loadings.shape}")
print(f"LRI factors shape: {lri_factors.shape}")

def plot_cells_per_patch(patch_metadata_df, save_path=None):
    """Plot distribution of cells per patch"""
    plt.figure(figsize=(8, 6))
    plt.hist(patch_metadata_df['n_cells'], bins=50, edgecolor='black')
    # plt.hist(patch_metadata_df['neighborhood_size'], bins=50, edgecolor='black')
    plt.xlabel('Number of cells per patch')
    plt.ylabel('Number of patches')
    plt.title('Distribution of cells per patch')
    plt.tight_layout()
    
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return plt.gcf()

fig = plot_cells_per_patch(
    results['patch_metadata_df'],
    save_path=f'{results_dir}/bptf_plots_{n_components}/cells_per_patch.pdf'
)

# Save BPTF results (includes automatic factor normalization)
al.save_bptf_results(
    model=model,
    patch_loadings=patch_loadings,
    lri_factors=lri_factors,
    column_names=results['column_names'],
    patch_lri_matrix=results['patch_lri_matrix'],
    output_dir=f"{results_dir}/bptf_{n_components}"
)

adata = sc.read_h5ad('/Users/jiayifan/tansey_lab/saha_ibd/SAHA_IBD_RNA.h5ad')
adata = adata[adata.obs['cell_type_general'].notna()].copy()
# Get unique cell types from your data
unique_ct = sorted(adata.obs['cell_type_general'].unique())

del adata

results = al.load_bptf_results(f"{results_dir}/bptf_{n_components}")
# patch_loadings = results['patch_loadings']
# lri_factors = results['lri_factors']
lri_motifs = results['lri_motifs']

# 2.3 Top LRI interactions (dot plot - RECOMMENDED)
print("3. Plotting top LRI interactions...")
fig = plot_top_lri_interactions_dot(
    lri_motifs,
    unique_ct,
    factor_col='factor_lrnorm',
    top_n=40,
    figsize_per_motif=(12, 12),
    save_path=f'{results_dir}/bptf_plots_{n_components}/top_lri_interactions_factor_lrnorm.pdf'
)

fig = plot_top_lri_interactions_by_pathway(
    lri_motifs,
    unique_ct,
    factor_col='factor_lrnorm',
    top_n=20,
    save_path=f'{results_dir}/bptf_plots_{n_components}/top_pathways_factor_lrnorm.pdf'
)

# 2.5 Network graphs (requires Graphviz)
print("5. Creating network graphs...")
fig=plot_lri_networks(
    lri_motifs,
    unique_ct,
    top_n=200,
    threshold=50,
    mode_filter=None,
    factor_col='factor_lrnorm',
    save_path=f'{results_dir}/bptf_plots_{n_components}/networks_factor_lrnorm.png'
)

fig = plot_celltype_communication_by_motif(
      lri_motifs,
      factor_col='factor_lrnorm',
      n_cols=5,
      save_path=f'{results_dir}/bptf_plots_{n_components}/celltype_communication_factor_lrnorm.pdf'
)



