import alarmist as al
import scanpy as sc
from bptf import load_bptf
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import os


results_dir = "/Users/jiayifan/tansey_lab/alarmist/results/SAHA_IBD_more"
lri_results = al.load_patch_lri_results(results_dir)
n_components = 20

adata = sc.read_h5ad('/Users/jiayifan/tansey_lab/saha_ibd/SAHA_IBD_RNA.h5ad')
adata.obs['cell_type'] = adata.obs['cell_type_general']
adata = adata[adata.obs['cell_type'].notna()].copy()
unique_ct = sorted(adata.obs['cell_type'].unique())
adata_dict = {batch: adata[adata.obs['section_ID'] == batch].copy() for batch in adata.obs['section_ID'].unique()}
del adata

print("Running neighborhood-based LRI analysis...")
patch_size_px = 50/0.12028

cell_analyzer = al.NeighborhoodLRIAnalyzer(
    neighborhood_size=patch_size_px,
    resource_name='cellchatdb'
)
cell_analyzer.cellchatdb_path = 'data/LRdatabase/CellChatDBv2.0.human.csv'

cell_results = cell_analyzer.run_neighborhood(
    adata_dict,  # {'P17_AIS': adata1, 'P17_LUAD': adata2, ...}
    output_dir=f'{results_dir}/single_cell',
    required_columns=lri_results['column_names']
)

model = load_bptf(Path(f"{results_dir}/bptf_{n_components}/bptf_model.npz/bptf_1.npz"))
np.random.seed(42)

# cell_loadings = al.project_cell_loadings(
#         model=model,
#         cell_lri_matrix=cell_results['cell_lri_matrix'],
#         max_iter=200,
#         chunk_size=5000,
#         verbose=True
#     )

import gc
import numpy as np

# 先读一次拿到 sample_info（然后立刻删掉大对象）
tmp = al.load_cell_lri_results(f'{results_dir}/single_cell')
sample_info = tmp['sample_info']
del tmp
gc.collect()

for i, row in sample_info.iterrows():
    print(i, row)
    cell_results = al.load_cell_lri_results(f'{results_dir}/single_cell')
    mat = cell_results['cell_lri_matrix']

    start = row['global_cell_idx_start']
    end   = row['global_cell_idx_end'] + 1

    subset = mat[start:end, :].copy()   # <- 关键：copy 断开对大矩阵的引用
    del mat
    del cell_results
    gc.collect()

    print("Projecting cell loadings from patch-level BPTF model...")

    cell_loadings = al.project_cell_loadings(
        model=model,
        cell_lri_matrix=subset,
        max_iter=200,
        chunk_size=5000,
        verbose=True
    )

    np.save(f'{results_dir}/cell_loadings_{i}.npy', cell_loadings)


# np.save(f'{results_dir}/single_cell/cell_loadings.npy', cell_loadings)
# print(f"Cell loadings saved to: {results_dir}/single_cell/cell_loadings.npy")

cell_results = al.load_cell_lri_results(f'{results_dir}/single_cell')
# cell_loadings = np.load(f'{results_dir}/single_cell/cell_loadings.npy')

adata = sc.read_h5ad('/Users/jiayifan/tansey_lab/saha_ibd/SAHA_IBD_RNA.h5ad')
adata.obs['cell_type'] = adata.obs['cell_type_general']
adata = adata[adata.obs['cell_type'].notna()].copy()
unique_ct = sorted(adata.obs['cell_type'].unique())
adata_dict = {batch: adata[adata.obs['section_ID'] == batch].copy() for batch in adata.obs['section_ID'].unique()}

# split up cell_loadings
cell_loadings_by_sample = {}
sample_info = cell_results['sample_info']
for i, info in sample_info.iterrows():
    sample_id = info['sample_id']

    cell_loadings = np.load(f'{results_dir}/cell_loadings_{i}.npy')
    
    cell_loadings_by_sample[sample_id] = cell_loadings
    print(f"{sample_id}: {cell_loadings.shape}")

ct_color_map = {
    ct: mcolors.to_hex(plt.get_cmap('tab20', len(unique_ct))(i))
    for i, ct in enumerate(unique_ct)
}


import alarmist as al
from alarmist.plotting import (
    plot_motif_celltype_composition,
    plot_motif_state_counts,
    plot_positive_motifs_distribution,
    plot_motif_spatial
)


# # ========== Step 1: Cell Type Composition Analysis ==========
# print("1. Analyzing cell type composition per motif...")
# # Collect cell metadata from all samples
# cell_meta_dfs = []
# for sample_id, adata in adata_dict.items():
#     df = adata.obs[['cell_type']].copy()
#     df['sample_id'] = sample_id
#     cell_meta_dfs.append(df)

# cell_meta_df = pd.concat(cell_meta_dfs, axis=0)

# print(len(cell_meta_df))
# print(cell_loadings.shape)
# # Compute weighted cell types for each motif
# tidy_df = al.weighted_celltypes_by_motif(
#     cell_loadings=cell_loadings,
#     metadata_df=cell_meta_df,
#     normalize=True,
#     top_n_per_motif=20,
#     other_label="Other"
# )

# # Visualize
# fig, ax = plot_motif_celltype_composition(
#     tidy_df,
#     title="Cell Type Composition per Motif (All Samples)",
#     color_map=ct_color_map,
#     save_path=f'{results_dir}/single_cell/motif_celltype_composition.png'
# )
# plt.show()


for sample_id, adata in adata_dict.items():
    print(f"\n--- Processing {sample_id} ---")
    gmm_summary = al.gmm_binarize_all_motifs(
        cell_loadings=cell_loadings_by_sample[sample_id],
        adata=adata,
        eps=1e-10,
        random_state=42
    )



# Results are now in adata.obs:
# - motif_0_state: 'negative' or 'positive'
# - motif_0_posprob: probability of being positive
# - motif_0_score_log: log-transformed loading

# ========== Step 3: ON/OFF Statistics ==========
print("\n3. Computing ON/OFF statistics...")

# 3a. Counts per motif
counts_df = al.compute_motif_state_counts(adata)
print("\nON/OFF counts per motif:")
print(counts_df)

# Visualize
fig, ax = plot_motif_state_counts(
    counts_df,
    title="Positive vs Negative Cells per Motif"
)
plt.savefig(f'{results_dir}/single_cell/motif_state_counts.png', dpi=300, bbox_inches='tight')
# plt.show()



import numpy as np
import matplotlib.lines as mlines

def plot_motif_spatial_grid(
    adata_dict: dict,
    motif_k: int,
    n_cols: int = 4,
    figsize_per_plot: tuple = (8, 8),
    point_size: float = 0.6,
    colors: dict = None,
    ct_color_map: dict = None,     # <- 新增：你的 cell type colormap
    celltype_col: str = "cell_type",
    save_path: str = None,
    show_celltype_legend: bool = False,  # 默认不显示（celltype太多会炸）
    legend_top_n: int = 12,              # 若要显示，只显示 top N 的 cell types
):
    if colors is None:
        colors = {'negative': '#d3d3d3'}
    if ct_color_map is None:
        raise ValueError("Please pass ct_color_map (dict: cell_type -> color hex).")

    n_samples = len(adata_dict)
    n_rows = (n_samples + n_cols - 1) // n_cols

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize_per_plot[0]*n_cols, figsize_per_plot[1]*n_rows)
    )
    axes = axes.flatten() if n_samples > 1 else [axes]

    state_col = f'motif_{motif_k}_state'

    for i, (sample_id, adata) in enumerate(adata_dict.items()):
        ax = axes[i]
        coords = adata.obsm['spatial'][:, :2]

        if state_col not in adata.obs.columns:
            ax.set_title(f'{sample_id}\n(no motif {motif_k} data)')
            ax.axis('off')
            continue

        states = adata.obs[state_col].astype(str)
        n_pos = (states == 'positive').sum()
        n_neg = (states == 'negative').sum()
        frac_pos = n_pos / len(states) * 100

        # --- background: negatives in grey ---
        mask_neg = (states == 'negative').values
        ax.scatter(
            coords[mask_neg, 0],
            coords[mask_neg, 1],
            c=colors['negative'],
            s=point_size,
            alpha=0.95,
            marker='o',
            edgecolors='none',
            linewidths=0,
            rasterized=True,
        )

        # --- positives: color by cell type ---
        mask_pos = (states == 'positive').values
        if mask_pos.any():
            ct = adata.obs[celltype_col].astype(str).values
            # 未知 celltype 给一个默认色（可选）
            pos_colors = np.array([ct_color_map.get(x, "#000000") for x in ct])[mask_pos]

            ax.scatter(
                coords[mask_pos, 0],
                coords[mask_pos, 1],
                c=pos_colors,
                s=point_size * 6,     # 适当放大，不然被灰底吞掉
                alpha=1,
                marker='o',
                edgecolors='none',
                linewidths=0,
                rasterized=True,
                zorder=3
            )

        ax.set_title(f'{sample_id}\nMotif {motif_k}: {frac_pos:.1f}% ON')
        ax.set_aspect('equal')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

        # --- legend：默认只给 ON/OFF 概览（不列 cell types）---
        off_handle = mlines.Line2D([], [], color=colors['negative'], marker='o',
                                  linestyle='None', markersize=6, label=f'OFF (n={n_neg:,})')
        on_handle = mlines.Line2D([], [], color='black', marker='o',
                                 linestyle='None', markersize=6, label=f'ON (n={n_pos:,}, {frac_pos:.1f}%)')

        handles = [off_handle, on_handle]

        # 可选：显示 top N 个 cell types（只对 ON 的 cell types 统计）
        if show_celltype_legend and mask_pos.any():
            ct_on = adata.obs.loc[mask_pos, celltype_col].astype(str)
            top_cts = ct_on.value_counts().head(legend_top_n).index.tolist()
            ct_handles = [
                mlines.Line2D([], [], color=ct_color_map.get(ct, "#000000"), marker='o',
                             linestyle='None', markersize=6, label=ct)
                for ct in top_cts
            ]
            handles += ct_handles

        ax.legend(handles=handles, loc='upper right', fontsize=8, frameon=False)

    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    plt.suptitle(f'Motif {motif_k} Spatial Distribution (ON colored by cell type)', fontsize=14, y=1.02)
    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig



for i in range(20):
    fig = plot_motif_spatial_grid(
        adata_dict,
        motif_k=i,
        n_cols=4,
        figsize_per_plot=(9, 9),
        point_size=0.2,
        ct_color_map=ct_color_map,
        celltype_col="cell_type",
        save_path=f'{results_dir}/motif_spatial_plots/motif_{i}_celltypeON.png',
        show_celltype_legend=False
    )
    print(f"Plotted motif {i}")