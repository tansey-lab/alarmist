"""
Differential Expression Analysis: Motif 24-only vs Double-positive vs Motif 10-only
===================================================================================

This script performs pairwise differential expression analysis between three cell groups
defined by their motif loading status:
- Motif 24-only: cells positive for motif 24 but not motif 10 (healthy vasculature)
- Double-positive: cells positive for both motif 24 and motif 10 (intermediate state)
- Motif 10-only: cells positive for motif 10 but not motif 24 (tumor-driven vasculature)

The analysis is performed for specific cell types of interest:
Tumor_epi, Macro, T, B, Neutro, Mast, pDC, cDC

Comparisons (Strategy A - pairwise):
1. Motif 24-only vs Motif 10-only (two extremes)
2. Motif 24-only vs Double-positive
3. Double-positive vs Motif 10-only
"""

import alarmist as al
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import os
from scipy import stats
import scipy.sparse as sp
from typing import Tuple, Optional, Dict, List
from alarmist.plotting import glm_plots
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Directories
RESULTS_DIR = '../results/AIS_LUAD'
OUTPUT_DIR = f'{RESULTS_DIR}/motif_de_analysis'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Motif indices (0-indexed)
MOTIF_24 = 24  # Healthy vasculature, pro-inflammatory pDC
MOTIF_10 = 10  # Tumor-driven vasculature, immune suppressive pDC

# Cell types of interest for DE analysis
CELL_TYPES_OF_INTEREST = ['Tumor_epi', 'Macro', 'T', 'B', 'Neutro', 'Mast', 'pDC', 'cDC']

# Comparison pairs for pairwise DE analysis
COMPARISONS = [
    ('M24_only', 'M10_only'),    # Two extremes - main comparison
    ('M24_only', 'Double'),      # Transition from healthy to intermediate
    ('Double', 'M10_only'),      # Transition from intermediate to tumor
]

# ============================================================================
# SECTION 1: DATA LOADING
# ============================================================================

print("=" * 70)
print("SECTION 1: Loading Data")
print("=" * 70)

# Load CellChat database
cellchatdb = pd.read_csv('../data/LRdatabase/CellChatDBv2.0.human.csv')

# Define input files
files = {
    'P17_AIS': "../data/linghua/P17_AIS_Xenium.h5ad",
    'P17_LUAD': "../data/linghua/P17_LUAD_Xenium.h5ad",
    'P21_AIS': "../data/linghua/P21_AIS_Xenium.h5ad",
    'P21_LUAD': "../data/linghua/P21_LUAD_Xenium.h5ad"
}

# Load and concatenate all samples
adata_list = []
for name, path in files.items():
    adata = sc.read_h5ad(path)
    adata.obs.rename(columns={"annotation_coarse": "cell_type"}, inplace=True)
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    adata.obs["sample_id"] = name
    mask = adata.obs['cell_type'].notna()
    print(f"{name}: {adata.n_obs} -> {mask.sum()} cells ({(~mask).sum()} nan removed)")
    adata_list.append(adata[mask].copy())

adata = sc.concat(adata_list, join="outer")
adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
print(f"\nCombined: {adata.n_obs} total cells")

# Set cell type colors
al.set_celltype_colors(adata, column='cell_type', palette='tab20')

# Add condition and patient metadata
adata.obs['condition'] = adata.obs['sample_id'].str.split('_').str[1]  # AIS or LUAD
adata.obs['patient'] = adata.obs['sample_id'].str.split('_').str[0]    # P17 or P21

# ============================================================================
# SECTION 2: MOTIF LOADING AND BINARIZATION
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 2: Loading Motif Data and Binarizing")
print("=" * 70)

# Load cell loadings from ALARMIST results
cell_loadings = np.load(f'{RESULTS_DIR}/single_cell/cell_loadings.npy')
print(f"Cell loadings shape: {cell_loadings.shape}")

# Binarize motif loadings using GMM (for multi-sample data)
gmm_summary = al.gmm_binarize_all_motifs(
    cell_loadings,
    adata,
    multi_sample=True,
    sample_column='sample_id'
)

# ============================================================================
# SECTION 3: DEFINE THREE CELL GROUPS BASED ON MOTIF STATUS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 3: Defining Cell Groups by Motif Status")
print("=" * 70)

# Get binary motif status for motif 24 and motif 10
# The gmm_summary should have added columns to adata.obs like 'motif_X_positive'
# If not, we need to extract from the binarization results

# Check if motif columns exist in adata.obs
motif_24_col = f'motif_{MOTIF_24}_positive'
motif_10_col = f'motif_{MOTIF_10}_positive'

if motif_24_col not in adata.obs.columns or motif_10_col not in adata.obs.columns:
    # If columns don't exist, create them from gmm_summary
    print("Creating motif positive columns from GMM binarization...")
    
    # Re-run binarization and store results
    for motif_idx in [MOTIF_24, MOTIF_10]:
        col_name = f'motif_{motif_idx}_positive'
        # Get the binary assignments from gmm_summary
        # This depends on how al.gmm_binarize_all_motifs returns results
        # Assuming it adds columns to adata.obs or returns a dict
        if col_name not in adata.obs.columns:
            # Fallback: use threshold-based binarization
            loadings = cell_loadings[:, motif_idx]
            threshold = np.percentile(loadings, 75)  # Top 25% as positive
            adata.obs[col_name] = loadings > threshold
            print(f"  {col_name}: {adata.obs[col_name].sum()} positive cells (threshold={threshold:.4f})")

# Extract boolean arrays for motif status
motif_24_positive = adata.obs[motif_24_col].values.astype(bool)
motif_10_positive = adata.obs[motif_10_col].values.astype(bool)

# Define three mutually exclusive groups
adata.obs['motif_group'] = 'Other'
adata.obs.loc[motif_24_positive & ~motif_10_positive, 'motif_group'] = 'M24_only'
adata.obs.loc[motif_24_positive & motif_10_positive, 'motif_group'] = 'Double'
adata.obs.loc[~motif_24_positive & motif_10_positive, 'motif_group'] = 'M10_only'

# Print group statistics
print("\nMotif group distribution:")
group_counts = adata.obs['motif_group'].value_counts()
for group, count in group_counts.items():
    print(f"  {group}: {count} cells ({100*count/adata.n_obs:.1f}%)")

# Print group distribution by cell type
print("\nMotif group distribution by cell type:")
ct_group_counts = pd.crosstab(adata.obs['cell_type'], adata.obs['motif_group'])
print(ct_group_counts.to_string())

# Save group assignments
adata.obs[['cell_type', 'motif_group', 'condition', 'patient']].to_csv(
    f'{OUTPUT_DIR}/cell_motif_groups.csv'
)

# ============================================================================
# SECTION 4: PREPROCESSING FOR DE ANALYSIS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 4: Preprocessing Data for DE Analysis")
print("=" * 70)

# Preprocess each sample separately to avoid batch effects in normalization
adatas = []
for sample_id in adata.obs['sample_id'].unique():
    adata_sub = adata[adata.obs['sample_id'] == sample_id].copy()
    
    # Filter genes with very low expression
    sc.pp.filter_genes(adata_sub, min_cells=20)
    
    # Normalize and log transform
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)
    
    adatas.append(adata_sub)
    print(f"Processed {sample_id}: {adata_sub.n_obs} cells, {adata_sub.n_vars} genes")

# Concatenate back (inner join to keep only shared genes)
adata_processed = sc.concat(adatas, join='inner')
print(f"\nAfter preprocessing: {adata_processed.n_obs} cells, {adata_processed.n_vars} genes")

# ============================================================================
# SECTION 5: DIFFERENTIAL EXPRESSION FUNCTIONS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 5: Setting Up DE Analysis Functions")
print("=" * 70)

def differential_expression(X, in_mask, out_mask=None, 
                           min_in_group_fraction=0.0001, 
                           min_out_group_fraction=0.0001):
    """
    Perform differential expression analysis using Mann-Whitney U test.
    
    Parameters
    ----------
    X : array-like
        Gene expression matrix (cells x genes)
    in_mask : array-like
        Boolean mask for the "in" group (positive logFC = higher in this group)
    out_mask : array-like, optional
        Boolean mask for the "out" group. If None, uses ~in_mask
    min_in_group_fraction : float
        Minimum fraction of cells expressing gene in the in-group
    min_out_group_fraction : float
        Minimum fraction of cells expressing gene in the out-group
        
    Returns
    -------
    dict
        Dictionary with p_values, p_adj (FDR), and logfoldchanges
    """
    if out_mask is None:
        out_mask = ~in_mask

    # Filter genes with very low expression in either group
    X_out = X[out_mask]
    X_in = X[in_mask]
    
    genes_mask = np.ones(X.shape[1], dtype=bool)
    
    if min_out_group_fraction > 0:
        control_genes_mask = (X_out > 0.5).mean(axis=0) >= min_out_group_fraction
        genes_mask = genes_mask & np.asarray(control_genes_mask).flatten()

    if min_in_group_fraction > 0:
        target_genes_mask = (X_in > 0.5).mean(axis=0) >= min_in_group_fraction
        genes_mask = genes_mask & np.asarray(target_genes_mask).flatten()

    genes_mask = np.array(genes_mask).flatten()
    X_out_filtered = X_out[:, genes_mask]
    X_in_filtered = X_in[:, genes_mask]

    # Calculate log-fold changes
    control_means = np.asarray(X_out_filtered.mean(axis=0)).ravel()
    target_means = np.asarray(X_in_filtered.mean(axis=0)).ravel()
    
    logfoldchanges = np.zeros(X.shape[1])
    eps = 1e-9
    logfoldchanges[genes_mask] = (np.log2(np.maximum(target_means, eps)) - 
                                   np.log2(np.maximum(control_means, eps)))

    # Mann-Whitney U test for each gene
    p_values = np.ones(X.shape[1])
    for j_idx, j in enumerate(np.where(genes_mask)[0]):
        x = X_in_filtered[:, j_idx]
        y = X_out_filtered[:, j_idx]
        
        # Handle sparse matrices
        if hasattr(x, 'toarray'):
            x = x.toarray().ravel()
        else:
            x = np.asarray(x).ravel()
        if hasattr(y, 'toarray'):
            y = y.toarray().ravel()
        else:
            y = np.asarray(y).ravel()
            
        try:
            p_values[j] = stats.mannwhitneyu(x, y, alternative='two-sided').pvalue
        except ValueError:
            p_values[j] = 1.0

    # FDR correction
    p_adj = np.ones(X.shape[1])
    p_adj[genes_mask] = stats.false_discovery_control(p_values[genes_mask])

    return {'p_values': p_values, 'p_adj': p_adj, 'logfoldchanges': logfoldchanges}


def run_pairwise_de(adata_ct, group1, group2, cell_type):
    """
    Run DE analysis between two motif groups for a specific cell type.
    
    Parameters
    ----------
    adata_ct : AnnData
        AnnData object subset to a specific cell type
    group1 : str
        Name of the first group (positive logFC = higher in group1)
    group2 : str
        Name of the second group
    cell_type : str
        Name of the cell type being analyzed
        
    Returns
    -------
    pd.DataFrame or None
        DataFrame with DE results, or None if insufficient cells
    """
    # Create masks for each group
    in_mask = (adata_ct.obs['motif_group'] == group1).values
    out_mask = (adata_ct.obs['motif_group'] == group2).values
    
    n_in = in_mask.sum()
    n_out = out_mask.sum()
    
    print(f"  {group1} vs {group2}: {n_in} vs {n_out} cells")
    
    # Check for sufficient cells
    if n_in < 5 or n_out < 5:
        print(f"    Skipping: too few cells in one group")
        return None
    
    # Run DE
    X = adata_ct.X
    result = differential_expression(X, in_mask, out_mask)
    
    # Create results dataframe
    de_df = pd.DataFrame({
        'gene': adata_ct.var_names,
        'logFC': result['logfoldchanges'],
        'pval': result['p_values'],
        'padj': result['p_adj'],
        'neglog10_padj': -np.log10(result['p_adj'] + 1e-300),
        'n_group1': n_in,
        'n_group2': n_out
    })
    de_df = de_df.sort_values('padj')
    
    n_sig = (de_df['padj'] < 0.1).sum()
    print(f"    {n_sig} significant genes (FDR < 0.1)")
    
    return de_df

# ============================================================================
# SECTION 6: RUN DE ANALYSIS FOR ALL CELL TYPES AND COMPARISONS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 6: Running Differential Expression Analysis")
print("=" * 70)

# Store all DE results
# Structure: de_results[cell_type][comparison] = DataFrame
de_results = {}

for ct in CELL_TYPES_OF_INTEREST:
    print(f"\n{'='*50}")
    print(f"Cell Type: {ct}")
    print(f"{'='*50}")
    
    # Subset to this cell type
    ct_mask = adata_processed.obs['cell_type'] == ct
    adata_ct = adata_processed[ct_mask].copy()
    
    if adata_ct.n_obs < 10:
        print(f"Skipping {ct}: too few cells ({adata_ct.n_obs})")
        continue
    
    # Print motif group distribution for this cell type
    group_dist = adata_ct.obs['motif_group'].value_counts()
    print(f"Motif group distribution:")
    for g, n in group_dist.items():
        print(f"  {g}: {n} cells")
    
    de_results[ct] = {}
    
    # Run DE for each comparison
    for group1, group2 in COMPARISONS:
        comparison_name = f"{group1}_vs_{group2}"
        print(f"\nComparison: {comparison_name}")
        
        de_df = run_pairwise_de(adata_ct, group1, group2, ct)
        
        if de_df is not None:
            de_results[ct][comparison_name] = de_df

# ============================================================================
# SECTION 7: SAVE DE RESULTS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 7: Saving DE Results")
print("=" * 70)

# Save individual DE results
for ct, comparisons in de_results.items():
    ct_dir = f'{OUTPUT_DIR}/{ct}'
    os.makedirs(ct_dir, exist_ok=True)
    
    for comp_name, de_df in comparisons.items():
        output_file = f'{ct_dir}/de_{comp_name}.csv'
        de_df.to_csv(output_file, index=False)
        print(f"Saved: {output_file}")

# Save combined summary
all_de = []
for ct, comparisons in de_results.items():
    for comp_name, de_df in comparisons.items():
        df = de_df.copy()
        df['cell_type'] = ct
        df['comparison'] = comp_name
        all_de.append(df)

if all_de:
    all_de_df = pd.concat(all_de, ignore_index=True)
    all_de_df.to_csv(f'{OUTPUT_DIR}/all_de_results.csv', index=False)
    print(f"\nSaved combined results: {OUTPUT_DIR}/all_de_results.csv")

# ============================================================================
# SECTION 8: VOLCANO PLOTS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 8: Generating Volcano Plots")
print("=" * 70)

# Create volcano plots for each cell type and comparison
for ct, comparisons in de_results.items():
    n_comps = len(comparisons)
    if n_comps == 0:
        continue
    
    fig, axes = plt.subplots(1, n_comps, figsize=(6*n_comps, 5))
    if n_comps == 1:
        axes = [axes]
    
    for idx, (comp_name, de_df) in enumerate(comparisons.items()):
        # Compute -log10(q) with jitter for very high values
        de_df = de_df.copy()
        de_df['neglog10_padj'] = -np.log10(de_df['padj'].clip(1e-300))
        m = de_df['neglog10_padj'] >= 300
        if m.any():
            de_df.loc[m, 'neglog10_padj'] = 300 + np.random.normal(0, 15, m.sum())
        
        ax = axes[idx]
        glm_plots.volcano_plot(
            de_df,
            x_col='logFC',
            y_col='neglog10_padj',
            label_col='gene',
            fdr=0.1,
            x_threshold=0.5,
            n_top=10,
            ax=ax
        )
        # Parse comparison name for title
        groups = comp_name.replace('_vs_', ' vs ')
        ax.set_title(f'{ct}\n{groups}', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/{ct}/volcano_plots.pdf', bbox_inches='tight', dpi=150)
    plt.savefig(f'{OUTPUT_DIR}/{ct}/volcano_plots.png', bbox_inches='tight', dpi=150)
    plt.close()
    print(f"Saved volcano plots for {ct}")

# Create combined figure with main comparison (M24_only vs M10_only) for all cell types
main_comparison = 'M24_only_vs_M10_only'
cts_with_main = [ct for ct in de_results if main_comparison in de_results[ct]]

if cts_with_main:
    n_cols = 4
    n_rows = (len(cts_with_main) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
    axes = axes.flatten()
    
    for idx, ct in enumerate(cts_with_main):
        de_df = de_results[ct][main_comparison].copy()
        
        # Compute -log10(q) with jitter
        de_df['neglog10_padj'] = -np.log10(de_df['padj'].clip(1e-300))
        m = de_df['neglog10_padj'] >= 300
        if m.any():
            de_df.loc[m, 'neglog10_padj'] = 300 + np.random.normal(0, 15, m.sum())
        
        ax = axes[idx]
        glm_plots.volcano_plot(
            de_df,
            x_col='logFC',
            y_col='neglog10_padj',
            label_col='gene',
            fdr=0.1,
            x_threshold=0.5,
            n_top=10,
            ax=ax
        )
        ax.set_title(f'{ct}\nM24-only vs M10-only', fontsize=12, fontweight='bold')
    
    # Hide empty subplots
    for idx in range(len(cts_with_main), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/volcano_M24_vs_M10_all_celltypes.pdf', bbox_inches='tight', dpi=150)
    plt.savefig(f'{OUTPUT_DIR}/volcano_M24_vs_M10_all_celltypes.png', bbox_inches='tight', dpi=150)
    plt.close()
    print(f"\nSaved combined volcano plot: {OUTPUT_DIR}/volcano_M24_vs_M10_all_celltypes.pdf")

# ============================================================================
# SECTION 9: GSEA ANALYSIS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 9: Running GSEA Analysis")
print("=" * 70)

import gseapy as gp
import time
import random

GSEAPY_AVAILABLE = True

def run_ranked_gsea_analysis(gene_ranking: pd.Series,
                             gene_sets: str = 'MSigDB_Hallmark_2020',
                             description: str = 'analysis',
                             max_retries: int = 8,
                             base_delay: float = 3.0) -> Optional[pd.DataFrame]:
    """
    Run GSEA analysis using ranked gene list with scores.
    
    Parameters
    ----------
    gene_ranking : pd.Series
        Gene names as index, ranking scores as values (e.g., logFC)
    gene_sets : str
        Gene set database to use
    description : str
        Description for logging
    max_retries : int
        Maximum number of retry attempts
    base_delay : float
        Base delay between retries (exponential backoff)
        
    Returns
    -------
    pd.DataFrame or None
        GSEA results dataframe
    """
    if len(gene_ranking) == 0:
        print("Empty gene ranking provided")
        return None
    
    print(f"  Running GSEA for {len(gene_ranking)} genes...")
    
    for attempt in range(max_retries):
        try:
            delay = base_delay * (2 ** attempt) + random.uniform(0, 1)
            if attempt > 0:
                print(f"    Retry {attempt}/{max_retries}, waiting {delay:.1f}s...")
                time.sleep(delay)
            
            gsea_res = gp.prerank(
                rnk=gene_ranking,
                gene_sets=gene_sets,
                processes=1,
                permutation_num=1000,
                outdir=None,
                no_plot=True,
                seed=42,
                verbose=False
            )
            
            if gsea_res.res2d is not None and len(gsea_res.res2d) > 0:
                print(f"    Success: {len(gsea_res.res2d)} pathways analyzed")
                return gsea_res.res2d
            else:
                print(f"    No results returned (attempt {attempt + 1})")
                
        except Exception as e:
            print(f"    Attempt {attempt + 1} failed: {str(e)[:100]}...")
            if attempt == max_retries - 1:
                print(f"    All {max_retries} attempts failed")
                return None
    
    return None


def gsea_barplot(gsea_df, title='GSEA Results', n_top=10, fdr_threshold=0.25, ax=None,
                 up_label='Up in Group1', down_label='Up in Group2'):
    """
    Plot GSEA results as a barplot.
    
    Parameters
    ----------
    gsea_df : pd.DataFrame
        GSEA results from gseapy prerank
    title : str
        Plot title
    n_top : int
        Number of top pathways to show (from each direction)
    fdr_threshold : float
        FDR threshold for significance
    ax : matplotlib axis
        Axis to plot on
    up_label : str
        Label for positive NES pathways
    down_label : str
        Label for negative NES pathways
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    
    df = gsea_df.copy()
    df['NES'] = df['NES'].astype(float)
    df['FDR q-val'] = df['FDR q-val'].astype(float)
    
    # Get top pathways from each direction
    df_pos = df[df['NES'] > 0].nsmallest(n_top, 'FDR q-val')
    df_neg = df[df['NES'] < 0].nsmallest(n_top, 'FDR q-val')
    
    df_plot = pd.concat([df_pos, df_neg]).drop_duplicates()
    df_plot = df_plot.sort_values('NES', ascending=True)
    
    if len(df_plot) == 0:
        ax.text(0.5, 0.5, 'No significant pathways', ha='center', va='center', 
                transform=ax.transAxes)
        ax.set_title(title)
        return ax
    
    # Colors based on NES direction and significance
    colors = []
    for _, row in df_plot.iterrows():
        if row['FDR q-val'] < fdr_threshold:
            colors.append('#E64B35' if row['NES'] > 0 else '#4DBBD5')
        else:
            colors.append('#CCCCCC')
    
    # Clean pathway names
    pathway_names = df_plot['Term'].apply(
        lambda x: x.replace('_', ' ')[:50] + '...' if len(x) > 50 else x.replace('_', ' ')
    )
    
    # Plot
    ax.barh(range(len(df_plot)), df_plot['NES'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(pathway_names, fontsize=9)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.axvline(0, color='black', linewidth=0.8)
    
    # Add FDR annotations
    for i, (idx, row) in enumerate(df_plot.iterrows()):
        fdr_text = f"q={row['FDR q-val']:.3f}"
        x_pos = row['NES'] + 0.05 if row['NES'] > 0 else row['NES'] - 0.05
        ha = 'left' if row['NES'] > 0 else 'right'
        ax.text(x_pos, i, fdr_text, va='center', ha=ha, fontsize=7, color='gray')
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E64B35', edgecolor='black', label=up_label),
        Patch(facecolor='#4DBBD5', edgecolor='black', label=down_label),
        Patch(facecolor='#CCCCCC', edgecolor='black', label=f'FDR ≥ {fdr_threshold}')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)
    
    return ax


# Run GSEA for each cell type and comparison
gsea_results = {}

# Gene set databases to use
GENE_SETS = ['MSigDB_Hallmark_2020']  # Can add 'GO_Biological_Process_2021', 'KEGG_2021_Human'

for ct, comparisons in de_results.items():
    print(f"\n{'='*50}")
    print(f"GSEA for {ct}")
    print(f"{'='*50}")
    
    gsea_results[ct] = {}
    
    for comp_name, de_df in comparisons.items():
        print(f"\nComparison: {comp_name}")
        
        # Create gene ranking from logFC
        de_df_clean = de_df[np.isfinite(de_df['logFC'])].copy()
        gene_ranking = de_df_clean.set_index('gene')['logFC'].sort_values(ascending=False)
        gene_ranking = gene_ranking[~gene_ranking.index.duplicated(keep='first')]
        
        print(f"  Gene ranking: {len(gene_ranking)} genes")
        
        for gene_set in GENE_SETS:
            print(f"\n  Gene set: {gene_set}")
            gsea_res = run_ranked_gsea_analysis(
                gene_ranking,
                gene_sets=gene_set,
                description=f'{ct}_{comp_name}'
            )
            
            if gsea_res is not None:
                key = f"{comp_name}_{gene_set}"
                gsea_results[ct][key] = gsea_res
                
                # Save GSEA results
                gsea_res.to_csv(f'{OUTPUT_DIR}/{ct}/gsea_{comp_name}_{gene_set}.csv', index=False)
                
                n_sig = (gsea_res['FDR q-val'].astype(float) < 0.25).sum()
                print(f"    {n_sig} significant pathways (FDR < 0.25)")

# ============================================================================
# SECTION 10: GSEA VISUALIZATION
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 10: Generating GSEA Plots")
print("=" * 70)

# Create GSEA barplots for main comparison
main_comparison = 'M24_only_vs_M10_only'
gene_set = 'MSigDB_Hallmark_2020'
main_gsea_key = f"{main_comparison}_{gene_set}"

cts_with_gsea = [ct for ct in gsea_results if main_gsea_key in gsea_results[ct]]

if cts_with_gsea:
    n_cols = 2
    n_rows = (len(cts_with_gsea) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 6*n_rows))
    axes = axes.flatten()
    
    for idx, ct in enumerate(cts_with_gsea):
        gsea_df = gsea_results[ct][main_gsea_key]
        title = f'{ct}: M24-only vs M10-only'
        gsea_barplot(
            gsea_df, title=title, ax=axes[idx],
            up_label='Up in M24-only (Healthy)',
            down_label='Up in M10-only (Tumor)'
        )
    
    # Hide empty subplots
    for idx in range(len(cts_with_gsea), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/gsea_M24_vs_M10_all_celltypes.pdf', bbox_inches='tight', dpi=150)
    plt.savefig(f'{OUTPUT_DIR}/gsea_M24_vs_M10_all_celltypes.png', bbox_inches='tight', dpi=150)
    plt.close()
    print(f"Saved combined GSEA plot: {OUTPUT_DIR}/gsea_M24_vs_M10_all_celltypes.pdf")

# ============================================================================
# SECTION 11: SUMMARY STATISTICS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 11: Summary Statistics")
print("=" * 70)

# Create summary table
summary_rows = []
for ct, comparisons in de_results.items():
    for comp_name, de_df in comparisons.items():
        n_up = ((de_df['padj'] < 0.1) & (de_df['logFC'] > 0.5)).sum()
        n_down = ((de_df['padj'] < 0.1) & (de_df['logFC'] < -0.5)).sum()
        n_total = de_df['n_group1'].iloc[0] + de_df['n_group2'].iloc[0]
        
        summary_rows.append({
            'cell_type': ct,
            'comparison': comp_name,
            'n_cells_total': n_total,
            'n_group1': de_df['n_group1'].iloc[0],
            'n_group2': de_df['n_group2'].iloc[0],
            'n_sig_up': n_up,
            'n_sig_down': n_down,
            'n_sig_total': n_up + n_down
        })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f'{OUTPUT_DIR}/de_summary.csv', index=False)
print("\nDE Analysis Summary:")
print(summary_df.to_string())

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
print(f"\nAll results saved to: {OUTPUT_DIR}")
print("\nKey output files:")
print(f"  - {OUTPUT_DIR}/all_de_results.csv (combined DE results)")
print(f"  - {OUTPUT_DIR}/de_summary.csv (summary statistics)")
print(f"  - {OUTPUT_DIR}/volcano_M24_vs_M10_all_celltypes.pdf")
print(f"  - {OUTPUT_DIR}/gsea_M24_vs_M10_all_celltypes.pdf")
print(f"  - {OUTPUT_DIR}/[cell_type]/ (individual cell type results)")