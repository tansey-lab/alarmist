"""
3-Point Module Score Analysis: M24-only vs Double vs M10-only
==============================================================

This script calculates gene module scores for key pathways and generates
boxplots comparing the three motif groups to visualize intermediate states.

Requirements:
- adata with motif_10_state and motif_24_state columns (values: 'positive'/'negative')
- Normalized, log-transformed expression data
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import gseapy as gp
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Paths - MODIFY THESE AS NEEDED
ADATA_PATH = '../results/AIS_LUAD/adata_with_motifs.h5ad'  # Your adata file
OUTPUT_DIR = '../results/AIS_LUAD/module_score_analysis'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Cell types of interest
CELL_TYPES_OF_INTEREST = ['Tumor_epi', 'Macro', 'T', 'B', 'Neutro', 'Mast', 'pDC', 'cDC']

# Pathways to analyze (from MSigDB Hallmark)
PATHWAYS_OF_INTEREST = [
    'Interferon Gamma Response',
    'Interferon Alpha Response', 
    'Epithelial Mesenchymal Transition',
    'Glycolysis',
    'G2-M Checkpoint',
    'Inflammatory Response',
    'TGF-beta Signaling',
    'KRAS Signaling Dn',
]

# Motif columns
MOTIF_24_COL = 'motif_24_state'
MOTIF_10_COL = 'motif_10_state'

# ============================================================================
# SECTION 1: LOAD DATA
# ============================================================================

print("=" * 70)
print("SECTION 1: Loading Data")
print("=" * 70)

# Load adata
adata = sc.read_h5ad(ADATA_PATH)
print(f"Loaded adata: {adata.n_obs} cells, {adata.n_vars} genes")

# ============================================================================
# SECTION 2: DEFINE MOTIF GROUPS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 2: Defining Motif Groups")
print("=" * 70)

# Create boolean masks
motif_24_positive = adata.obs[MOTIF_24_COL] == 'positive'
motif_10_positive = adata.obs[MOTIF_10_COL] == 'positive'

# Define three mutually exclusive groups
adata.obs['motif_group'] = 'Other'
adata.obs.loc[motif_24_positive & ~motif_10_positive, 'motif_group'] = 'M24_only'
adata.obs.loc[motif_24_positive & motif_10_positive, 'motif_group'] = 'Double'
adata.obs.loc[~motif_24_positive & motif_10_positive, 'motif_group'] = 'M10_only'

# Print distribution
print("\nOverall motif group distribution:")
print(adata.obs['motif_group'].value_counts())

print("\nMotif group distribution by cell type:")
for ct in CELL_TYPES_OF_INTEREST:
    ct_mask = adata.obs['cell_type'] == ct
    if ct_mask.sum() > 0:
        print(f"\n{ct}:")
        print(adata.obs.loc[ct_mask, 'motif_group'].value_counts().to_string())

# ============================================================================
# SECTION 3: GET GENE SETS FROM MSIGDB
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 3: Fetching Gene Sets from MSigDB")
print("=" * 70)

# Download Hallmark gene sets
try:
    hallmark_genesets = gp.get_library('MSigDB_Hallmark_2020')
    print(f"Downloaded {len(hallmark_genesets)} Hallmark gene sets")
except Exception as e:
    print(f"Error downloading gene sets: {e}")
    print("Trying alternative method...")
    # Alternative: use enrichr libraries
    hallmark_genesets = gp.get_library_name()
    hallmark_genesets = gp.get_library('MSigDB_Hallmark_2020')

# Print available gene sets for reference
print(f"\nAvailable gene sets:")
for k in hallmark_genesets.keys():
    print(f"  - {k}")

# Extract gene lists for pathways of interest
# gseapy uses the same names as PATHWAYS_OF_INTEREST (e.g., 'Interferon Gamma Response')
pathway_genes = {}
for pathway in PATHWAYS_OF_INTEREST:
    if pathway in hallmark_genesets:
        genes = hallmark_genesets[pathway]
        # Filter to genes in adata
        genes_in_adata = [g for g in genes if g in adata.var_names]
        pathway_genes[pathway] = genes_in_adata
        print(f"{pathway}: {len(genes_in_adata)}/{len(genes)} genes found in adata")
    else:
        # Try case-insensitive match
        matched = False
        for key in hallmark_genesets.keys():
            if pathway.lower() == key.lower():
                genes = hallmark_genesets[key]
                genes_in_adata = [g for g in genes if g in adata.var_names]
                pathway_genes[pathway] = genes_in_adata
                print(f"{pathway}: {len(genes_in_adata)}/{len(genes)} genes found in adata")
                matched = True
                break
        if not matched:
            print(f"Warning: {pathway} not found in gene sets")

# ============================================================================
# SECTION 4: CALCULATE MODULE SCORES
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 4: Calculating Module Scores")
print("=" * 70)

# Check if data is normalized
if adata.X.max() > 100:
    print("Data appears to be raw counts. Normalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

# Calculate module scores for each pathway
for pathway, genes in pathway_genes.items():
    score_name = pathway.replace(' ', '_')
    print(f"Calculating score for: {pathway} ({len(genes)} genes)")
    
    sc.tl.score_genes(
        adata, 
        gene_list=genes, 
        score_name=score_name,
        use_raw=False
    )

print("\nModule scores added to adata.obs")

# ============================================================================
# SECTION 5: VIOLIN PLOT FUNCTIONS
# ============================================================================

def plot_module_boxplot(adata, score_col, cell_type, ax, title=None):
    """
    Plot violin plot of module scores for three motif groups.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with motif_group and score columns
    score_col : str
        Column name for module score in adata.obs
    cell_type : str
        Cell type to subset
    ax : matplotlib axis
        Axis to plot on
    title : str, optional
        Plot title
    """
    # Subset to cell type and motif groups of interest
    mask = (adata.obs['cell_type'] == cell_type) & \
           (adata.obs['motif_group'].isin(['M24_only', 'Double', 'M10_only']))
    
    df = adata.obs.loc[mask, ['motif_group', score_col]].copy()
    
    if len(df) < 10:
        ax.text(0.5, 0.5, f'Too few cells\n(n={len(df)})', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title or score_col)
        return
    
    # Order groups for progression visualization
    group_order = ['M24_only', 'Double', 'M10_only']
    df['motif_group'] = pd.Categorical(df['motif_group'], categories=group_order, ordered=True)
    
    # Colors: RED for M24 (healthy), PURPLE for Double, BLUE for M10 (tumor)
    colors = {'M24_only': '#E64B35', 'Double': '#9B59B6', 'M10_only': '#4DBBD5'}
    
    # Create violin plot
    sns.violinplot(
        data=df, 
        x='motif_group', 
        y=score_col,
        order=group_order,
        palette=colors,
        ax=ax,
        width=0.8,
        inner='box',
        cut=0,
    )
    
    # Add sample sizes
    y_min, y_max = df[score_col].min(), df[score_col].max()
    y_range = y_max - y_min
    for i, group in enumerate(group_order):
        n = (df['motif_group'] == group).sum()
        ax.text(i, y_min - 0.22 * y_range, f'n={n}', ha='center', va='top', fontsize=8)
    
    # Labels - separate rows to avoid overlap
    ax.set_xlabel('')
    ax.set_ylabel('Module Score')
    ax.set_title(title or score_col, fontsize=11, fontweight='bold')
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['M24-only', 'Double', 'M10-only'], fontsize=9)
    
    # Add second row of labels using text
    ax.text(0, y_min - 0.12 * y_range, '(Healthy)', ha='center', fontsize=8, color='gray')
    ax.text(1, y_min - 0.12 * y_range, '(Interm.)', ha='center', fontsize=8, color='gray')
    ax.text(2, y_min - 0.12 * y_range, '(Tumor)', ha='center', fontsize=8, color='gray')
    
    # Adjust y-axis
    ax.set_ylim(y_min - 0.25 * y_range, y_max + 0.1 * y_range)
    
    # Add trend line (connecting medians)
    medians = [df[df['motif_group'] == g][score_col].median() for g in group_order]
    ax.plot(range(3), medians, 'k--', alpha=0.5, linewidth=1.5)


def plot_module_boxplot_with_stats(adata, score_col, cell_type, ax, title=None):
    """
    Plot violin plot with statistical comparisons (all three pairwise tests).
    """
    from scipy import stats
    
    # Subset to cell type and motif groups of interest
    mask = (adata.obs['cell_type'] == cell_type) & \
           (adata.obs['motif_group'].isin(['M24_only', 'Double', 'M10_only']))
    
    df = adata.obs.loc[mask, ['motif_group', score_col]].copy()
    
    if len(df) < 10:
        ax.text(0.5, 0.5, f'Too few cells\n(n={len(df)})', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title or score_col)
        return
    
    # Order groups
    group_order = ['M24_only', 'Double', 'M10_only']
    df['motif_group'] = pd.Categorical(df['motif_group'], categories=group_order, ordered=True)
    
    # Colors: RED for M24 (healthy), PURPLE for Double, BLUE for M10 (tumor)
    colors = {'M24_only': '#E64B35', 'Double': '#9B59B6', 'M10_only': '#4DBBD5'}
    
    # Create violin plot
    sns.violinplot(
        data=df, 
        x='motif_group', 
        y=score_col,
        order=group_order,
        palette=colors,
        ax=ax,
        width=0.8,
        inner='box',
        cut=0,
    )
    
    # Statistical tests (Mann-Whitney U)
    groups_data = {g: df[df['motif_group'] == g][score_col].values for g in group_order}
    
    # Test M24 vs Double
    if len(groups_data['M24_only']) > 5 and len(groups_data['Double']) > 5:
        _, p1 = stats.mannwhitneyu(groups_data['M24_only'], groups_data['Double'])
    else:
        p1 = 1.0
    
    # Test Double vs M10
    if len(groups_data['Double']) > 5 and len(groups_data['M10_only']) > 5:
        _, p2 = stats.mannwhitneyu(groups_data['Double'], groups_data['M10_only'])
    else:
        p2 = 1.0
    
    # Test M24 vs M10 (overall trend)
    if len(groups_data['M24_only']) > 5 and len(groups_data['M10_only']) > 5:
        _, p3 = stats.mannwhitneyu(groups_data['M24_only'], groups_data['M10_only'])
    else:
        p3 = 1.0
    
    def p_to_stars(p):
        if p < 0.001: return '***'
        elif p < 0.01: return '**'
        elif p < 0.05: return '*'
        else: return 'ns'
    
    # Get y range for annotation positioning
    y_min, y_max = df[score_col].min(), df[score_col].max()
    y_range = y_max - y_min
    
    # Option B: Show all three comparisons
    # Level 1: M24 vs Double (lowest bar)
    y1 = y_max + 0.08 * y_range
    ax.plot([0, 1], [y1, y1], 'k-', linewidth=1)
    ax.text(0.5, y1 + 0.02 * y_range, p_to_stars(p1), ha='center', fontsize=9)
    
    # Level 2: Double vs M10 (same level as level 1)
    ax.plot([1, 2], [y1, y1], 'k-', linewidth=1)
    ax.text(1.5, y1 + 0.02 * y_range, p_to_stars(p2), ha='center', fontsize=9)
    
    # Level 3: M24 vs M10 (higher bar spanning all)
    y2 = y_max + 0.22 * y_range
    ax.plot([0, 2], [y2, y2], 'k-', linewidth=1)
    ax.text(1, y2 + 0.02 * y_range, p_to_stars(p3), ha='center', fontsize=9)
    
    # Adjust y-axis to fit annotations and n labels
    ax.set_ylim(y_min - 0.25 * y_range, y_max + 0.35 * y_range)
    
    # Add sample sizes below x-axis - move further down
    for i, group in enumerate(group_order):
        n = (df['motif_group'] == group).sum()
        ax.text(i, y_min - 0.22 * y_range, f'n={n}', ha='center', va='top', fontsize=8)
    
    # Labels - remove line breaks to avoid overlap
    ax.set_xlabel('')
    ax.set_ylabel('Module Score')
    ax.set_title(title or score_col, fontsize=11, fontweight='bold')
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['M24-only', 'Double', 'M10-only'], fontsize=9)
    
    # Add second row of labels using text
    ax.text(0, y_min - 0.12 * y_range, '(Healthy)', ha='center', fontsize=8, color='gray')
    ax.text(1, y_min - 0.12 * y_range, '(Interm.)', ha='center', fontsize=8, color='gray')
    ax.text(2, y_min - 0.12 * y_range, '(Tumor)', ha='center', fontsize=8, color='gray')
    
    # Add trend line connecting medians
    medians = [df[df['motif_group'] == g][score_col].median() for g in group_order]
    ax.plot(range(3), medians, 'k--', alpha=0.5, linewidth=1.5)

# ============================================================================
# SECTION 6: GENERATE PLOTS - BY PATHWAY (ALL CELL TYPES)
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 6: Generating Violin Plots by Pathway")
print("=" * 70)

for pathway in PATHWAYS_OF_INTEREST:
    score_col = pathway.replace(' ', '_')
    
    if score_col not in adata.obs.columns:
        print(f"Skipping {pathway}: score not calculated")
        continue
    
    print(f"\nPlotting: {pathway}")
    
    # Get cell types that have enough cells in all three groups
    valid_cts = []
    for ct in CELL_TYPES_OF_INTEREST:
        mask = (adata.obs['cell_type'] == ct) & \
               (adata.obs['motif_group'].isin(['M24_only', 'Double', 'M10_only']))
        group_counts = adata.obs.loc[mask, 'motif_group'].value_counts()
        if all(group_counts.get(g, 0) >= 5 for g in ['M24_only', 'Double', 'M10_only']):
            valid_cts.append(ct)
    
    if not valid_cts:
        print(f"  No cell types with sufficient cells in all groups")
        continue
    
    # Create figure - taller and thinner subplots
    n_cols = 4
    n_rows = (len(valid_cts) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.5*n_cols, 5*n_rows))
    axes = axes.flatten()
    
    for idx, ct in enumerate(valid_cts):
        plot_module_boxplot_with_stats(
            adata, 
            score_col, 
            ct, 
            axes[idx],
            title=ct
        )
    
    # Hide empty subplots
    for idx in range(len(valid_cts), len(axes)):
        axes[idx].set_visible(False)
    
    fig.suptitle(f'{pathway}', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    # Save
    filename = pathway.replace(' ', '_')
    plt.savefig(f'{OUTPUT_DIR}/violin_{filename}_all_celltypes.pdf', bbox_inches='tight', dpi=150)
    plt.savefig(f'{OUTPUT_DIR}/violin_{filename}_all_celltypes.png', bbox_inches='tight', dpi=150)
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/violin_{filename}_all_celltypes.pdf")

# ============================================================================
# SECTION 7: GENERATE PLOTS - BY CELL TYPE (ALL PATHWAYS)
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 7: Generating Violin Plots by Cell Type")
print("=" * 70)

for ct in CELL_TYPES_OF_INTEREST:
    print(f"\nPlotting: {ct}")
    
    # Check if cell type has enough cells
    mask = (adata.obs['cell_type'] == ct) & \
           (adata.obs['motif_group'].isin(['M24_only', 'Double', 'M10_only']))
    
    if mask.sum() < 30:
        print(f"  Skipping: too few cells ({mask.sum()})")
        continue
    
    # Get valid pathways (those with calculated scores)
    valid_pathways = [p for p in PATHWAYS_OF_INTEREST 
                      if p.replace(' ', '_') in adata.obs.columns]
    
    if not valid_pathways:
        print(f"  No valid pathways")
        continue
    
    # Create figure - taller and thinner subplots
    n_cols = 3
    n_rows = (len(valid_pathways) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.5*n_cols, 5*n_rows))
    axes = axes.flatten()
    
    for idx, pathway in enumerate(valid_pathways):
        score_col = pathway.replace(' ', '_')
        plot_module_boxplot_with_stats(
            adata,
            score_col,
            ct,
            axes[idx],
            title=pathway
        )
    
    # Hide empty subplots
    for idx in range(len(valid_pathways), len(axes)):
        axes[idx].set_visible(False)
    
    fig.suptitle(f'{ct}', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    # Save
    plt.savefig(f'{OUTPUT_DIR}/violin_{ct}_all_pathways.pdf', bbox_inches='tight', dpi=150)
    plt.savefig(f'{OUTPUT_DIR}/violin_{ct}_all_pathways.png', bbox_inches='tight', dpi=150)
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/violin_{ct}_all_pathways.pdf")

# ============================================================================
# SECTION 8: SUMMARY FIGURE - KEY PATHWAYS x KEY CELL TYPES
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 8: Generating Summary Violin Plot")
print("=" * 70)

# Key combinations to highlight
key_combos = [
    ('Interferon Gamma Response', 'T'),
    ('Interferon Gamma Response', 'pDC'),
    ('Interferon Alpha Response', 'cDC'),
    ('Epithelial Mesenchymal Transition', 'Macro'),
    ('Epithelial Mesenchymal Transition', 'Mast'),
    ('Glycolysis', 'Tumor_epi'),
    ('G2-M Checkpoint', 'Tumor_epi'),
]

# Filter to valid combinations
valid_combos = []
for pathway, ct in key_combos:
    score_col = pathway.replace(' ', '_')
    if score_col in adata.obs.columns:
        mask = (adata.obs['cell_type'] == ct) & \
               (adata.obs['motif_group'].isin(['M24_only', 'Double', 'M10_only']))
        if mask.sum() >= 30:
            valid_combos.append((pathway, ct))

if valid_combos:
    n_cols = 4
    n_rows = (len(valid_combos) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.5*n_cols, 5*n_rows))
    axes = axes.flatten()
    
    for idx, (pathway, ct) in enumerate(valid_combos):
        score_col = pathway.replace(' ', '_')
        plot_module_boxplot_with_stats(
            adata,
            score_col,
            ct,
            axes[idx],
            title=f'{ct}\n{pathway}'
        )
    
    # Hide empty subplots
    for idx in range(len(valid_combos), len(axes)):
        axes[idx].set_visible(False)
    
    fig.suptitle('Key Pathway-Cell Type Combinations', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    plt.savefig(f'{OUTPUT_DIR}/violin_summary_key_combos.pdf', bbox_inches='tight', dpi=150)
    plt.savefig(f'{OUTPUT_DIR}/violin_summary_key_combos.png', bbox_inches='tight', dpi=150)
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/violin_summary_key_combos.pdf")

# ============================================================================
# SECTION 9: SAVE MODULE SCORES
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 9: Saving Module Scores")
print("=" * 70)

# Save module scores to CSV
score_cols = [p.replace(' ', '_') for p in PATHWAYS_OF_INTEREST if p.replace(' ', '_') in adata.obs.columns]
scores_df = adata.obs[['cell_type', 'motif_group', 'condition', 'patient'] + score_cols].copy()
scores_df.to_csv(f'{OUTPUT_DIR}/module_scores.csv')
print(f"Saved: {OUTPUT_DIR}/module_scores.csv")

# ============================================================================
# COMPLETE
# ============================================================================

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
print(f"\nOutput directory: {OUTPUT_DIR}")
print("\nGenerated files:")
print("  - violin_[Pathway]_all_celltypes.pdf  (one per pathway)")
print("  - violin_[CellType]_all_pathways.pdf  (one per cell type)")
print("  - violin_summary_key_combos.pdf       (key combinations)")
print("  - module_scores.csv                    (all scores)")