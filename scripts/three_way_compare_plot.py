"""
Plot All Three Comparisons from Saved DE and GSEA Results
==========================================================

This script reads the saved DE and GSEA results and generates combined plots
for all three comparisons:
1. M24_only vs M10_only
2. M24_only vs Double
3. Double vs M10_only
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Patch
from alarmist.plotting import glm_plots

# ============================================================================
# CONFIGURATION
# ============================================================================

RESULTS_DIR = '../results/AIS_LUAD'
OUTPUT_DIR = f'{RESULTS_DIR}/motif_de_analysis'

CELL_TYPES_OF_INTEREST = ['Tumor_epi', 'Macro', 'T', 'B', 'Neutro', 'Mast', 'pDC', 'cDC']

COMPARISONS = [
    'M24_only_vs_M10_only',
    'M24_only_vs_Double',
    'Double_vs_M10_only'
]

# Labels for plots
COMPARISON_LABELS = {
    'M24_only_vs_M10_only': ('M24-only (Healthy)', 'M10-only (Tumor)'),
    'M24_only_vs_Double': ('M24-only (Healthy)', 'Double (Intermediate)'),
    'Double_vs_M10_only': ('Double (Intermediate)', 'M10-only (Tumor)')
}

GENE_SET = 'MSigDB_Hallmark_2020'

# ============================================================================
# SECTION 1: LOAD SAVED DE RESULTS
# ============================================================================

print("=" * 70)
print("SECTION 1: Loading Saved DE Results")
print("=" * 70)

de_results = {}

for ct in CELL_TYPES_OF_INTEREST:
    ct_dir = f'{OUTPUT_DIR}/{ct}'
    if not os.path.exists(ct_dir):
        print(f"Skipping {ct}: directory not found")
        continue
    
    de_results[ct] = {}
    
    for comp in COMPARISONS:
        de_file = f'{ct_dir}/de_{comp}.csv'
        if os.path.exists(de_file):
            de_df = pd.read_csv(de_file)
            de_results[ct][comp] = de_df
            print(f"Loaded: {ct} - {comp} ({len(de_df)} genes)")
        else:
            print(f"Not found: {de_file}")

# ============================================================================
# SECTION 2: LOAD SAVED GSEA RESULTS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 2: Loading Saved GSEA Results")
print("=" * 70)

gsea_results = {}

for ct in CELL_TYPES_OF_INTEREST:
    ct_dir = f'{OUTPUT_DIR}/{ct}'
    if not os.path.exists(ct_dir):
        continue
    
    gsea_results[ct] = {}
    
    for comp in COMPARISONS:
        gsea_file = f'{ct_dir}/gsea_{comp}_{GENE_SET}.csv'
        if os.path.exists(gsea_file):
            gsea_df = pd.read_csv(gsea_file)
            gsea_results[ct][comp] = gsea_df
            print(f"Loaded: {ct} - {comp} ({len(gsea_df)} pathways)")
        else:
            print(f"Not found: {gsea_file}")

# ============================================================================
# SECTION 3: GSEA BARPLOT FUNCTION
# ============================================================================

def gsea_barplot(gsea_df, title='GSEA Results', n_top=10, fdr_threshold=0.25, ax=None,
                 up_label='Up in Group1', down_label='Up in Group2'):
    """
    Plot GSEA results as a barplot.
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
    legend_elements = [
        Patch(facecolor='#E64B35', edgecolor='black', label=up_label),
        Patch(facecolor='#4DBBD5', edgecolor='black', label=down_label),
        Patch(facecolor='#CCCCCC', edgecolor='black', label=f'FDR ≥ {fdr_threshold}')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)
    
    return ax

# ============================================================================
# SECTION 4: COMBINED VOLCANO PLOTS FOR ALL COMPARISONS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 4: Generating Combined Volcano Plots")
print("=" * 70)

for comp in COMPARISONS:
    print(f"\nProcessing: {comp}")
    
    # Get cell types that have this comparison
    cts_with_comp = [ct for ct in de_results if comp in de_results[ct]]
    
    if not cts_with_comp:
        print(f"  No cell types have results for {comp}")
        continue
    
    # Create figure
    n_cols = 4
    n_rows = (len(cts_with_comp) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
    axes = axes.flatten()
    
    for idx, ct in enumerate(cts_with_comp):
        de_df = de_results[ct][comp].copy()
        
        # Compute -log10(q) with jitter for very high values
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
        
        # Title with comparison name
        comp_display = comp.replace('_vs_', ' vs ').replace('_', '-')
        ax.set_title(f'{ct}\n{comp_display}', fontsize=12, fontweight='bold')
    
    # Hide empty subplots
    for idx in range(len(cts_with_comp), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    
    # Save
    output_file = f'{OUTPUT_DIR}/volcano_{comp}_all_celltypes'
    plt.savefig(f'{output_file}.pdf', bbox_inches='tight', dpi=150)
    plt.savefig(f'{output_file}.png', bbox_inches='tight', dpi=150)
    plt.close()
    print(f"  Saved: {output_file}.pdf")

# ============================================================================
# SECTION 5: COMBINED GSEA PLOTS FOR ALL COMPARISONS
# ============================================================================

print("\n" + "=" * 70)
print("SECTION 5: Generating Combined GSEA Plots")
print("=" * 70)

for comp in COMPARISONS:
    print(f"\nProcessing: {comp}")
    
    # Get cell types that have GSEA results for this comparison
    cts_with_gsea = [ct for ct in gsea_results if comp in gsea_results[ct]]
    
    if not cts_with_gsea:
        print(f"  No cell types have GSEA results for {comp}")
        continue
    
    # Create figure
    n_cols = 2
    n_rows = (len(cts_with_gsea) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 6*n_rows))
    axes = axes.flatten()
    
    # Get labels for this comparison
    up_label, down_label = COMPARISON_LABELS[comp]
    
    for idx, ct in enumerate(cts_with_gsea):
        gsea_df = gsea_results[ct][comp]
        
        comp_display = comp.replace('_vs_', ' vs ').replace('_', '-')
        title = f'{ct}: {comp_display}'
        
        gsea_barplot(
            gsea_df, 
            title=title, 
            ax=axes[idx],
            up_label=f'Up in {up_label}',
            down_label=f'Up in {down_label}'
        )
    
    # Hide empty subplots
    for idx in range(len(cts_with_gsea), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    
    # Save
    output_file = f'{OUTPUT_DIR}/gsea_{comp}_all_celltypes'
    plt.savefig(f'{output_file}.pdf', bbox_inches='tight', dpi=150)
    plt.savefig(f'{output_file}.png', bbox_inches='tight', dpi=150)
    plt.close()
    print(f"  Saved: {output_file}.pdf")

# ============================================================================
# SECTION 6: SUMMARY
# ============================================================================

print("\n" + "=" * 70)
print("PLOTTING COMPLETE")
print("=" * 70)

print("\nGenerated files:")
for comp in COMPARISONS:
    print(f"\n{comp}:")
    print(f"  - {OUTPUT_DIR}/volcano_{comp}_all_celltypes.pdf")
    print(f"  - {OUTPUT_DIR}/gsea_{comp}_all_celltypes.pdf")