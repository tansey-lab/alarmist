#!/usr/bin/env python3
"""
BPTF Visualization Utilities

Helper functions for creating all BPTF analysis visualizations.
Adapted from epithelioid_spatial_new/notebooks/bptf_analysis_50 copy.ipynb
"""

import io
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
from PIL import Image
from pathlib import Path
from typing import List, Optional, Dict, Tuple
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch, Rectangle
import networkx as nx
from graphviz import Digraph


def parse_lri_full(lri_name):
    """Parse LRI name into components"""
    parts = lri_name.split('|')
    if len(parts) >= 5:
        return parts[0], parts[1], parts[2], parts[3], parts[4]
    elif len(parts) == 4:
        c1, c2, ligand, receptor = parts
        mode = 'autocrine' if c1 == c2 else 'paracrine'
        return c1, c2, ligand, receptor, mode
    elif len(parts) == 2:
        return 'unknown', 'unknown', parts[0], parts[1], 'unknown'
    else:
        return 'unknown', 'unknown', lri_name, lri_name, 'unknown'


def get_cell_type_colors(unique_ct):
    """Generate color map for cell types"""
    ct_cmap = plt.get_cmap('tab20', len(unique_ct))
    return {ct: ct_cmap(i) for i, ct in enumerate(unique_ct)}


def plot_cells_per_patch(patch_metadata_df, save_path=None):
    """Plot distribution of cells per patch"""
    plt.figure(figsize=(8, 6))
    plt.hist(patch_metadata_df['n_cells'], bins=50, edgecolor='black')
    plt.xlabel('Number of cells per patch')
    plt.ylabel('Number of patches')
    plt.title('Distribution of cells per patch')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return plt.gcf()


def plot_factor_distributions(patch_loadings, lri_factors, save_path=None):
    """Plot LRI participation and patch loading distributions"""
    # Compute data
    lri_participation = lri_factors.sum(axis=0)
    patch_loading_totals = patch_loadings.sum(axis=1)

    # Create side-by-side subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left subplot: LRI participation distribution
    axes[0].hist(lri_participation, bins=100, alpha=0.7, edgecolor='black')
    axes[0].set_title('LRI Factor Distribution')
    axes[0].set_xlabel('Total Factor')
    axes[0].set_ylabel('Number of LRIs (log)')
    axes[0].set_yscale('log')
    axes[0].grid(True)

    # Right subplot: Patch loading distribution
    axes[1].hist(patch_loading_totals, bins=100, alpha=0.7, edgecolor='black')
    axes[1].set_title('Patch Loading Distribution')
    axes[1].set_xlabel('Total Loading')
    axes[1].set_ylabel('Number of Patches')
    axes[1].set_yscale('log')
    axes[1].grid(True)

    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return fig


def plot_factor_sparsity(patch_loadings, lri_factors, save_path=None):
    """Plot factor matrix sparsity"""
    plt.figure(figsize=(8, 4))
    sparsity_threshold = 1e-6
    patch_sparsity = (patch_loadings < sparsity_threshold).mean()
    lri_sparsity = (lri_factors < sparsity_threshold).mean()
    plt.bar(['Patch Loadings', 'LRI Factors'], [patch_sparsity, lri_sparsity])
    plt.title('Factor Matrix Sparsity')
    plt.ylabel('Fraction of < 1e-6')
    plt.ylim(0, 1)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return plt.gcf()


def plot_motif_activities(patch_loadings, save_path=None):
    """Plot motif activity rankings"""
    motif_activities = patch_loadings.sum(axis=0)
    ranked_indices = np.argsort(motif_activities)[::-1]
    ranked_motifs = {
        'motif_indices': ranked_indices.tolist(),
        'activities': motif_activities[ranked_indices].tolist(),
        'activity_fractions': (motif_activities[ranked_indices] / motif_activities.sum()).tolist()
    }
    
    print("Motif Activity Rankings:")
    for i, (motif_idx, activity, fraction) in enumerate(zip(
            ranked_motifs['motif_indices'],
            ranked_motifs['activities'],
            ranked_motifs['activity_fractions']
        )):
        print(f"  {i+1}. Motif {motif_idx}: {fraction:.1%} of total activity")
    
    # Create bar plot
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(ranked_motifs['motif_indices'])), ranked_motifs['activity_fractions'])
    plt.xlabel('Motif Rank')
    plt.ylabel('Activity Fraction')
    plt.title('Motif Activity Distribution')
    plt.xticks(range(len(ranked_motifs['motif_indices'])), 
               [f"M{i}" for i in ranked_motifs['motif_indices']], rotation=45)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return plt.gcf()


def plot_lri_clustermap(lri_factors, column_names, save_path=None):
    """Create clustermap of LRI cell type pairs across motifs"""
    # Construct LRI x motif DataFrame
    H_obs = lri_factors  # shape (n_motifs, n_lris)
    n_motifs = H_obs.shape[0]
    
    df_all = pd.DataFrame(
        H_obs.T,
        index=column_names,
        columns=range(n_motifs)
    ).reset_index().rename(columns={'index':'lri'})

    # Extract cell_pair (first two parts: celltype1-celltype2)
    df_all['cell_pair'] = df_all['lri'].str.split('|').str[:2].str.join('|')

    # Aggregate by cell_pair
    pivot = df_all.groupby('cell_pair')[list(range(n_motifs))].sum()
    
    # Filter out granulocyte pairs
    pivot_filled = pivot.fillna(0)
    pivot_filled = pivot_filled[~pivot_filled.index.str.contains('granulocyte', case=False, na=False)]

    # Create clustermap
    g = sns.clustermap(
        pivot_filled,
        cmap='Blues',
        norm=LogNorm(vmin=pivot_filled[pivot_filled>0].min().min(),
                     vmax=pivot_filled.max().max()),
        linewidths=0.5,
        cbar_kws={'label': 'Activity Value (log scale)'},
        figsize=(20, 40),
        cbar_pos=(1.02, 0.2, 0.03, 0.6),
        xticklabels=True,
        yticklabels=True
    )

    g.ax_heatmap.set_xlabel('Motif')
    g.ax_heatmap.set_ylabel('Ligand-Receptor Cell Type Pair')
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    plt.suptitle('Clustermap of LRI Cell Type Pairs Across Motifs', y=0.98)
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"Saved: {save_path}")
    
    return g


def plot_celltype_communication_heatmap(lri_factors, column_names, save_path=None):
    """Create cell-cell communication activity heatmap"""
    # Same setup as clustermap
    H_obs = lri_factors
    n_motifs = H_obs.shape[0]
    
    df_all = pd.DataFrame(
        H_obs.T,
        index=column_names,
        columns=range(n_motifs)
    ).reset_index().rename(columns={'index':'lri'})

    df_all['cell_pair'] = df_all['lri'].str.split('|').str[:2].str.join('|')
    pivot = df_all.groupby('cell_pair')[list(range(n_motifs))].sum()
    pivot_filled = pivot.fillna(0)
    pivot_filled = pivot_filled[~pivot_filled.index.str.contains('granulocyte', case=False, na=False)]

    # Calculate total activity per cell pair
    cellpair_activity_sum = pivot_filled.sum(axis=1)

    # Split cell pairs into ligand and receptor cell types
    ligand_cts = []
    receptor_cts = []
    activities = []

    for cellpair, activity in cellpair_activity_sum.items():
        if '|' in cellpair:
            parts = cellpair.split('|')
        else:
            parts = [cellpair[:len(cellpair)//2], cellpair[len(cellpair)//2:]]
        
        if len(parts) >= 2:
            ligand_ct = parts[0].strip()
            receptor_ct = parts[1].strip()
            
            ligand_cts.append(ligand_ct)
            receptor_cts.append(receptor_ct)
            activities.append(activity)

    # Create DataFrame and pivot table
    df_split = pd.DataFrame({
        'ligand_ct': ligand_cts,
        'receptor_ct': receptor_cts,
        'activity': activities
    })

    heatmap_data = df_split.pivot_table(
        index='ligand_ct', 
        columns='receptor_ct', 
        values='activity', 
        fill_value=0,
        aggfunc='sum'
    )

    # Create heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        heatmap_data,
        cmap='Blues',
        norm=LogNorm(vmin=pivot_filled[pivot_filled>0].min().min(),
                    vmax=pivot_filled.max().max()),
        cbar_kws={'label': 'Total Activity Across Motifs (log scale)'},
        linewidths=0.5,
        square=True,
        xticklabels=True,
        yticklabels=True
    )

    plt.title('Cell-Cell Communication Activity Heatmap\n(Ligand Cell Type x Receptor Cell Type)', 
              fontsize=14, fontweight='bold')
    plt.xlabel('Receptor Cell Type', fontsize=12)
    plt.ylabel('Ligand Cell Type', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return plt.gcf()


def plot_celltype_communication_by_motif(lri_factors, column_names, suffix="", save_path=None):
    """Create cell-cell communication heatmaps for each motif"""
    H_obs = lri_factors
    n_motifs = H_obs.shape[0]
    
    df_all = pd.DataFrame(
        H_obs.T,
        index=column_names,
        columns=range(n_motifs)
    ).reset_index().rename(columns={'index':'lri'})

    df_all['cell_pair'] = df_all['lri'].str.split('|').str[:2].str.join('|')
    pivot = df_all.groupby('cell_pair')[list(range(n_motifs))].sum()
    pivot_filled = pivot.fillna(0)
    filtered_pivot = pivot_filled[~pivot_filled.index.str.contains('granulocyte', case=False, na=False)]

    def create_cellpair_heatmap(motif_data, motif_name):
        """Create heatmap for single motif"""
        ligand_cts = []
        receptor_cts = []
        activities = []
        
        for cellpair, activity in motif_data.items():
            if activity == 0 or pd.isna(activity):
                continue
                
            if '|' in cellpair:
                parts = cellpair.split('|')
            else:
                parts = [cellpair[:len(cellpair)//2], cellpair[len(cellpair)//2:]]
            
            if len(parts) >= 2:
                ligand_ct = parts[0].strip()
                receptor_ct = parts[1].strip()
                
                ligand_cts.append(ligand_ct)
                receptor_cts.append(receptor_ct)
                activities.append(activity)
        
        if not ligand_cts:
            return None
        
        df_split = pd.DataFrame({
            'ligand_ct': ligand_cts,
            'receptor_ct': receptor_cts,
            'activity': activities
        })
        
        heatmap_data = df_split.pivot_table(
            index='ligand_ct', 
            columns='receptor_ct', 
            values='activity', 
            fill_value=0,
            aggfunc='sum'
        )
        
        return heatmap_data

    # Get motifs and setup layout
    motifs = filtered_pivot.columns
    n_motifs = len(motifs)
    n_cols = 5
    n_rows = (n_motifs + n_cols - 1) // n_cols

    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    # Global scale for consistent colors
    all_positive_values = filtered_pivot[filtered_pivot > 0].values.flatten()
    all_positive_values = all_positive_values[~np.isnan(all_positive_values)]
    global_vmin = all_positive_values.min() if len(all_positive_values) > 0 else 1
    global_vmax = filtered_pivot.max().max()

    # Plot each motif
    for idx, motif in enumerate(motifs):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        
        motif_data = filtered_pivot[motif]
        heatmap_data = create_cellpair_heatmap(motif_data, motif)
        
        if heatmap_data is not None and not heatmap_data.empty:
            sns.heatmap(
                heatmap_data,
                cmap='Blues',
                cbar_kws={'label': 'Activity'},
                linewidths=0.3,
                square=True,
                xticklabels=True,
                yticklabels=True,
                ax=ax
            )
        else:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center', transform=ax.transAxes)
        
        ax.set_title(f'{motif}', fontsize=10, fontweight='bold')
        ax.set_xlabel('Receptor Cell Type', fontsize=8)
        ax.set_ylabel('Ligand Cell Type', fontsize=8)
        ax.tick_params(axis='x', rotation=45, labelsize=5)
        ax.tick_params(axis='y', rotation=0, labelsize=6)

    # Hide extra subplots
    for idx in range(n_motifs, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)

    plt.suptitle('Cell-Cell Communication Activity by Motif\n(Excluding Granulocyte pairs)', 
                 fontsize=16, fontweight='bold', y=1)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return fig


def plot_top_lri_interactions(lri_motifs, unique_ct, suffix="", save_path=None):
    """Plot top LRI interactions per motif with three signaling types:
    - Autocrine: diagonal stripes (///)
    - Paracrine: no pattern (solid)
    - Juxtacrine: opposite diagonal stripes (\\\\)
    """
    # Setup colors
    ct_color_map = get_cell_type_colors(unique_ct)
    
    motifs = sorted(lri_motifs['motif_idx'].unique())
    n_prog = len(motifs)
    top_n = 25

    # Layout
    cols = 3
    rows = int(np.ceil(n_prog / cols))
    fig, axes = plt.subplots(rows, cols,
                             figsize=(cols * 10, rows * 6),
                             constrained_layout=False)
    axes = axes.flatten()

    fig.subplots_adjust(right=0.85, left=0.15, top=0.94, bottom=0.06, wspace=1.35, hspace=0.4)

    for ax, prog in zip(axes, motifs):
        dfp = lri_motifs[lri_motifs['motif_idx'] == prog].copy()
        if dfp.empty:
            ax.set_visible(False)
            continue

        # Parse interaction names
        dfp[['celltype1','celltype2','ligand','receptor','signaling_type']] = (
            dfp['lri_name'].apply(lambda x: pd.Series(parse_lri_full(x)))
        )

        # Filter to unique receptors per signaling type
        dfp_sorted = dfp.sort_values('factor', ascending=False)
        dfp_filtered = dfp_sorted.drop_duplicates(
            subset=['signaling_type','receptor'],
            keep='first'
        )

        # Get top N
        top_df = dfp_filtered.nlargest(top_n, 'factor').reset_index(drop=True)
        y = np.arange(len(top_df))

        # Draw dual-segment bars with patterns for each signaling type
        for yi, row in top_df.iterrows():
            total = row['factor']
            half = total / 2.0
            
            # Get signaling type (normalize names)
            sig_type = row['signaling_type'].lower()
            if sig_type in ['auto', 'autocrine']:
                sig_type = 'autocrine'
                hatch_pattern = '///'
            elif sig_type in ['juxta', 'juxtacrine']:
                sig_type = 'juxtacrine'
                hatch_pattern = '\\\\\\\\'
            else:  # paracrine or default
                sig_type = 'paracrine'
                hatch_pattern = None

            # Colors for sender and receiver
            c1 = row['celltype1']
            col1 = ct_color_map.get(c1, 'gray')
            c2 = row['celltype2']
            col2 = ct_color_map.get(c2, 'gray')

            # Draw bars
            ax.barh(yi, half, color=col1, height=0.8, edgecolor='none')
            ax.barh(yi, half, left=half, color=col2, height=0.8, edgecolor='none')
            
            # Add patterns based on signaling type
            if hatch_pattern:
                ax.barh(yi, half, color='none', height=0.8, 
                       edgecolor='black', linewidth=1, hatch=hatch_pattern)
                ax.barh(yi, half, left=half, color='none', height=0.8,
                       edgecolor='black', linewidth=1, hatch=hatch_pattern)

        # Create labels
        labels = []
        for _, row in top_df.iterrows():
            sender = row['celltype1'].ljust(12)[:12]
            receiver = row['celltype2'].ljust(12)[:12]
            ligand = row['ligand'].ljust(10)[:10]
            receptor = row['receptor'].ljust(10)[:10]
            label = f"{sender} → {receiver} | {ligand} → {receptor}"
            labels.append(label)
        
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9, fontfamily='monospace')
        
        ax.invert_yaxis()
        ax.set_xlabel('Factor', fontsize=10)
        ax.set_title(f'Motif {prog}', fontsize=12, fontweight='bold')

    # Hide unused subplots
    for ax in axes[n_prog:]:
        ax.set_visible(False)

    # Create legend with all three signaling types
    legend_handles = []
    for ct, col in ct_color_map.items():
        legend_handles.append(Patch(facecolor=col, label=ct, edgecolor='black', linewidth=1))

    legend_handles.append(Patch(facecolor='none', edgecolor='none', label=''))
    
    # Add pattern legends for all three signaling types
    autocrine_patch = Patch(facecolor='lightgray', edgecolor='black', 
                           linewidth=1, hatch='///', label='Autocrine')
    paracrine_patch = Patch(facecolor='lightgray', edgecolor='black', 
                           linewidth=1, hatch='', label='Paracrine')
    juxtacrine_patch = Patch(facecolor='lightgray', edgecolor='black', 
                            linewidth=1, hatch='\\\\\\\\', label='Juxtacrine')
    
    legend_handles.extend([autocrine_patch, paracrine_patch, juxtacrine_patch])

    fig.legend(handles=legend_handles,
               title='Cell Types & Signaling',
               loc='center right',
               bbox_to_anchor=(0.92, 0.5),
               fontsize=8,
               title_fontsize=10,
               frameon=True,
               fancybox=True,
               shadow=False)

    fig.suptitle(
        'Top LR Interactions per Motif\n' + 
        'Format: Sender → Receiver | Ligand → Receptor\n' +
        '(Patterns indicate signaling type)',
        fontsize=14, y=1, fontweight='bold'
    )

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return fig

def plot_top_lri_interactions_dot(lri_motifs, unique_ct, suffix="", save_path=None):
    """Plot top LRI interactions per motif with line bars and colored dots at endpoints.
    Signaling types are indicated by letters: A (Autocrine), P (Paracrine), J (Juxtacrine)
    """
    # Setup colors
    ct_color_map = get_cell_type_colors(unique_ct)
    
    motifs = sorted(lri_motifs['motif_idx'].unique())
    n_prog = len(motifs)
    top_n = 25

    # Layout
    cols = 3
    rows = int(np.ceil(n_prog / cols))
    fig, axes = plt.subplots(rows, cols,
                             figsize=(cols * 10, rows * 6),
                             constrained_layout=False)
    axes = axes.flatten()

    fig.subplots_adjust(right=0.85, left=0.15, top=0.94, bottom=0.06, wspace=1.35, hspace=0.4)

    for ax, prog in zip(axes, motifs):
        dfp = lri_motifs[lri_motifs['motif_idx'] == prog].copy()
        if dfp.empty:
            ax.set_visible(False)
            continue

        # Parse interaction names
        dfp[['celltype1','celltype2','ligand','receptor','signaling_type']] = (
            dfp['lri_name'].apply(lambda x: pd.Series(parse_lri_full(x)))
        )

        # Filter to unique receptors per signaling type
        dfp_sorted = dfp.sort_values('factor', ascending=False)
        dfp_filtered = dfp_sorted.drop_duplicates(
            subset=['signaling_type','receptor'],
            keep='first'
        )

        # Get top N
        top_df = dfp_filtered.nlargest(top_n, 'factor').reset_index(drop=True)
        y = np.arange(len(top_df))

        # Draw line bars with dots at endpoints
        for yi, row in top_df.iterrows():
            total = row['factor']
            
            # Get signaling type and determine marker shape
            sig_type = row['signaling_type'].lower()
            if sig_type in ['auto', 'autocrine']:
                marker_shape = 's'  # square for autocrine
            elif sig_type in ['juxta', 'juxtacrine']:
                marker_shape = '^'  # triangle for juxtacrine
            else:  # paracrine or default
                marker_shape = 'o'  # circle for paracrine

            # Colors for sender and receiver
            c1 = row['celltype1']
            col1 = ct_color_map.get(c1, 'gray')
            c2 = row['celltype2']
            col2 = ct_color_map.get(c2, 'gray')

            # Draw thin horizontal line (the bar)
            ax.plot([0, total], [yi, yi], color='gray', linewidth=1.5, zorder=1)
            
            # Draw shaped markers at start and end
            marker_size = 80  # Adjust size as needed
            ax.scatter(0, yi, color=col1, s=marker_size, marker=marker_shape, 
                      zorder=2, edgecolors='black', linewidth=0.5)
            ax.scatter(total, yi, color=col2, s=marker_size, marker=marker_shape,
                      zorder=2, edgecolors='black', linewidth=0.5)

        # Create labels
        labels = []
        for _, row in top_df.iterrows():
            sender = row['celltype1'].ljust(12)[:12]
            receiver = row['celltype2'].ljust(12)[:12]
            ligand = row['ligand'].ljust(10)[:10]
            receptor = row['receptor'].ljust(10)[:10]
            label = f"{sender} → {receiver} | {ligand} → {receptor}"
            labels.append(label)
        
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9, fontfamily='monospace')
        
        ax.invert_yaxis()
        ax.set_xlabel('Factor', fontsize=10)
        ax.set_title(f'Motif {prog}', fontsize=12, fontweight='bold')
        
        # Set x-axis limits with some padding
        if len(top_df) > 0:
            max_val = top_df['factor'].max()
            ax.set_xlim(-max_val*0.02, max_val*1.05)

    # Hide unused subplots
    for ax in axes[n_prog:]:
        ax.set_visible(False)

    # Create legend
    legend_handles = []
    
    # Cell type legend (colored dots)
    for ct, col in ct_color_map.items():
        legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor=col, markersize=8, 
                                        markeredgecolor='black', markeredgewidth=0.5,
                                        label=ct))

    # Add separator
    legend_handles.append(plt.Line2D([0], [0], color='none', label=''))
    
    # Add signaling type legend with marker shapes
    legend_handles.append(plt.Line2D([0], [0], color='none', label='Signaling Types:'))
    legend_handles.append(plt.Line2D([0], [0], marker='s', color='w', 
                                    markerfacecolor='lightgray', markersize=8,
                                    markeredgecolor='black', markeredgewidth=0.5,
                                    label='■ = Autocrine'))
    legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', 
                                    markerfacecolor='lightgray', markersize=8,
                                    markeredgecolor='black', markeredgewidth=0.5,
                                    label='● = Paracrine'))
    legend_handles.append(plt.Line2D([0], [0], marker='^', color='w', 
                                    markerfacecolor='lightgray', markersize=8,
                                    markeredgecolor='black', markeredgewidth=0.5,
                                    label='▲ = Juxtacrine'))

    fig.legend(handles=legend_handles,
               title='Cell Types & Signaling',
               loc='center right',
               bbox_to_anchor=(0.92, 0.5),
               fontsize=8,
               title_fontsize=10,
               frameon=True,
               fancybox=True,
               shadow=False)

    fig.suptitle(
        'Top LR Interactions per Motif\n' + 
        'Format: Sender → Receiver | Ligand → Receptor\n' +
        '(Start marker = Sender, End marker = Receiver, Shape = Signaling type)',
        fontsize=14, y=1, fontweight='bold'
    )

    if save_path:
        # plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return fig


def plot_lri_networks(lri_motifs, unique_ct, suffix="", threshold=500, top_n=200, 
                     annotate_edges=False, save_path=None):
    """Plot LRI networks for motifs using Graphviz with three signaling types:
    - Autocrine: red edges
    - Paracrine: black edges
    - Juxtacrine: blue edges
    """
    # Setup colors
    ct_color_map = {ct: mcolors.to_hex(plt.get_cmap('tab20', len(unique_ct))(i)) 
                   for i, ct in enumerate(unique_ct)}
    
    motifs = sorted(lri_motifs['motif_idx'].unique())
    n_cols = 5
    n_rows = (len(motifs) + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols*3, n_rows*3),
                             constrained_layout=False)
    axes = axes.flatten()
    fig.subplots_adjust(right=0.90, bottom=0.05)

    # Plot each motif
    for i, motif in enumerate(motifs):
        ax = axes[i]
        df = lri_motifs[lri_motifs['motif_idx'] == motif].copy()
        if df.empty:
            ax.axis('off')
            continue

        # Parse and aggregate
        parsed = df['lri_name'].apply(parse_lri_full).tolist()
        df[['cell1','cell2','ligand','receptor','mode']] = pd.DataFrame(parsed, index=df.index)
        
        # Normalize signaling type names
        def normalize_signaling(mode):
            mode_lower = str(mode).lower()
            if mode_lower in ['auto', 'autocrine']:
                return 'autocrine'
            elif mode_lower in ['juxta', 'juxtacrine']:
                return 'juxtacrine'
            else:
                return 'paracrine'
        
        df['mode'] = df['mode'].apply(normalize_signaling)
        
        df2 = df.nlargest(top_n, 'factor')
        
        if annotate_edges:
            # Get top 20 LR pairs for edge annotation
            df_sorted = df.sort_values('factor', ascending=False)
            df_filtered = df_sorted.drop_duplicates(subset=['mode', 'receptor'], keep='first')
            top20_df = df_filtered.nlargest(20, 'factor')[['cell1', 'cell2', 'ligand', 'receptor', 'mode']]
            top20_set = set([(row['cell1'], row['cell2'], row['ligand'], row['receptor'], row['mode']) 
                            for _, row in top20_df.iterrows()])
        
        # Aggregate by cell pairs
        agg = df2.groupby(['cell1','cell2','mode'])['factor'].sum().reset_index(name='weight')
        agg = agg[agg['weight'] > threshold]
        
        if agg.empty:
            ax.axis('off')
            continue

        # Edge annotations if requested
        edge_annotations = {}
        if annotate_edges:
            for _, row in agg.iterrows():
                key = (row['cell1'], row['cell2'], row['mode'])
                mask = (df2['cell1'] == row['cell1']) & \
                       (df2['cell2'] == row['cell2']) & \
                       (df2['mode'] == row['mode'])
                lr_pairs = df2[mask].copy()
                
                if not lr_pairs.empty:
                    lr_pairs['in_top20'] = lr_pairs.apply(
                        lambda x: (x['cell1'], x['cell2'], x['ligand'], x['receptor'], x['mode']) in top20_set, 
                        axis=1
                    )
                    top20_pairs = lr_pairs[lr_pairs['in_top20']].nlargest(3, 'factor')
                    
                    lr_names = []
                    if len(top20_pairs) > 0:
                        for _, lr_row in top20_pairs.iterrows():
                            lr_names.append(f"{lr_row['ligand']}-{lr_row['receptor']}")
                    elif len(lr_pairs) > 0:
                        lr_row = lr_pairs.iloc[0]
                        lr_names.append(f"{lr_row['ligand']}-{lr_row['receptor']}")
                    
                    if lr_names:
                        edge_annotations[key] = '\n'.join(lr_names)

        # Build Graphviz graph
        dot = Digraph(engine='neato', format='png')
        dot.graph_attr['dpi'] = '300'
        dot.attr(splines='curved')
        dot.attr('node',
                 shape='circle',
                 style='filled',
                 fontsize='8',
                 width='0.5', height='0.5',
                 fixedsize='true')

        # Add nodes
        cells = sorted(set(agg['cell1']) | set(agg['cell2']))
        for c in cells:
            fill = ct_color_map.get(c, '#CCCCCC')
            dot.node(c, fillcolor=fill)

        # Add edges with colors based on signaling type
        max_w = agg['weight'].max()
        for _, row in agg.iterrows():
            pen = str(max(0.4, (2.5 * row['weight'] / max_w)))
            
            # Set edge color based on signaling type
            if row['mode'] == 'autocrine':
                ec = 'red'
            elif row['mode'] == 'juxtacrine':
                ec = 'blue'
            else:  # paracrine
                ec = 'black'
            
            if annotate_edges:
                key = (row['cell1'], row['cell2'], row['mode'])
                label = edge_annotations.get(key, '')
                dot.edge(row['cell1'], row['cell2'], 
                        penwidth=pen, color=ec, label=label,
                        fontsize='6', fontcolor='darkblue')
            else:
                dot.edge(row['cell1'], row['cell2'], penwidth=pen, color=ec)

        # Render and plot
        try:
            png = dot.pipe()
            img = Image.open(io.BytesIO(png))
            ax.imshow(img)
            ax.set_title(f'Motif {motif}', fontsize=8)
            ax.axis('off')
        except Exception as e:
            ax.text(0.5, 0.5, f'Error rendering\nMotif {motif}', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')

    # Hide unused subplots
    for ax in axes[len(motifs):]:
        ax.axis('off')

    # Legends
    node_handles = [Patch(facecolor=col, label=ct) for ct, col in ct_color_map.items()]
    
    fig.legend(handles=node_handles,
               title='Cell Types',
               loc='center right',
               bbox_to_anchor=(0.98, 0.5),
               fontsize=6,
               title_fontsize=8,
               ncol=1)

    # Add signaling type legend with all three types
    par_line = mlines.Line2D([], [], color='black', linewidth=2, label='Paracrine')
    auto_line = mlines.Line2D([], [], color='red', linewidth=2, label='Autocrine')
    juxta_line = mlines.Line2D([], [], color='blue', linewidth=2, label='Juxtacrine')
    
    fig.legend(handles=[auto_line, par_line, juxta_line],
               loc='lower center',
               ncol=3, frameon=False,
               fontsize=6)

    edge_label = " (with LR annotations)" if annotate_edges else ""
    fig.suptitle(f'LRI Networks Across Motifs (weight>{threshold}){edge_label}',
                 fontsize=12, y=0.98)

    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        # plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return fig

def plot_all_punches_by_cell_type(adata, cell_type_column='cell_type', 
                                  cell_types_to_show=None, n_cols=8,
                                  spot_size=1, figsize_per_subplot=(6, 6),
                                  title='Cell Type Distribution Across TMAs',
                                  color_palette=None,
                                  save_path=None):
    """Plot spatial distribution across all TMA punches"""
    print("🎨 PLOTTING ALL TMA PUNCHES BY CELL TYPE")
    print("=" * 50)
    
    # Get unique TMA IDs
    tma_ids = sorted(adata.obs['tma_id'].unique())
    n_punches = len(tma_ids)
    print(f"   Found {n_punches} TMA punches")
    
    # Setup cell types
    if cell_types_to_show is None:
        if hasattr(adata.obs[cell_type_column], 'cat'):
            cell_types_to_show = list(adata.obs[cell_type_column].cat.categories)
        else:
            cell_types_to_show = sorted(adata.obs[cell_type_column].unique())
        print(f"   Showing all {len(cell_types_to_show)} cell types")
    else:
        print(f"   Showing {len(cell_types_to_show)} specified cell types")
    
    # Setup grid
    n_rows = int(np.ceil(n_punches / n_cols))
    fig_width = n_cols * figsize_per_subplot[0]
    fig_height = n_rows * figsize_per_subplot[1]
    
    print(f"   Grid layout: {n_rows} rows × {n_cols} columns")
    print(f"   Figure size: {fig_width:.1f} × {fig_height:.1f} inches")
    
    # Create colors
    if not color_palette:
        if len(cell_types_to_show) <= 10:
            colors = sns.color_palette("tab10", len(cell_types_to_show))
        elif len(cell_types_to_show) <= 20:
            colors = sns.color_palette("tab20", len(cell_types_to_show))
        else:
            colors1 = sns.color_palette("tab20", 20)
            colors2 = sns.color_palette("husl", len(cell_types_to_show) - 20)
            colors = list(colors1) + list(colors2)
    else:
        colors = color_palette
        if len(colors) < len(cell_types_to_show):
            raise ValueError("Provided color_palette has fewer colors than cell_types_to_show")
    
    color_dict = {ct: colors[i] for i, ct in enumerate(sorted(cell_types_to_show))}
    
    # Create figure
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    print("   Plotting individual punches...")
    
    # Plot each TMA
    for idx, tma_id in enumerate(tma_ids):
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        
        adata_tma = adata[adata.obs['tma_id'] == tma_id].copy()
        
        if adata_tma.n_obs == 0:
            ax.text(0.5, 0.5, f'TMA {tma_id}\n(No cells)', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
        else:
            # Color cells
            cell_colors = []
            for cell_type in adata_tma.obs[cell_type_column]:
                if cell_type in color_dict:
                    cell_colors.append(color_dict[cell_type])
                else:
                    cell_colors.append('lightgray')
            
            # Plot
            spatial_coords = adata_tma.obsm['spatial']
            ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1], 
                      c=cell_colors, s=spot_size, alpha=0.7)
        
        ax.set_title(f'TMA {tma_id}', fontsize=10, pad=5)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')
        
        if (idx + 1) % 10 == 0:
            print(f"     Completed {idx + 1}/{n_punches} punches")
    
    # Hide empty subplots
    for idx in range(n_punches, n_rows * n_cols):
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        ax.axis('off')
    
    # Create legend
    print("   Creating legend...")
    legend_elements = []
    for cell_type in sorted(cell_types_to_show):
        if cell_type in color_dict:
            legend_elements.append(
                plt.Line2D([0], [0], marker='o', color='w', 
                          markerfacecolor=color_dict[cell_type], 
                          markersize=8, label=cell_type)
            )
    
    if legend_elements:
        fig.legend(handles=legend_elements, 
                  loc='center right', 
                  bbox_to_anchor=(1.08, 0.5),
                  fontsize=9,
                  title='Cell Types',
                  title_fontsize=10)
    
    plt.suptitle(title, fontsize=14, y=0.995)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, format='pdf', dpi=300, bbox_inches='tight')
        print(f"✅ Plot saved: {save_path}")
    
    return fig