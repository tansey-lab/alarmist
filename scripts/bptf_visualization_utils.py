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
    # plt.hist(patch_metadata_df['neighborhood_size'], bins=50, edgecolor='black')
    plt.xlabel('Number of cells per patch')
    plt.ylabel('Number of patches')
    plt.title('Distribution of cells per patch')
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

    df_all = df_all[~df_all['lri'].str.startswith('GENE')]

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


def plot_celltype_communication_by_motif(lri_factors, column_names, suffix="", save_path=None):
    """Create cell-cell communication heatmaps for each motif"""
    H_obs = lri_factors
    n_motifs = H_obs.shape[0]
    
    df_all = pd.DataFrame(
        H_obs.T,
        index=column_names,
        columns=range(n_motifs)
    ).reset_index().rename(columns={'index':'lri'})

    df_all = df_all[~df_all['lri'].str.startswith('GENE')]
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
        dfp = dfp[~dfp['lri_name'].str.startswith('GENE')]
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

def plot_top_lri_interactions_dot(lri_motifs, unique_ct, use_normalized=True, suffix="", save_path=None):
    """Plot top LRI interactions per motif with line bars and colored dots at endpoints.
    Signaling types are indicated by letters: A (Autocrine), P (Paracrine), J (Juxtacrine)
    
    Parameters:
    - use_normalized: If True, use factor_norm column; if False, use factor column
    """
    # Setup colors
    ct_color_map = get_cell_type_colors(unique_ct)
    
    # Choose which column to use
    factor_col = 'factor_norm' if use_normalized else 'factor'
    print(f"Using column: {factor_col}")
    
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

    for ax, prog in zip(axes, motifs):
        dfp = lri_motifs[lri_motifs['motif_idx'] == prog].copy()
        print(dfp.head())
        if dfp.empty:
            ax.set_visible(False)
            continue

        # Parse interaction names
        dfp = dfp[~dfp['lri_name'].str.startswith('GENE')]
        dfp[['celltype1','celltype2','ligand','receptor','signaling_type']] = (
            dfp['lri_name'].apply(lambda x: pd.Series(parse_lri_full(x)))
        )

        # Filter to unique receptors per signaling type
        dfp_sorted = dfp.sort_values(factor_col, ascending=False) # Changed to use factor_col 
        dfp_filtered = dfp_sorted.drop_duplicates(
            subset=['signaling_type','receptor','celltype1','celltype2'], 
            keep='first'
            )

        # Get top N
        top_df = dfp_filtered.nlargest(top_n, factor_col).reset_index(drop=True)  # Changed to use factor_col
        y = np.arange(len(top_df))

        # Draw line bars with dots at endpoints
        for yi, row in top_df.iterrows():
            total = row[factor_col]  # Changed to use factor_col
            
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
            ligand = row['ligand'].ljust(16)[:16]
            receptor = row['receptor'].ljust(16)[:16]
            label = f"{sender} → {receiver} | {ligand} → {receptor}"
            labels.append(label)
        
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9, fontfamily='monospace')
        
        ax.invert_yaxis()
        ax.set_xlabel('Normalized Factor' if use_normalized else 'Factor', fontsize=10)  # Changed label
        ax.set_title(f'Motif {prog}', fontsize=12, fontweight='bold')
        
        # Set x-axis limits with some padding
        if len(top_df) > 0:
            max_val = top_df[factor_col].max()  # Changed to use factor_col
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

    title_suffix = ' (Normalized by Mean)' if use_normalized else ''
    fig.suptitle(
        f'Top LR Interactions per Motif{title_suffix}\n' + 
        'Format: Sender → Receiver | Ligand → Receptor\n' +
        '(Start marker = Sender, End marker = Receiver, Shape = Signaling type)',
        fontsize=14, y=1, fontweight='bold'
    )

    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return fig

def build_master_edge_gate(lri_motifs, top_n=200, threshold=1500):
    """
    For each motif, aggregate across ALL signaling types,
    keep edges (cell1,cell2) whose total weight > threshold.
    Return: dict[motif_idx] -> set of (cell1, cell2)
    """
    gates = {}
    motifs = sorted(lri_motifs['motif_idx'].unique())
    for motif in motifs:
        df = lri_motifs[lri_motifs['motif_idx'] == motif].copy()
        if df.empty:
            gates[motif] = set(); continue
        df = df[~df['lri_name'].str.startswith('GENE')]
        parsed = df['lri_name'].apply(parse_lri_full).tolist()
        df[['cell1','cell2','ligand','receptor','mode']] = pd.DataFrame(parsed, index=df.index)
        df2 = df.nlargest(top_n, 'factor')  # reduce clutter before aggregation
        agg_all = df2.groupby(['cell1','cell2'])['factor'].sum().reset_index(name='weight')
        keep = agg_all.loc[agg_all['weight'] > threshold, ['cell1','cell2']]
        gates[motif] = set(map(tuple, keep.values))
    return gates


def plot_lri_networks(lri_motifs, unique_ct, suffix="", threshold=2000, top_n=200, 
                      annotate_edges=False, save_path=None, mode_filter=None, use_normalized=True,
                      edge_gate=None):
    """
    Plot LRI networks for each motif using Graphviz.

    mode_filter:
        - None: aggregate across ALL signaling types (default)
        - "paracrine": keep only paracrine edges
        - "juxtacrine": keep only juxtacrine edges
        - (supports "autocrine" too if ever needed)
    """
    # Setup node colors by cell type
    ct_color_map = {ct: mcolors.to_hex(plt.get_cmap('tab20', len(unique_ct))(i)) 
                    for i, ct in enumerate(unique_ct)}

    factor_col = 'factor_norm' if use_normalized else 'factor'
    print(f"Using column: {factor_col}")

    def _norm_mode(x):
        t = str(x).strip().lower()
        if t in {"auto","autocrine","a"}: return "autocrine"
        if t in {"juxta","juxtacrine","j"}: return "juxtacrine"
        return "paracrine"

    motifs = sorted(lri_motifs['motif_idx'].unique())
    n_cols = 5
    n_rows = (len(motifs) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*3, n_rows*3), constrained_layout=False)
    axes = axes.flatten()
    fig.subplots_adjust(right=0.90, bottom=0.05)

    for i, motif in enumerate(motifs):
        ax = axes[i]
        df = lri_motifs[lri_motifs['motif_idx'] == motif].copy()
        if df.empty: ax.axis('off'); continue

        df = df[~df['lri_name'].str.startswith('GENE')]
        parsed = df['lri_name'].apply(parse_lri_full).tolist()
        df[['cell1','cell2','ligand','receptor','mode_raw']] = pd.DataFrame(parsed, index=df.index)
        df['mode'] = df['mode_raw'].apply(_norm_mode)

        if mode_filter is not None:
            mf = _norm_mode(mode_filter)
            df = df[df['mode'] == mf]
        if df.empty: ax.axis('off'); continue

        df2 = df.nlargest(top_n, factor_col)

        # aggregate (within current filter)
        agg = df2.groupby(['cell1','cell2'])[factor_col].sum().reset_index(name='weight')

        # --- apply global "gate from ALL" if provided ---
        if edge_gate is not None:
            keep = edge_gate.get(motif, set())
            if keep:
                mask = [(r['cell1'], r['cell2']) in keep for _, r in agg.iterrows()]
                agg = agg.loc[mask]
        else:
            # fallback to per-figure thresholding if no gate provided
            agg = agg[agg['weight'] > threshold]

        if agg.empty: ax.axis('off'); continue

        # optional annotations (不变，略)
        edge_annotations = {}
        if annotate_edges:
            df_sorted = df.sort_values(factor_col, ascending=False)
            df_filtered = df_sorted.drop_duplicates(subset=['receptor'], keep='first')
            top20_df = df_filtered.nlargest(20, factor_col)[['cell1','cell2','ligand','receptor']]
            top20_set = set((r['cell1'], r['cell2'], r['ligand'], r['receptor']) 
                            for _, r in top20_df.iterrows())
            for _, row in agg.iterrows():
                key = (row['cell1'], row['cell2'])
                lr_pairs = df2[(df2['cell1']==row['cell1']) & (df2['cell2']==row['cell2'])] \
                            .sort_values(factor_col, ascending=False)
                if not lr_pairs.empty:
                    lr_pairs['in_top20'] = lr_pairs.apply(
                        lambda x: (x['cell1'], x['cell2'], x['ligand'], x['receptor']) in top20_set,
                        axis=1
                    )
                    top_sel = lr_pairs[lr_pairs['in_top20']].head(3)
                    if top_sel.empty: top_sel = lr_pairs.head(3)
                    edge_annotations[key] = '\n'.join(f"{r['ligand']}-{r['receptor']}" for _, r in top_sel.iterrows())

        # graphviz (同前)
        dot = Digraph(engine='neato', format='png')
        dot.graph_attr['dpi'] = '300'
        dot.attr(splines='curved')
        dot.attr('node', shape='circle', style='filled', fontsize='8',
                 width='0.5', height='0.5', fixedsize='true')

        cells = sorted(set(agg['cell1']) | set(agg['cell2']))
        for c in cells:
            dot.node(c, fillcolor=ct_color_map.get(c, '#CCCCCC'))

        max_w = agg['weight'].max()
        for _, row in agg.iterrows():
            pen = str(max(0.4, (2.5 * row['weight'] / max_w)))
            if annotate_edges:
                label = edge_annotations.get((row['cell1'], row['cell2']), '')
                dot.edge(row['cell1'], row['cell2'], penwidth=pen, color='black',
                         label=label, fontsize='6', fontcolor='darkblue')
            else:
                dot.edge(row['cell1'], row['cell2'], penwidth=pen, color='black')

        try:
            png = dot.pipe()
            img = Image.open(io.BytesIO(png))
            ax.imshow(img); ax.set_title(f'Motif {motif}', fontsize=8); ax.axis('off')
        except Exception:
            ax.text(0.5, 0.5, f'Error rendering\nMotif {motif}', ha='center', va='center',
                    transform=ax.transAxes); ax.axis('off')

    for ax in axes[len(motifs):]:
        ax.axis('off')

    node_handles = [Patch(facecolor=col, label=ct) for ct, col in ct_color_map.items()]
    fig.legend(handles=node_handles, title='Cell Types',
               loc='center right', bbox_to_anchor=(0.98, 0.5),
               fontsize=6, title_fontsize=8, ncol=1)

    mf_text = "" if mode_filter is None else f" | only { _norm_mode(mode_filter) }"
    edge_text = " (with LR annotations)" if annotate_edges else ""
    gate_text = "" if edge_gate is None else " | gated by ALL>threshold"
    fig.suptitle(f'LRI Networks Across Motifs{mf_text}{gate_text}{edge_text}',
                 fontsize=12, y=0.98)

    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight'); print(f"Saved: {save_path}")
    return fig