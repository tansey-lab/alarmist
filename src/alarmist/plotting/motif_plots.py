"""
BPTF Motif Visualization Functions

Functions for visualizing BPTF analysis results including:
- LRI clustering and communication heatmaps
- Top LRI interactions per motif
- LRI network graphs
- Pathway-aggregated interactions

Adapted from scripts/bptf_visualization_utils.py for notebook-friendly usage.
"""

import io
import warnings
from typing import List, Optional, Dict, Tuple, Set
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.colors import LogNorm

# Optional dependencies
try:
    from PIL import Image
    from graphviz import Digraph
    GRAPHVIZ_AVAILABLE = True
except ImportError:
    GRAPHVIZ_AVAILABLE = False


# ==============================================================================
# Utility Functions
# ==============================================================================

def parse_lri_full(lri_name: str) -> Tuple[str, str, str, str, str]:
    """Parse LRI name into components

    Parameters
    ----------
    lri_name : str
        LRI name in format: celltype1|celltype2|ligand|receptor|mode

    Returns
    -------
    tuple
        (celltype1, celltype2, ligand, receptor, signaling_type)
    """
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


def get_cell_type_colors(unique_ct: List[str]) -> Dict[str, tuple]:
    """Generate color map for cell types

    Parameters
    ----------
    unique_ct : list of str
        List of unique cell type names

    Returns
    -------
    dict
        Mapping from cell type to RGB color tuple
    """
    ct_cmap = plt.get_cmap('tab20', len(unique_ct))
    return {ct: ct_cmap(i) for i, ct in enumerate(unique_ct)}


# ==============================================================================
# Data Preprocessing Functions
# ==============================================================================

def add_lri_components(lri_motifs_df: pd.DataFrame,
                       lri_col: str = 'lri_name') -> pd.DataFrame:
    """Add parsed LRI components to DataFrame

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with LRI names
    lri_col : str, default 'lri_name'
        Column containing LRI names

    Returns
    -------
    pd.DataFrame
        DataFrame with added columns: celltype1, celltype2, ligand, receptor, signaling_type
    """
    df = lri_motifs_df.copy()

    # Filter out GENE entries
    df = df[~df[lri_col].str.startswith('GENE')]

    # Parse LRI names
    df[['celltype1', 'celltype2', 'ligand', 'receptor', 'signaling_type']] = (
        df[lri_col].apply(lambda x: pd.Series(parse_lri_full(x)))
    )

    return df


def annotate_pathways(lri_motifs_df: pd.DataFrame,
                      cellchatdb_df: pd.DataFrame,
                      ligand_col: str = 'ligand',
                      receptor_col: str = 'receptor') -> pd.DataFrame:
    """Annotate LRI interactions with pathway information from CellChatDB

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with ligand and receptor columns
    cellchatdb_df : pd.DataFrame
        CellChatDB DataFrame with 'ligand', 'receptor', 'pathway' columns
    ligand_col : str, default 'ligand'
        Column name for ligand
    receptor_col : str, default 'receptor'
        Column name for receptor

    Returns
    -------
    pd.DataFrame
        DataFrame with added 'pathway' column
    """
    df = lri_motifs_df.copy()

    # Build ligand-receptor to pathway mapping
    lr_to_pathway = {}
    for _, row in cellchatdb_df.iterrows():
        ligands = row['ligand'].split('_') if pd.notna(row['ligand']) else []
        receptors = row['receptor'].split('_') if pd.notna(row['receptor']) else []
        pathway = row['pathway']

        for lig in ligands:
            for rec in receptors:
                lr_to_pathway[(lig, rec)] = pathway
        lr_to_pathway[(row['ligand'], row['receptor'])] = pathway

    # Map pathways
    def get_pathway(row):
        pathway = lr_to_pathway.get((row[ligand_col], row[receptor_col]))
        if pathway:
            return pathway

        # Fuzzy matching
        for (lig, rec), path in lr_to_pathway.items():
            if (row[ligand_col] in lig or lig in row[ligand_col]) and \
               (row[receptor_col] in rec or rec in row[receptor_col]):
                return path
        return 'Unknown'

    df['pathway'] = df.apply(get_pathway, axis=1)

    return df


# normalize_by_column_mean() has been removed
# Normalization is now done automatically in save_bptf_results()
# The saved lri_motifs.csv already contains 'factor_norm' and 'mean' columns


# ==============================================================================
# Visualization Functions
# ==============================================================================

def plot_lri_clustermap(lri_factors: np.ndarray,
                       column_names: List[str],
                       save_path: Optional[str] = None,
                       figsize: Tuple[int, int] = (20, 40)) -> object:
    """Create clustermap of LRI cell type pairs across motifs

    Parameters
    ----------
    lri_factors : np.ndarray
        LRI factors matrix (n_motifs × n_lris)
    column_names : list of str
        LRI column names
    save_path : str, optional
        Path to save figure. If None, displays plot.
    figsize : tuple, default (20, 40)
        Figure size

    Returns
    -------
    seaborn.ClusterGrid
        Clustermap object
    """
    H_obs = lri_factors
    n_motifs = H_obs.shape[0]

    # Create DataFrame
    df_all = pd.DataFrame(
        H_obs.T,
        index=column_names,
        columns=range(n_motifs)
    ).reset_index().rename(columns={'index':'lri'})

    df_all = df_all[~df_all['lri'].str.startswith('GENE')]

    # Extract cell_pair
    df_all['cell_pair'] = df_all['lri'].str.split('|').str[:2].str.join('|')

    # Aggregate by cell_pair
    pivot = df_all.groupby('cell_pair')[list(range(n_motifs))].sum()
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
        figsize=figsize,
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


def plot_celltype_communication_by_motif(lri_factors: np.ndarray,
                                        column_names: List[str],
                                        save_path: Optional[str] = None,
                                        n_cols: int = 5,
                                        figsize_per_motif: Tuple[float, float] = (5, 4)) -> plt.Figure:
    """Create cell-cell communication heatmaps for each motif

    Parameters
    ----------
    lri_factors : np.ndarray
        LRI factors matrix (n_motifs × n_lris)
    column_names : list of str
        LRI column names
    save_path : str, optional
        Path to save figure. If None, displays plot.
    n_cols : int, default 5
        Number of columns in subplot grid
    figsize_per_motif : tuple, default (5, 4)
        Size per motif subplot

    Returns
    -------
    matplotlib.figure.Figure
        Figure object
    """
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

    # Setup layout
    motifs = filtered_pivot.columns
    n_motifs = len(motifs)
    n_rows = (n_motifs + n_cols - 1) // n_cols

    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols,
                            figsize=(figsize_per_motif[0]*n_cols, figsize_per_motif[1]*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

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


def plot_top_lri_interactions_dot(
    lri_motifs_df: pd.DataFrame,
    unique_ct: List[str],
    factor_col: str = "factor",
    top_n: int = 35,
    save_path: Optional[str] = None,
    n_cols: int = 3,
    figsize_per_motif: Tuple[float, float] = (12, 9),
) -> plt.Figure:
    """
    Plot top LRI interactions per motif with lollipop lines + colored dots at endpoints.

    Changes vs original:
      - No long y-axis labels
      - Ligand label placed further LEFT of the x=0 axis line
      - Receptor label placed further RIGHT of the right endpoint
      - Remove top/right subplot spines

    Signaling types by marker:
      ■ autocrine, ● paracrine, ▲ juxtacrine
    """
    required_cols = {"motif_idx", "celltype1", "celltype2", "ligand", "receptor", "signaling_type", factor_col}
    missing = required_cols - set(lri_motifs_df.columns)
    if missing:
        raise ValueError(f"Missing columns in lri_motifs_df: {sorted(missing)}")

    ct_color_map = get_cell_type_colors(unique_ct)

    motifs = sorted(lri_motifs_df["motif_idx"].unique())
    n_prog = len(motifs)

    n_rows = int(np.ceil(n_prog / n_cols))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * figsize_per_motif[0], n_rows * figsize_per_motif[1]),
        constrained_layout=False,
    )
    axes = np.array(axes).reshape(-1)

    # Give some room for text on both sides
    fig.subplots_adjust(right=0.85, left=0.10, top=0.94, bottom=0.06, wspace=0.9, hspace=0.4)

    # --- label spacing knobs (in "points") ---
    # More negative => ligand further left; more positive => receptor further right
    ligand_dx_pts = -20
    receptor_dx_pts = 12

    for ax, prog in zip(axes, motifs):
        dfp = lri_motifs_df[lri_motifs_df["motif_idx"] == prog].copy()
        if dfp.empty:
            ax.set_visible(False)
            continue

        # Filter: keep top per (signaling_type, receptor, direction)
        dfp_sorted = dfp.sort_values(factor_col, ascending=False)
        dfp_filtered = dfp_sorted.drop_duplicates(
            subset=["signaling_type", "receptor", "celltype1", "celltype2"],
            keep="first",
        )

        top_df = dfp_filtered.nlargest(top_n, factor_col).reset_index(drop=True)
        y = np.arange(len(top_df))

        # No y tick labels (avoid huge block)
        ax.set_yticks(y)
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)

        for yi, row in top_df.iterrows():
            total = float(row[factor_col])

            sig_type = str(row["signaling_type"]).lower()
            if sig_type in ["auto", "autocrine"]:
                marker_shape = "s"
            elif sig_type in ["juxta", "juxtacrine"]:
                marker_shape = "^"
            else:
                marker_shape = "o"

            col1 = ct_color_map.get(row["celltype1"], "gray")
            col2 = ct_color_map.get(row["celltype2"], "gray")

            # lollipop line
            ax.plot([0, total], [yi, yi], color="gray", linewidth=1.5, zorder=1)

            marker_size = 80
            ax.scatter(
                0, yi,
                color=col1, s=marker_size, marker=marker_shape,
                zorder=2, edgecolors="black", linewidth=0.5
            )
            ax.scatter(
                total, yi,
                color=col2, s=marker_size, marker=marker_shape,
                zorder=2, edgecolors="black", linewidth=0.5
            )

            # Ligand: further LEFT of the left endpoint (x=0), so it sits left of the axis line
            ax.annotate(
                str(row["ligand"]),
                xy=(0, yi),
                xytext=(ligand_dx_pts, 0),
                textcoords="offset points",
                ha="right",
                va="center",
                fontsize=9,
                clip_on=False,
            )

            # Receptor: further RIGHT of the right endpoint
            ax.annotate(
                str(row["receptor"]),
                xy=(total, yi),
                xytext=(receptor_dx_pts, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=9,
                clip_on=False,
            )

        ax.invert_yaxis()
        ax.set_xlabel(f"{factor_col}", fontsize=10)
        ax.set_title(f"Motif {prog}", fontsize=12, fontweight="bold")

        # More padding so text doesn't feel cramped
        if len(top_df) > 0:
            max_val = float(top_df[factor_col].max())
            ax.set_xlim(-max_val * 0.04, max_val * 1.12)

        # Remove top/right spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Hide unused subplots
    for ax in axes[n_prog:]:
        ax.set_visible(False)

    # Legend
    legend_handles = []
    for ct, col in ct_color_map.items():
        legend_handles.append(
            plt.Line2D(
                [0], [0], marker="o", color="w",
                markerfacecolor=col, markersize=8,
                markeredgecolor="black", markeredgewidth=0.5,
                label=ct
            )
        )

    legend_handles.append(plt.Line2D([0], [0], color="none", label=""))
    legend_handles.append(plt.Line2D([0], [0], color="none", label="Signaling Types:"))
    legend_handles.append(
        plt.Line2D([0], [0], marker="s", color="w",
                   markerfacecolor="lightgray", markersize=8,
                   markeredgecolor="black", markeredgewidth=0.5,
                   label="■ = Autocrine")
    )
    legend_handles.append(
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="lightgray", markersize=8,
                   markeredgecolor="black", markeredgewidth=0.5,
                   label="● = Paracrine")
    )
    legend_handles.append(
        plt.Line2D([0], [0], marker="^", color="w",
                   markerfacecolor="lightgray", markersize=8,
                   markeredgecolor="black", markeredgewidth=0.5,
                   label="▲ = Juxtacrine")
    )

    fig.legend(
        handles=legend_handles,
        title="Cell Types & Signaling",
        loc="center right",
        bbox_to_anchor=(0.92, 0.5),
        fontsize=8,
        title_fontsize=10,
        frameon=True,
        fancybox=True,
        shadow=False,
    )

    fig.suptitle(
        "Top LR Interactions per Motif\n"
        "(Left text = Ligand, Right text = Receptor; Shape = Signaling type)",
        fontsize=14,
        y=1,
        fontweight="bold",
    )

    if save_path:
        plt.savefig(save_path, bbox_inches="tight", dpi=300)
        print(f"Saved: {save_path}")

    return fig


def plot_top_lri_interactions_by_pathway(lri_motifs_df: pd.DataFrame,
                                         unique_ct: List[str],
                                         top_n: int = 35,
                                         save_path: Optional[str] = None,
                                         n_cols: int = 3,
                                         factor_col: str = "factor",
                                         figsize_per_motif: Tuple[float, float] = (12, 9)) -> plt.Figure:
    """Plot top LRI interactions per motif aggregated by pathway

    Same sender, receiver, and pathway combinations are merged (factors summed).

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with parsed components and pathway column
    unique_ct : list of str
        List of unique cell types
    top_n : int, default 35
        Number of top pathways to show per motif
    save_path : str, optional
        Path to save figure. If None, displays plot.
    n_cols : int, default 3
        Number of columns in subplot grid
    figsize_per_motif : tuple, default (12, 9)
        Size per motif subplot

    Returns
    -------
    matplotlib.figure.Figure
        Figure object
    """
    if 'pathway' not in lri_motifs_df.columns:
        raise ValueError("DataFrame must have 'pathway' column. Use annotate_pathways() first.")

    # Setup colors
    ct_color_map = get_cell_type_colors(unique_ct)

    motifs = sorted(lri_motifs_df['motif_idx'].unique())
    n_prog = len(motifs)

    # Layout
    n_rows = int(np.ceil(n_prog / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * figsize_per_motif[0], n_rows * figsize_per_motif[1]),
                             constrained_layout=False)
    axes = axes.flatten()

    fig.subplots_adjust(right=0.85, left=0.15, top=0.94, bottom=0.06, wspace=1.35, hspace=0.4)

    for ax, prog in zip(axes, motifs):
        dfp = lri_motifs_df[lri_motifs_df['motif_idx'] == prog].copy()
        if dfp.empty:
            ax.set_visible(False)
            continue

        # Aggregate by sender, receiver, and pathway
        pathway_agg = dfp.groupby(['celltype1', 'celltype2', 'pathway']).agg({
            factor_col: 'sum',
            'ligand': lambda x: ', '.join(x.unique()[:2]),
            'receptor': lambda x: ', '.join(x.unique()[:2])
        }).reset_index()

        # Sort and get top N
        pathway_agg = pathway_agg.sort_values(factor_col, ascending=False)
        top_df = pathway_agg.head(top_n).reset_index(drop=True)
        y = np.arange(len(top_df))

        # Draw line bars with dots
        for yi, row in top_df.iterrows():
            total = row[factor_col]

            c1 = row['celltype1']
            col1 = ct_color_map.get(c1, 'gray')
            c2 = row['celltype2']
            col2 = ct_color_map.get(c2, 'gray')

            ax.plot([0, total], [yi, yi], color='gray', linewidth=1.5, zorder=1)

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

        if len(top_df) > 0:
            max_val = top_df[factor_col].max()
            ax.set_xlim(-max_val*0.02, max_val*1.05)

    # Hide unused subplots
    for ax in axes[n_prog:]:
        ax.set_visible(False)

    # Create legend
    legend_handles = []
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
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"Saved: {save_path}")

    return fig


def build_master_edge_gate(lri_motifs_df: pd.DataFrame,
                           top_n: int = 200,
                            factor_col: str = "factor",
                           threshold: float = 1500) -> Dict[int, Set[Tuple[str, str]]]:
    """Build edge gates for network filtering

    For each motif, aggregate across all signaling types and keep edges
    (cell1, cell2) whose total weight > threshold.

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with parsed LRI components
    top_n : int, default 200
        Number of top interactions to consider
    threshold : float, default 1500
        Threshold for edge weight

    Returns
    -------
    dict
        Mapping from motif_idx to set of (cell1, cell2) tuples
    """
    gates = {}
    motifs = sorted(lri_motifs_df['motif_idx'].unique())

    for motif in motifs:
        df = lri_motifs_df[lri_motifs_df['motif_idx'] == motif].copy()
        if df.empty:
            gates[motif] = set()
            continue

        # Take top N to reduce clutter
        df2 = df.nlargest(top_n, factor_col)

        # Aggregate by cell pair
        agg_all = df2.groupby(['celltype1', 'celltype2'])[factor_col].sum().reset_index(name='weight')
        keep = agg_all.loc[agg_all['weight'] > threshold, ['celltype1', 'celltype2']]
        gates[motif] = set(map(tuple, keep.values))

    return gates


def plot_lri_networks(lri_motifs_df: pd.DataFrame,
                     unique_ct: List[str],
                     threshold: float = 2000,
                     top_n: int = 200,
                     factor_col: str = "factor",
                     annotate_edges: bool = False,
                     mode_filter: Optional[str] = None,
                     edge_gate: Optional[Dict[int, Set[Tuple[str, str]]]] = None,
                     save_path: Optional[str] = None,
                     n_cols: int = 5,
                     figsize_per_motif: Tuple[float, float] = (3, 3)) -> plt.Figure:
    """Plot LRI networks for each motif using Graphviz

    Parameters
    ----------
    lri_motifs_df : pd.DataFrame
        DataFrame with parsed LRI components
    unique_ct : list of str
        List of unique cell types
    threshold : float, default 2000
        Threshold for edge weight filtering
    top_n : int, default 200
        Number of top interactions to consider
    annotate_edges : bool, default False
        Whether to annotate edges with L-R pairs
    mode_filter : str, optional
        Filter by signaling mode: None (all), 'paracrine', 'juxtacrine', 'autocrine'
    edge_gate : dict, optional
        Pre-computed edge gates from build_master_edge_gate()
    save_path : str, optional
        Path to save figure. If None, displays plot.
    n_cols : int, default 5
        Number of columns in subplot grid
    figsize_per_motif : tuple, default (3, 3)
        Size per motif subplot

    Returns
    -------
    matplotlib.figure.Figure
        Figure object

    Raises
    ------
    ImportError
        If Graphviz is not installed
    """
    if not GRAPHVIZ_AVAILABLE:
        raise ImportError("Graphviz required for network plots. Install: pip install graphviz pillow")

    import matplotlib.colors as mcolors

    # Setup colors
    ct_color_map = {ct: mcolors.to_hex(plt.get_cmap('tab20', len(unique_ct))(i))
                    for i, ct in enumerate(unique_ct)}

    def _norm_mode(x):
        t = str(x).strip().lower()
        if t in {"auto", "autocrine", "a"}:
            return "autocrine"
        if t in {"juxta", "juxtacrine", "j"}:
            return "juxtacrine"
        return "paracrine"

    motifs = sorted(lri_motifs_df['motif_idx'].unique())
    n_motifs = len(motifs)
    n_rows = (n_motifs + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols,
                            figsize=(n_cols*figsize_per_motif[0], n_rows*figsize_per_motif[1]),
                            constrained_layout=False)
    axes = axes.flatten()
    fig.subplots_adjust(right=0.90, bottom=0.05)

    for i, motif in enumerate(motifs):
        ax = axes[i]
        df = lri_motifs_df[lri_motifs_df['motif_idx'] == motif].copy()
        if df.empty:
            ax.axis('off')
            continue

        # Normalize mode
        df['mode'] = df['signaling_type'].apply(_norm_mode)

        # Filter by mode if specified
        if mode_filter is not None:
            mf = _norm_mode(mode_filter)
            df = df[df['mode'] == mf]
        if df.empty:
            ax.axis('off')
            continue

        df2 = df.nlargest(top_n, factor_col)

        # Aggregate by cell pair
        agg = df2.groupby(['celltype1', 'celltype2'])[factor_col].sum().reset_index(name='weight')

        # Apply edge gate or threshold
        if edge_gate is not None:
            keep = edge_gate.get(motif, set())
            if keep:
                mask = [(r['celltype1'], r['celltype2']) in keep for _, r in agg.iterrows()]
                agg = agg.loc[mask]
        else:
            agg = agg[agg['weight'] > threshold]

        if agg.empty:
            ax.axis('off')
            continue

        # Edge annotations
        edge_annotations = {}
        if annotate_edges:
            df_sorted = df.sort_values(factor_col, ascending=False)
            df_filtered = df_sorted.drop_duplicates(subset=['receptor'], keep='first')
            top20_df = df_filtered.nlargest(20, factor_col)[['celltype1', 'celltype2', 'ligand', 'receptor']]
            top20_set = set((r['celltype1'], r['celltype2'], r['ligand'], r['receptor'])
                           for _, r in top20_df.iterrows())

            for _, row in agg.iterrows():
                key = (row['celltype1'], row['celltype2'])
                lr_pairs = df2[(df2['celltype1']==row['celltype1']) &
                              (df2['celltype2']==row['celltype2'])].sort_values(factor_col, ascending=False)
                if not lr_pairs.empty:
                    lr_pairs['in_top20'] = lr_pairs.apply(
                        lambda x: (x['celltype1'], x['celltype2'], x['ligand'], x['receptor']) in top20_set,
                        axis=1
                    )
                    top_sel = lr_pairs[lr_pairs['in_top20']].head(3)
                    if top_sel.empty:
                        top_sel = lr_pairs.head(3)
                    edge_annotations[key] = '\n'.join(f"{r['ligand']}-{r['receptor']}" for _, r in top_sel.iterrows())

        # Build Graphviz graph
        dot = Digraph(engine='neato', format='png')
        dot.graph_attr['dpi'] = '300'
        dot.attr(splines='curved')
        dot.attr('node', shape='circle', style='filled', fontsize='8',
                 width='0.5', height='0.5', fixedsize='true')

        cells = sorted(set(agg['celltype1']) | set(agg['celltype2']))
        for c in cells:
            dot.node(c, fillcolor=ct_color_map.get(c, '#CCCCCC'))

        max_w = agg['weight'].max()
        for _, row in agg.iterrows():
            pen = str(max(0.4, (2.5 * row['weight'] / max_w)))
            if annotate_edges:
                label = edge_annotations.get((row['celltype1'], row['celltype2']), '')
                dot.edge(row['celltype1'], row['celltype2'], penwidth=pen, color='black',
                        label=label, fontsize='6', fontcolor='darkblue')
            else:
                dot.edge(row['celltype1'], row['celltype2'], penwidth=pen, color='black')

        # Render
        try:
            png = dot.pipe()
            img = Image.open(io.BytesIO(png))
            ax.imshow(img)
            ax.set_title(f'Motif {motif}', fontsize=8)
            ax.axis('off')
        except Exception as e:
            ax.text(0.5, 0.5, f'Error rendering\nMotif {motif}', ha='center', va='center',
                   transform=ax.transAxes)
            ax.axis('off')

    # Hide unused subplots
    for ax in axes[len(motifs):]:
        ax.axis('off')

    # Legend
    node_handles = [Patch(facecolor=col, label=ct) for ct, col in ct_color_map.items()]
    fig.legend(handles=node_handles, title='Cell Types',
               loc='center right', bbox_to_anchor=(0.98, 0.5),
               fontsize=6, title_fontsize=8, ncol=1)

    mf_text = "" if mode_filter is None else f" | only {_norm_mode(mode_filter)}"
    edge_text = " (with LR annotations)" if annotate_edges else ""
    gate_text = "" if edge_gate is None else " | gated by ALL>threshold"
    fig.suptitle(f'LRI Networks Across Motifs{mf_text}{gate_text}{edge_text}',
                 fontsize=12, y=0.98)

    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig
