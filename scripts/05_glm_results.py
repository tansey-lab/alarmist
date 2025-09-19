#!/usr/bin/env python3
"""
04 - GLM Results Analysis and Visualization

Performs downstream analysis of GLM results including volcano plots, forest plots,
and marker gene filtering to remove bleeding effects between cell types.

Usage:
    python scripts/04_glm_results.py --data-file data.h5ad --results-dir results/bptf_de --output-dir results/plots
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import scipy.sparse as sp
from adjustText import adjust_text
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc


def differential_expression(X, in_mask, out_mask=None, min_in_group_fraction=0.0001, min_out_group_fraction=0.0001):
    """Perform differential expression analysis using Mann-Whitney U test"""
    if out_mask is None:
        out_mask = ~in_mask

    # Filter genes that are basically never expressed in the control
    X_out = X[out_mask]
    genes_mask = np.ones(X.shape[1], dtype=bool)
    if min_out_group_fraction > 0:
        control_genes_mask = (X_out > 0.5).mean(axis=0) >= min_out_group_fraction
        genes_mask = genes_mask & control_genes_mask

    # Filter genes that are basically never expressed in the target
    X_in = X[in_mask]
    if min_out_group_fraction > 0:
        target_genes_mask = (X_in > 0.5).mean() >= min_out_group_fraction
        genes_mask = genes_mask & target_genes_mask

    genes_mask = np.array(genes_mask).flatten()
    X_out = X_out[:, genes_mask]
    X_in = X_in[:, genes_mask]

    # Calculate the log-fold changes
    control_means = np.asarray(X_out.mean(axis=0)).ravel()
    target_means = np.asarray(X_in.mean(axis=0)).ravel()
    logfoldchanges = np.zeros(X.shape[1])
    logfoldchanges[genes_mask] = np.log2(target_means.clip(1e-300)) - np.log2(control_means.clip(1e-300))

    p_values = np.ones(X.shape[1])
    for j_idx, j in enumerate(np.where(genes_mask)[0]):
        x = X_in[:, j_idx].toarray().ravel()
        y = X_out[:, j_idx].toarray().ravel()
        p_values[j] = stats.mannwhitneyu(x, y).pvalue

    p_adj = np.ones(X.shape[1])
    p_adj[genes_mask] = stats.false_discovery_control(p_values[genes_mask])

    return {'p_values': p_values, 'p_adj': p_adj, 'logfoldchanges': logfoldchanges}


def compute_exclusion_mask(adata, marker_lfc=1, marker_pvalue=1e-5, marker_subsample=50000):
    """
    Identify marker genes for each cell type to exclude later.
    Subsamples large groups to marker_subsample size for speed.
    Returns cell_types, gene list, and exclusion mask matrix.
    """
    genes = np.array(adata.var_names)
    cell_types = np.unique(adata.obs['cell_type'])
    n_types = len(cell_types)
    exclusion_mask = np.zeros((n_types, len(genes)), dtype=bool)

    print(f"Computing marker genes for {len(cell_types)} cell types...")

    for cidx, cell_type in enumerate(cell_types):
        print(f"Processing {cell_type} ({cidx+1}/{len(cell_types)})...")
        
        # Create boolean masks for in-group and out-group cells
        in_mask = adata.obs['cell_type'] == cell_type
        out_mask = ~in_mask

        # Subsample in-group if too large
        if in_mask.sum() > marker_subsample:
            tmp = np.zeros(len(in_mask), dtype=bool)
            tmp[np.random.choice(
                np.where(in_mask)[0],
                size=marker_subsample, replace=False
            )] = True
            in_mask = tmp

        # Subsample out-group if too large
        if out_mask.sum() > marker_subsample:
            tmp = np.zeros(len(out_mask), dtype=bool)
            tmp[np.random.choice(
                np.where(out_mask)[0],
                size=marker_subsample, replace=False
            )] = True
            out_mask = tmp

        print(f"  {in_mask.sum()} in-group cells, {out_mask.sum()} out-group cells")
        
        # Perform differential expression
        deg = differential_expression(
            adata.X, in_mask=in_mask, out_mask=out_mask
        )

        # Identify markers by p_adj and logfoldchanges thresholds
        marker_mask = (
            (deg['p_adj'] <= marker_pvalue) &
            (deg['logfoldchanges'] >= marker_lfc)
        )
        exclusion_mask[cidx, marker_mask] = True
        print(f"  Found {marker_mask.sum()} marker genes")

    return cell_types, genes, exclusion_mask


# def volcano_plot(df, x_col, y_col, label_col=None, fdr=0.1,
#                  x_threshold=1, marker='o', n_top=10,
#                  figsize=(10, 10), fontsize=8, ax=None):
#     """Draw a volcano plot with customized coloring and labeling"""
#     if ax is None:
#         fig, ax = plt.subplots(figsize=figsize)
    
#     plt.rcParams["font.family"] = "Arial"
#     colors = ["#E64B35FF", "#3C5488FF", "#00A087FF", "#4DBBD5FF",
#               "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
#               "#7E6148FF", "#B09C85FF"]
#     sns.set_palette(sns.color_palette(colors))
#     sns.set(rc={'figure.figsize': figsize, "font.size": fontsize})
#     sns.set_style("white")
    
#     threshold_fdr = -np.log10(fdr)
#     x = df[x_col].values
#     y = df[y_col].values
    
#     sigmask = y >= threshold_fdr
#     effmask = np.abs(x) >= x_threshold
#     rightmask = x >= x_threshold
#     leftmask = x <= -x_threshold
    
#     xr = np.abs(x).max() if np.abs(x).max() > 0 else 1
#     yr = y.max() if y.max() > 0 else 1
#     upper_right = x/xr + y/yr
#     upper_left = -x/xr + y/yr
#     top_right = np.argsort(upper_right)[::-1]
#     top_left = np.argsort(upper_left)[::-1]

#     indices = np.arange(len(y))
#     top_right = top_right[sigmask[top_right] & rightmask[top_right]]
#     top_left = top_left[sigmask[top_left] & leftmask[top_left]]
#     bottom_right = indices[(~sigmask) & rightmask]
#     bottom_left = indices[(~sigmask) & leftmask]
#     middle = indices[~effmask]

#     ax.scatter(x[middle], y[middle], color='gray', alpha=0.2, marker=marker)
#     ax.scatter(x[bottom_right], y[bottom_right], color='black', alpha=0.7, marker=marker)
#     ax.scatter(x[bottom_left], y[bottom_left], color='black', alpha=0.7, marker=marker)
#     ax.scatter(x[top_left], y[top_left], color='blue', alpha=0.7, marker=marker)
#     ax.scatter(x[top_right], y[top_right], color='red', alpha=0.7, marker=marker)
#     ax.set_ylabel('-$\\log_{10}(q)$', fontsize=fontsize+4)

#     # eff_y = df.loc[effmask, y_col]
#     # if eff_y.notnull().any() and eff_y.max() > 0:
#     #     ax.set_ylim([0, eff_y.max()*1.1])
#     # else:
#     #     ax.set_ylim([0, threshold_fdr*1.2])
#     # ---- Y 轴：按所有点 ----
#     all_y = df[y_col].to_numpy()
#     valid_y = np.isfinite(all_y)
#     if valid_y.any():
#         y_max = np.nanmax(all_y[valid_y])
#         y_max = max(y_max, threshold_fdr)  # 至少能看到阈值线
#         y_max = max(y_max, 1e-8)           # 防 0
#         ax.set_ylim(0, 1.1 * y_max)
#     else:
#         ax.set_ylim(0, threshold_fdr * 1.2)

#     # ---- X 轴：用 |x|<10 的点决定范围，空则退回所有有效点 ----
#     valid_x = np.isfinite(x)
#     if valid_x.any():
#         in10 = np.abs(x) < 10
#         candidates = valid_x & in10
#         if candidates.any():
#             max_x = np.nanmax(np.abs(x[candidates]))
#         else:
#             max_x = np.nanmax(np.abs(x[valid_x]))  # 退回所有有效点
#         max_x = max(max_x, 1e-8)  # 防 0
#         ax.set_xlim(-1.1 * max_x, 1.1 * max_x)
#     else:
#         ax.set_xlim(-1.1 * x_threshold, 1.1 * x_threshold)

#     head_idx = np.concatenate([top_right[:n_top], top_left[:n_top]])
#     head_df = df.iloc[head_idx]
#     texts = []
#     if label_col:
#         for xi, yi, lbl in zip(head_df[x_col], head_df[y_col], head_df[label_col]):
#             texts.append(ax.text(xi, yi, lbl, ha='center', va='center', fontsize=fontsize))
#         adjust_text(texts, ax=ax,
#                     arrowprops=dict(arrowstyle='-', lw=0.5, color='gray'),
#                     expand_points=(3, 3), expand_text=(4, 4),
#                     force_text=(0.75, 0.75), force_points=(0.3, 0.3), lim=100)
    
#     ax.axvline(x_threshold, ls='--', color='black', alpha=0.5)
#     ax.axvline(-x_threshold, ls='--', color='black', alpha=0.5)
#     ax.axhline(threshold_fdr, ls='--', color='black', alpha=0.5)
#     return ax

def volcano_plot(df, x_col, y_col, label_col=None, fdr=0.1,
                 x_threshold=1, marker='o', n_top=30,
                 figsize=(10, 10), fontsize=8, ax=None):
    """Draw a volcano plot with customized coloring and labeling."""
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # ---- styles ----
    plt.rcParams["font.family"] = "Arial"
    colors = ["#E64B35FF", "#3C5488FF", "#00A087FF", "#4DBBD5FF",
              "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
              "#7E6148FF", "#B09C85FF"]
    sns.set_palette(sns.color_palette(colors))
    sns.set(rc={'figure.figsize': figsize, "font.size": fontsize})
    sns.set_style("white")

    # ---- data & masks ----
    threshold_fdr = -np.log10(fdr)
    x = df[x_col].to_numpy()
    y = df[y_col].to_numpy()
    valid = np.isfinite(x) & np.isfinite(y)

    # basic masks (always AND with valid)
    sigmask   = (y >= threshold_fdr) & valid
    effmask   = (np.abs(x) >= x_threshold) & valid
    rightmask = (x >= x_threshold) & valid
    leftmask  = (x <= -x_threshold) & valid

    # ---- scoring for "top" candidate selection (independent of limits) ----
    xr = np.nanmax(np.abs(x[valid])) if np.any(valid) else 1.0
    yr = np.nanmax(y[valid]) if np.any(valid) else 1.0
    xr = xr if xr > 0 else 1.0
    yr = yr if yr > 0 else 1.0

    upper_right = (x / xr) + (y / yr)
    upper_left  = (-x / xr) + (y / yr)
    top_right = np.argsort(upper_right)[::-1]
    top_left  = np.argsort(upper_left)[::-1]

    # keep only significant & on each side
    top_right = top_right[sigmask[top_right] & rightmask[top_right]]
    top_left  = top_left[sigmask[top_left] & leftmask[top_left]]

    # ---- set axis limits BEFORE labeling ----
    # Y: based on all finite points (ensure at least threshold visible)
    if np.any(valid):
        y_max = np.nanmax(y[valid])
        y_max = max(y_max, threshold_fdr, 1e-8)
        ax.set_ylim(0, 1.1 * y_max)
    else:
        ax.set_ylim(0, threshold_fdr * 1.2)

    # X: use points with |x| < 10 if available; else all finite
    if np.any(valid):
        in10 = (np.abs(x) < 10) & valid
        if np.any(in10):
            max_x = np.nanmax(np.abs(x[in10]))
        else:
            max_x = np.nanmax(np.abs(x[valid]))
        max_x = max(max_x, 1e-8)
        ax.set_xlim(-1.1 * max_x, 1.1 * max_x)
    else:
        ax.set_xlim(-1.1 * x_threshold, 1.1 * x_threshold)

    # ---- plot points (finite only) ----
    indices = np.arange(len(y))
    bottom_right = indices[(~sigmask) & rightmask]
    bottom_left  = indices[(~sigmask) & leftmask]
    middle       = indices[valid & (~effmask)]  # non-effective (|x|<thr), finite

    ax.scatter(x[middle],        y[middle],        color='gray',  alpha=0.2, marker=marker, clip_on=True)
    ax.scatter(x[bottom_right],  y[bottom_right],  color='black', alpha=0.7, marker=marker, clip_on=True)
    ax.scatter(x[bottom_left],   y[bottom_left],   color='black', alpha=0.7, marker=marker, clip_on=True)
    ax.scatter(x[top_left],      y[top_left],      color='blue',  alpha=0.7, marker=marker, clip_on=True)
    ax.scatter(x[top_right],     y[top_right],     color='red',   alpha=0.7, marker=marker, clip_on=True)
    ax.set_ylabel('-$\\log_{10}(q)$', fontsize=fontsize+4)

    # ---- label only points that are within the current axes window ----
    if label_col:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        in_view = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax) & valid

        # filter the top candidates to those in view
        top_right_in = top_right[in_view[top_right]]
        top_left_in  = top_left[in_view[top_left]]

        # take up to n_top per side that are in view
        head_idx = np.concatenate([top_right_in[:n_top], top_left_in[:n_top]])
        if head_idx.size > 0:
            head_df = df.iloc[head_idx]
            texts = []
            for xi, yi, lbl in zip(head_df[x_col], head_df[y_col], head_df[label_col]):
                texts.append(ax.text(xi, yi, lbl, ha='center', va='center',
                                     fontsize=fontsize, clip_on=True))
            adjust_text(
                texts, ax=ax,
                arrowprops=dict(arrowstyle='-', lw=0.5, color='gray'),
                expand_points=(3, 3), expand_text=(4, 4),
                force_text=(0.75, 0.75), force_points=(0.3, 0.3),
                lim=100
            )

    # ---- thresholds ----
    ax.axvline(x_threshold,  ls='--', color='black', alpha=0.5)
    ax.axvline(-x_threshold, ls='--', color='black', alpha=0.5)
    ax.axhline(threshold_fdr, ls='--', color='black', alpha=0.5)

    return ax


def forest_plot(df, ax=None, n_top=20,
                effect_col='logFC', q_col='qval', gene_col='gene'):
    """Draw a horizontal forest plot of the strongest absolute logFC signals"""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 8))

    # rank and pick top
    df['_score'] = np.abs(df[effect_col]) * -np.log10(df[q_col].clip(1e-300))
    df = df[np.abs(df[effect_col]) <= 10].copy()
    top = (df.sort_values('_score', ascending=False)
             .head(n_top)
             .sort_values(effect_col))  # negatives on top
    
    # top = top[np.abs(top[effect_col]) <= 10].copy()

    y = np.arange(len(top))
    colors = top[effect_col].map(lambda x: 'red' if x > 0 else 'blue')

    # 95% Wald CI: effect ± 1.96 * SE
    se_x = top['se']
    lo = top[effect_col] - 1.96 * se_x
    hi = top[effect_col] + 1.96 * se_x
    xerr = np.vstack([top[effect_col] - lo, hi - top[effect_col]])

    # thin CI bars (no SE available)
    # ax.hlines(y, xmin=0, xmax=top[effect_col], color='lightgray', lw=2, zorder=1)
    ax.errorbar(x=top[effect_col], y=y, xerr=xerr, fmt='none', elinewidth=1, capsize=3, zorder=2)
    ax.scatter(top[effect_col], y, c=colors, s=25, zorder=2)

    # gene labels & q-values on the right
    xmax = ax.get_xlim()[1]
    for yi, gene, q in zip(y, top[gene_col], top[q_col]):
        ax.text(xmax * 1.02, yi, f'{gene}  (q={q:.1e})',
                va='center', fontsize=7)

    ax.axvline(0, color='black', lw=1)
    ax.set_yticks([])
    ax.set_xlabel('log$_2$ Fold-Change')
    sns.despine(ax=ax, left=True)
    ax.set_title(f'Top {n_top} DE genes', fontsize=9)
    plt.tight_layout()

# def forest_plot(
#     df,
#     ax=None,
#     n_top=30,
#     effect_col='logFC',      # effect column in your CSV
#     q_col='qval',
#     gene_col='gene',
#     se_col='se',             # SE column you added
#     effect_is_log2=False     # False if effect_col is natural-log; True if already log2
# ):
#     """
#     Draw a horizontal forest plot with Wald 95% CIs computed on the fly.

#     If effect_is_log2 is False, converts effects/SE from natural-log to log2.
#     Requires a 'se' column giving the SE on the same scale as effect_col.
#     """
#     import numpy as np
#     import matplotlib.pyplot as plt
#     import seaborn as sns

#     if se_col not in df.columns:
#         raise ValueError(f"Missing '{se_col}' column with standard errors.")

#     if ax is None:
#         _, ax = plt.subplots(figsize=(6, 8))

#     # --- prepare effects on log2 scale ---
#     LOG2 = np.log(2.0)
#     eff = df[effect_col].astype(float).to_numpy()
#     se  = df[se_col].astype(float).to_numpy()

#     if not effect_is_log2:
#         eff_log2 = eff / LOG2
#         se_log2  = se  / LOG2
#     else:
#         eff_log2 = eff
#         se_log2  = se

#     # score for ranking (robust to p=0 by clipping)
#     score = np.abs(eff_log2) * -np.log10(df[q_col].clip(1e-300).astype(float))
#     tmp = df.copy()
#     tmp['_eff_log2'] = eff_log2
#     tmp['_se_log2']  = se_log2
#     tmp['_score']    = score

#     top = (tmp.sort_values('_score', ascending=False)
#               .head(n_top)
#               .sort_values('_eff_log2'))  # negatives on top

#     print(top)
#     y = np.arange(len(top))
#     x = top['_eff_log2'].to_numpy()
#     se_x = top['_se_log2'].to_numpy()

#     # 95% Wald CI: effect ± 1.96 * SE
#     lo = x - 1.96 * se_x
#     hi = x + 1.96 * se_x
#     xerr = np.vstack([x - lo, hi - x])

#     # points + error bars
#     colors = ['red' if xi > 0 else 'blue' for xi in x]
#     ax.errorbar(x=x, y=y, xerr=xerr, fmt='o', elinewidth=1, capsize=3, zorder=2)
#     ax.scatter(x, y, c=colors, s=50, zorder=3)

#     # right-side labels
#     xmax = ax.get_xlim()[1]
#     for yi, gene, q in zip(y, top[gene_col], top[q_col]):
#         ax.text(xmax * 1.02, yi, f'{gene}  (q={float(q):.1e})',
#                 va='center', fontsize=7)

#     ax.axvline(0, color='black', lw=1)
#     ax.set_yticks([])
#     ax.set_xlabel('log$_2$ Fold-Change')
#     sns.despine(ax=ax, left=True)
#     ax.set_title(f'Top {n_top} DE genes', fontsize=9)
#     plt.tight_layout()
#     # return ax



def main():
    """Main GLM results analysis pipeline"""
    parser = argparse.ArgumentParser(description='GLM results analysis and visualization')
    parser.add_argument('--data-file', required=True,
                       help='Input h5ad file with expression data')
    parser.add_argument('--results-dir', required=True,
                       help='Directory containing GLM results CSVs')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for plots')
    parser.add_argument('--n-motifs', type=int, default=20,
                       help='Number of motifs to analyze')
    parser.add_argument('--min-expression-frac', type=float, default=0.02,
                       help='Minimum expression fraction in cell type')
    parser.add_argument('--marker-lfc', type=float, default=1.0,
                       help='Log fold change threshold for marker genes')
    parser.add_argument('--marker-pvalue', type=float, default=1e-5,
                       help='P-value threshold for marker genes')
    parser.add_argument('--marker-subsample', type=int, default=50000,
                       help='Maximum cells to use for marker gene detection')
    parser.add_argument('--fdr-threshold', type=float, default=0.05,
                       help='FDR threshold for volcano plots')
    parser.add_argument('--lfc-threshold', type=float, default=0.5,
                       help='Log fold change threshold for volcano plots')
    parser.add_argument('--n-top-genes', type=int, default=10,
                       help='Number of top genes to label in volcano plots')
    parser.add_argument('--random-state', type=int, default=42,
                       help='Random seed for reproducibility')
    parser.add_argument('--force-recompute-markers', action='store_true',
                       help='Force recompute marker genes even if cached files exist')
    
    args = parser.parse_args()
    
    print("="*60)
    print("04 - GLM RESULTS ANALYSIS AND VISUALIZATION")
    print("="*60)
    print(f"Data file: {args.data_file}")
    print(f"Results directory: {args.results_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Number of motifs: {args.n_motifs}")
    print("="*60)
    
    # Set random seed
    np.random.seed(args.random_state)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Check if data file exists
    if not os.path.exists(args.data_file):
        print(f"Error: Data file not found: {args.data_file}")
        return
        
    # Check if results directory exists
    if not os.path.exists(args.results_dir):
        print(f"Error: Results directory not found: {args.results_dir}")
        return
    
    # Load data
    print("Loading expression data...")
    adata = sc.read_h5ad(args.data_file)
    print(f"Data shape: {adata.shape}")
    print(f"Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")
    
    # Load or calculate marker exclusion mask
    marker_dir = os.path.join(args.output_dir, "marker_genes")
    exclusion_file = os.path.join(marker_dir, "exclusion_matrix.csv")
    
    if os.path.exists(exclusion_file) and not args.force_recompute_markers:
        print(f"\nFound existing marker genes at: {marker_dir}")
        print("Loading pre-computed marker gene exclusion mask...")
        
        # Load exclusion matrix
        exclusion_df = pd.read_csv(exclusion_file, index_col=0)
        exclusion_mask = exclusion_df.values.T  # Transpose back to (cell_types, genes)
        all_genes = np.array(exclusion_df.index)
        cell_types = np.array(exclusion_df.columns)
        
        print(f"Loaded exclusion mask for {len(cell_types)} cell types and {len(all_genes)} genes")
        
        # Quick validation
        if len(all_genes) != adata.n_vars:
            print(f"Warning: Exclusion mask has {len(all_genes)} genes, but adata has {adata.n_vars} genes")
            print("Gene count mismatch - will recompute exclusion mask")
            raise ValueError("Gene count mismatch")
            
    else:
        if args.force_recompute_markers:
            print("\nForced recompute: Computing marker gene exclusion mask...")
        else:
            print("\nNo existing marker genes found. Computing marker gene exclusion mask...")
        cell_types, all_genes, exclusion_mask = compute_exclusion_mask(
            adata,
            marker_lfc=args.marker_lfc,
            marker_pvalue=args.marker_pvalue,
            marker_subsample=args.marker_subsample
        )
        print(f"Exclusion mask computed for {len(cell_types)} cell types and {len(all_genes)} genes")
    
    # Save marker genes information (only if newly computed)
    if not os.path.exists(exclusion_file):
        print("\nSaving marker genes information...")
        os.makedirs(marker_dir, exist_ok=True)
        
        # Save individual marker gene lists for each cell type
        marker_summary = []
        for i, cell_type in enumerate(cell_types):
            marker_genes = all_genes[exclusion_mask[i, :]]
            
            # Save individual cell type markers
            marker_df = pd.DataFrame({
                'gene': marker_genes,
                'cell_type': cell_type
            })
            marker_file = os.path.join(marker_dir, f"{cell_type.replace('/', '_').replace(' ', '_')}_markers.csv")
            marker_df.to_csv(marker_file, index=False)
            
            # Add to summary
            marker_summary.append({
                'cell_type': cell_type,
                'n_marker_genes': len(marker_genes),
                'marker_percentage': len(marker_genes) / len(all_genes) * 100
            })
            print(f"  {cell_type}: {len(marker_genes)} marker genes ({len(marker_genes)/len(all_genes)*100:.1f}%)")
        
        # Save marker summary
        summary_df = pd.DataFrame(marker_summary)
        summary_file = os.path.join(marker_dir, "marker_genes_summary.csv")
        summary_df.to_csv(summary_file, index=False)
        
        # Save complete exclusion matrix
        exclusion_df = pd.DataFrame(
            exclusion_mask.T,  # Transpose so genes are rows, cell types are columns
            index=all_genes,
            columns=cell_types
        )
        exclusion_df.to_csv(exclusion_file)
        
        # Create a detailed marker genes report
        all_markers_data = []
        for i, cell_type in enumerate(cell_types):
            marker_mask = exclusion_mask[i, :]
            marker_genes = all_genes[marker_mask]
            for gene in marker_genes:
                all_markers_data.append({
                    'cell_type': cell_type,
                    'gene': gene,
                    'is_marker': True
                })
        
        all_markers_df = pd.DataFrame(all_markers_data)
        all_markers_file = os.path.join(marker_dir, "all_marker_genes.csv")
        all_markers_df.to_csv(all_markers_file, index=False)
        
        print(f"Marker genes information saved to: {marker_dir}")
        print(f"  - Individual cell type files: {len(cell_types)} CSV files")
        print(f"  - Summary: marker_genes_summary.csv")
        print(f"  - Complete matrix: exclusion_matrix.csv")
        print(f"  - All markers: all_marker_genes.csv")
    else:
        print(f"Using existing marker genes from: {marker_dir}")
        
        # Print summary from existing files
        summary_file = os.path.join(marker_dir, "marker_genes_summary.csv")
        if os.path.exists(summary_file):
            summary_df = pd.read_csv(summary_file)
            print(f"Marker gene summary:")
            for _, row in summary_df.iterrows():
                print(f"  {row['cell_type']}: {row['n_marker_genes']} marker genes ({row['marker_percentage']:.1f}%)")
    
    # Generate volcano plots
    print("\nGenerating volcano plots...")
    out_volcano = os.path.join(args.output_dir, "volcano_plots_filtered.pdf")
    
    with PdfPages(out_volcano) as pdf:
        for motif_id in range(args.n_motifs):
            files = sorted(glob.glob(os.path.join(
                args.results_dir, f"motif_{motif_id}_celltype_*_de_results.csv"
            )))
            if not files:
                print(f"No results found for motif {motif_id}")
                continue

            fig, axes = plt.subplots(4, 4, figsize=(16, 16))
            axes = axes.flatten()

            for ax, fp in zip(axes, files[:16]):
                df = pd.read_csv(fp)
                ct = os.path.basename(fp)\
                        .split(f"motif_{motif_id}_celltype_")[1]\
                        .replace("_de_results.csv", "")
                
                print(f"  Motif {motif_id}, {ct}: Original genes = {len(df)}")
                
                try:
                    cidx = list(cell_types).index(ct)
                except ValueError:
                    print(f"Warning: Cell type {ct} not found in exclusion mask")
                    continue

                # Keep only genes with sufficient expression in this cell type
                # IMPORTANT: Only filter genes that are actually in the GLM results
                glm_genes_set = set(df['gene'])
                
                subX = adata[adata.obs['cell_type'] == ct].X
                if sp.issparse(subX):
                    expr_frac = np.asarray((subX > 0).sum(axis=0)).ravel() / subX.shape[0]
                else:
                    expr_frac = np.sum(subX > 0, axis=0) / subX.shape[0]
                
                # Only check expression for genes that exist in GLM results
                gene_to_idx = {gene: i for i, gene in enumerate(all_genes)}
                genes_to_keep = []
                for gene in glm_genes_set:
                    if gene in gene_to_idx:
                        gene_idx = gene_to_idx[gene]
                        if expr_frac[gene_idx] >= args.min_expression_frac:
                            genes_to_keep.append(gene)
                
                df = df[df['gene'].isin(genes_to_keep)].copy()
                print(f"    After expression filter: {len(df)} genes")

                # Exclude markers of other cell types
                # Again, only consider genes that exist in GLM results
                other = [i for i in range(len(cell_types)) if i != cidx]
                drop_mask = exclusion_mask[other, :].any(axis=0)
                marker_genes_to_drop = set(all_genes[drop_mask]) & glm_genes_set
                df = df[~df['gene'].isin(marker_genes_to_drop)]
                print(f"    After marker filter: {len(df)} genes")

                # Compute –log10(q) and jitter saturated values
                df['neg_log10_q'] = -np.log10(df['qval'].clip(1e-300))
                m = df['neg_log10_q'] >= 300
                if m.any():
                    df.loc[m, 'neg_log10_q'] = 300 + np.random.normal(0, 15, m.sum())

                if (motif_id == 1) & (ct == 'fibroblast'):
                    print(len(df), df.head())
                # Plot volcano
                volcano_plot(
                    df, 'logFC', 'neg_log10_q', label_col='gene',
                    fdr=args.fdr_threshold, x_threshold=0.2, 
                    marker='o', n_top=30, fontsize=5, ax=ax
                )
                ax.set_title(ct, fontsize=10, pad=5)

            # Remove unused subplots
            for j in range(len(files[:16]), 16):
                fig.delaxes(axes[j])

            plt.suptitle(f"Motif {motif_id} Volcano Plots (marker genes filtered)", fontsize=14)
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig(fig)
            plt.close(fig)
            print(f"Motif {motif_id} volcano plots completed")

    print(f"Volcano plots saved to: {out_volcano}")

    # Generate forest plots
    print("\nGenerating forest plots...")
    out_forest = os.path.join(args.output_dir, 'forest_plots_filtered.pdf')
    
    with PdfPages(out_forest) as pdf:
        for motif_id in range(args.n_motifs):
            csv_paths = sorted(glob.glob(
                os.path.join(args.results_dir, f"motif_{motif_id}_celltype_*_de_results.csv")
            ))
            if not csv_paths:
                continue

            fig, axes = plt.subplots(4, 4, figsize=(16, 24))
            axes = axes.flatten()

            for ax, csv_fp in zip(axes, csv_paths[:16]):
                df = pd.read_csv(csv_fp)
                ct = os.path.basename(csv_fp)\
                       .split(f"motif_{motif_id}_celltype_")[1]\
                       .replace("_de_results.csv", "")
                
                try:
                    ct_idx = list(cell_types).index(ct)
                except ValueError:
                    continue

                # Expression filter
                # IMPORTANT: Only filter genes that are actually in the GLM results
                glm_genes_set = set(df['gene'])
                
                subX = adata[adata.obs['cell_type'] == ct].X
                if sp.issparse(subX):
                    expr_frac = np.asarray((subX > 0).sum(axis=0)).ravel() / subX.shape[0]
                else:
                    expr_frac = (subX > 0).sum(axis=0) / subX.shape[0]
                
                # Only check expression for genes that exist in GLM results
                gene_to_idx = {gene: i for i, gene in enumerate(all_genes)}
                genes_to_keep = []
                for gene in glm_genes_set:
                    if gene in gene_to_idx:
                        gene_idx = gene_to_idx[gene]
                        if expr_frac[gene_idx] >= args.min_expression_frac:
                            genes_to_keep.append(gene)
                
                df = df[df['gene'].isin(genes_to_keep)].copy()

                # Remove markers of other cell types
                # Again, only consider genes that exist in GLM results
                others = [i for i in range(len(cell_types)) if i != ct_idx]
                drop_mask = exclusion_mask[others, :].any(axis=0)
                marker_genes_to_drop = set(all_genes[drop_mask]) & glm_genes_set
                df = df[~df['gene'].isin(marker_genes_to_drop)]

                forest_plot(df, ax=ax, n_top=30)
                # forest_plot(df, effect_col='logFC', se_col='se', effect_is_log2=True)
                ax.set_title(ct, fontsize=9)

            # Remove unused subplots
            for j in range(len(csv_paths[:16]), 16):
                fig.delaxes(axes[j])

            plt.suptitle(f"Motif {motif_id} – Forest plots (marker genes filtered)",
                         fontsize=14)
            plt.tight_layout(rect=[0, 0, 1, 0.965])
            pdf.savefig(fig)
            plt.close(fig)
            print(f"Motif {motif_id} forest plots completed")

    print(f"Forest plots saved to: {out_forest}")
    
    print("\n" + "="*60)
    print("GLM RESULTS ANALYSIS COMPLETED SUCCESSFULLY!")
    print("="*60)
    print(f"Results saved to: {args.output_dir}")
    print("Files created:")
    # print(f"  - {os.path.basename(out_volcano)} (volcano plots)")
    print(f"  - {os.path.basename(out_forest)} (forest plots)")


if __name__ == '__main__':
    main()