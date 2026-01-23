#!/usr/bin/env python3
"""
BPTF-based Differential Expression Analysis

Performs differential expression analysis using BPTF factor loadings as continuous 
covariates in Poisson GLM regression. Supports both univariate and multivariate modes.

Usage:
    python bptf_de_analysis.py --mode univariate --results_dir ../data/bptf_results
    python bptf_de_analysis.py --mode multivariate --alpha 0.01
"""

import pandas as pd
import numpy as np
import scipy.sparse as sp
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from joblib import Parallel, delayed
import warnings
from sklearn.linear_model import PoissonRegressor
from scipy.stats import norm
warnings.filterwarnings('ignore')


def load_bptf_results(results_dir: str) -> Dict:
    """Load BPTF results from directory"""
    print(f"Loading BPTF results from {results_dir}...")
    
    # Adjust these paths based on your actual file structure
    results = {}
    
    # Load patch loadings (W matrix)
    results['cell_loadings'] = np.load(os.path.join(results_dir, 'cell_loadings.npy'))
    
    # Load LRI motifs info
    # results['lri_motifs'] = pd.read_csv(os.path.join(results_dir, 'cell_lri_columns.csv'))
    
    print(f"✓ Loaded results with {results['cell_loadings'].shape[0]} cells, "
          f"{results['cell_loadings'].shape[1]} motifs")
    
    return results

def extract_lri_genes(lri_names: pd.Series, splitter: str = '|', gene_splitter: str = '_') -> set:
    """
    Extract unique ligand and receptor genes from LRI names.

    Expected LRI name format:
        sender_celltype | receiver_celltype | ligand | receptor | mode

    - Ligand and receptor fields may contain multiple genes joined by `gene_splitter`
      (e.g., "TGFBR2_TGFBR1").
    
    Args
    ----
    lri_names : pd.Series
        A series of LRI name strings.
    splitter : str, default '|'
        Delimiter separating fields in each LRI name.
    gene_splitter : str, default '_'
        Delimiter separating multiple genes within ligand/receptor fields.
    
    Returns
    -------
    set
        A set of unique gene symbols extracted from ligand and receptor fields.
    """
    
    print("Extracting LRI genes...")

    lri_genes = set()

    # Show sample names for debugging
    sample_names = lri_names.head(5).tolist()
    print(f"Sample LRI names: {sample_names}")

    for name in lri_names:
        if pd.isna(name):
            continue
        
        parts = str(name).split(splitter)
        # Expect at least 4 fields: ct1, ct2, ligand, receptor
        if len(parts) < 4:
            continue
        
        ligand_field = parts[2]
        receptor_field = parts[3]

        # Split multi-gene ligand/receptor entries
        ligand_genes = [g.strip() for g in ligand_field.split(gene_splitter) if g.strip()]
        receptor_genes = [g.strip() for g in receptor_field.split(gene_splitter) if g.strip()]

        # Add to final set
        lri_genes.update(ligand_genes)
        lri_genes.update(receptor_genes)

    print(f"✓ Extracted {len(lri_genes)} unique LRI genes")
    return lri_genes

def run_univariate_de_sklearn_by_celltype(
    cell_df: pd.DataFrame,
    counts: np.ndarray,
    gene_names: List[str],
    non_lri_genes: List[str],
    n_motifs: int,
    output_dir: Optional[str] = None,
    alpha: float = 0.05,
    print_every: int = 1000
) -> Dict[str, pd.DataFrame]:
    """
    Run univariate DE by motif AND cell type using scikit-learn's PoissonRegressor.
    Now uses direct cell-level loadings instead of patch-inherited loadings.
    """
    print("Running univariate DE with scikit-learn PoissonRegressor (by cell type)...")
    idx_map = {g: i for i, g in enumerate(gene_names)}
    de_results = {}

    # iterate motifs
    for k in range(n_motifs):
        for ct in cell_df['cell_type'].unique():
            if ct == 'granulocyte':
                print(f"Skipping cell type '{ct}' for motif {k} (granulocytes not included)")
                continue
            
            subset_mask = (cell_df['cell_type'] == ct)
            X_all = cell_df.loc[subset_mask, f'prog_{k}_loading'].values
            counts_all = counts[subset_mask, :]

            # if negative loadings
            valid = ~np.isnan(X_all) & (X_all > 0)
            
            if valid.sum() == 0:
                print(f"motif {k}, cell_type '{ct}': no valid positive loadings, skipping")
                continue

            X_log = np.log(X_all[valid])
            # Z-score within this cell type
            mu = X_log.mean()
            sigma = X_log.std(ddof=0) if X_log.std(ddof=0) > 0 else 1.0
            X = ((X_log - mu) / sigma).reshape(-1, 1)

            Y = counts_all[valid, :]

            print(f"\n🚀 motif {k}, cell_type '{ct}': {X.shape[0]} cells")
            genes, coefs, pvals, ses = [], [], [], []
            total = len(non_lri_genes)

            for idx, gene in enumerate(non_lri_genes, 1):
                gi = idx_map.get(gene)
                if gi is None:
                    continue
                    
                y = Y[:, gi]
                
                model = PoissonRegressor(alpha=0.0, fit_intercept=True,
                                        max_iter=2000, tol=1e-6)
                model.fit(X, y)
                beta1 = model.coef_[0]
                mu = model.predict(X)
                I = np.sum(mu * (X.ravel()**2))
                se = np.sqrt(1.0 / I) if I > 0 else np.nan
                z = beta1 / se if se and se > 0 else 0.0
                pval = 2 * (1 - norm.cdf(abs(z)))
                
                genes.append(gene)
                coefs.append(beta1)
                pvals.append(pval)
                ses.append(se)

                if idx % print_every == 0 or idx == total:
                    print(f"   ✔ {idx}/{total} genes processed")

            df = pd.DataFrame({
                'gene': genes,
                'logFC': coefs,
                'se': ses,
                'pval': pvals
            })

            if not df.empty:
                reject, qvals, _, _ = multipletests(df['pval'], alpha=alpha, method='fdr_bh')
                df['qval'] = qvals
                df['significant'] = reject
                df = df.sort_values('qval')

            key = f'motif_{k}_celltype_{ct}'
            de_results[key] = df

            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                fname = f"{key}_de_results.csv"
                path = os.path.join(output_dir, fname)
                df.to_csv(path, index=False)
                print(f"✓ Saved {fname}")

    return de_results

def save_cell_type_result(df: pd.DataFrame, motif_idx: int, cell_type: str, output_dir: str):
    """Save results for a single motif-celltype combination immediately"""
    
    # Create subdirectory for this motif
    motif_dir = os.path.join(output_dir, f"motif_{motif_idx}")
    os.makedirs(motif_dir, exist_ok=True)
    
    # Save with cell type in filename (replace special characters)
    safe_cell_type = cell_type.replace('/', '_').replace(' ', '_').replace('+', 'pos')
    filename = f"motif_{motif_idx}_{safe_cell_type}_de_results.csv"
    filepath = os.path.join(motif_dir, filename)
    
    df.to_csv(filepath, index=False)
    print(f"      💾 Saved: {filename}")


def generate_analysis_summary(de_results: Dict, output_dir: str):
    """Generate summary statistics across all motifs and cell types"""
    
    print("\n📋 Generating analysis summary...")
    
    summary_data = []
    volcano_summary = []
    
    for motif_key, motif_data in de_results.items():
        motif_idx = int(motif_key.split('_')[1])
        
        for cell_type, df in motif_data.items():
            if len(df) == 0:
                continue
                
            # Basic statistics
            summary_data.append({
                'motif': motif_idx,
                'cell_type': cell_type,
                'total_genes': len(df),
                'significant_genes': df['significant'].sum(),
                'sig_percentage': df['significant'].sum() / len(df) * 100,
                'min_qval': df['qval'].min(),
                'max_abs_logFC': df['abs_logFC'].max(),
                'upregulated': ((df['logFC'] > 0) & df['significant']).sum(),
                'downregulated': ((df['logFC'] < 0) & df['significant']).sum()
            })
            
            # Volcano plot summary (top hits)
            sig_df = df[df['significant']].copy()
            if len(sig_df) > 0:
                # Top upregulated
                top_up = sig_df[sig_df['logFC'] > 0].nlargest(5, 'abs_logFC')
                # Top downregulated  
                top_down = sig_df[sig_df['logFC'] < 0].nlargest(5, 'abs_logFC')
                
                for _, row in pd.concat([top_up, top_down]).iterrows():
                    volcano_summary.append({
                        'motif': motif_idx,
                        'cell_type': cell_type,
                        'gene': row['gene'],
                        'logFC': row['logFC'],
                        'qval': row['qval'],
                        'neg_log10_qval': row['neg_log10_qval'],
                        'direction': 'up' if row['logFC'] > 0 else 'down'
                    })
    
    # Save summary files
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, "analysis_summary.csv"), index=False)
    
    volcano_df = pd.DataFrame(volcano_summary)
    volcano_df.to_csv(os.path.join(output_dir, "volcano_plot_highlights.csv"), index=False)
    
    print(f"   💾 Saved analysis_summary.csv ({len(summary_df)} motif-celltype combinations)")
    print(f"   💾 Saved volcano_plot_highlights.csv ({len(volcano_df)} top hits)")
    
    # Print top-level summary
    total_sig = summary_df['significant_genes'].sum()
    total_tested = summary_df['total_genes'].sum()
    
    print(f"\n📈 Final Summary:")
    print(f"   • Total tests performed: {total_tested:,}")
    print(f"   • Total significant genes: {total_sig:,}")
    print(f"   • Overall significance rate: {total_sig/total_tested*100:.2f}%")
    print(f"   • motifs analyzed: {summary_df['motif'].nunique()}")
    print(f"   • Cell types analyzed: {summary_df['cell_type'].nunique()}")



def save_results(results, output_dir: str, mode: str):
    """Save DE results"""
    os.makedirs(output_dir, exist_ok=True)
    
    if mode == 'univariate':
        for motif, df in results.items():
            filename = f"{motif}_de_results.csv"
            filepath = os.path.join(output_dir, filename)
            # Skip if already saved (e.g., motif_0)
            if not os.path.exists(filepath):
                df.to_csv(filepath, index=False)
                print(f"✓ Saved {filename}")
            else:
                print(f"✓ {filename} already exists, skipping")
        
        long_data = []
        for motif, df in results.items():
            motif_num = int(motif.split('_')[1])  # extract number
            df_copy = df.copy()
            df_copy['motif'] = motif_num
            long_data.append(df_copy)
        
        long_df = pd.concat(long_data, ignore_index=True)
        # 重排列顺序
        long_df = long_df[['gene', 'motif', 'logFC', 'pval', 'qval', 'significant']]
        long_df = long_df.sort_values(['motif', 'qval'])
        
        long_df.to_csv(os.path.join(output_dir, "univariate_de_long.csv"), index=False)
        print("✓ Saved univariate_de_long.csv")
    
    elif mode == 'multivariate':
        filename = "multivariate_de_results.csv"
        results.to_csv(os.path.join(output_dir, filename), index=False)
        print(f"✓ Saved {filename}")


def prepare_cell_data_memory_efficient(cell_loadings: np.ndarray, 
                                       adata_file: str,
                                       count_layer: str = 'X',
                                       cell_metadata_file: str = None,
                                       keep_sparse: bool = True) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    """
    Memory-efficient cell data preparation with sparse matrix support
    Now directly uses cell-level loadings instead of patch-to-cell mapping
    """
    import scanpy as sc
    import gc
    
    print("=== Memory-Efficient Data Loading ===")
    
    # Check initial memory
    print("Initial memory check:")
    check_memory_usage()
    
    # Load AnnData with memory monitoring
    print(f"Loading AnnData from {adata_file}...")
    print("This may take a while for large files...")
    
    adata = sc.read_h5ad(adata_file)
    print(f"Full AnnData shape: {adata.shape}")
    print(f"Data type: {type(adata.X)}")
    print(f"Is sparse: {sp.issparse(adata.X)}")
    
    print("Memory after loading full AnnData:")
    check_memory_usage()
    
    print("Verifying cell loadings alignment...")
    if cell_loadings.shape[0] != adata.shape[0]:
        raise ValueError(f"Mismatch: cell_loadings has {cell_loadings.shape[0]} cells, "
                        f"adata has {adata.shape[0]} cells")
    
    print(f"✓ Cell loadings shape: {cell_loadings.shape}")
    print(f"✓ AnnData shape: {adata.shape}")
    
    # create cell_df, directly use cell_loadings
    print("Creating cell dataframe with loadings...")
    n_cells, n_motifs = cell_loadings.shape
    
    # Initialize cell_df with adata metadata
    cell_df = pd.DataFrame({
        'cell_id': adata.obs.index.astype(str),
        'cell_type': adata.obs['cell_type'] if 'cell_type' in adata.obs.columns else 'unknown'
    })
    
    # If an additional cell-metadata file is provided 
    # (for example, the original cell_patch_correspondence), 
    # it can be merged
    if cell_metadata_file is not None:
        print(f"Loading additional cell metadata from {cell_metadata_file}...")
        cell_meta = pd.read_csv(cell_metadata_file)
        # avoid repetition
        if 'cell_type' in cell_meta.columns and 'cell_type' not in adata.obs.columns:
            cell_df['cell_type'] = cell_meta['cell_type'].values
    
    # Directly add the cell-level loadings to the dataframe
    print("Adding cell-level motif loadings...")
    for k in range(n_motifs):
        cell_df[f'prog_{k}_loading'] = cell_loadings[:, k]
    
    print(f"✓ Added {n_motifs} motif loadings for {n_cells} cells")
    print("Memory after creating cell dataframe:")
    check_memory_usage()
    
    # Handle count matrix according to count_layer
    print(f"Processing count matrix from '{count_layer}'...")
    
    if count_layer == 'X':
        counts = adata.X
        print("Using adata.X for analysis")
    elif count_layer == 'raw':
        if adata.raw is None:
            raise ValueError("adata.raw is None but count_layer='raw' was specified")
        counts = adata.raw.X
        print("Using adata.raw.X (raw counts) for analysis")
    elif count_layer.startswith('layer'):
        layer_name = count_layer.split(':', 1)[1]
        if layer_name not in adata.layers:
            raise ValueError(f"Layer '{layer_name}' not found in adata.layers. "
                           f"Available layers: {list(adata.layers.keys())}")
        counts = adata.layers[layer_name]
        print(f"Using adata.layers['{layer_name}'] for analysis")
    else:
        raise ValueError(f"Invalid count_layer value: '{count_layer}'. "
                       f"Use 'X', 'raw', or 'layer:NAME'")
    
    print(f"Count matrix shape: {counts.shape}")
    print(f"Count matrix type: {type(counts)}")
    print(f"Is sparse: {sp.issparse(counts)}")

    # Get gene names
    if count_layer == 'raw' and adata.raw is not None:
        gene_names = list(adata.raw.var_names)
        print(f"Using gene names from adata.raw.var_names ({len(gene_names)} genes)")
    else:
        gene_names = list(adata.var_names)
        print(f"Using gene names from adata.var_names ({len(gene_names)} genes)")
    
    if sp.issparse(counts):
        if keep_sparse:
            print("Keeping data in sparse format")
            counts_sparse = counts
            gene_names = list(adata.var_names)
            
            # Clean up adata
            del adata
            gc.collect()
            
            print("Memory after cleanup (sparse):")
            check_memory_usage()
            
            return cell_df, counts_sparse, gene_names
        else:
            print("Converting sparse to dense (warning: high memory usage)")
            print("Memory before sparse->dense conversion:")
            check_memory_usage()
            
            counts = counts.toarray()
            
            print("Memory after sparse->dense conversion:")
            check_memory_usage()
    
    # Clean up
    del adata
    gc.collect()
    
    print("Final memory after cleanup:")
    check_memory_usage()
    
    return cell_df, counts, gene_names


def check_memory_usage():
    """Check current memory usage"""
    try:
        import psutil
        import os
        
        process = psutil.Process(os.getpid())
        memory_mb = process.memory_info().rss / 1024 / 1024
        available_mb = psutil.virtual_memory().available / 1024 / 1024
        
        print(f"   Memory usage: {memory_mb:.1f} MB")
        print(f"   Available: {available_mb:.1f} MB")
        
        return memory_mb, available_mb
    except ImportError:
        print("   psutil not available for memory monitoring")
        return None, None


def main():
    parser = argparse.ArgumentParser(description='BPTF Differential Expression Analysis')
    parser.add_argument('--mode', choices=['univariate', 'multivariate'], 
                       default='univariate', help='Analysis mode')
    parser.add_argument('--results-dir', default='results/bptf_patch_lri',
                       help='BPTF results directory')
    parser.add_argument('--patch-lri-dir', default=None,
                       help='Patch-LRI results directory (if different from BPTF dir parent)')
    parser.add_argument('--data-file', default='data/processed/preprocessed_xenium_data_subset.h5ad',
                        help='AnnData file path (h5ad format)')
    parser.add_argument('--output-dir', default='results/bptf_de_results',
                       help='Output directory')
    parser.add_argument('--alpha', type=float, default=0.05,
                       help='FDR significance threshold')
    parser.add_argument('--spliter', type=str,
                       default='|', help='cell-gene or cell|gene ...')
    parser.add_argument('--random-state', type=int, default=42,
                       help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    np.random.seed(args.random_state)
    
    print(f"Starting BPTF DE analysis in {args.mode} mode...")
    # Fix argument attribute names (argparse converts - to _)
    results_dir = getattr(args, 'results_dir', None)
    patch_lri_dir = getattr(args, 'patch_lri_dir', None)
    output_dir = getattr(args, 'output_dir', None) 
    data_file = getattr(args, 'data_file', None)
    
    # Determine patch-LRI directory
    if patch_lri_dir is None:
        # Try to infer from BPTF directory structure: results/project/bptf -> results/project/patch_lri
        results_path = Path(results_dir)
        if results_path.name == 'bptf':
            patch_lri_dir = str(results_path.parent / 'patch_lri')
        else:
            patch_lri_dir = results_dir  # Fallback to same directory
    
    print(f"BPTF results directory: {results_dir}")
    print(f"Patch-LRI directory: {patch_lri_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Random seed: {args.random_state}")
    
    # Load BPTF results
    bptf_results = load_bptf_results(results_dir)
    lri_columns = pd.read_csv(os.path.join(patch_lri_dir, 'patch_lri_columns.csv'))
    
    # Load AnnData (you'll need to implement this part based on your data)
    print("Loading AnnData...")
    if data_file:
        import scanpy as sc
        adata = sc.read_h5ad(data_file)
    else:
        # Placeholder - replace with your actual adata loading
        print("Warning: No AnnData file specified. Please modify the script to load your data.")
        return
    
    # Extract LRI genes
    lri_genes = extract_lri_genes(lri_columns['column_name'], 
                                 args.spliter)
    
    # Get all non-LRI genes
    all_genes = list(adata.var_names)
    non_lri_genes = [g for g in all_genes if g not in lri_genes]
    print(f"Analyzing {len(non_lri_genes)} non-LRI genes")

    del adata  # Free memory after loading AnnData
    
    # Prepare cell data
    # patch_map = bptf_results['patch_motifs'][['patch_idx', 'patch_id']].drop_duplicates()
    
    # CRITICAL: Validate patch indices alignment
    print("Validating patch indices...")
    print(f"BPTF cell_loadings shape: {bptf_results['cell_loadings'].shape}")
    # print(f"Patch map contains {len(patch_map)} unique patches")
    # print(f"patch_idx range: {patch_map['patch_idx'].min()} to {patch_map['patch_idx'].max()}")
    
    # Check if patch indices are within bounds
    # max_patch_idx = patch_map['patch_idx'].max()
    # n_patches_bptf = bptf_results['patch_loadings'].shape[0]
    
    # if max_patch_idx >= n_patches_bptf:
    #     raise ValueError(f"patch_idx max ({max_patch_idx}) exceeds BPTF patch_loadings rows ({n_patches_bptf})")
    
    # print("✓ Patch indices validation passed")
    
    # cell_patch_file = os.path.join(patch_lri_dir, "cell_patch_correspondence.csv")
    cell_df, counts, gene_names = prepare_cell_data_memory_efficient(
        bptf_results['cell_loadings'],  # 直接传入cell_loadings
        data_file,
        count_layer = 'layers:counts',
        # cell_metadata_file=cell_patch_file,  # 可选：如果需要额外的metadata
        keep_sparse=False
    )
    n_motifs = bptf_results['cell_loadings'].shape[1]

    print(cell_df.head())
    
    # Run DE analysis  
    results = run_univariate_de_sklearn_by_celltype(
        cell_df,
        counts, 
        gene_names, 
        non_lri_genes, 
        n_motifs,
        output_dir=output_dir,
        print_every=1000
    )
    
    # Save results
    save_results(results, output_dir, args.mode)
    
    print("✓ Analysis completed!")


if __name__ == "__main__":
    main()