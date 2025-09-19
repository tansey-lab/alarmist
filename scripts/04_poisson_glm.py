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
    results['patch_loadings'] = np.load(os.path.join(results_dir, 'patch_loadings.npy'))
    
    # Load patch motifs info
    results['patch_motifs'] = pd.read_csv(os.path.join(results_dir, 'patch_motifs.csv'))
    
    # Load LRI motifs info
    results['lri_motifs'] = pd.read_csv(os.path.join(results_dir, 'lri_motifs.csv'))
    
    print(f"✓ Loaded results with {results['patch_loadings'].shape[0]} patches, "
          f"{results['patch_loadings'].shape[1]} motifs")
    
    return results


def extract_lri_genes(lri_names: pd.Series, spliter = '|') -> set:
    """
    Extract LRI genes from LRI names
    
    Args:
        lri_names: Series of LRI names
    """
    print("Extracting LRI genes...")
    
    lri_genes = set()
    sample_names = lri_names.head(5).tolist()
    print(f"Sample LRI names: {sample_names}")
    
    # Try to detect format from first few names
    sample_splits = [name.split(spliter) for name in sample_names]
    avg_parts = np.mean([len(parts) for parts in sample_splits])

    for name in lri_names:
        if pd.isna(name):
            continue
            
        parts = name.split(spliter)
        lri_genes.update(parts[2:4])
    
    print(f"✓ Extracted {len(lri_genes)} unique LRI genes")
    return lri_genes


def run_univariate_de_sklearn_by_celltype(
    cell_patch: pd.DataFrame,
    counts: np.ndarray,
    gene_names: List[str],
    non_lri_genes: List[str],
    n_motifs: int,
    output_dir: Optional[str] = None,
    alpha: float = 0.05,
    print_every: int = 20
) -> Dict[str, pd.DataFrame]:
    """
    Run univariate DE by motif AND cell type using scikit-learn's PoissonRegressor.
    Saves a CSV for each motif and cell_type combo immediately.
    """
    print("Running univariate DE with scikit-learn PoissonRegressor (by cell type)...")
    idx_map = {g: i for i, g in enumerate(gene_names)}
    de_results = {}

    # iterate motifs (example limited to motif 7 if desired)
    for k in range(n_motifs):
        for ct in cell_patch['cell_type'].unique():
            if ct == 'granulocyte':
                print(f"Skipping cell type '{ct}' for motif {k} (granulocytes not included)")
                continue
            subset_mask = (cell_patch['cell_type'] == ct)
            X_all = cell_patch.loc[subset_mask, f'prog_{k}_loading'].values
            counts_all = counts[subset_mask, :]

            valid = ~np.isnan(X_all)
            # X = np.log(X_all[valid]).reshape(-1, 1)

            X_log = np.log(X_all[valid])
            # Z-score within this cell type
            mu = X_log.mean()
            sigma = X_log.std(ddof=0) if X_log.std(ddof=0) > 0 else 1.0
            X = ((X_log - mu) / sigma).reshape(-1, 1)

            Y = counts_all[valid, :]

            if X.size == 0:
                print(f"motif {k}, cell_type '{ct}': no valid data, skipping")
                continue

            print(f"\n🚀 motif {k}, cell_type '{ct}': {X.shape[0]} cells")
            genes, coefs, pvals, ses = [], [], [], []
            total = len(non_lri_genes)

            for idx, gene in enumerate(non_lri_genes, 1):
                gi = idx_map.get(gene)
                if gi is None:
                    print(gi,'is none')
                    continue
                y = Y[:, gi]
                # try:
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
                # except:
                #     pass

                if idx % print_every == 0 or idx == total:
                    print(f"   ✔ {idx}/{total} genes processed")

            df = pd.DataFrame({
                'gene': genes,
                'logFC': coefs,
                'se': ses,              # natural-log SE
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


def prepare_cell_data_memory_efficient(cell_patch_file: str, 
                                       patch_map: pd.DataFrame, 
                                       patch_loadings: np.ndarray, 
                                       adata_file: str,
                                       keep_sparse: bool = True) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    """
    Memory-efficient cell data preparation with sparse matrix support
    
    IMPORTANT: The 'cell_id' column in cell_patch_correspondence.csv now contains
    actual adata.obs['cell_id'] values. The file maintains the same row order as 
    adata.obs when possible, enabling direct 1:1 alignment, with (cell_id, tma_id) 
    pairs used for robust mapping verification.
    """
    import scanpy as sc
    import gc
    
    print("=== Memory-Efficient Data Loading ===")
    
    # Check initial memory
    print("Initial memory check:")
    check_memory_usage()
    
    print("Loading cell-patch correspondence...")
    cell_patch = pd.read_csv(cell_patch_file)
    
    print("Mapping patch IDs to indices...")
    cell_patch = cell_patch.merge(patch_map, on='patch_id', how='left')
    
    # Check for missing mappings
    missing_patches = cell_patch['patch_idx'].isna().sum()
    if missing_patches > 0:
        print(f"Warning: {missing_patches} cells have missing patch mappings")
        cell_patch = cell_patch.dropna(subset=['patch_idx'])
    
    print(f"Valid cells after filtering: {len(cell_patch)}")
    
    # Add motif loadings for each cell
    print("Adding motif loadings to cells...")
    n_motifs = patch_loadings.shape[1]
    
    def get_loading(patch_idx, k):
        if pd.isna(patch_idx):
            return np.nan
        patch_idx = int(patch_idx)
        if 0 <= patch_idx < patch_loadings.shape[0]:
            return patch_loadings[patch_idx, k]
        else:
            return np.nan
    
    for k in range(n_motifs):
        cell_patch[f'prog_{k}_loading'] = cell_patch['patch_idx'].apply(
            lambda i: get_loading(i, k)
        )

    print("Memory after loading cell metadata:")
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
    
    # CRITICAL: Verify cell alignment using actual cell_id values
    print("Verifying cell-patch correspondence alignment...")
    
    # Now cell_patch_correspondence.csv contains actual adata.obs['cell_id'] values
    # We need to create a proper mapping using (cell_id, tma_id) pairs
    print(f"AnnData shape: {adata.shape}")
    print(f"Cell-patch correspondence rows: {len(cell_patch)}")
    
    if len(adata) != len(cell_patch):
        raise ValueError(f"Mismatch: AnnData has {len(adata)} cells, cell_patch has {len(cell_patch)} rows")
    
    # Verify that cell_id matches adata.obs['cell_id'] and maintain row order compatibility
    adata_cell_ids = adata.obs['cell_id'].astype(str)
    adata_tma_ids = adata.obs['tma_id'].astype(str)
    patch_cell_ids = cell_patch['cell_id'].astype(str)
    patch_tma_ids = cell_patch['tma_id'].astype(str)
    
    # Create mapping keys for verification
    adata_keys = [f"{cid}_{tid}" for cid, tid in zip(adata_cell_ids, adata_tma_ids)]
    patch_keys = [f"{cid}_{tid}" for cid, tid in zip(patch_cell_ids, patch_tma_ids)]
    
    # Check if they match exactly (in order) - this preserves row alignment
    if adata_keys == patch_keys:
        print("✓ Cell-patch correspondence verified - perfect 1:1 alignment with (cell_id, tma_id) keys")
    else:
        print("Warning: Row order doesn't match exactly")
        print("Checking if all cells are present with valid (cell_id, tma_id) mappings...")
        
        # Check if all mapping keys exist in both
        adata_key_set = set(adata_keys)
        patch_key_set = set(patch_keys)
        
        missing_in_adata = patch_key_set - adata_key_set
        missing_in_patch = adata_key_set - patch_key_set
        
        if missing_in_adata:
            print(f"Error: {len(missing_in_adata)} (cell_id, tma_id) pairs in patch file not found in adata")
            raise ValueError("Cell ID mapping mismatch detected!")
        
        if missing_in_patch:
            print(f"Error: {len(missing_in_patch)} adata (cell_id, tma_id) pairs not found in patch file") 
            raise ValueError("Cell ID mapping mismatch detected!")
        
        print("All (cell_id, tma_id) pairs present but order differs")
        print("This maintains correct mapping but order preservation is recommended for the pipeline")
        
    print("✓ Cell-patch mapping verified using (cell_id, tma_id) keys")
    print("Memory after verification:")
    check_memory_usage()
    
    # Handle count matrix
    print("Processing count matrix...")
    counts = adata.layers["counts"]
    print("Using 'counts' layer (raw counts) for analysis")
    
    if sp.issparse(counts):
        if keep_sparse:
            print("Keeping data in sparse format")
            # Keep sparse - we'll handle this differently
            counts_sparse = counts
            gene_names = list(adata.var_names)
            
            # Clean up adata
            del adata
            gc.collect()
            
            print("Memory after cleanup (sparse):")
            check_memory_usage()
            
            return cell_patch, counts_sparse, gene_names
        else:
            print("Converting sparse to dense (warning: high memory usage)")
            print("Memory before sparse->dense conversion:")
            check_memory_usage()
            
            counts = counts.toarray()
            
            print("Memory after sparse->dense conversion:")
            check_memory_usage()
    
    # Get gene names
    gene_names = list(adata.var_names)
    
    # Clean up
    del adata
    gc.collect()
    
    print("Final memory after cleanup:")
    check_memory_usage()
    
    return cell_patch, counts, gene_names



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
    lri_genes = extract_lri_genes(bptf_results['lri_motifs']['lri_name'], 
                                 args.spliter)
    
    # Get all non-LRI genes
    all_genes = list(adata.var_names)
    non_lri_genes = [g for g in all_genes if g not in lri_genes]
    print(f"Analyzing {len(non_lri_genes)} non-LRI genes")

    del adata  # Free memory after loading AnnData
    
    # Prepare cell data
    patch_map = bptf_results['patch_motifs'][['patch_idx', 'patch_id']].drop_duplicates()
    
    # CRITICAL: Validate patch indices alignment
    print("Validating patch indices...")
    print(f"BPTF patch_loadings shape: {bptf_results['patch_loadings'].shape}")
    print(f"Patch map contains {len(patch_map)} unique patches")
    print(f"patch_idx range: {patch_map['patch_idx'].min()} to {patch_map['patch_idx'].max()}")
    
    # Check if patch indices are within bounds
    max_patch_idx = patch_map['patch_idx'].max()
    n_patches_bptf = bptf_results['patch_loadings'].shape[0]
    
    if max_patch_idx >= n_patches_bptf:
        raise ValueError(f"patch_idx max ({max_patch_idx}) exceeds BPTF patch_loadings rows ({n_patches_bptf})")
    
    print("✓ Patch indices validation passed")
    
    cell_patch_file = os.path.join(patch_lri_dir, "cell_patch_correspondence.csv")
    cell_patch, counts, gene_names = prepare_cell_data_memory_efficient(
        cell_patch_file, patch_map, 
        bptf_results['patch_loadings'], data_file,
        keep_sparse=False
    )
    
    n_motifs = bptf_results['patch_loadings'].shape[1]

    print(cell_patch.head())
    
    # Run DE analysis
    if args.mode == 'univariate':
        results = run_univariate_de_sklearn_by_celltype(
                    cell_patch, counts, gene_names, non_lri_genes, n_motifs,
                    output_dir=output_dir,
                    print_every=100
                )
    else:
        print("Multivariate analysis mode is not yet implemented. " \
        "Please use '--mode univariate' for now.")
    
    # Save results
    save_results(results, output_dir, args.mode)
    
    print("✓ Analysis completed!")


if __name__ == "__main__":
    main()