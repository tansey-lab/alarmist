"""
Poisson GLM differential expression analysis

Based on scripts/04_poisson_glm.py
"""

import anndata
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.stats as stats
import scanpy as sc
import gc
import os
from typing import Dict, List, Tuple, Optional
from sklearn.linear_model import PoissonRegressor
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
from pathlib import Path
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Import plotting functions
from alarmist.plotting import glm_plots


def load_bptf_results(results_dir: str) -> Dict:
    """
    Load BPTF results from directory

    Parameters
    ----------
    results_dir : str
        Directory containing BPTF results

    Returns
    -------
    dict
        Dictionary with cell_loadings array
    """
    import os

    print(f"Loading BPTF results from {results_dir}...")

    results = {}

    # Load cell loadings (W matrix)
    results['cell_loadings'] = np.load(os.path.join(results_dir, 'cell_loadings.npy'))

    print(f"✓ Loaded results with {results['cell_loadings'].shape[0]} cells, "
          f"{results['cell_loadings'].shape[1]} motifs")

    return results


def extract_lri_genes(lri_names: pd.Series, splitter: str = '|',
                     gene_splitter: str = '_') -> set:
    """
    Extract unique ligand and receptor genes from LRI names

    Expected LRI name format:
        sender_celltype | receiver_celltype | ligand | receptor | mode

    Parameters
    ----------
    lri_names : pd.Series
        Series of LRI name strings
    splitter : str, default '|'
        Delimiter separating fields in each LRI name
    gene_splitter : str, default '_'
        Delimiter separating multiple genes within ligand/receptor fields

    Returns
    -------
    set
        Set of unique gene symbols from ligand and receptor fields
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
    alpha: float = 0.05
) -> Dict[str, pd.DataFrame]:
    """
    Run univariate DE by motif AND cell type using scikit-learn's PoissonRegressor

    Uses direct cell-level loadings instead of patch-inherited loadings

    Parameters
    ----------
    cell_df : pd.DataFrame
        DataFrame with cell metadata and motif loadings
    counts : np.ndarray
        Expression count matrix (cells × genes)
    gene_names : List[str]
        List of gene names
    non_lri_genes : List[str]
        List of non-LRI genes to analyze
    n_motifs : int
        Number of motifs
    output_dir : str, optional
        Output directory for results
    alpha : float, default 0.05
        FDR significance threshold

    Returns
    -------
    dict
        Dictionary mapping motif-celltype combinations to DE results DataFrames
    """
    import os

    print("Running univariate DE with scikit-learn PoissonRegressor (by cell type)...")
    idx_map = {g: i for i, g in enumerate(gene_names)}
    de_results = {}

    # Iterate motifs
    for k in range(n_motifs):
        for ct in cell_df['cell_type'].unique():
            if ct == 'granulocyte':
                print(f"Skipping cell type '{ct}' for motif {k} (granulocytes not included)")
                continue

            subset_mask = (cell_df['cell_type'] == ct)
            X_all = cell_df.loc[subset_mask, f'prog_{k}_loading'].values
            counts_all = counts[subset_mask, :]

            # Filter negative loadings
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

            for gene in tqdm(non_lri_genes, desc=f"  Motif {k}, {ct}", unit="gene", leave=False):
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


def prepare_cell_data_from_adata(
    cell_loadings: np.ndarray,
    adata: 'anndata.AnnData',
    count_layer: str = 'X',
    keep_sparse: bool = True
) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    """
    Prepare cell data from already-loaded AnnData object

    Parameters
    ----------
    cell_loadings : np.ndarray
        Cell-level BPTF loadings (n_cells × n_motifs)
    adata : anndata.AnnData
        Annotated data object
    count_layer : str, default 'X'
        Which layer to use: 'X', 'raw', or 'layers:NAME'
    keep_sparse : bool, default True
        Keep counts in sparse format

    Returns
    -------
    tuple
        (cell_df, counts, gene_names)
    """
    import gc

    print("=== Preparing Cell Data ===")

    # Validate alignment
    if cell_loadings.shape[0] != adata.shape[0]:
        raise ValueError(f"Mismatch: cell_loadings has {cell_loadings.shape[0]} cells, "
                        f"adata has {adata.shape[0]} cells")

    print(f"✓ Cell loadings shape: {cell_loadings.shape}")
    print(f"✓ AnnData shape: {adata.shape}")

    # Create cell_df with direct cell loadings
    n_cells, n_motifs = cell_loadings.shape

    # Initialize cell_df with adata metadata
    cell_df = pd.DataFrame({
        'cell_id': adata.obs.index.astype(str),
        'cell_type': adata.obs['cell_type'] if 'cell_type' in adata.obs.columns else 'unknown'
    })

    # Add cell-level motif loadings
    print(f"Adding {n_motifs} motif loadings...")
    for k in range(n_motifs):
        cell_df[f'prog_{k}_loading'] = cell_loadings[:, k]

    print(f"✓ Cell dataframe created: {cell_df.shape}")

    # Handle count matrix
    print(f"Processing count matrix from '{count_layer}'...")

    if count_layer == 'X':
        counts = adata.X
        print("Using adata.X for analysis")
    elif count_layer == 'raw':
        if adata.raw is None:
            raise ValueError("adata.raw is None but count_layer='raw' was specified")
        counts = adata.raw.X
        print("Using adata.raw.X (raw counts) for analysis")
    elif count_layer.startswith('layers:'):
        layer_name = count_layer.split(':', 1)[1]
        if layer_name not in adata.layers:
            raise ValueError(f"Layer '{layer_name}' not found in adata.layers. "
                           f"Available layers: {list(adata.layers.keys())}")
        counts = adata.layers[layer_name]
        print(f"Using adata.layers['{layer_name}'] for analysis")
    else:
        raise ValueError(f"Invalid count_layer value: '{count_layer}'. "
                       f"Use 'X', 'raw', or 'layers:NAME'")

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

    if sp.issparse(counts) and not keep_sparse:
        print("Converting sparse to dense (warning: high memory usage)")
        counts = counts.toarray()
        print("Memory after sparse->dense conversion:")
        check_memory_usage()

    print("✓ Cell data preparation complete")

    return cell_df, counts, gene_names


def prepare_cell_data_memory_efficient(
    cell_loadings: np.ndarray,
    adata_file: str,
    count_layer: str = 'X',
    cell_metadata_file: str = None,
    keep_sparse: bool = True
) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    """
    Memory-efficient cell data preparation with sparse matrix support

    Directly uses cell-level loadings instead of patch-to-cell mapping

    Parameters
    ----------
    cell_loadings : np.ndarray
        Cell-level BPTF loadings (n_cells × n_motifs)
    adata_file : str
        Path to AnnData h5ad file
    count_layer : str, default 'X'
        Which layer to use: 'X', 'raw', or 'layers:NAME'
    cell_metadata_file : str, optional
        Additional cell metadata CSV file
    keep_sparse : bool, default True
        Keep counts in sparse format

    Returns
    -------
    tuple
        (cell_df, counts, gene_names)
    """
    print("=== Memory-Efficient Data Loading ===")

    # Check initial memory
    print("Initial memory check:")
    check_memory_usage()

    # Load AnnData
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

    # Create cell_df with direct cell loadings
    print("Creating cell dataframe with loadings...")
    n_cells, n_motifs = cell_loadings.shape

    # Initialize cell_df with adata metadata
    cell_df = pd.DataFrame({
        'cell_id': adata.obs.index.astype(str),
        'cell_type': adata.obs['cell_type'] if 'cell_type' in adata.obs.columns else 'unknown'
    })

    # Add additional cell metadata if provided
    if cell_metadata_file is not None:
        print(f"Loading additional cell metadata from {cell_metadata_file}...")
        cell_meta = pd.read_csv(cell_metadata_file)
        if 'cell_type' in cell_meta.columns and 'cell_type' not in adata.obs.columns:
            cell_df['cell_type'] = cell_meta['cell_type'].values

    # Add cell-level motif loadings
    print("Adding cell-level motif loadings...")
    for k in range(n_motifs):
        cell_df[f'prog_{k}_loading'] = cell_loadings[:, k]

    print(f"✓ Added {n_motifs} motif loadings for {n_cells} cells")
    print("Memory after creating cell dataframe:")
    check_memory_usage()

    # Handle count matrix
    print(f"Processing count matrix from '{count_layer}'...")

    if count_layer == 'X':
        counts = adata.X
        print("Using adata.X for analysis")
    elif count_layer == 'raw':
        if adata.raw is None:
            raise ValueError("adata.raw is None but count_layer='raw' was specified")
        counts = adata.raw.X
        print("Using adata.raw.X (raw counts) for analysis")
    elif count_layer.startswith('layers:'):
        layer_name = count_layer.split(':', 1)[1]
        if layer_name not in adata.layers:
            raise ValueError(f"Layer '{layer_name}' not found in adata.layers. "
                           f"Available layers: {list(adata.layers.keys())}")
        counts = adata.layers[layer_name]
        print(f"Using adata.layers['{layer_name}'] for analysis")
    else:
        raise ValueError(f"Invalid count_layer value: '{count_layer}'. "
                       f"Use 'X', 'raw', or 'layers:NAME'")

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


def save_de_results(results: Dict, output_dir: str, mode: str = 'univariate'):
    """
    Save DE results to files

    Parameters
    ----------
    results : dict
        Dictionary mapping motif keys to DataFrames
    output_dir : str
        Output directory
    mode : str, default 'univariate'
        Analysis mode
    """
    import os

    os.makedirs(output_dir, exist_ok=True)

    if mode == 'univariate':
        for motif, df in results.items():
            filename = f"{motif}_de_results.csv"
            filepath = os.path.join(output_dir, filename)
            if not os.path.exists(filepath):
                df.to_csv(filepath, index=False)
                print(f"✓ Saved {filename}")
            else:
                print(f"✓ {filename} already exists, skipping")

        # Create long-format file
        long_data = []
        for motif, df in results.items():
            motif_num = int(motif.split('_')[1])
            df_copy = df.copy()
            df_copy['motif'] = motif_num
            long_data.append(df_copy)

        long_df = pd.concat(long_data, ignore_index=True)
        long_df = long_df[['gene', 'motif', 'logFC', 'pval', 'qval', 'significant']]
        long_df = long_df.sort_values(['motif', 'qval'])

        long_df.to_csv(os.path.join(output_dir, "univariate_de_long.csv"), index=False)
        print("✓ Saved univariate_de_long.csv")

    elif mode == 'multivariate':
        filename = "multivariate_de_results.csv"
        results.to_csv(os.path.join(output_dir, filename), index=False)
        print(f"✓ Saved {filename}")


def run_poisson_glm_analysis(
    cell_loadings: np.ndarray,
    adata: 'anndata.AnnData',
    lri_column_names: pd.Series,
    output_dir: Optional[str] = None,
    count_layer: str = 'X',
    splitter: str = '|',
    alpha: float = 0.05,
    random_state: int = 42,
    keep_sparse: bool = False
):
    """
    Run Poisson GLM differential expression analysis

    Parameters
    ----------
    cell_loadings : np.ndarray
        Cell-level BPTF loadings (n_cells × n_motifs)
    adata : anndata.AnnData
        Annotated data object with expression counts
    lri_column_names : pd.Series
        Series of LRI column names (for extracting LRI genes to exclude)
    output_dir : str, optional
        Output directory for results. If None, results are not saved to disk.
    count_layer : str, default 'X'
        Which layer to use: 'X', 'raw', or 'layers:NAME'
    splitter : str, default '|'
        Separator for LRI names
    alpha : float, default 0.05
        FDR significance threshold
    random_state : int, default 42
        Random seed
    keep_sparse : bool, default False
        Keep sparse format for counts

    Returns
    -------
    dict
        Dictionary of DE results

    Example
    -------
    >>> # After getting cell loadings from projection
    >>> results = al.run_poisson_glm_analysis(
    ...     cell_loadings=cell_loadings,
    ...     adata=adata,
    ...     lri_column_names=lri_columns['column_name'],
    ...     output_dir="results/glm"  # Optional
    ... )
    """
    print("="*60)
    print("BPTF DIFFERENTIAL EXPRESSION ANALYSIS")
    print("="*60)
    if output_dir:
        print(f"Output directory: {output_dir}")
    else:
        print("Output directory: None (results not saved)")
    print(f"Random seed: {random_state}")
    print("="*60)

    # Set random seed
    np.random.seed(random_state)

    # Validate inputs
    if cell_loadings.shape[0] != adata.shape[0]:
        raise ValueError(f"Shape mismatch: cell_loadings has {cell_loadings.shape[0]} cells, "
                        f"adata has {adata.shape[0]} cells")

    print(f"Cell loadings shape: {cell_loadings.shape}")
    print(f"AnnData shape: {adata.shape}")

    # Extract LRI genes
    lri_genes = extract_lri_genes(lri_column_names, splitter)

    # Get non-LRI genes
    if count_layer == 'raw' and adata.raw is not None:
        all_genes = list(adata.raw.var_names)
    else:
        all_genes = list(adata.var_names)

    non_lri_genes = [g for g in all_genes if g not in lri_genes]
    print(f"Analyzing {len(non_lri_genes)} non-LRI genes (excluding {len(lri_genes)} LRI genes)")

    # Prepare cell data
    cell_df, counts, gene_names = prepare_cell_data_from_adata(
        cell_loadings,
        adata,
        count_layer=count_layer,
        keep_sparse=keep_sparse
    )
    n_motifs = cell_loadings.shape[1]

    print(cell_df.head())

    # Run DE analysis
    results = run_univariate_de_sklearn_by_celltype(
        cell_df,
        counts,
        gene_names,
        non_lri_genes,
        n_motifs,
        output_dir=output_dir,
        alpha=alpha
    )

    # Save results if output_dir provided
    if output_dir:
        save_de_results(results, output_dir, mode='univariate')

    print("✓ Analysis completed!")

    return results


# ============================================================================
# GLM Results Analysis Functions
# (From glm/analysis.py)
# ============================================================================

def differential_expression(X, in_mask, out_mask=None, min_in_group_fraction=0.0001,
                           min_out_group_fraction=0.0001):
    """
    Perform differential expression analysis using Mann-Whitney U test

    Parameters
    ----------
    X : array-like
        Expression matrix (cells × genes)
    in_mask : array-like of bool
        Boolean mask for target cells
    out_mask : array-like of bool, optional
        Boolean mask for control cells (default: ~in_mask)
    min_in_group_fraction : float, default 0.0001
        Minimum expression fraction in target group
    min_out_group_fraction : float, default 0.0001
        Minimum expression fraction in control group

    Returns
    -------
    dict
        Dictionary with 'p_values', 'p_adj', and 'logfoldchanges'
    """
    if out_mask is None:
        out_mask = ~in_mask

    # Filter genes basically never expressed in control
    X_out = X[out_mask]
    genes_mask = np.ones(X.shape[1], dtype=bool)
    if min_out_group_fraction > 0:
        control_genes_mask = (X_out > 0.5).mean(axis=0) >= min_out_group_fraction
        genes_mask = genes_mask & control_genes_mask

    # Filter genes basically never expressed in target
    X_in = X[in_mask]
    if min_out_group_fraction > 0:
        target_genes_mask = (X_in > 0.5).mean() >= min_out_group_fraction
        genes_mask = genes_mask & target_genes_mask

    genes_mask = np.array(genes_mask).flatten()
    X_out = X_out[:, genes_mask]
    X_in = X_in[:, genes_mask]

    # Calculate log-fold changes
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



def load_exclusion_mask(csv_path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load previously saved exclusion mask from CSV

    Parameters
    ----------
    csv_path : str
        Path to exclusion_matrix.csv file

    Returns
    -------
    tuple
        (cell_types, genes, exclusion_mask)

    Example
    -------
    >>> ct, genes, mask = al.load_exclusion_mask("results/markers/exclusion_matrix.csv")
    """
    print(f"Loading exclusion mask from: {csv_path}")
    exclusion_df = pd.read_csv(csv_path, index_col=0)
    genes = np.array(exclusion_df.index)
    cell_types = np.array(exclusion_df.columns)
    exclusion_mask = exclusion_df.values.T  # Transpose to (n_cell_types, n_genes)

    print(f"✓ Loaded exclusion mask for {len(cell_types)} cell types, {len(genes)} genes")
    return cell_types, genes, exclusion_mask


def compute_exclusion_mask(
    adata,
    marker_lfc=1,
    marker_pvalue=1e-5,
    marker_subsample=50000,
    output_dir: Optional[str] = None
):
    """
    Identify marker genes for each cell type to exclude from plots

    Subsamples large groups to marker_subsample size for speed.
    This can be slow for large datasets - consider saving results.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    marker_lfc : float, default 1
        Log fold change threshold for markers
    marker_pvalue : float, default 1e-5
        P-value threshold for markers
    marker_subsample : int, default 50000
        Maximum cells to use for marker gene detection
    output_dir : str, optional
        Directory to save exclusion_matrix.csv. If None, results are not saved.

    Returns
    -------
    tuple
        (cell_types, genes, exclusion_mask)

    Example
    -------
    >>> # Without saving
    >>> ct, genes, mask = al.compute_exclusion_mask(adata)
    >>>
    >>> # With saving (recommended for reuse)
    >>> ct, genes, mask = al.compute_exclusion_mask(
    ...     adata,
    ...     output_dir="results/markers"
    ... )
    >>> # Saved to: results/markers/exclusion_matrix.csv
    """
    # Compute exclusion mask
    genes = np.array(adata.var_names)
    cell_types = np.unique(adata.obs['cell_type'])
    n_types = len(cell_types)
    exclusion_mask = np.zeros((n_types, len(genes)), dtype=bool)

    print(f"Computing marker genes for {len(cell_types)} cell types...")
    print("This may take a while for large datasets...")

    for cidx, cell_type in enumerate(cell_types):
        print(f"Processing {cell_type} ({cidx+1}/{len(cell_types)})...")

        # Create boolean masks
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

        # Identify markers
        marker_mask = (
            (deg['p_adj'] <= marker_pvalue) &
            (deg['logfoldchanges'] >= marker_lfc)
        )
        exclusion_mask[cidx, marker_mask] = True
        print(f"  Found {marker_mask.sum()} marker genes")

    print(f"✓ Computed exclusion mask for {len(cell_types)} cell types, {len(genes)} genes")

    # Save if output_dir provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        csv_path = os.path.join(output_dir, 'exclusion_matrix.csv')

        exclusion_df = pd.DataFrame(
            exclusion_mask.T,  # Transpose to (n_genes, n_cell_types)
            index=genes,
            columns=cell_types
        )
        exclusion_df.to_csv(csv_path)
        print(f"✓ Saved exclusion matrix to: {csv_path}")
    else:
        print("\nℹ️  Tip: Consider saving the exclusion mask for reuse:")
        print("   ct, genes, mask = al.compute_exclusion_mask(adata, output_dir='results/markers')")
        print("   Then load later: pd.read_csv('results/markers/exclusion_matrix.csv', index_col=0)")

    return cell_types, genes, exclusion_mask



def analyze_glm_results(
    adata: 'anndata.AnnData',
    de_results: Dict[str, pd.DataFrame],
    n_motifs: int,
    output_dir: Optional[str] = None,
    min_expression_frac: float = 0.02,
    marker_lfc: float = 1.0,
    marker_pvalue: float = 1e-5,
    marker_subsample: int = 50000,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 0.5,
    n_top_genes: int = 10,
    random_state: int = 42,
    exclusion_mask: Optional[np.ndarray] = None,
    cell_types: Optional[np.ndarray] = None,
    all_genes: Optional[np.ndarray] = None
):
    """
    Analyze and visualize GLM differential expression results

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data object with expression data
    de_results : dict
        Dictionary mapping motif-celltype keys to DE results DataFrames
        (output from run_poisson_glm_analysis)
    n_motifs : int
        Number of motifs analyzed
    output_dir : str, optional
        Output directory for plots. If None, plots are not saved.
    min_expression_frac : float, default 0.02
        Minimum expression fraction in cell type
    marker_lfc : float, default 1.0
        Log fold change threshold for marker genes
    marker_pvalue : float, default 1e-5
        P-value threshold for marker genes
    marker_subsample : int, default 50000
        Maximum cells for marker detection
    fdr_threshold : float, default 0.05
        FDR threshold for volcano plots
    lfc_threshold : float, default 0.5
        Log fold change threshold for volcano plots
    n_top_genes : int, default 10
        Number of top genes to label in plots
    random_state : int, default 42
        Random seed
    exclusion_mask : np.ndarray, optional
        Pre-computed marker gene exclusion mask. If None, will compute.
    cell_types : np.ndarray, optional
        Cell type names (required if exclusion_mask provided)
    all_genes : np.ndarray, optional
        Gene names (required if exclusion_mask provided)

    Returns
    -------
    tuple
        (cell_types, all_genes, exclusion_mask)

    Example
    -------
    >>> # After running GLM analysis
    >>> cell_types, genes, mask = al.analyze_glm_results(
    ...     adata=adata,
    ...     de_results=glm_results,
    ...     n_motifs=15,
    ...     output_dir="results/glm_viz"  # Optional
    ... )
    """
    import os

    print("="*60)
    print("GLM RESULTS ANALYSIS AND VISUALIZATION")
    print("="*60)
    if output_dir:
        print(f"Output directory: {output_dir}")
    else:
        print("Output directory: None (plots not saved)")
    print(f"Number of motifs: {n_motifs}")
    print("="*60)

    # Set random seed
    np.random.seed(random_state)

    print(f"Data shape: {adata.shape}")
    print(f"Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

    # Compute or validate marker exclusion mask
    if exclusion_mask is None:
        print("\nComputing marker gene exclusion mask...")
        cell_types, all_genes, exclusion_mask = compute_exclusion_mask(
            adata,
            marker_lfc=marker_lfc,
            marker_pvalue=marker_pvalue,
            marker_subsample=marker_subsample
        )
        print(f"Exclusion mask computed for {len(cell_types)} cell types and {len(all_genes)} genes")

        # Save marker genes if output_dir provided
        if output_dir:
            marker_dir = os.path.join(output_dir, "marker_genes")
            _save_marker_genes(marker_dir, cell_types, all_genes, exclusion_mask)
    else:
        print("\nUsing provided exclusion mask...")
        if cell_types is None or all_genes is None:
            raise ValueError("If exclusion_mask is provided, cell_types and all_genes must also be provided")
        print(f"Exclusion mask shape: {exclusion_mask.shape}")

    # Convert de_results dict to files structure for plotting functions
    # (The plotting functions currently expect files in a directory)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

        # Save DE results as CSV files temporarily for plotting
        results_temp_dir = os.path.join(output_dir, "_temp_de_results")
        os.makedirs(results_temp_dir, exist_ok=True)

        for key, df in de_results.items():
            csv_path = os.path.join(results_temp_dir, f"{key}_de_results.csv")
            df.to_csv(csv_path, index=False)

        # Generate volcano plots
        print("\nGenerating volcano plots...")
        glm_plots.generate_volcano_plots(
            results_temp_dir, output_dir, n_motifs, adata, cell_types, all_genes,
            exclusion_mask, min_expression_frac, fdr_threshold, lfc_threshold, n_top_genes
        )

        # Generate forest plots
        print("\nGenerating forest plots...")
        glm_plots.generate_forest_plots(
            results_temp_dir, output_dir, n_motifs, adata, cell_types, all_genes,
            exclusion_mask, min_expression_frac
        )

        # Clean up temporary files
        import shutil
        shutil.rmtree(results_temp_dir)

        print("\n" + "="*60)
        print("GLM RESULTS ANALYSIS COMPLETED SUCCESSFULLY!")
        print("="*60)
        print(f"Results saved to: {output_dir}")
    else:
        print("\nNo output directory specified - skipping plot generation")

    return cell_types, all_genes, exclusion_mask



def glm_volcano(
    adata: 'anndata.AnnData',
    de_results: Dict[str, pd.DataFrame],
    cell_types: np.ndarray,
    all_genes: np.ndarray,
    exclusion_mask: np.ndarray,
    motif_id: Optional[List[int]] = None,
    output_dir: Optional[str] = None,
    min_expression_frac: float = 0.02,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 0.2,
    n_top_genes: int = 30,
    figsize: tuple = (16, 16)
):
    """
    Generate volcano plots for GLM differential expression results

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data object
    de_results : dict
        Dictionary mapping motif-celltype keys to DE results DataFrames
    cell_types : np.ndarray
        Cell type names
    all_genes : np.ndarray
        Gene names
    exclusion_mask : np.ndarray
        Marker gene exclusion mask (n_cell_types × n_genes)
    motif_id : list of int, optional
        List of specific motif IDs to plot. If None, plots all motifs.
    output_dir : str, optional
        Directory to save plots. Required if motif_id is None.
    min_expression_frac : float, default 0.02
        Minimum expression fraction in cell type
    fdr_threshold : float, default 0.05
        FDR threshold for significance
    lfc_threshold : float, default 0.2
        Log fold change threshold
    n_top_genes : int, default 30
        Number of top genes to label
    figsize : tuple, default (16, 16)
        Figure size for multi-panel plots

    Returns
    -------
    matplotlib.figure.Figure or None
        Figure object if motif_id is specified, None if saving PDF

    Example
    -------
    >>> # Plot specific motifs
    >>> fig = al.glm_volcano(
    ...     adata, de_results, cell_types, genes, mask,
    ...     motif_id=[0, 1, 2]
    ... )
    >>> plt.show()
    >>>
    >>> # Generate PDF for all motifs
    >>> al.glm_volcano(
    ...     adata, de_results, cell_types, genes, mask,
    ...     output_dir="results/volcano"  # Required for all motifs
    ... )
    """
    import os
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    # Validation
    if motif_id is None and output_dir is None:
        raise ValueError("output_dir is required when motif_id is None (generating all motifs as PDF)")

    # Determine which motifs to plot
    if motif_id is not None:
        motifs_to_plot = motif_id
        save_as_pdf = False
    else:
        # Extract all unique motif IDs from de_results keys
        motifs_to_plot = sorted(set(
            int(key.split('_')[1]) for key in de_results.keys() if key.startswith('motif_')
        ))
        save_as_pdf = True

    print(f"Generating volcano plots for motifs: {motifs_to_plot}")

    if save_as_pdf:
        # Multi-page PDF for all motifs
        os.makedirs(output_dir, exist_ok=True)
        pdf_path = os.path.join(output_dir, "volcano_plots_filtered.pdf")

        with PdfPages(pdf_path) as pdf:
            for mid in motifs_to_plot:
                print(f"  Processing motif {mid}...")
                fig = _create_volcano_figure_for_motif(
                    mid, de_results, adata, cell_types, all_genes, exclusion_mask,
                    min_expression_frac, fdr_threshold, lfc_threshold, n_top_genes, figsize
                )
                pdf.savefig(fig)
                plt.close(fig)

        print(f"Volcano plots saved to: {pdf_path}")
        return None

    else:
        # Create figures for specified motifs
        if len(motifs_to_plot) == 1:
            # Single motif - create one figure
            fig = _create_volcano_figure_for_motif(
                motifs_to_plot[0], de_results, adata, cell_types, all_genes, exclusion_mask,
                min_expression_frac, fdr_threshold, lfc_threshold, n_top_genes, figsize
            )

            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                save_path = os.path.join(output_dir, f"volcano_motif_{motifs_to_plot[0]}.pdf")
                fig.savefig(save_path, bbox_inches='tight')
                print(f"Volcano plot saved to: {save_path}")

            return fig
        else:
            # Multiple motifs - create separate figures and save as multi-page PDF
            figs = []
            for mid in motifs_to_plot:
                print(f"  Creating volcano plot for motif {mid}...")
                fig = _create_volcano_figure_for_motif(
                    mid, de_results, adata, cell_types, all_genes, exclusion_mask,
                    min_expression_frac, fdr_threshold, lfc_threshold, n_top_genes, figsize
                )
                figs.append(fig)

            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                pdf_path = os.path.join(output_dir, f"volcano_motifs_{'_'.join(map(str, motifs_to_plot))}.pdf")

                with PdfPages(pdf_path) as pdf:
                    for fig in figs:
                        pdf.savefig(fig)
                        plt.close(fig)

                print(f"Volcano plots saved to: {pdf_path}")
                return None
            else:
                # Return list of figures
                return figs


def _create_volcano_figure_for_motif(
    motif_id, de_results, adata, cell_types, all_genes, exclusion_mask,
    min_expression_frac, fdr_threshold, lfc_threshold, n_top_genes, figsize
):
    """Helper to create volcano figure for a single motif"""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 4, figsize=figsize)
    axes = axes.flatten()

    # Find all results for this motif
    motif_files = {k: v for k, v in de_results.items() if k.startswith(f'motif_{motif_id}_celltype_')}

    for ax_idx, (key, df) in enumerate(sorted(motif_files.items())[:16]):
        ct = key.split(f'motif_{motif_id}_celltype_')[1].replace('_de_results', '')

        try:
            cidx = list(cell_types).index(ct)
        except ValueError:
            print(f"Warning: Cell type {ct} not found")
            continue

        # Filter genes
        df = _filter_genes_for_volcano(df, ct, adata, all_genes, exclusion_mask, cidx, min_expression_frac)

        # Compute -log10(q) with jitter
        df['neg_log10_q'] = -np.log10(df['qval'].clip(1e-300))
        m = df['neg_log10_q'] >= 300
        if m.any():
            df.loc[m, 'neg_log10_q'] = 300 + np.random.normal(0, 15, m.sum())

        # Plot
        glm_plots.volcano_plot(
            df, 'logFC', 'neg_log10_q', label_col='gene',
            fdr=fdr_threshold, x_threshold=lfc_threshold,
            marker='o', n_top=n_top_genes, fontsize=5, ax=axes[ax_idx]
        )
        axes[ax_idx].set_title(ct, fontsize=10, pad=5)

    # Remove unused subplots
    for j in range(len(motif_files), 16):
        fig.delaxes(axes[j])

    plt.suptitle(f"Motif {motif_id} Volcano Plots (marker genes filtered)", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    return fig


def _filter_genes_for_volcano(df, ct, adata, all_genes, exclusion_mask, cidx, min_expression_frac):
    """Filter genes by expression and marker status"""
    glm_genes_set = set(df['gene'])

    # Expression filter
    subX = adata[adata.obs['cell_type'] == ct].X
    if sp.issparse(subX):
        expr_frac = np.asarray((subX > 0).sum(axis=0)).ravel() / subX.shape[0]
    else:
        expr_frac = np.sum(subX > 0, axis=0) / subX.shape[0]

    gene_to_idx = {gene: i for i, gene in enumerate(all_genes)}
    genes_to_keep = []
    for gene in glm_genes_set:
        if gene in gene_to_idx:
            gene_idx = gene_to_idx[gene]
            if expr_frac[gene_idx] >= min_expression_frac:
                genes_to_keep.append(gene)

    df = df[df['gene'].isin(genes_to_keep)].copy()

    # Marker filter
    other = [i for i in range(len(exclusion_mask)) if i != cidx]
    drop_mask = exclusion_mask[other, :].any(axis=0)
    marker_genes_to_drop = set(all_genes[drop_mask]) & glm_genes_set
    df = df[~df['gene'].isin(marker_genes_to_drop)]

    return df


def glm_forest(
    adata: 'anndata.AnnData',
    de_results: Dict[str, pd.DataFrame],
    cell_types: np.ndarray,
    all_genes: np.ndarray,
    exclusion_mask: np.ndarray,
    motif_id: Optional[List[int]] = None,
    output_dir: Optional[str] = None,
    min_expression_frac: float = 0.02,
    n_top: int = 30,
    figsize: tuple = (16, 24)
):
    """
    Generate forest plots for GLM differential expression results

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data object
    de_results : dict
        Dictionary mapping motif-celltype keys to DE results DataFrames
    cell_types : np.ndarray
        Cell type names
    all_genes : np.ndarray
        Gene names
    exclusion_mask : np.ndarray
        Marker gene exclusion mask (n_cell_types × n_genes)
    motif_id : list of int, optional
        List of specific motif IDs to plot. If None, plots all motifs.
    output_dir : str, optional
        Directory to save plots. Required if motif_id is None.
    min_expression_frac : float, default 0.02
        Minimum expression fraction in cell type
    n_top : int, default 30
        Number of top genes to show per plot
    figsize : tuple, default (16, 24)
        Figure size for multi-panel plots

    Returns
    -------
    matplotlib.figure.Figure or None
        Figure object if motif_id is specified, None if saving PDF

    Example
    -------
    >>> # Plot specific motifs
    >>> fig = al.glm_forest(
    ...     adata, de_results, cell_types, genes, mask,
    ...     motif_id=[0, 1]
    ... )
    >>> plt.show()
    >>>
    >>> # Generate PDF for all motifs
    >>> al.glm_forest(
    ...     adata, de_results, cell_types, genes, mask,
    ...     output_dir="results/forest"  # Required
    ... )
    """
    import os
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    # Validation
    if motif_id is None and output_dir is None:
        raise ValueError("output_dir is required when motif_id is None (generating all motifs as PDF)")

    # Determine which motifs to plot
    if motif_id is not None:
        motifs_to_plot = motif_id
        save_as_pdf = False
    else:
        # Extract all unique motif IDs from de_results keys
        motifs_to_plot = sorted(set(
            int(key.split('_')[1]) for key in de_results.keys() if key.startswith('motif_')
        ))
        save_as_pdf = True

    print(f"Generating forest plots for motifs: {motifs_to_plot}")

    if save_as_pdf:
        # Multi-page PDF for all motifs
        os.makedirs(output_dir, exist_ok=True)
        pdf_path = os.path.join(output_dir, "forest_plots_filtered.pdf")

        with PdfPages(pdf_path) as pdf:
            for mid in motifs_to_plot:
                print(f"  Processing motif {mid}...")
                fig = _create_forest_figure_for_motif(
                    mid, de_results, adata, cell_types, all_genes, exclusion_mask,
                    min_expression_frac, n_top, figsize
                )
                pdf.savefig(fig)
                plt.close(fig)

        print(f"Forest plots saved to: {pdf_path}")
        return None

    else:
        # Create figures for specified motifs
        if len(motifs_to_plot) == 1:
            # Single motif - create one figure
            fig = _create_forest_figure_for_motif(
                motifs_to_plot[0], de_results, adata, cell_types, all_genes, exclusion_mask,
                min_expression_frac, n_top, figsize
            )

            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                save_path = os.path.join(output_dir, f"forest_motif_{motifs_to_plot[0]}.pdf")
                fig.savefig(save_path, bbox_inches='tight')
                print(f"Forest plot saved to: {save_path}")

            return fig
        else:
            # Multiple motifs - create separate figures and save as multi-page PDF
            figs = []
            for mid in motifs_to_plot:
                print(f"  Creating forest plot for motif {mid}...")
                fig = _create_forest_figure_for_motif(
                    mid, de_results, adata, cell_types, all_genes, exclusion_mask,
                    min_expression_frac, n_top, figsize
                )
                figs.append(fig)

            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                pdf_path = os.path.join(output_dir, f"forest_motifs_{'_'.join(map(str, motifs_to_plot))}.pdf")

                with PdfPages(pdf_path) as pdf:
                    for fig in figs:
                        pdf.savefig(fig)
                        plt.close(fig)

                print(f"Forest plots saved to: {pdf_path}")
                return None
            else:
                # Return list of figures
                return figs


def _create_forest_figure_for_motif(
    motif_id, de_results, adata, cell_types, all_genes, exclusion_mask,
    min_expression_frac, n_top, figsize
):
    """Helper to create forest figure for a single motif"""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 4, figsize=figsize)
    axes = axes.flatten()

    # Find all results for this motif
    motif_files = {k: v for k, v in de_results.items() if k.startswith(f'motif_{motif_id}_celltype_')}

    for ax_idx, (key, df) in enumerate(sorted(motif_files.items())[:16]):
        ct = key.split(f'motif_{motif_id}_celltype_')[1].replace('_de_results', '')

        try:
            ct_idx = list(cell_types).index(ct)
        except ValueError:
            continue

        # Filter genes
        df = _filter_genes_for_volcano(df, ct, adata, all_genes, exclusion_mask, ct_idx, min_expression_frac)

        glm_plots.forest_plot(df, ax=axes[ax_idx], n_top=n_top)
        axes[ax_idx].set_title(ct, fontsize=9)

    # Remove unused subplots
    for j in range(len(motif_files), 16):
        fig.delaxes(axes[j])

    plt.suptitle(f"Motif {motif_id} – Forest plots (marker genes filtered)", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.965])

    return fig


def _save_marker_genes(marker_dir, cell_types, all_genes, exclusion_mask):
    """Save marker gene information"""
    print("\nSaving marker genes information...")
    os.makedirs(marker_dir, exist_ok=True)

    marker_summary = []
    for i, cell_type in enumerate(cell_types):
        marker_genes = all_genes[exclusion_mask[i, :]]

        # Save individual cell type markers
        marker_df = pd.DataFrame({
            'gene': marker_genes,
            'cell_type': cell_type
        })
        marker_file = os.path.join(marker_dir,
                                   f"{cell_type.replace('/', '_').replace(' ', '_')}_markers.csv")
        marker_df.to_csv(marker_file, index=False)

        marker_summary.append({
            'cell_type': cell_type,
            'n_marker_genes': len(marker_genes),
            'marker_percentage': len(marker_genes) / len(all_genes) * 100
        })
        print(f"  {cell_type}: {len(marker_genes)} marker genes ({len(marker_genes)/len(all_genes)*100:.1f}%)")

    # Save summary
    summary_df = pd.DataFrame(marker_summary)
    summary_df.to_csv(os.path.join(marker_dir, "marker_genes_summary.csv"), index=False)

    # Save exclusion matrix
    exclusion_df = pd.DataFrame(
        exclusion_mask.T,
        index=all_genes,
        columns=cell_types
    )
    exclusion_df.to_csv(os.path.join(marker_dir, "exclusion_matrix.csv"))

    print(f"Marker genes information saved to: {marker_dir}")



