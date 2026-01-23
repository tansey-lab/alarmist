"""
Data loading functions for alarmist results
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata
import os
from typing import Dict, Optional


def load_patch_lri_results(input_dir: str,
                           sparse_matrix_name: str = 'patch_lri_matrix.npz',
                           column_df_name: str = 'patch_lri_columns.csv') -> Dict:
    """
    Load patch-based LRI analysis results

    Parameters
    ----------
    input_dir : str
        Directory containing results
    sparse_matrix_name : str
        Name of sparse matrix file
    column_df_name : str
        Name of column names file

    Returns
    -------
    dict
        Dictionary with:
        - patch_lri_matrix: sparse matrix
        - column_names: list
        - parameters: DataFrame
        - sample_info: DataFrame (only if multi-sample, i.e., sample_info.csv exists)
    """
    print(f"Loading patch-LRI results from: {input_dir}")

    # Load sparse matrix
    mat_path = os.path.join(input_dir, sparse_matrix_name)
    patch_lri_matrix = sp.load_npz(mat_path)

    # Load column names
    cols_path = os.path.join(input_dir, column_df_name)
    cols_df = pd.read_csv(cols_path)
    column_names = cols_df['column_name'].tolist()

    # Load parameters
    params_file = os.path.join(input_dir, 'analysis_parameters.csv')
    params_df = pd.read_csv(params_file)

    print(f"Loaded matrix shape: {patch_lri_matrix.shape}")
    print(f"Matrix sparsity: {params_df[params_df['parameter'] == 'matrix_sparsity']['value'].iloc[0]}")

    results = {
        'patch_lri_matrix': patch_lri_matrix,
        'column_names': column_names,
        'parameters': params_df
    }

    # Auto-detect multi-sample: load sample_info if exists
    sample_info_file = os.path.join(input_dir, 'sample_info.csv')
    if os.path.exists(sample_info_file):
        sample_info_df = pd.read_csv(sample_info_file)
        # Convert to dict format matching run_patchify return
        sample_info = {}
        for _, row in sample_info_df.iterrows():
            sample_id = row['sample_id']
            sample_info[sample_id] = {
                'n_cells': row['n_cells'],
                'n_patches': row['n_patches'],
                'global_patch_idx_start': row['global_patch_idx_start'],
                'global_patch_idx_end': row['global_patch_idx_end']
            }
        results['sample_info'] = sample_info
        print(f"Multi-sample detected: {len(sample_info)} samples")

    return results

def load_cell_lri_results(output_dir: str) -> Dict:
    """
    Load previously saved cell-LRI analysis results.

    Parameters
    ----------
    output_dir : str
        Directory containing saved results

    Returns
    -------
    results : dict
        Dictionary containing:
        - cell_lri_matrix: sparse matrix (n_cells × n_lris)
        - column_names: list of LRI column names
        - cell_metadata_df: cell metadata
        - parameters: analysis parameters
    """
    print(f"Loading cell-LRI results from: {output_dir}")

    # Load sparse matrix
    matrix_file = os.path.join(output_dir, 'cell_lri_matrix.npz')
    cell_lri_matrix = sp.load_npz(matrix_file)

    # Load column names
    columns_file = os.path.join(output_dir, 'cell_lri_columns.csv')
    column_names = pd.read_csv(columns_file)['column_name'].tolist()

    # Load metadata
    metadata_file = os.path.join(output_dir, 'cell_metadata.csv')
    cell_metadata_df = pd.read_csv(metadata_file)

    # Load parameters
    params_file = os.path.join(output_dir, 'analysis_parameters.csv')
    params_df = pd.read_csv(params_file)

    # Load sample info if available
    sample_info_file = os.path.join(output_dir, 'sample_info.csv')
    if os.path.exists(sample_info_file):
        sample_info_df = pd.read_csv(sample_info_file)

    print(f"Loaded matrix shape: {cell_lri_matrix.shape}")
    print(f"Matrix sparsity: {params_df[params_df['parameter'] == 'matrix_sparsity']['value'].iloc[0]}")

    return {
        'cell_lri_matrix': cell_lri_matrix,
        'column_names': column_names,
        'cell_metadata_df': cell_metadata_df,
        'parameters': params_df,
        'sample_info': sample_info_df if os.path.exists(sample_info_file) else None
    }


def load_bptf_results(results_dir: str, load_rescaled: bool = False) -> Dict:
    """
    Load BPTF factorization results

    Parameters
    ----------
    results_dir : str
        Directory containing BPTF results
    load_rescaled : bool, default False
        Whether to also load rescaled matrices (patch_loadings_rescaled, lri_factors_rescaled)

    Returns
    -------
    dict
        Dictionary with:
        - patch_loadings: np.ndarray
        - lri_factors: np.ndarray
        - lri_motifs: pd.DataFrame
        - patch_loadings_rescaled: np.ndarray (only if load_rescaled=True)
        - lri_factors_rescaled: np.ndarray (only if load_rescaled=True)
    """
    print(f"Loading BPTF results from: {results_dir}")

    patch_loadings = np.load(os.path.join(results_dir, 'patch_loadings.npy'))
    lri_factors = np.load(os.path.join(results_dir, 'lri_factors.npy'))
    lri_motifs = pd.read_csv(os.path.join(results_dir, 'lri_motifs.csv'))

    print(f"Loaded results:")
    print(f"  - Patch loadings: {patch_loadings.shape}")
    print(f"  - LRI factors: {lri_factors.shape}")

    results = {
        'patch_loadings': patch_loadings,
        'lri_factors': lri_factors,
        'lri_motifs': lri_motifs
    }

    if load_rescaled:
        patch_loadings_rescaled = np.load(os.path.join(results_dir, 'patch_loadings_rescaled.npy'))
        lri_factors_rescaled = np.load(os.path.join(results_dir, 'lri_factors_rescaled.npy'))
        results['patch_loadings_rescaled'] = patch_loadings_rescaled
        results['lri_factors_rescaled'] = lri_factors_rescaled
        print(f"  - Patch loadings rescaled: {patch_loadings_rescaled.shape}")
        print(f"  - LRI factors rescaled: {lri_factors_rescaled.shape}")

    return results
