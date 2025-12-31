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
                           column_df_name: str = 'patch_lri_columns.csv',
                           meta_df_name: str = 'cell_patch_correspondence.csv',
                           neighborhood: bool = False,
                           single_cell: bool = False) -> Dict:
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
    meta_df_name : str
        Name of metadata file
    neighborhood : bool
        Whether this is neighborhood-based analysis
    single_cell : bool
        Whether this is single-cell analysis

    Returns
    -------
    dict
        Dictionary with:
        - patch_lri_matrix: sparse matrix
        - column_names: list
        - cell_patch_df: DataFrame
        - patch_tma_df: DataFrame (if not single_cell)
    """
    print(f"Loading patch-LRI results from: {input_dir}")

    # Load sparse matrix
    mat_path = os.path.join(input_dir, sparse_matrix_name)
    patch_lri_matrix = sp.load_npz(mat_path)

    # Load column names
    cols_path = os.path.join(input_dir, column_df_name)
    cols_df = pd.read_csv(cols_path)
    column_names = cols_df['column_name'].tolist()

    # Load metadata dataframes
    if not single_cell:
        cell_patch_file = os.path.join(input_dir, 'cell_patch_correspondence.csv')
        cell_patch_df = pd.read_csv(cell_patch_file)

    # Load parameters
    params_file = os.path.join(input_dir, 'analysis_parameters.csv')
    params_df = pd.read_csv(params_file)

    print(f"Loaded matrix shape: {patch_lri_matrix.shape}")
    print(f"Matrix sparsity: {params_df[params_df['parameter'] == 'matrix_sparsity']['value'].iloc[0]}")

    results = {
        'patch_lri_matrix': patch_lri_matrix,
        'column_names': column_names,
    }

    if neighborhood:
        results['cell_patch_df'] = cell_patch_df
        results['parameters'] = params_df
    elif not single_cell:
        if meta_df_name == 'cell_patch_correspondence.csv':
            patch_tma_file = os.path.join(input_dir, 'patch_tma_correspondence.csv')
        else:
            patch_tma_file = os.path.join(input_dir, meta_df_name)
        patch_tma_df = pd.read_csv(patch_tma_file)
        results['patch_tma_df'] = patch_tma_df
        results['cell_patch_df'] = cell_patch_df
        results['parameters'] = params_df
    else:
        patch_tma_file = os.path.join(input_dir, meta_df_name)
        patch_tma_df = pd.read_csv(patch_tma_file)
        results['patch_tma_df'] = patch_tma_df
        results['parameters'] = params_df

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

    print(f"Loaded matrix shape: {cell_lri_matrix.shape}")
    print(f"Matrix sparsity: {params_df[params_df['parameter'] == 'matrix_sparsity']['value'].iloc[0]}")

    return {
        'cell_lri_matrix': cell_lri_matrix,
        'column_names': column_names,
        'cell_metadata_df': cell_metadata_df,
        'parameters': params_df
    }


def load_bptf_results(results_dir: str) -> Dict:
    """
    Load BPTF factorization results

    Parameters
    ----------
    results_dir : str
        Directory containing BPTF results

    Returns
    -------
    dict
        Dictionary with patch_loadings, lri_factors, patch_motifs, lri_motifs
    """
    print(f"Loading BPTF results from: {results_dir}")

    patch_loadings = np.load(os.path.join(results_dir, 'patch_loadings.npy'))
    lri_factors = np.load(os.path.join(results_dir, 'lri_factors.npy'))
    patch_motifs = pd.read_csv(os.path.join(results_dir, 'patch_motifs.csv'))
    lri_motifs = pd.read_csv(os.path.join(results_dir, 'lri_motifs.csv'))

    print(f"Loaded results:")
    print(f"  - Patch loadings: {patch_loadings.shape}")
    print(f"  - LRI factors: {lri_factors.shape}")

    return {
        'patch_loadings': patch_loadings,
        'lri_factors': lri_factors,
        'patch_motifs': patch_motifs,
        'lri_motifs': lri_motifs
    }
