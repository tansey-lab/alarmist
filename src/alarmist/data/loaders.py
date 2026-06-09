"""
Data loading functions for alarmist results
"""

import logging
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp

from alarmist.constants import (
    COLUMN_NAME_AVG_NEIGHBORHOOD_SIZE,
    COLUMN_NAME_COLUMN_NAME,
    COLUMN_NAME_GLOBAL_CELL_IDX_END,
    COLUMN_NAME_GLOBAL_CELL_IDX_START,
    COLUMN_NAME_GLOBAL_PATCH_IDX_END,
    COLUMN_NAME_GLOBAL_PATCH_IDX_START,
    COLUMN_NAME_N_CELLS,
    COLUMN_NAME_N_PATCHES,
    COLUMN_NAME_PARAMETER,
    COLUMN_NAME_SAMPLE_ID,
    COLUMN_NAME_VALUE,
)

logger = logging.getLogger(__name__)


def load_patch_lri_results(
    input_dir: str,
    sparse_matrix_name: str = "patch_lri_matrix.npz",
    column_df_name: str = "patch_lri_columns.csv",
) -> dict:
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
    logger.debug(f"Loading patch-LRI results from: {input_dir}")

    # Load sparse matrix
    mat_path = os.path.join(input_dir, sparse_matrix_name)
    patch_lri_matrix = sp.load_npz(mat_path)

    # Load column names
    cols_path = os.path.join(input_dir, column_df_name)
    cols_df = pd.read_csv(cols_path)
    column_names = cols_df[COLUMN_NAME_COLUMN_NAME].tolist()

    # Load parameters
    params_file = os.path.join(input_dir, "analysis_parameters.csv")
    params_df = pd.read_csv(params_file)

    logger.debug(f"Loaded matrix shape: {patch_lri_matrix.shape}")
    logger.debug(
        f"Matrix sparsity: {params_df[params_df[COLUMN_NAME_PARAMETER] == 'matrix_sparsity'][COLUMN_NAME_VALUE].iloc[0]}"
    )

    results = {
        "patch_lri_matrix": patch_lri_matrix,
        "column_names": column_names,
        "parameters": params_df,
    }

    # Auto-detect multi-sample: load sample_info if exists
    sample_info_file = os.path.join(input_dir, "sample_info.csv")
    if os.path.exists(sample_info_file):
        sample_info_df = pd.read_csv(sample_info_file)
        # Convert to dict format matching run_patchify return
        sample_info = {}
        for _, row in sample_info_df.iterrows():
            sample_id = row[COLUMN_NAME_SAMPLE_ID]
            sample_info[sample_id] = {
                COLUMN_NAME_N_CELLS: row[COLUMN_NAME_N_CELLS],
                COLUMN_NAME_N_PATCHES: row[COLUMN_NAME_N_PATCHES],
                COLUMN_NAME_GLOBAL_PATCH_IDX_START: row[
                    COLUMN_NAME_GLOBAL_PATCH_IDX_START
                ],
                COLUMN_NAME_GLOBAL_PATCH_IDX_END: row[COLUMN_NAME_GLOBAL_PATCH_IDX_END],
            }
        results["sample_info"] = sample_info
        logger.debug(f"Multi-sample detected: {len(sample_info)} samples")

    return results


def load_cell_lri_results(output_dir: str) -> dict:
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
        - parameters: analysis parameters
        - sample_info: dict (only if multi-sample, i.e., sample_info.csv exists)
    """
    logger.debug(f"Loading cell-LRI results from: {output_dir}")

    # Load sparse matrix
    matrix_file = os.path.join(output_dir, "cell_lri_matrix.npz")
    cell_lri_matrix = sp.load_npz(matrix_file)

    # Load column names
    columns_file = os.path.join(output_dir, "cell_lri_columns.csv")
    column_names = pd.read_csv(columns_file)[COLUMN_NAME_COLUMN_NAME].tolist()

    # Load parameters
    params_file = os.path.join(output_dir, "analysis_parameters.csv")
    params_df = pd.read_csv(params_file)

    logger.debug(f"Loaded matrix shape: {cell_lri_matrix.shape}")
    logger.debug(
        f"Matrix sparsity: {params_df[params_df[COLUMN_NAME_PARAMETER] == 'matrix_sparsity'][COLUMN_NAME_VALUE].iloc[0]}"
    )

    results = {
        "cell_lri_matrix": cell_lri_matrix,
        "column_names": column_names,
        "parameters": params_df,
    }

    # Auto-detect multi-sample: load sample_info if exists
    sample_info_file = os.path.join(output_dir, "sample_info.csv")
    if os.path.exists(sample_info_file):
        sample_info_df = pd.read_csv(sample_info_file)
        # Convert to dict format matching run_neighborhood return
        sample_info = {}
        for _, row in sample_info_df.iterrows():
            sample_id = row[COLUMN_NAME_SAMPLE_ID]
            sample_info[sample_id] = {
                COLUMN_NAME_N_CELLS: row[COLUMN_NAME_N_CELLS],
                COLUMN_NAME_GLOBAL_CELL_IDX_START: row[
                    COLUMN_NAME_GLOBAL_CELL_IDX_START
                ],
                COLUMN_NAME_GLOBAL_CELL_IDX_END: row[COLUMN_NAME_GLOBAL_CELL_IDX_END],
                COLUMN_NAME_AVG_NEIGHBORHOOD_SIZE: row[
                    COLUMN_NAME_AVG_NEIGHBORHOOD_SIZE
                ],
            }
        results["sample_info"] = sample_info
        logger.debug(f"Multi-sample detected: {len(sample_info)} samples")

    return results


def load_bptf_results(results_dir: str, load_rescaled: bool = False) -> dict:
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
    logger.debug(f"Loading BPTF results from: {results_dir}")

    patch_loadings = np.load(os.path.join(results_dir, "patch_loadings.npy"))
    lri_factors = np.load(os.path.join(results_dir, "lri_factors.npy"))
    lri_motifs = pd.read_csv(os.path.join(results_dir, "lri_motifs.csv"))

    logger.debug("Loaded results:")
    logger.debug(f"  - Patch loadings: {patch_loadings.shape}")
    logger.debug(f"  - LRI factors: {lri_factors.shape}")

    results = {
        "patch_loadings": patch_loadings,
        "lri_factors": lri_factors,
        "lri_motifs": lri_motifs,
    }

    if load_rescaled:
        patch_loadings_rescaled = np.load(
            os.path.join(results_dir, "patch_loadings_rescaled.npy")
        )
        lri_factors_rescaled = np.load(
            os.path.join(results_dir, "lri_factors_rescaled.npy")
        )
        results["patch_loadings_rescaled"] = patch_loadings_rescaled
        results["lri_factors_rescaled"] = lri_factors_rescaled
        logger.debug(f"  - Patch loadings rescaled: {patch_loadings_rescaled.shape}")
        logger.debug(f"  - LRI factors rescaled: {lri_factors_rescaled.shape}")

    return results
