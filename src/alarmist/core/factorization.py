"""
BPTF (Bayesian Poisson Tensor Factorization) wrapper functions
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
import warnings
import os
from typing import Optional, Tuple, Dict

try:
    import sparse
    from bptf import BPTF, save_bptf, load_bptf
    BPTF_AVAILABLE = True
except ImportError:
    BPTF_AVAILABLE = False
    warnings.warn("BPTF not available. Install from: https://github.com/aschein/bptf")


def run_bptf(mat: sp.spmatrix,
             n_components: int = 15,
             max_iter: int = 10000,
             verbose: bool = True,
             random_state: int = 0) -> object:
    """
    Run BPTF matrix factorization

    Parameters
    ----------
    mat : scipy.sparse matrix
        Input sparse matrix
    n_components : int, default 15
        Number of latent factors/motifs
    max_iter : int, default 10000
        Maximum iterations
    verbose : bool, default True
        Print progress
    random_state : int, default 0
        Random seed for reproducibility

    Returns
    -------
    model : BPTF model object
    """
    if not BPTF_AVAILABLE:
        raise ImportError("BPTF not available. Install: pip install git+https://github.com/aschein/bptf.git")

    # Convert to sparse.COO format
    if isinstance(mat, sp.spmatrix):
        data = sparse.COO(mat)
    else:
        data = mat

    np.random.seed(random_state)

    model = BPTF(data_shape=data.shape, n_components=n_components)
    model.fit(data, max_iter=max_iter, verbose=verbose)

    return model


def extract_factors(model: object) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract factor matrices from BPTF model

    Returns
    -------
    patch_loadings : np.ndarray
        Patch loadings (patches × motifs)
    lri_factors : np.ndarray
        LRI factors (motifs × LRIs)
    """
    G_DK_M = model.G_DK_M
    patch_loadings = G_DK_M[0]
    lri_factors = G_DK_M[1].T  # Transpose to (motifs × LRIs)

    return patch_loadings, lri_factors


def get_top_motifs(patch_loadings: np.ndarray, top_k: int = 10) -> Dict:
    """Get top motifs by total activity"""
    motif_activities = patch_loadings.sum(axis=0)
    top_indices = np.argsort(motif_activities)[::-1][:top_k]

    return {
        'motif_indices': top_indices.tolist(),
        'activities': motif_activities[top_indices].tolist(),
        'activity_fractions': (motif_activities[top_indices] / motif_activities.sum()).tolist()
    }


def project_cell_loadings(
    model: object,
    cell_lri_matrix: sp.spmatrix,
    model_lri_columns: Optional[np.ndarray] = None,
    cell_lri_columns: Optional[np.ndarray] = None,
    max_iter: int = 200,
    chunk_size: int = 50000,
    verbose: bool = True
) -> np.ndarray:
    """
    Project cell-LRI matrix to cell loadings using a fitted BPTF model.

    This fixes the LRI factors (mode=1) from the patch-level model and only
    updates the cell loadings (mode=0), allowing us to get cell-level motif
    activities from the patch-level factorization.

    Parameters
    ----------
    model : BPTF model object
        Fitted patch-level BPTF model with learned LRI factors
    cell_lri_matrix : scipy.sparse matrix
        Cell-LRI interaction matrix (n_cells × n_lris)
    model_lri_columns : np.ndarray, optional
        LRI column names from the model (for alignment)
    cell_lri_columns : np.ndarray, optional
        LRI column names from the cell matrix (for alignment)
    max_iter : int, default 200
        Maximum iterations for projection
    chunk_size : int, default 50000
        Number of cells to process per chunk (for memory efficiency)
    verbose : bool, default True
        Print progress

    Returns
    -------
    cell_loadings : np.ndarray
        Cell loadings matrix (n_cells × n_motifs)

    Notes
    -----
    - If model_lri_columns and cell_lri_columns are provided, will align
      the cell matrix columns to match the model's LRI order
    - Uses chunked processing to handle large cell matrices
    - Only mode=0 (cell loadings) is updated; mode=1 (LRI factors) stays fixed
    """
    if not BPTF_AVAILABLE:
        raise ImportError("BPTF not available. Install: pip install git+https://github.com/aschein/bptf.git")

    import gc

    # Align LRI columns if names provided
    if model_lri_columns is not None and cell_lri_columns is not None:
        # Quick check: if columns are already aligned, skip alignment
        if (len(model_lri_columns) == len(cell_lri_columns) and
            all(a == b for a, b in zip(model_lri_columns, cell_lri_columns))):
            if verbose:
                print(f"✓ Columns already aligned ({len(model_lri_columns)} LRIs)")
        else:
            # Columns need alignment
            if verbose:
                print("Aligning cell-LRI matrix columns to model...")

            # Find intersection and preserve model's order
            from collections import defaultdict, deque

            a = cell_lri_columns
            b = model_lri_columns

            # Which b elements appear in a (preserving b's order and duplicates)
            mask = np.isin(b, a)
            idx_in_model = np.flatnonzero(mask)

            # Build queue for one-to-one matching
            pos = defaultdict(deque)
            for i, v in enumerate(a):
                pos[v].append(i)

            # Match each b[mask] value to next unused position in a
            idx_in_mat = np.array([pos[v].popleft() if pos[v] else -1 for v in b[mask]], dtype=int)

            if (idx_in_mat == -1).any():
                raise ValueError("Some model LRIs not found in cell matrix. "
                               "Use run_neighborhood(..., required_columns=patch_columns) to ensure alignment.")

            # Reorder cell matrix columns to match model
            cell_lri_matrix = cell_lri_matrix[:, idx_in_mat]

            if cell_lri_matrix.shape[1] != len(model_lri_columns):
                raise ValueError(f"Column mismatch: cell matrix has {cell_lri_matrix.shape[1]} "
                               f"columns after alignment, model expects {len(model_lri_columns)}")

            if verbose:
                print(f"✓ Aligned to {cell_lri_matrix.shape[1]} LRI columns")

    # Convert to CSR if needed
    if not sp.isspmatrix_csr(cell_lri_matrix):
        if verbose:
            print("Converting to CSR format...")
        cell_lri_matrix = cell_lri_matrix.tocsr(copy=False)

    # Ensure data types are efficient
    cell_lri_matrix.data = cell_lri_matrix.data.astype(np.float64, copy=False)
    cell_lri_matrix.indices = cell_lri_matrix.indices.astype(np.int32, copy=False)
    cell_lri_matrix.indptr = cell_lri_matrix.indptr.astype(np.int64, copy=False)

    n_cells, n_lri = cell_lri_matrix.shape
    K = model.n_components
    alpha = model.alpha

    if verbose:
        print(f"Projecting {n_cells:,} cells to {K} motifs...")
        print(f"Cell-LRI matrix: {n_cells:,} cells × {n_lri:,} LRIs")
        print(f"BPTF model: K={K}, alpha={alpha}")

    # Extract fixed LRI factors (mode=1)
    shp1 = model.shp_DK_M[1]
    rte1 = model.rte_DK_M[1]
    bet1 = model.beta_M[1]

    # Initialize cell loadings array
    cell_loadings = np.zeros((n_cells, K), dtype=np.float64)

    # Process in chunks to manage memory
    n_chunks = (n_cells + chunk_size - 1) // chunk_size

    if verbose:
        print(f"Processing in {n_chunks} chunks of up to {chunk_size:,} cells...")

    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, n_cells)

        if verbose:
            print(f"  Chunk {chunk_idx + 1}/{n_chunks}: cells {start:,}-{end:,}")

        # Extract chunk
        mat_chunk_csr = cell_lri_matrix[start:end, :]

        # Convert to sparse.COO
        coo = mat_chunk_csr.tocoo(copy=False)
        coo.data = coo.data.astype(np.float64, copy=False)
        coo.row = coo.row.astype(np.int32, copy=False)
        coo.col = coo.col.astype(np.int32, copy=False)
        X_chunk = sparse.COO.from_scipy_sparse(coo)

        # Ensure consistent dtypes
        X_chunk.coords = X_chunk.coords.astype(np.int32, copy=False)
        if getattr(X_chunk, "has_duplicates", False):
            X_chunk = X_chunk.sum_duplicates()

        # Initialize small BPTF for this chunk
        proj = BPTF(data_shape=(end - start, n_lri), n_components=K, alpha=alpha)
        proj._init(modes=[0, 1])

        # Copy fixed column parameters (LRI factors)
        proj.shp_DK_M[1][...] = shp1
        proj.rte_DK_M[1][...] = rte1
        proj.beta_M[1] = bet1
        proj._update_cache(1)
        proj._clamp_component(1)

        # Update only mode=0 (cell loadings) for this chunk
        proj._update(X_chunk, modes=[0], max_iter=max_iter, verbose=False)

        # Store results
        cell_loadings[start:end, :] = proj.G_DK_M[0]

        # Clean up
        del proj, X_chunk, coo, mat_chunk_csr
        gc.collect()

    if verbose:
        print(f"✓ Projection complete. Cell loadings shape: {cell_loadings.shape}")

    return cell_loadings


def save_bptf_results(model: object,
                      patch_loadings: np.ndarray,
                      lri_factors: np.ndarray,
                      column_names: list,
                      patch_metadata_df: pd.DataFrame,
                      patch_lri_matrix: sp.spmatrix,
                      output_dir: str,
                      elbo_hist: Optional[list] = None,
                      delta_hist: Optional[list] = None):
    """
    Save all BPTF results to files

    Parameters
    ----------
    model : BPTF model
    patch_loadings : np.ndarray
        Patch loadings (n_patches × K)
    lri_factors : np.ndarray
        LRI factors (K × n_lris)
    column_names : list
        LRI column names
    patch_metadata_df : pd.DataFrame
        Patch metadata
    patch_lri_matrix : scipy.sparse matrix
        Original patch-LRI matrix for computing column means
    output_dir : str
        Output directory
    elbo_hist : list, optional
        ELBO history
    delta_hist : list, optional
        Delta history
    """
    os.makedirs(output_dir, exist_ok=True)

    # Save model
    model_path = Path(output_dir) / 'bptf_model.npz'
    save_bptf(model, model_path)
    print(f"BPTF model saved to: {model_path}")

    # Save factor matrices
    np.save(os.path.join(output_dir, 'patch_loadings.npy'), patch_loadings)
    np.save(os.path.join(output_dir, 'lri_factors.npy'), lri_factors)

    # Save convergence history
    if elbo_hist is not None and delta_hist is not None:
        history_df = pd.DataFrame({
            'iteration': range(len(elbo_hist)),
            'elbo': elbo_hist,
            'delta': delta_hist
        })
        history_df.to_csv(os.path.join(output_dir, 'iteration_history.csv'), index=False)

    # Analyze and save patch-motif relationships
    _save_patch_motif_analysis(patch_loadings, patch_metadata_df, output_dir)

    # Analyze and save LRI-motif relationships (with normalization)
    _save_lri_motif_analysis(lri_factors, column_names, patch_lri_matrix, output_dir)

    # Save model parameters
    params_df = pd.DataFrame({
        'parameter': ['n_components', 'n_patches', 'n_lris', 'method'],
        'value': [
            model.n_components,
            patch_loadings.shape[0],
            lri_factors.shape[1],
            'BPTF'
        ]
    })
    params_df.to_csv(os.path.join(output_dir, 'factorization_parameters.csv'), index=False)

    print(f"All results saved to: {output_dir}")


def _save_patch_motif_analysis(patch_loadings, patch_metadata_df, output_dir):
    """Save patch-motif relationship analysis"""
    patch_motifs = []
    if 'patch_id' in patch_metadata_df.columns:
        patch_ids = sorted(patch_metadata_df['patch_id'].astype(str).unique().tolist())
    else:
        patch_ids = [f"patch_{i}" for i in range(len(patch_loadings))]

    for i, patch_id in enumerate(patch_ids):
        total_loading = patch_loadings[i].sum()
        for k in range(patch_loadings.shape[1]):
            loading = patch_loadings[i, k]
            fraction = loading / total_loading if total_loading > 0 else 0

            patch_motifs.append({
                'patch_idx': i,
                'patch_id': patch_id,
                'motif': k,
                'loading': loading,
                'loading_normalized': fraction
            })

    patch_motifs_df = pd.DataFrame(patch_motifs)
    patch_motifs_df.to_csv(os.path.join(output_dir, 'patch_motifs.csv'), index=False)


def _save_lri_motif_analysis(lri_factors, column_names, patch_lri_matrix, output_dir):
    """Save LRI-motif relationship analysis with normalization

    Parameters
    ----------
    lri_factors : np.ndarray
        LRI factors (K × n_lris)
    column_names : list
        LRI column names
    patch_lri_matrix : scipy.sparse matrix
        Original patch-LRI matrix for computing column means
    output_dir : str
        Output directory
    """
    # Compute column means for normalization
    column_means = np.array(patch_lri_matrix.mean(axis=0)).flatten()
    lri_to_mean = dict(zip(column_names, column_means))

    lri_motifs = []
    for j, column_name in enumerate(column_names):
        motif_factors = lri_factors[:, j]  # factors for this LRI across motifs
        mean_expr = lri_to_mean.get(column_name, 0)

        for k in range(len(motif_factors)):
            factor = motif_factors[k]

            # Normalize by column mean
            if mean_expr > 0:
                factor_norm = factor / mean_expr
            else:
                factor_norm = factor

            lri_motifs.append({
                'lri_idx': j,
                'motif_idx': k,
                'lri_name': column_name,
                'factor': factor,
                'factor_norm': factor_norm,
                'mean': mean_expr
            })

    lri_motifs_df = pd.DataFrame(lri_motifs)
    lri_motifs_df.to_csv(os.path.join(output_dir, 'lri_motifs.csv'), index=False)

    # Save motif summary
    motif_activities = lri_factors.sum(axis=1)
    top_motifs = get_top_motifs(lri_factors.T, top_k=len(motif_activities))

    top_motifs_df = pd.DataFrame({
        'motif_idx': top_motifs['motif_indices'],
        'total_activity': top_motifs['activities'],
        'activity_fraction': top_motifs['activity_fractions']
    })
    top_motifs_df.to_csv(os.path.join(output_dir, 'top_motifs.csv'), index=False)
