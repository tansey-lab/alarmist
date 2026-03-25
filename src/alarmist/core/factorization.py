"""
BPTF (Bayesian Poisson Tensor Factorization) wrapper functions
"""

import logging
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import sparse
    from bptf import BPTF, save_bptf

    BPTF_AVAILABLE = True
except ImportError:
    BPTF_AVAILABLE = False
    warnings.warn("BPTF not available. Install from: https://github.com/aschein/bptf")

from alarmist.plotting import add_lri_components, annotate_pathways

logger = logging.getLogger(__name__)


def run_bptf(
    mat: sp.spmatrix,
    n_components: int = 15,
    max_iter: int = 10000,
    verbose: bool = True,
    random_state: int = 0,
) -> object:
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
        raise ImportError(
            "BPTF not available. Install: pip install git+https://github.com/aschein/bptf.git"
        )

    # Convert to sparse.COO format
    if isinstance(mat, sp.spmatrix):
        data = sparse.COO(mat)
    else:
        data = mat

    np.random.seed(random_state)

    model = BPTF(data_shape=data.shape, n_components=n_components)
    model.fit(data, max_iter=max_iter, verbose=verbose)

    return model


def extract_factors(model: object) -> tuple[np.ndarray, np.ndarray]:
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


def get_top_motifs(patch_loadings: np.ndarray, top_k: int = 10) -> dict:
    """Get top motifs by total activity"""
    motif_activities = patch_loadings.sum(axis=0)
    top_indices = np.argsort(motif_activities)[::-1][:top_k]

    return {
        "motif_indices": top_indices.tolist(),
        "activities": motif_activities[top_indices].tolist(),
        "activity_fractions": (
            motif_activities[top_indices] / motif_activities.sum()
        ).tolist(),
    }


def project_cell_loadings(
    model: object,
    cell_lri_matrix: sp.spmatrix,
    model_lri_columns: np.ndarray | None = None,
    cell_lri_columns: np.ndarray | None = None,
    max_iter: int = 200,
    chunk_size: int = 50000,
    verbose: bool = True,
    output_dir: str | None = None,
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
    output_dir : str, optional
        Directory to save cell_loadings.npy. If None, results are not saved.

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
        raise ImportError(
            "BPTF not available. Install: pip install git+https://github.com/aschein/bptf.git"
        )

    import gc

    # Align LRI columns if names provided
    if model_lri_columns is not None and cell_lri_columns is not None:
        # Quick check: if columns are already aligned, skip alignment
        if len(model_lri_columns) == len(cell_lri_columns) and all(
            a == b for a, b in zip(model_lri_columns, cell_lri_columns)
        ):
            if verbose:
                logger.debug(f"Columns already aligned ({len(model_lri_columns)} LRIs)")
        else:
            # Columns need alignment
            if verbose:
                logger.debug("Aligning cell-LRI matrix columns to model...")

            # Find intersection and preserve model's order
            from collections import defaultdict, deque

            a = cell_lri_columns
            b = model_lri_columns

            # Which b elements appear in a (preserving b's order and duplicates)
            mask = np.isin(b, a)
            _idx_in_model = np.flatnonzero(mask)

            # Build queue for one-to-one matching
            pos = defaultdict(deque)
            for i, v in enumerate(a):
                pos[v].append(i)

            # Match each b[mask] value to next unused position in a
            idx_in_mat = np.array(
                [pos[v].popleft() if pos[v] else -1 for v in b[mask]], dtype=int
            )

            if (idx_in_mat == -1).any():
                raise ValueError(
                    "Some model LRIs not found in cell matrix. "
                    "Use run_neighborhood(..., required_columns=patch_columns) to ensure alignment."
                )

            # Reorder cell matrix columns to match model
            cell_lri_matrix = cell_lri_matrix[:, idx_in_mat]

            if cell_lri_matrix.shape[1] != len(model_lri_columns):
                raise ValueError(
                    f"Column mismatch: cell matrix has {cell_lri_matrix.shape[1]} "
                    f"columns after alignment, model expects {len(model_lri_columns)}"
                )

            if verbose:
                logger.debug(f"Aligned to {cell_lri_matrix.shape[1]} LRI columns")

    # Convert to CSR if needed
    if not sp.isspmatrix_csr(cell_lri_matrix):
        if verbose:
            logger.debug("Converting to CSR format...")
        cell_lri_matrix = cell_lri_matrix.tocsr(copy=False)

    # Ensure data types are efficient
    cell_lri_matrix.data = cell_lri_matrix.data.astype(np.float64, copy=False)
    cell_lri_matrix.indices = cell_lri_matrix.indices.astype(np.int32, copy=False)
    cell_lri_matrix.indptr = cell_lri_matrix.indptr.astype(np.int64, copy=False)

    n_cells, n_lri = cell_lri_matrix.shape
    K = model.n_components
    alpha = model.alpha

    if verbose:
        logger.debug(f"Projecting {n_cells:,} cells to {K} motifs...")
        logger.debug(f"Cell-LRI matrix: {n_cells:,} cells × {n_lri:,} LRIs")
        logger.debug(f"BPTF model: K={K}, alpha={alpha}")

    # Extract fixed LRI factors (mode=1)
    shp1 = model.shp_DK_M[1]
    rte1 = model.rte_DK_M[1]
    bet1 = model.beta_M[1]

    # Initialize cell loadings array
    cell_loadings = np.zeros((n_cells, K), dtype=np.float64)

    # Process in chunks to manage memory
    n_chunks = (n_cells + chunk_size - 1) // chunk_size

    if verbose:
        logger.debug(
            f"Processing in {n_chunks} chunks of up to {chunk_size:,} cells..."
        )

    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, n_cells)

        if verbose:
            logger.debug(f"  Chunk {chunk_idx + 1}/{n_chunks}: cells {start:,}-{end:,}")

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
        logger.debug(f"Projection complete. Cell loadings shape: {cell_loadings.shape}")

    # Save results if output_dir provided
    if output_dir is not None:
        import os

        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, "cell_loadings.npy")
        np.save(save_path, cell_loadings)
        if verbose:
            logger.debug(f"Cell loadings saved to: {save_path}")

    return cell_loadings


def process_bptf_results(
    model: object,
    results: dict,
    output_dir: str | None = None,
    cellchatdb_path: str = "data/LRdatabase/CellChatDBv2.0.human.csv",
) -> dict:
    """
    Process BPTF model results and optionally save to disk.

    Parameters
    ----------
    model : BPTF model
        Fitted BPTF model object
    results : dict
        Results from run_patchify or load_patch_lri_results, containing:
        - patch_lri_matrix: sparse matrix
        - column_names: list of column names
    output_dir : str, optional
        Output directory. If None, results are not saved to disk.
    cellchatdb_path : str
        Path to CellChatDB annotation file for pathway annotation

    Returns
    -------
    dict
        Dictionary with:
        - patch_loadings: np.ndarray (n_patches × K)
        - lri_factors: np.ndarray (K × n_lris)
        - patch_loadings_rescaled: np.ndarray (n_patches × K), max=1 per motif
        - lri_factors_rescaled: np.ndarray (K × n_lris)
        - lri_motifs: pd.DataFrame with all scores and annotations

    Examples
    --------
    >>> model = al.run_bptf(results['patch_lri_matrix'], n_components=20)
    >>> bptf_results = al.process_bptf_results(model, results, output_dir='bptf_results/')
    >>> # Or without saving:
    >>> bptf_results = al.process_bptf_results(model, results)
    """
    # Extract factors from model
    patch_loadings, lri_factors = extract_factors(model)

    # Get data from results
    column_names = results["column_names"]
    patch_lri_matrix = results["patch_lri_matrix"]

    # Compute rescaled matrices
    logger.debug("Rescaling motif matrices...")
    patch_loadings_rescaled, lri_factors_rescaled, motif_scales = (
        rescale_motif_matrices(patch_loadings, lri_factors, verify=True)
    )

    # Create LRI motifs DataFrame
    logger.debug("Creating LRI motifs DataFrame...")
    column_means = np.array(patch_lri_matrix.mean(axis=0)).flatten()
    lri_to_mean = dict(zip(column_names, column_means))

    lri_motifs_list = []
    for j, column_name in enumerate(column_names):
        motif_factors = lri_factors[:, j]
        mean_expr = lri_to_mean.get(column_name, 0)

        for k in range(len(motif_factors)):
            lri_motifs_list.append(
                {
                    "lri_idx": j,
                    "motif_idx": k,
                    "lri_name": column_name,
                    "factor": motif_factors[k],
                    "mean": mean_expr,
                }
            )

    lri_motifs = pd.DataFrame(lri_motifs_list)
    lri_motifs = lri_motifs[lri_motifs["mean"] > 0]
    logger.debug(f"Total entries: {len(lri_motifs)}")

    # Parse LRI components
    logger.debug("Parsing LRI components...")
    lri_motifs = add_lri_components(lri_motifs)

    # Annotate pathways
    logger.debug("Annotating pathways...")
    if os.path.exists(cellchatdb_path):
        cellchatdb = pd.read_csv(cellchatdb_path)
        lri_motifs = annotate_pathways(lri_motifs, cellchatdb)
    else:
        logger.debug(
            f"Warning: CellChatDB not found at {cellchatdb_path}, skipping pathway annotation"
        )

    # Add normalized scores
    logger.debug("Computing normalized scores...")
    lri_motifs = add_normalized_scores(lri_motifs, motif_scales)

    # Build return dict
    bptf_results = {
        "patch_loadings": patch_loadings,
        "lri_factors": lri_factors,
        "patch_loadings_rescaled": patch_loadings_rescaled,
        "lri_factors_rescaled": lri_factors_rescaled,
        "lri_motifs": lri_motifs,
    }

    # Save if output_dir provided
    if output_dir is not None:
        logger.debug(f"Saving to {output_dir}...")
        os.makedirs(output_dir, exist_ok=True)

        # Save model
        model_path = Path(output_dir)
        save_bptf(model, model_path)
        logger.debug("  BPTF model saved")

        # Save factor matrices
        np.save(os.path.join(output_dir, "patch_loadings.npy"), patch_loadings)
        np.save(os.path.join(output_dir, "lri_factors.npy"), lri_factors)
        np.save(
            os.path.join(output_dir, "patch_loadings_rescaled.npy"),
            patch_loadings_rescaled,
        )
        np.save(
            os.path.join(output_dir, "lri_factors_rescaled.npy"), lri_factors_rescaled
        )

        # Save lri_motifs
        lri_motifs.to_csv(os.path.join(output_dir, "lri_motifs.csv"), index=False)

        # Save convergence history if available
        if hasattr(model, "elbo_hist") and hasattr(model, "delta_hist"):
            elbo_hist = model.elbo_hist
            delta_hist = model.delta_hist
            if elbo_hist is not None and delta_hist is not None:
                history_df = pd.DataFrame(
                    {
                        "iteration": range(len(elbo_hist)),
                        "elbo": elbo_hist,
                        "delta": delta_hist,
                    }
                )
                history_df.to_csv(
                    os.path.join(output_dir, "iteration_history.csv"), index=False
                )

        # Save model parameters
        params_df = pd.DataFrame(
            {
                "parameter": ["n_components", "n_patches", "n_lris", "method"],
                "value": [
                    model.n_components,
                    patch_loadings.shape[0],
                    lri_factors.shape[1],
                    "BPTF",
                ],
            }
        )
        params_df.to_csv(
            os.path.join(output_dir, "factorization_parameters.csv"), index=False
        )

        logger.debug(f"All results saved to: {output_dir}")

    return bptf_results


# Keep save_bptf_results as alias for backwards compatibility
def save_bptf_results(model: object, results: dict, output_dir: str) -> dict:
    """Alias for process_bptf_results with output_dir. See process_bptf_results for docs."""
    return process_bptf_results(model, results, output_dir=output_dir)


def rescale_motif_matrices(
    W: np.ndarray, V: np.ndarray, verify: bool = True, eps: float = 1e-10
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Rescale W and V so each motif's patch loading has max=1.

    Transformation preserves reconstruction: W @ V = W_tilde @ V_tilde

    Parameters
    ----------
    W : (n_patches, n_motifs) patch loadings
    V : (n_motifs, n_lris) LRI factors
    verify : whether to check reconstruction invariance
    eps : small constant for numerical stability

    Returns
    -------
    W_tilde, V_tilde, c_safe (scaling factors per motif)
    """
    c = W.max(axis=0)  # (n_motifs,)
    c_safe = np.where(c > 0, c, 1.0)

    W_tilde = W / c_safe[None, :]
    V_tilde = V * c_safe[:, None]

    if verify:
        _verify_reconstruction(W, V, W_tilde, V_tilde, eps=eps)

    return W_tilde, V_tilde, c_safe


def _verify_reconstruction(
    W: np.ndarray,
    V: np.ndarray,
    W_tilde: np.ndarray,
    V_tilde: np.ndarray,
    n_samples: int = 200,
    seed: int = 0,
    eps: float = 1e-10,
) -> None:
    """Sanity check that rescaling preserves reconstruction."""
    rng = np.random.default_rng(seed)
    patch_idx = rng.choice(W.shape[0], size=min(n_samples, W.shape[0]), replace=False)
    lri_idx = rng.choice(V.shape[1], size=min(n_samples, V.shape[1]), replace=False)

    X1 = W[patch_idx, :] @ V[:, lri_idx]
    X2 = W_tilde[patch_idx, :] @ V_tilde[:, lri_idx]
    rel_err = np.linalg.norm(X1 - X2) / (np.linalg.norm(X1) + eps)

    logger.debug(f"Reconstruction error: {rel_err:.2e} (should be ~0)")
    assert rel_err < 1e-10, f"Reconstruction failed: {rel_err}"


def compute_lr_global_prevalence(df: pd.DataFrame) -> pd.Series:
    """
    Compute global prevalence for each (ligand, receptor) pair.

    Uses unique lri_idx only to avoid overcounting repeated entries across motifs.
    """
    unique_lris = df.drop_duplicates(subset=["lri_idx"])
    return unique_lris.groupby(["ligand", "receptor"])["mean"].sum()


def add_normalized_scores(
    df: pd.DataFrame, motif_scales: np.ndarray, eps: float = 1
) -> pd.DataFrame:
    """
    Add rescaled and normalized score columns to lri_motifs.

    Columns added:
    - W_max: motif scaling factor (max patch loading for that motif)
    - factor_rescaled: factor × W_max
    - lr_global_mean: background prevalence of (ligand, receptor) pair
    - score: factor_rescaled / lr_global_mean
    """
    df = df.copy()

    # Motif-wise rescaling
    df["W_max"] = df["motif_idx"].map(dict(enumerate(motif_scales))).astype(float)
    df["factor_rescaled"] = df["factor"] * df["W_max"]

    # LR global prevalence normalization
    lr_prevalence = compute_lr_global_prevalence(df)
    lr_keys = list(zip(df["ligand"], df["receptor"]))
    df["lr_global_mean"] = (
        pd.Series(lr_keys, index=df.index).map(lr_prevalence).astype(float)
    )

    df["factor_lrnorm"] = df["factor"] / (df["lr_global_mean"] + eps)

    # Final scores
    df["score"] = df["factor_rescaled"] / (df["lr_global_mean"] + eps)

    return df


def process_and_save_lri_motif_analysis(
    lri_factors: np.ndarray,
    patch_loadings: np.ndarray,
    column_names: list,
    patch_lri_matrix,
    output_dir: str,
    cellchatdb_path: str = "data/LRdatabase/CellChatDBv2.0.human.csv",
    save: bool = True,
    eps_rescale: float = 1e-10,
    eps_norm: float = 1,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """
    Full pipeline: create LRI motifs DataFrame, process, and save.

    Parameters
    ----------
    lri_factors : np.ndarray
        LRI factors (K × n_lris)
    patch_loadings : np.ndarray
        Patch loadings (n_patches × K)
    column_names : list
        LRI column names
    patch_lri_matrix : scipy.sparse matrix
        Original patch-LRI matrix for computing column means
    output_dir : str
        Output directory
    cellchatdb_path : str
        Path to CellChatDB annotation file
    save : bool
        Whether to save outputs to disk
    eps_rescale : float
        Small constant for numerical stability in rescaling
    eps_norm : float
        Small constant for numerical stability in normalization

    Returns
    -------
    lri_motifs : pd.DataFrame
        Fully processed LRI motifs with all scores
    W_tilde : np.ndarray
        Rescaled patch loadings
    V_tilde : np.ndarray
        Rescaled LRI factors
    """
    # Step 1: Create initial DataFrame
    logger.debug("Creating LRI motifs DataFrame...")
    column_means = np.array(patch_lri_matrix.mean(axis=0)).flatten()
    lri_to_mean = dict(zip(column_names, column_means))

    lri_motifs = []
    for j, column_name in enumerate(column_names):
        motif_factors = lri_factors[:, j]
        mean_expr = lri_to_mean.get(column_name, 0)

        for k in range(len(motif_factors)):
            lri_motifs.append(
                {
                    "lri_idx": j,
                    "motif_idx": k,
                    "lri_name": column_name,
                    "factor": motif_factors[k],
                    "mean": mean_expr,
                }
            )

    lri_motifs_df = pd.DataFrame(lri_motifs)
    logger.debug(f"Total entries: {len(lri_motifs_df)}")

    # Step 2: Parse LRI components
    logger.debug("Parsing LRI components...")
    lri_motifs_df = add_lri_components(lri_motifs_df)

    # Step 3: Annotate pathways
    logger.debug("Annotating pathways...")
    cellchatdb = pd.read_csv(cellchatdb_path)
    lri_motifs_df = annotate_pathways(lri_motifs_df, cellchatdb)

    # Step 4: Rescale matrices
    logger.debug("Rescaling motif matrices...")
    W_tilde, V_tilde, motif_scales = rescale_motif_matrices(
        patch_loadings, lri_factors, verify=True, eps=eps_rescale
    )

    # Step 5: Add normalized scores
    logger.debug("Computing normalized scores...")
    lri_motifs_df = add_normalized_scores(lri_motifs_df, motif_scales, eps=eps_norm)

    # Step 6: Save
    if save:
        logger.debug(f"Saving to {output_dir}...")
        os.makedirs(output_dir, exist_ok=True)
        lri_motifs_df.to_csv(os.path.join(output_dir, "lri_motifs.csv"), index=False)
        np.save(os.path.join(output_dir, "W_tilde.npy"), W_tilde)
        np.save(os.path.join(output_dir, "V_tilde.npy"), V_tilde)
        logger.debug("Done!")


# def _save_lri_motif_analysis(lri_factors, column_names, patch_lri_matrix, output_dir):
#     """Save LRI-motif relationship analysis with normalization

#     Parameters
#     ----------
#     lri_factors : np.ndarray
#         LRI factors (K × n_lris)
#     column_names : list
#         LRI column names
#     patch_lri_matrix : scipy.sparse matrix
#         Original patch-LRI matrix for computing column means
#     output_dir : str
#         Output directory
#     """
#     # Compute column means for normalization
#     column_means = np.array(patch_lri_matrix.mean(axis=0)).flatten()
#     lri_to_mean = dict(zip(column_names, column_means))

#     lri_motifs = []
#     for j, column_name in enumerate(column_names):
#         motif_factors = lri_factors[:, j]  # factors for this LRI across motifs
#         mean_expr = lri_to_mean.get(column_name, 0)

#         for k in range(len(motif_factors)):
#             factor = motif_factors[k]

#             lri_motifs.append({
#                 'lri_idx': j,
#                 'motif_idx': k,
#                 'lri_name': column_name,
#                 'factor': factor,
#                 'mean': mean_expr
#             })

#     lri_motifs_df = pd.DataFrame(lri_motifs)
#     lri_motifs_df.to_csv(os.path.join(output_dir, 'lri_motifs.csv'), index=False)
