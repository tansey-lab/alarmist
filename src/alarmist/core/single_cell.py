"""
Single cell analysis functions for motif calling

Includes GMM-based binarization and cell type composition analysis.
"""

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from typing import Optional, Tuple, Dict, Union, List
import anndata


def weighted_celltypes_by_motif(
    cell_loadings: np.ndarray,
    metadata_df: pd.DataFrame,
    normalize: bool = True,
    top_n_per_motif: Optional[int] = None,
    other_label: str = "Other"
) -> pd.DataFrame:
    """
    Compute cell type composition for each motif based on cell loadings.

    Parameters
    ----------
    cell_loadings : np.ndarray
        Cell loadings matrix (n_cells × n_motifs)
    metadata_df : pd.DataFrame
        Cell metadata, must contain 'cell_type' column
    normalize : bool, default True
        If True, per-motif weights sum to 1 (for percent stacked bar)
    top_n_per_motif : int, optional
        If set, keep top-N cell types per motif; rest collapsed to 'Other'
    other_label : str, default "Other"
        Label for collapsed cell types

    Returns
    -------
    pd.DataFrame
        Tidy dataframe with columns: [motif, cell_type, weight]

    Examples
    --------
    >>> tidy = weighted_celltypes_by_motif(
    ...     cell_loadings,
    ...     cell_meta_df,
    ...     normalize=True,
    ...     top_n_per_motif=8
    ... )
    """
    n_cells, n_motifs = cell_loadings.shape

    if len(metadata_df) != n_cells:
        raise ValueError("metadata_df length must match cell_loadings rows (n_cells).")
    if "cell_type" not in metadata_df.columns:
        raise ValueError("metadata_df must contain a 'cell_type' column.")

    # Build tidy long table of weights
    records = []
    cell_types = metadata_df["cell_type"].to_numpy()

    for k in range(n_motifs):
        w = cell_loadings[:, k]
        df_k = pd.DataFrame({"cell_type": cell_types, "weight": w})
        agg = df_k.groupby("cell_type", as_index=False)["weight"].sum()
        agg["motif"] = k

        # Optional: keep only top N cell types per motif
        if top_n_per_motif is not None and top_n_per_motif > 0 and len(agg) > top_n_per_motif:
            agg = agg.sort_values("weight", ascending=False)
            keep = agg.head(top_n_per_motif).copy()
            other_w = agg["weight"].iloc[top_n_per_motif:].sum()
            if other_w > 0:
                keep = pd.concat([
                    keep,
                    pd.DataFrame({
                        "cell_type": [other_label],
                        "weight": [other_w],
                        "motif": [k]
                    })
                ], ignore_index=True)
            agg = keep

        # Optional: normalize to sum to 1
        if normalize:
            s = agg["weight"].sum()
            if s > 0:
                agg["weight"] = agg["weight"] / s

        records.append(agg)

    result = pd.concat(records, ignore_index=True)
    return result


def gmm_binarize_all_motifs(
    cell_loadings: Union[np.ndarray, Dict[str, np.ndarray]],
    adata: Union[anndata.AnnData, Dict[str, anndata.AnnData], None] = None,
    eps: float = 1e-10,
    random_state: int = 0,
    return_arrays: bool = False
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Dict[str, np.ndarray], Dict[str, np.ndarray]]]:
    """
    Use Gaussian Mixture Model to binarize each motif into ON/OFF states.

    Supports both single-sample and multi-sample modes. In multi-sample mode,
    GMM is fit on all samples combined for consistent thresholds across samples.

    Parameters
    ----------
    cell_loadings : np.ndarray or Dict[str, np.ndarray]
        Single sample: Cell loadings matrix (n_cells × n_motifs)
        Multi-sample: Dict mapping sample_id -> cell loadings matrix
    adata : AnnData, Dict[str, AnnData], or None
        Single sample: AnnData object. Results added to adata.obs.
        Multi-sample: Dict mapping sample_id -> AnnData. Results added to each.
        If None, results are only returned (not attached to any adata).
    eps : float, default 1e-10
        Small constant added before log transform
    random_state : int, default 0
        Random seed for GMM
    return_arrays : bool, default False
        If True, also return state and posprob arrays (useful when adata is None)

    Returns
    -------
    pd.DataFrame
        Summary dataframe with GMM component statistics for each motif.
        Columns: [motif, mean0, mean1, weight0, weight1, positive_component]

    If return_arrays=True, returns tuple:
        (summary_df, states_dict, posprob_dict)
        - states_dict: Dict[str, np.ndarray] mapping sample_id -> state array (n_cells, n_motifs)
        - posprob_dict: Dict[str, np.ndarray] mapping sample_id -> posprob array (n_cells, n_motifs)
        For single sample mode, dict key is 'default'.

    Examples
    --------
    Single sample:
    >>> gmm_summary = gmm_binarize_all_motifs(cell_loadings, adata)

    Multi-sample (fit GMM on all samples combined):
    >>> cell_loadings_dict = {'P17_AIS': loadings1, 'P17_LUAD': loadings2}
    >>> adata_dict = {'P17_AIS': adata1, 'P17_LUAD': adata2}
    >>> gmm_summary = gmm_binarize_all_motifs(cell_loadings_dict, adata_dict)

    Return arrays without modifying adata:
    >>> summary, states, posprobs = gmm_binarize_all_motifs(
    ...     cell_loadings_dict, adata=None, return_arrays=True
    ... )
    """
    # Determine if multi-sample mode
    is_multi = isinstance(cell_loadings, dict)

    if is_multi:
        # Multi-sample mode: concatenate all loadings
        sample_ids = list(cell_loadings.keys())
        sample_sizes = {sid: cell_loadings[sid].shape[0] for sid in sample_ids}

        # Stack all loadings
        all_loadings = np.vstack([cell_loadings[sid] for sid in sample_ids])
        n_cells_total, K = all_loadings.shape

        # Validate adata_dict if provided
        if adata is not None:
            if not isinstance(adata, dict):
                raise ValueError("When cell_loadings is a dict, adata must also be a dict or None")
            for sid in sample_ids:
                if sid not in adata:
                    raise ValueError(f"Sample '{sid}' in cell_loadings but not in adata")
                if adata[sid].n_obs != sample_sizes[sid]:
                    raise ValueError(f"Sample '{sid}': adata has {adata[sid].n_obs} cells but loadings has {sample_sizes[sid]}")
    else:
        # Single-sample mode
        sample_ids = ['default']
        sample_sizes = {'default': cell_loadings.shape[0]}
        all_loadings = cell_loadings
        n_cells_total, K = all_loadings.shape

        # Wrap single adata in dict for uniform processing
        if adata is not None:
            adata = {'default': adata}

    # Initialize storage for results
    all_states = np.empty((n_cells_total, K), dtype=object)
    all_posprobs = np.zeros((n_cells_total, K))
    all_scores_log = np.zeros((n_cells_total, K))

    summary = []

    # Fit GMM for each motif on ALL cells combined
    for k in range(K):
        raw = all_loadings[:, k]
        scores_log = np.log(raw + eps)
        all_scores_log[:, k] = scores_log

        # Handle degenerate motifs (near-constant)
        if np.isfinite(scores_log).sum() < 2 or np.nanstd(scores_log) < 1e-6:
            state = np.array(["negative"] * n_cells_total)
            posprob = np.zeros(n_cells_total)
            means = [float(np.nanmean(scores_log)), float(np.nanmean(scores_log))]
            weights = [1.0, 0.0]
            pos_comp = 0
        else:
            gmm = GaussianMixture(n_components=2, random_state=random_state)
            gmm.fit(scores_log.reshape(-1, 1))
            means = gmm.means_.flatten().tolist()
            weights = gmm.weights_.flatten().tolist()

            labels = gmm.predict(scores_log.reshape(-1, 1))
            probs = gmm.predict_proba(scores_log.reshape(-1, 1))

            # Positive component is the one with higher mean
            pos_comp = int(np.argmax(gmm.means_.flatten()))
            state = np.where(labels == pos_comp, "positive", "negative")
            posprob = probs[:, pos_comp]

        all_states[:, k] = state
        all_posprobs[:, k] = posprob

        summary.append({
            "motif": k,
            "mean0": means[0],
            "mean1": means[1],
            "weight0": weights[0],
            "weight1": weights[1],
            "positive_component": pos_comp
        })

    # Split results back to samples and attach to adata
    states_dict = {}
    posprob_dict = {}
    scores_log_dict = {}

    start_idx = 0
    for sid in sample_ids:
        n_cells = sample_sizes[sid]
        end_idx = start_idx + n_cells

        sample_states = all_states[start_idx:end_idx, :]
        sample_posprobs = all_posprobs[start_idx:end_idx, :]
        sample_scores_log = all_scores_log[start_idx:end_idx, :]

        states_dict[sid] = sample_states
        posprob_dict[sid] = sample_posprobs
        scores_log_dict[sid] = sample_scores_log

        # Attach to adata if provided
        if adata is not None and sid in adata:
            for k in range(K):
                adata[sid].obs[f"motif_{k}_score_log"] = sample_scores_log[:, k]
                adata[sid].obs[f"motif_{k}_state"] = pd.Categorical(
                    sample_states[:, k],
                    categories=["negative", "positive"]
                )
                adata[sid].obs[f"motif_{k}_posprob"] = sample_posprobs[:, k]

        start_idx = end_idx

    # For single-sample mode, unwrap from dict
    if not is_multi:
        states_dict = states_dict['default']
        posprob_dict = posprob_dict['default']
        scores_log_dict = scores_log_dict['default']

    summary_df = pd.DataFrame(summary)

    if return_arrays:
        return summary_df, states_dict, posprob_dict
    else:
        return summary_df


def compute_motif_state_counts(adata) -> pd.DataFrame:
    """
    Count positive/negative cells for each motif.

    Parameters
    ----------
    adata : anndata.AnnData
        Must contain motif_{k}_state columns in obs

    Returns
    -------
    pd.DataFrame
        Counts dataframe with columns [positive, negative] and motif names as index

    Examples
    --------
    >>> counts_df = compute_motif_state_counts(adata)
    >>> print(counts_df)
    """
    # Extract only motif_state columns
    motif_state_cols = [col for col in adata.obs.columns if col.endswith("_state")]

    # Count positive/negative per motif
    counts = {}
    for col in motif_state_cols:
        motif_name = col.replace("_state", "")
        counts[motif_name] = adata.obs[col].value_counts()

    # Convert to DataFrame and fill missing values with 0
    df_counts = pd.DataFrame(counts).fillna(0).astype(int).T
    if 'positive' in df_counts.columns and 'negative' in df_counts.columns:
        df_counts = df_counts[['positive', 'negative']]

    return df_counts


def compute_positive_motifs_per_cell(adata) -> pd.Series:
    """
    Compute distribution of how many motifs are positive per cell.

    Parameters
    ----------
    adata : anndata.AnnData
        Must contain motif_{k}_state columns in obs

    Returns
    -------
    pd.Series
        Counts of cells for each number of positive motifs (0 to K)

    Examples
    --------
    >>> dist = compute_positive_motifs_per_cell(adata)
    >>> print(dist)
    """
    # Find all state columns
    state_cols = [c for c in adata.obs.columns if c.endswith('_state')]
    K = len(state_cols)

    # Count positive motifs per cell
    n_pos_per_cell = adata.obs[state_cols].eq('positive').sum(axis=1)

    # Full count from 0 to K
    counts = n_pos_per_cell.value_counts().reindex(range(K + 1), fill_value=0).sort_index()

    return counts
