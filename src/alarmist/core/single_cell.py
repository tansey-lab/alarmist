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
    multi_sample: bool = False,
    sample_column: Optional[str] = None,
    eps: float = 1e-10,
    random_state: int = 42,
    return_arrays: bool = False
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Dict[str, np.ndarray], Dict[str, np.ndarray]]]:
    """
    Use Gaussian Mixture Model to binarize each motif into ON/OFF states.

    Supports three input modes:
    1. Single AnnData: standard single-sample analysis
    2. Dict of AnnData objects: multi-sample with dict
    3. Single merged AnnData with sample_column: multi-sample merged

    In multi-sample mode, GMM is fit on all samples combined for consistent
    thresholds across samples.

    Parameters
    ----------
    cell_loadings : np.ndarray or Dict[str, np.ndarray]
        Single sample: Cell loadings matrix (n_cells × n_motifs)
        Multi-sample dict: Dict mapping sample_id -> cell loadings matrix
        Multi-sample merged: Cell loadings matrix for merged samples
    adata : AnnData, Dict[str, AnnData], or None
        Single sample: AnnData object. Results added to adata.obs.
        Multi-sample dict: Dict mapping sample_id -> AnnData. Results added to each.
        Multi-sample merged: Single AnnData with sample_column in obs.
        If None, results are only returned (not attached to any adata).
        When provided, adds to adata.obs:
        - motif_{k}_loading: original cell loading value
        - motif_{k}_state: 'positive' or 'negative'
    multi_sample : bool, default False
        Whether to treat single AnnData as containing multiple samples.
        Only used when cell_loadings is np.ndarray and adata is single AnnData.
    sample_column : str, optional
        Column name in adata.obs identifying samples (for merged AnnData).
        Required when multi_sample=True.
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

    Multi-sample with dict:
    >>> cell_loadings_dict = {'P17_AIS': loadings1, 'P17_LUAD': loadings2}
    >>> adata_dict = {'P17_AIS': adata1, 'P17_LUAD': adata2}
    >>> gmm_summary = gmm_binarize_all_motifs(cell_loadings_dict, adata_dict)

    Multi-sample with merged AnnData:
    >>> gmm_summary = gmm_binarize_all_motifs(
    ...     cell_loadings, adata_merged,
    ...     multi_sample=True, sample_column='patient_id'
    ... )

    Return arrays without modifying adata:
    >>> summary, states, posprobs = gmm_binarize_all_motifs(
    ...     cell_loadings_dict, adata=None, return_arrays=True
    ... )
    """
    # Determine input mode
    is_dict_mode = isinstance(cell_loadings, dict)
    is_merged_mode = (not is_dict_mode and multi_sample and sample_column is not None)

    if is_dict_mode:
        # Mode 2: Multi-sample dict mode
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

        # Track original adata for merged mode processing
        adata_merged_original = None

    elif is_merged_mode:
        # Mode 3: Multi-sample merged mode
        if adata is None:
            raise ValueError("adata must be provided when multi_sample=True")
        if sample_column not in adata.obs.columns:
            raise ValueError(f"sample_column '{sample_column}' not found in adata.obs")

        # Validate cell count match
        if adata.n_obs != cell_loadings.shape[0]:
            raise ValueError(f"adata has {adata.n_obs} cells but cell_loadings has {cell_loadings.shape[0]} rows")

        all_loadings = cell_loadings
        n_cells_total, K = all_loadings.shape

        # Get sample info from adata
        sample_labels = adata.obs[sample_column].values
        unique_samples = adata.obs[sample_column].unique()
        sample_ids = list(unique_samples)

        # Compute sample sizes and indices
        sample_sizes = {}
        sample_indices = {}
        for sid in sample_ids:
            mask = sample_labels == sid
            sample_sizes[sid] = mask.sum()
            sample_indices[sid] = np.where(mask)[0]

        # Store original merged adata for later
        adata_merged_original = adata
        adata = None  # Will handle separately

    else:
        # Mode 1: Single-sample mode
        sample_ids = ['default']
        sample_sizes = {'default': cell_loadings.shape[0]}
        all_loadings = cell_loadings
        n_cells_total, K = all_loadings.shape

        # Wrap single adata in dict for uniform processing
        if adata is not None:
            adata = {'default': adata}

        # Track original adata for merged mode processing
        adata_merged_original = None

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

    # Handle merged mode: attach results directly to merged adata
    if is_merged_mode and adata_merged_original is not None:
        for k in range(K):
            adata_merged_original.obs[f"motif_{k}_loading"] = all_loadings[:, k]
            adata_merged_original.obs[f"motif_{k}_state"] = pd.Categorical(
                all_states[:, k],
                categories=["negative", "positive"]
            )

        # For merged mode, return arrays directly (not split by sample)
        summary_df = pd.DataFrame(summary)

        if return_arrays:
            return summary_df, all_states, all_posprobs
        else:
            return summary_df

    # Split results back to samples and attach to adata (dict mode and single mode)
    states_dict = {}
    posprob_dict = {}

    start_idx = 0
    for sid in sample_ids:
        n_cells = sample_sizes[sid]
        end_idx = start_idx + n_cells

        sample_states = all_states[start_idx:end_idx, :]
        sample_posprobs = all_posprobs[start_idx:end_idx, :]
        sample_loadings = all_loadings[start_idx:end_idx, :]

        states_dict[sid] = sample_states
        posprob_dict[sid] = sample_posprobs

        # Attach to adata if provided
        if adata is not None and sid in adata:
            for k in range(K):
                adata[sid].obs[f"motif_{k}_loading"] = sample_loadings[:, k]
                adata[sid].obs[f"motif_{k}_state"] = pd.Categorical(
                    sample_states[:, k],
                    categories=["negative", "positive"]
                )

        start_idx = end_idx

    # For single-sample mode, unwrap from dict
    if not is_dict_mode:
        states_dict = states_dict['default']
        posprob_dict = posprob_dict['default']

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
