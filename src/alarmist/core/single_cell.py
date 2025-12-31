"""
Single cell analysis functions for motif calling

Includes GMM-based binarization and cell type composition analysis.
"""

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from typing import Optional, Tuple


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
    cell_loadings: np.ndarray,
    adata,
    eps: float = 1e-10,
    random_state: int = 0
) -> pd.DataFrame:
    """
    Use Gaussian Mixture Model to binarize each motif into ON/OFF states.

    For each motif (column in cell_loadings), fit a 2-component GMM on log(loading+eps),
    assign positive/negative states, and attach results to adata.obs.

    Parameters
    ----------
    cell_loadings : np.ndarray
        Cell loadings matrix (n_cells × n_motifs)
    adata : anndata.AnnData
        Annotated data object. Results will be added to adata.obs:
        - motif_{k}_score_log: log-transformed loading
        - motif_{k}_state: categorical ('negative', 'positive')
        - motif_{k}_posprob: probability of positive state
    eps : float, default 1e-10
        Small constant added before log transform
    random_state : int, default 0
        Random seed for GMM

    Returns
    -------
    pd.DataFrame
        Summary dataframe with GMM component statistics for each motif.
        Columns: [motif, mean0, mean1, weight0, weight1, positive_component]

    Examples
    --------
    >>> gmm_summary = gmm_binarize_all_motifs(cell_loadings, adata)
    >>> print(gmm_summary)
    >>> # Check states
    >>> print(adata.obs['motif_0_state'].value_counts())
    """
    n_cells, K = cell_loadings.shape
    summary = []

    for k in range(K):
        raw = cell_loadings[:, k]
        scores_log = np.log(raw + eps)

        # Handle degenerate motifs (near-constant)
        if np.isfinite(scores_log).sum() < 2 or np.nanstd(scores_log) < 1e-6:
            state = np.array(["negative"] * n_cells)
            posprob = np.zeros(n_cells)
            means = [scores_log.mean(), scores_log.mean()]
            weights = [1.0, 0.0]
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

        # Attach to adata.obs
        adata.obs[f"motif_{k}_score_log"] = scores_log
        adata.obs[f"motif_{k}_state"] = pd.Categorical(
            state,
            categories=["negative", "positive"]
        )
        adata.obs[f"motif_{k}_posprob"] = posprob

        summary.append({
            "motif": k,
            "mean0": means[0],
            "mean1": means[1],
            "weight0": weights[0],
            "weight1": weights[1],
            "positive_component": int(np.argmax(means))
        })

    return pd.DataFrame(summary)


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
