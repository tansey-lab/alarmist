"""
Similarity between the motifs of two ALARMIST runs.

Motif indices are arbitrary across runs, so comparing two fits (two platforms,
two datasets, patch vs neighbourhood) means comparing every motif of one run with
every motif of the other on a shared feature axis (e.g. LRI loadings), then
pairing them up.
"""

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import rankdata


def _row_center_normalize(X: np.ndarray) -> np.ndarray:
    X = X - X.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(X, axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    return X / norm


def row_similarity(
    A: np.ndarray,
    B: np.ndarray,
    method: str = "cosine",
) -> np.ndarray:
    """
    Pairwise similarity between the rows of two matrices.

    Parameters
    ----------
    A : np.ndarray
        (n_motifs_a, n_features), e.g. LRI loadings of run A.
    B : np.ndarray
        (n_motifs_b, n_features) on the same, identically ordered features.
    method : {'cosine', 'pearson', 'spearman'}, default 'cosine'
        Cosine on the raw rows, Pearson on mean-centred rows, or Spearman
        (Pearson on per-row ranks). Constant rows give 0 similarity.

    Returns
    -------
    np.ndarray
        (n_motifs_a, n_motifs_b) similarity matrix.
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    if A.ndim != 2 or B.ndim != 2 or A.shape[1] != B.shape[1]:
        raise ValueError(
            f"A and B must be 2-D with the same number of columns, got {A.shape} "
            f"and {B.shape}"
        )

    if method == "cosine":
        na = np.linalg.norm(A, axis=1, keepdims=True)
        nb = np.linalg.norm(B, axis=1, keepdims=True)
        na[na == 0] = 1.0
        nb[nb == 0] = 1.0
        return (A / na) @ (B / nb).T
    if method == "spearman":
        A = np.apply_along_axis(rankdata, 1, A)
        B = np.apply_along_axis(rankdata, 1, B)
    elif method != "pearson":
        raise ValueError(f"Unknown method '{method}'")
    return _row_center_normalize(A) @ _row_center_normalize(B).T


def match_motifs(
    similarity: np.ndarray,
    min_similarity: float | None = None,
) -> list[tuple[int, int]]:
    """
    One-to-one motif matching that maximises total similarity (Hungarian).

    Parameters
    ----------
    similarity : np.ndarray
        (n_a, n_b) matrix, e.g. from :func:`row_similarity`.
    min_similarity : float, optional
        Drop matched pairs below this similarity (e.g. 0.4 for cosine).

    Returns
    -------
    list of (int, int)
        Matched (row, column) index pairs, sorted by row.
    """
    S = np.asarray(similarity, dtype=float)
    rows, cols = linear_sum_assignment(-np.nan_to_num(S, nan=-np.inf))
    pairs = [(int(r), int(c)) for r, c in zip(rows, cols, strict=True)]
    if min_similarity is not None:
        pairs = [(r, c) for r, c in pairs if S[r, c] >= min_similarity]
    return pairs
