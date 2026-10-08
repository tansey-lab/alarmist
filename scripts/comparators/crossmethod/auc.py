"""One-dimensional-classifier scoring for the cross-method benchmark.

WHY AUC AND NOT A FITTED LOGISTIC COEFFICIENT
---------------------------------------------
Wesley asked for a one-dimensional logistic regression per feature.  With a
single predictor the fitted model is monotone in that predictor, so the ROC
curve -- and therefore the AUC -- of the fitted probabilities is IDENTICAL to
the rank AUC of the raw feature.  The fit adds nothing to the ranking and only
costs an optimiser per feature (there are ~10k features per section).  We
therefore score by rank AUC and *prove* the identity numerically against
sklearn's `LogisticRegression` in `validate_logistic_identity`, rather than
asserting it.  That check is part of the run and its result is recorded.

  AUC = U1 / (n1 * n0)   with U1 the Mann-Whitney U of the positive class.

TWO IMPLEMENTATIONS, ON PURPOSE
-------------------------------
CLAUDE.md: "Grep the package before writing a statistical primitive... When
both exist, prefer the package function and record the numerical comparison."
`alarmist.core.glm.mann_whitney_u_sparse_nonneg` already computes exactly this
U column-wise on nonnegative sparse CSC with average-rank tie correction, so
sparse nonnegative score matrices (LIANA inflow) go through it.  SpatialDM's
local Moran z takes negative values, which that function's implicit-zero
accounting cannot represent, so those go through a dense `scipy.stats.rankdata`
path.  `compare_implementations` runs both on the same data and reports the
maximum discrepancy.

DIRECTION
---------
`auc` is signed: > 0.5 means the feature is HIGHER in the positive (tumour)
class.  `auc_abs = max(auc, 1 - auc)` is the direction-agnostic "ability to
predict this label", which is what the distribution comparison uses -- an LR
pair that is strongly DOWN in tumour is just as predictive, and giving every
comparator feature that credit makes the test harder for ALARMIST, not easier.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse, stats

__all__ = ['auc_dense', 'auc_sparse_nonneg', 'auc_abs', 'validate_logistic_identity',
           'compare_implementations']


def auc_abs(auc: np.ndarray) -> np.ndarray:
    """Direction-agnostic predictive ability."""
    return np.maximum(auc, 1.0 - auc)


def auc_dense(X, y, *, chunk: int = 256) -> np.ndarray:
    """Rank AUC per column of a dense 2-D array.

    Parameters
    ----------
    X : (n_units, n_features) array; may contain negatives.  NaNs are not
        allowed -- callers must have dropped or imputed them, because a NaN
        would silently sort to the top of `rankdata`.
    y : (n_units,) bool -- the positive class.
    chunk : columns per `rankdata` call, to bound the float64 rank buffer.
    """
    X = np.asarray(X)
    if X.ndim == 1:
        X = X[:, None]
    y = np.asarray(y, dtype=bool)
    if X.shape[0] != y.shape[0]:
        raise ValueError(f'X has {X.shape[0]} rows, y has {y.shape[0]}')
    if not np.isfinite(X).all():
        raise ValueError('auc_dense: X contains non-finite values')

    n1 = int(y.sum())
    n0 = int((~y).sum())
    if n1 == 0 or n0 == 0:
        raise ValueError(f'degenerate label: n_pos={n1}, n_neg={n0}')

    out = np.empty(X.shape[1], dtype=np.float64)
    for a in range(0, X.shape[1], chunk):
        b = min(a + chunk, X.shape[1])
        R = stats.rankdata(X[:, a:b], method='average', axis=0)
        r1 = R[y].sum(axis=0)
        out[a:b] = (r1 - n1 * (n1 + 1) / 2.0) / (n1 * n0)
    return out


def auc_sparse_nonneg(X, y) -> np.ndarray:
    """Rank AUC per column via the PACKAGE's tie-corrected sparse U.

    X must be nonnegative CSC; zeros are the implicit minimum, which is what
    `mann_whitney_u_sparse_nonneg` assumes.
    """
    import sys
    from pathlib import Path
    src = Path(__file__).resolve().parents[3] / 'src'
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    from alarmist.core.glm import mann_whitney_u_sparse_nonneg

    X = sparse.csc_matrix(X)
    if X.shape[0] != y.shape[0]:
        raise ValueError(f'X has {X.shape[0]} rows, y has {y.shape[0]}')
    if X.data.size and X.data.min() < 0:
        raise ValueError('auc_sparse_nonneg: negative values present; use auc_dense')
    y = np.asarray(y, dtype=bool)
    n1, n0 = int(y.sum()), int((~y).sum())
    # X1 = positives -> returned u_stats is U for group 1 (see glm.py:1091).
    u1, _ = mann_whitney_u_sparse_nonneg(X[y], X[~y], alternative='two-sided')
    return u1 / (n1 * n0)


def validate_logistic_identity(X, y, *, n_features: int = 12, seed: int = 0) -> dict:
    """Prove AUC(1-D logistic fit) == rank AUC of the raw feature.

    Fits sklearn `LogisticRegression` on one standardised feature at a time and
    compares `roc_auc_score` of the fitted probabilities with `auc_dense`.
    Returns the max absolute difference; the caller records it.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    X = np.asarray(X)
    y = np.asarray(y, dtype=bool)
    rng = np.random.default_rng(seed)
    # Only features with some spread can be fitted at all.
    spread = X.max(axis=0) - X.min(axis=0)
    usable = np.flatnonzero(spread > 0)
    pick = rng.choice(usable, size=min(n_features, usable.size), replace=False)

    diffs, rows = [], []
    for j in pick:
        x = X[:, j].astype(np.float64)
        s = x.std()
        xs = ((x - x.mean()) / s if s > 0 else x)[:, None]
        lr = LogisticRegression(max_iter=1000).fit(xs, y)
        auc_fit = roc_auc_score(y, lr.predict_proba(xs)[:, 1])
        auc_rank = float(auc_dense(x[:, None], y)[0])
        diffs.append(abs(auc_fit - auc_rank))
        rows.append({'feature_index': int(j), 'auc_logistic': float(auc_fit),
                     'auc_rank': auc_rank, 'abs_diff': float(diffs[-1])})
    return {'n_features_checked': len(rows), 'max_abs_diff': float(max(diffs)),
            'per_feature': rows}


def compare_implementations(X_sparse, y, *, n_features: int = 200, seed: int = 0) -> dict:
    """Run the package sparse U and the dense rankdata path on the same columns."""
    rng = np.random.default_rng(seed)
    X_sparse = sparse.csc_matrix(X_sparse)
    pick = rng.choice(X_sparse.shape[1], size=min(n_features, X_sparse.shape[1]),
                      replace=False)
    sub = X_sparse[:, pick]
    a_pkg = auc_sparse_nonneg(sub, y)
    a_dns = auc_dense(np.asarray(sub.todense()), y)
    d = np.abs(a_pkg - a_dns)
    return {'n_features_compared': int(pick.size), 'max_abs_diff': float(d.max()),
            'median_abs_diff': float(np.median(d))}
