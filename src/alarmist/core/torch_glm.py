"""Batched Poisson GLM via torch (device-agnostic: CPU or CUDA).

Within a single (motif x cell-type) DE comparison every gene is regressed on the
*same* design matrix ``X = [1, z(log motif-loading)]`` (n_cells x 2); only the
response ``y`` (that gene's counts) differs. That makes the whole set of genes a
batch of independent 2-parameter Poisson GLMs sharing one ``X`` — ideal for a
vectorized IRLS/Newton solve on the GPU.

The math (per Newton step, all genes at once):
    eta = X @ B.T               # (N, G)
    mu  = exp(eta)
    grad = X.T @ (Y - mu)       # (2, G)
    H    = X.T diag(mu) X        # (G, 2, 2), built from column sums Sigma mu*{1,x,x^2}
    B  += solve(H, grad)         # closed-form 2x2

Standard errors use the same Fisher-information expression as the sklearn path
(``fisher = sum(mu * x^2)``, ``se = 1/sqrt(fisher)``) so results are directly
comparable; validated to ~1e-6 on beta and ~1e-9 on se against sklearn.
"""

from __future__ import annotations

import logging

import numpy as np
import scipy.sparse as sp

logger = logging.getLogger(__name__)


def _resolve_device(device: str):
    import torch

    if device in (None, "auto"):
        if torch.cuda.is_available():
            return torch.device("cuda")
        if (
            getattr(torch.backends, "mps", None) is not None
            and torch.backends.mps.is_available()
        ):
            return torch.device("mps")
        return torch.device("cpu")
    return torch.device(device)


def batched_poisson_glm(
    x: np.ndarray,
    Y,
    device: str = "auto",
    max_iter: int = 100,
    tol: float = 1e-8,
    gene_tile: int = 2048,
    dtype: str = "float64",
):
    """Fit ``y_g ~ Poisson(exp(b0_g + b1_g * x))`` for every column g of ``Y``.

    Parameters
    ----------
    x : np.ndarray, shape (n_cells,)
        Single predictor (e.g. z-scored log motif loading). Intercept added here.
    Y : np.ndarray or scipy.sparse, shape (n_cells, n_genes)
        Count matrix for the genes to test. Densified one tile at a time.
    device : {"auto", "cpu", "cuda", ...}
        Torch device. "auto" -> cuda if available else cpu.
    max_iter : int
        Max Newton iterations per tile.
    tol : float
        Convergence tol on max abs parameter update.
    gene_tile : int
        Number of genes per tile (bounds memory: ~ n_cells * gene_tile * 8 bytes).
    dtype : {"float64", "float32"}
        float64 recommended for parity with sklearn; float32 ~2x faster on GPU.

    Returns
    -------
    beta1 : np.ndarray (n_genes,)
        Slope estimate per gene.
    se : np.ndarray (n_genes,)
        Standard error of beta1 (Fisher-info form matching the sklearn path).
    n_iter : np.ndarray (n_genes,)
        Newton iterations used (per tile, broadcast to genes).
    """
    import torch

    dev = _resolve_device(device)
    tdt = torch.float64 if dtype == "float64" else torch.float32
    logger.info(
        "torch GLM backend: device=%s (requested=%s), dtype=%s", dev, device, dtype
    )
    if dev.type == "mps" and tdt == torch.float64:
        # MPS does not support float64; fall back to float32 and warn.
        logger.warning("MPS does not support float64; falling back to float32.")
        tdt = torch.float32
    # Floor tol at ~10x dtype epsilon — Newton steps can't shrink below precision,
    # so a stricter tol would falsely flag converged fits as non-convergent.
    eps_floor = 10.0 * float(torch.finfo(tdt).eps)
    if tol < eps_floor:
        logger.info(
            "Raising tol from %.0e to %.0e (10x %s epsilon).", tol, eps_floor, tdt
        )
        tol = eps_floor

    n_cells = x.shape[0]
    n_genes = Y.shape[1]
    issparse = sp.issparse(Y)
    if issparse:
        Y = Y.tocsc()

    x_t = torch.as_tensor(np.asarray(x, dtype=np.float64), dtype=tdt, device=dev)
    ones = torch.ones(n_cells, dtype=tdt, device=dev)
    X = torch.stack([ones, x_t], dim=1)  # (N, 2)
    x2 = x_t * x_t

    beta1 = np.empty(n_genes, dtype=np.float64)
    se = np.empty(n_genes, dtype=np.float64)
    n_iter_out = np.zeros(n_genes, dtype=np.int32)
    not_converged_total = 0

    for start in range(0, n_genes, gene_tile):
        end = min(start + gene_tile, n_genes)
        block = Y[:, start:end]
        if issparse:
            block = block.toarray()
        Yt = torch.as_tensor(np.asarray(block, dtype=np.float64), dtype=tdt, device=dev)
        G = Yt.shape[1]

        B = torch.zeros((2, G), dtype=tdt, device=dev)  # [b0; b1]
        last_iter = max_iter
        for it in range(max_iter):
            eta = X @ B  # (N, G)
            mu = torch.exp(torch.clamp(eta, -30.0, 30.0))
            r = Yt - mu
            g0 = r.sum(0)
            g1 = (x_t[:, None] * r).sum(0)
            S0 = mu.sum(0)
            S1 = (x_t[:, None] * mu).sum(0)
            S2 = (x2[:, None] * mu).sum(0)
            det = S0 * S2 - S1 * S1
            det = torch.where(det.abs() < 1e-12, torch.full_like(det, 1e-12), det)
            d0 = (S2 * g0 - S1 * g1) / det
            d1 = (-S1 * g0 + S0 * g1) / det
            B[0] += d0
            B[1] += d1
            # Per-gene step magnitude; loop exits when *every* gene has converged.
            per_gene_step = torch.maximum(d0.abs(), d1.abs())
            if float(per_gene_step.max()) < tol:
                last_iter = it + 1
                break
        else:
            # max_iter exhausted without all genes converging — count and warn.
            n_bad = int((per_gene_step >= tol).sum().item())
            if n_bad:
                worst = float(per_gene_step.max().item())
                logger.warning(
                    "torch GLM tile [%d:%d]: %d/%d genes did not converge in "
                    "max_iter=%d (worst step=%.3e, tol=%.0e). Results for those "
                    "genes may be unreliable; consider raising max_iter.",
                    start,
                    end,
                    n_bad,
                    G,
                    max_iter,
                    worst,
                    tol,
                )
                not_converged_total += n_bad

        eta = X @ B
        mu = torch.exp(torch.clamp(eta, -30.0, 30.0))
        fisher = (x2[:, None] * mu).sum(0)
        se_blk = torch.where(
            fisher > 0, 1.0 / torch.sqrt(fisher), torch.full_like(fisher, float("nan"))
        )

        beta1[start:end] = B[1].detach().cpu().numpy().astype(np.float64)
        se[start:end] = se_blk.detach().cpu().numpy().astype(np.float64)
        n_iter_out[start:end] = last_iter

    if not_converged_total:
        raise RuntimeError(
            f"torch GLM: {not_converged_total}/{n_genes} genes failed to converge "
            f"within max_iter={max_iter}. Consider raising max_iter."
        )

    return beta1, se, n_iter_out
