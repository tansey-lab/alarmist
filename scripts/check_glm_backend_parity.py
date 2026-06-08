"""End-to-end parity check: torch vs sklearn Poisson GLM backends.

Simulates `y_g ~ Poisson(exp(b0_g + b1_g * x))` for many genes on a shared
design matrix `X = [1, z(log loading)]`, fits both backends, and reports the
worst-case differences in beta1, SE, and two-sided Wald p-value.

Run:
    python scripts/check_glm_backend_parity.py
"""

from __future__ import annotations

import argparse
import logging
import sys
import time

import numpy as np
from scipy.stats import norm
from sklearn.linear_model import PoissonRegressor

from alarmist.core.torch_glm import batched_poisson_glm


def simulate(n_cells: int, n_genes: int, seed: int):
    rng = np.random.default_rng(seed)
    loading = rng.lognormal(mean=0.0, sigma=1.0, size=n_cells)
    x_log = np.log(loading)
    x = (x_log - x_log.mean()) / (x_log.std() or 1.0)

    b0 = rng.normal(0.5, 0.5, size=n_genes)
    b1 = rng.normal(0.0, 0.4, size=n_genes)
    eta = b0[None, :] + x[:, None] * b1[None, :]
    mu = np.exp(np.clip(eta, -20, 20))
    Y = rng.poisson(mu).astype(np.float64)
    return x, Y, b1


def fit_sklearn(x: np.ndarray, Y: np.ndarray):
    X = x.reshape(-1, 1)
    n_genes = Y.shape[1]
    beta1 = np.empty(n_genes)
    se = np.empty(n_genes)
    for g in range(n_genes):
        m = PoissonRegressor(alpha=0.0, fit_intercept=True, max_iter=2000, tol=1e-6)
        m.fit(X, Y[:, g])
        beta1[g] = m.coef_[0]
        mu = m.predict(X)
        fisher = np.sum(mu * (x**2))
        se[g] = np.sqrt(1.0 / fisher) if fisher > 0 else np.nan
    return beta1, se


def pvals(beta1: np.ndarray, se: np.ndarray) -> np.ndarray:
    z = np.divide(beta1, se, out=np.zeros_like(beta1), where=(se > 0) & np.isfinite(se))
    return 2.0 * (1.0 - norm.cdf(np.abs(z)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-cells", type=int, default=500)
    ap.add_argument("--n-genes", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", default="auto")
    ap.add_argument("--dtype", default="float64")
    ap.add_argument("--beta-tol", type=float, default=1e-4)
    ap.add_argument("--se-tol", type=float, default=1e-4)
    ap.add_argument("--p-tol", type=float, default=1e-4)
    ap.add_argument("--max-iter", type=int, default=25)
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s %(name)s: %(message)s"
    )

    print(f"Simulating n_cells={args.n_cells} n_genes={args.n_genes} seed={args.seed}")
    x, Y, _ = simulate(args.n_cells, args.n_genes, args.seed)

    print("Fitting sklearn (one gene at a time)...")
    t0 = time.perf_counter()
    b_sk, se_sk = fit_sklearn(x, Y)
    t_sk = time.perf_counter() - t0
    p_sk = pvals(b_sk, se_sk)

    print(f"Fitting torch (batched, device={args.device}, dtype={args.dtype})...")
    t0 = time.perf_counter()
    b_t, se_t, _ = batched_poisson_glm(
        x,
        Y,
        device=args.device,
        dtype=args.dtype,
        gene_tile=2048,
        max_iter=args.max_iter,
    )
    t_torch = time.perf_counter() - t0
    p_t = pvals(b_t, se_t)

    print(
        f"Wall time: sklearn={t_sk:.2f}s  torch={t_torch:.2f}s  "
        f"speedup={t_sk / t_torch:.1f}x  (delta={t_sk - t_torch:+.2f}s)"
    )

    db = np.abs(b_sk - b_t)
    ds = np.abs(se_sk - se_t)
    dp = np.abs(p_sk - p_t)

    def report(name, d, tol):
        ok = np.nanmax(d) <= tol
        print(
            f"  {name}: max={np.nanmax(d):.3e}  mean={np.nanmean(d):.3e}  "
            f"tol={tol:.0e}  {'OK' if ok else 'FAIL'}"
        )
        return ok

    print("Differences (sklearn vs torch):")
    ok = all(
        [
            report("beta1", db, args.beta_tol),
            report("se   ", ds, args.se_tol),
            report("pval ", dp, args.p_tol),
        ]
    )
    # Spot-print the worst gene.
    worst = int(np.nanargmax(db))
    print(
        f"Worst-beta gene idx={worst}: "
        f"sklearn(b={b_sk[worst]:.6f}, se={se_sk[worst]:.6f}, p={p_sk[worst]:.3e})  "
        f"torch  (b={b_t[worst]:.6f}, se={se_t[worst]:.6f}, p={p_t[worst]:.3e})"
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
