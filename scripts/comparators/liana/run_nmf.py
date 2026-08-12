#!/usr/bin/env python
"""NMF on LIANA+ CCC output, following the NMF section of bivariate.ipynb.

    li.multi.nmf(lrdata, n_components=None, inplace=True, random_state=0, max_iter=200)
    li.ut.get_variable_loadings(lrdata, varm_key='NMF_H')
    li.ut.get_factor_scores(lrdata, obsm_key='NMF_W')

Run on BOTH LIANA outputs and kept side by side:
  * bivariate  -- the tutorial-sanctioned composition (NMF is demonstrated only here)
  * inflow     -- the resolution-appropriate input; NMF on it is OUR composition. The authors'
                  own unsupervised route for single-cell spatial is Inflow + MOFA-Flex
                  (inflow_mofaflex.ipynb), a DIFFERENT factorisation. See DEVIATIONS.md.

This is the only comparator output sharing ALARMIST's shape: locations x factors (W) and
features x factors (H).

Usage: python run_nmf.py --input {bivariate,inflow} --in-dir DIR --out-dir DIR [--k-max N]
"""
import argparse, json, os, time
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import liana as li

p = argparse.ArgumentParser()
p.add_argument("--input", required=True, choices=["bivariate", "inflow"])
p.add_argument("--in-dir", required=True)
p.add_argument("--out-dir", required=True)
p.add_argument("--k-min", type=int, default=1)
p.add_argument("--k-max", type=int, default=11, help="exclusive upper bound, as range()")
p.add_argument("--k-step", type=int, default=1)
p.add_argument("--max-iter", type=int, default=200)
p.add_argument("--top-n", type=int, default=10)
p.add_argument("--seed", type=int, default=0)
p.add_argument("--min-punches", type=int, default=5,
               help="cross-punch reproducibility filter applied BEFORE NMF. The LIANA paper "
                    "keeps interactions present in >=10 of 28 sections; scaled to 13 TMA "
                    "punches that is >=5/13. A feature counts as present in a punch if it is "
                    "non-zero in at least one cell of that punch. 0 disables.")
a = p.parse_args()

OUT = a.out_dir
PLOTS, DATA = os.path.join(OUT, "plots"), os.path.join(OUT, "data")
for d in (PLOTS, DATA): os.makedirs(d, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)

# ------------------------------------------------------------------- load
meta = pd.read_csv(os.path.join(a.in_dir, "data", "cell_meta.csv"))
if a.input == "inflow":
    lrdata = sc.read_h5ad(os.path.join(a.in_dir, "data", "inflow_lrdata.h5ad"))
else:
    z = np.load(os.path.join(a.in_dir, "data", "local_scores.npz"), allow_pickle=True)
    lrdata = ad.AnnData(X=np.asarray(z["values"], dtype=np.float32),
                        obs=pd.DataFrame(index=z["cells"].astype(str)),
                        var=pd.DataFrame(index=z["pairs"].astype(str)))
if "spatial" not in lrdata.obsm:
    lrdata.obsm["spatial"] = meta[["x", "y"]].values.astype(float)
lrdata.obs["cell_type"] = pd.Categorical(meta["cell_type"].astype(str).values)
X = lrdata.X.toarray() if hasattr(lrdata.X, "toarray") else np.asarray(lrdata.X)
log(f"{a.input}: {lrdata.shape[0]} locations x {lrdata.shape[1]} features | "
    f"min {X.min():.4g} | negatives {(X < 0).sum()}")

# ---- cross-punch reproducibility filter (LIANA paper: >=10/28 sections -> >=5/13 punches) ----
n_before = lrdata.shape[1]
if a.min_punches > 0:
    if "core" not in meta.columns:
        raise SystemExit("STOP: cell_meta.csv has no 'core' column; re-run the upstream script "
                         "so the punch ID is written, or pass --min-punches 0.")
    # cell_meta['cell'] may be int64 while obs_names are str -> cast BOTH before aligning
    _m = meta.copy(); _m["cell"] = _m["cell"].astype(str)
    cores = _m.set_index("cell").loc[lrdata.obs_names.astype(str), "core"].astype(str).values
    uniq = np.unique(cores)
    present = np.zeros((len(uniq), lrdata.shape[1]), dtype=bool)
    for i, c in enumerate(uniq):
        present[i] = (X[cores == c] > 0).any(axis=0)
    n_punch = present.sum(axis=0)
    keep = n_punch >= a.min_punches
    log(f"  cross-punch filter (>={a.min_punches} of {len(uniq)} punches): "
        f"{n_before} -> {int(keep.sum())} features "
        f"({100*keep.mean():.1f}% kept); median punches/feature = {int(np.median(n_punch))}")
    pd.DataFrame({"feature": lrdata.var_names.astype(str), "n_punches": n_punch,
                  "kept": keep}).to_csv(os.path.join(DATA, "punch_presence.csv"), index=False)
    if keep.sum() < 2:
        raise SystemExit("STOP: fewer than 2 features survive the cross-punch filter.")
    lrdata = lrdata[:, keep].copy()
    X = lrdata.X.toarray() if hasattr(lrdata.X, "toarray") else np.asarray(lrdata.X)
n_after = lrdata.shape[1]
if (X < 0).sum():
    raise SystemExit("STOP: negative values -- NMF is not applicable. Do not clip; report.")

# --------------------------------------------------------------- run NMF
k_range = range(a.k_min, a.k_max, a.k_step)
log(f"li.multi.nmf(n_components=None, k_range={k_range}, max_iter={a.max_iter}, "
    f"random_state={a.seed})")
plt.close("all")
li.multi.nmf(lrdata, n_components=None, k_range=k_range, inplace=True,
             random_state=a.seed, max_iter=a.max_iter, verbose=True)
plt.close("all")

rank = int(lrdata.uns["nmf_rank"])
errors = lrdata.uns.get("nmf_errors")

# li.multi.nmf calls _plot_elbow() internally, but that draws with PLOTNINE, not matplotlib --
# plt.gcf() is empty afterwards, so the figure cannot be captured. Redraw it here from the
# errors that li.multi.nmf itself stored, with the same content (error vs k, dashed line at the
# chosen rank). Done inside the runner so it survives a re-run.
if errors is not None:
    _e = pd.DataFrame(errors)
    _fig, _ax = plt.subplots(figsize=(7, 4))
    _ax.plot(_e["k"], _e["error"], "-o", ms=3.5, lw=1.2, color="#4C72B0")
    _ax.axvline(rank, ls="--", c="red", lw=1.2)
    _ax.annotate(f"rank = {rank}", (rank, _e["error"].max()), color="red", fontsize=9,
                 ha="left", va="top", xytext=(4, -2), textcoords="offset points")
    _ax.set_xlabel("Component number (k)"); _ax.set_ylabel("Reconstruction error")
    _ax.set_title(f"NMF elbow — {a.input} (rank={rank}, k_range "
                  f"{a.k_min}..{a.k_max - 1})", fontsize=10)
    _ax.grid(alpha=.3); _fig.tight_layout()
    _fig.savefig(os.path.join(PLOTS, "elbow.png"), dpi=200, bbox_inches="tight")
    plt.close("all")
    log("  saved elbow.png (redrawn from nmf_errors; _plot_elbow uses plotnine)")
log(f"  ELBOW RANK = {rank}  (k_range {a.k_min}..{a.k_max - 1} step {a.k_step})")
if rank in (a.k_min, a.k_max - 1):
    log("  *** WARNING: rank sits at a k_range boundary -- the elbow may be truncated, "
        "not detected. Widen k_range. ***")

W = lrdata.obsm["NMF_W"]; H = lrdata.varm["NMF_H"]
log(f"  W (locations x factors) {W.shape} | H (features x factors) {H.shape}")

loadings = li.ut.get_variable_loadings(lrdata, varm_key="NMF_H").set_index("index")
scores = li.ut.get_factor_scores(lrdata, obsm_key="NMF_W")
loadings.to_csv(os.path.join(DATA, "NMF_H_loadings.csv"))
scores.to_csv(os.path.join(DATA, "NMF_W_factor_scores.csv"), index=False)
np.savez_compressed(os.path.join(DATA, "nmf_WH.npz"), W=np.asarray(W, np.float32),
                    H=np.asarray(H, np.float32),
                    features=np.asarray(lrdata.var_names).astype(str),
                    cells=np.asarray(lrdata.obs_names).astype(str))
if errors is not None:
    pd.DataFrame(errors).to_csv(os.path.join(DATA, "nmf_errors.csv"), index=False)

# top-N loadings per factor
rows = []
for f in loadings.columns:
    for r, (feat, val) in enumerate(loadings[f].sort_values(ascending=False).head(a.top_n).items(), 1):
        rows.append(dict(factor=f, rank=r, feature=feat, loading=float(val)))
top = pd.DataFrame(rows)
top.to_csv(os.path.join(DATA, f"top{a.top_n}_loadings_per_factor.csv"), index=False)
log(f"  saved top-{a.top_n} loadings for {loadings.shape[1]} factors")

# --------------------------------------------------------- spatial maps
nmf = ad.AnnData(X=np.asarray(W), obs=lrdata.obs.copy(),
                 var=pd.DataFrame(index=loadings.columns), obsm=dict(spatial=lrdata.obsm["spatial"]))
factors = list(nmf.var.index)
try:
    ncol = 3
    nrow = int(np.ceil(len(factors) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow))
    axes = np.atleast_1d(axes).ravel()
    xy = nmf.obsm["spatial"]
    for i, f in enumerate(factors):
        v = np.asarray(nmf[:, f].X).ravel()
        axes[i].scatter(xy[:, 0], xy[:, 1], c=v, s=0.6, cmap="magma",
                        vmax=np.quantile(v, 0.99))
        axes[i].set_title(f, fontsize=9); axes[i].set_aspect("equal"); axes[i].axis("off")
    for j in range(len(factors), len(axes)): axes[j].axis("off")
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "factor_maps.png"), dpi=180,
                                    bbox_inches="tight"); plt.close("all")
    log(f"  saved factor_maps.png ({len(factors)} factors)")
except Exception as e:
    plt.close("all"); log(f"  factor_maps FAILED: {type(e).__name__}: {e}")

# ------- fit quality, so "rank 6" / "rank 7" is never reported bare -------
# li.multi.nmf picks its elbow on MAE. On a matrix that is ~99% zeros the MAE of the zero
# predictor is a meaningful floor, so record it next to the achieved error, plus the relative
# Frobenius error, which is what a reader actually wants when comparing to a BPTF fit.
_Wf, _Hf = np.asarray(W, float), np.asarray(H, float)
_rec = _Wf @ _Hf.T
_rel_frob = float(np.linalg.norm(X - _rec) / np.linalg.norm(X))
_mae = float(np.abs(X - _rec).mean())
_mae_zero = float(np.abs(X).mean())
log(f"  fit quality: rel_frobenius {_rel_frob:.4f} "
    f"({100*(1-_rel_frob**2):.1f}% of SS captured) | MAE {_mae:.6f} "
    f"vs zero-predictor MAE {_mae_zero:.6f}")
if _mae >= _mae_zero:
    log("  *** WARNING: the fit's MAE does not beat the zero predictor. The elbow is a property "
        "of an L1 criterion on a near-all-zero matrix, not evidence of a good rank. ***")

try:
    import subprocess
    _sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True).stdout.strip() or None
except Exception:
    _sha = None
json.dump({"script": "run_nmf.py", "input": a.input, "in_dir": a.in_dir,
           "dataset": "GBM", "tier": "cellchatdb2",
           "n_locations": int(lrdata.shape[0]),
           "n_features_before_punch_filter": int(n_before),
           "n_features": int(n_after), "min_punches": a.min_punches,
           "n_punches": int(len(uniq)) if a.min_punches > 0 else None,
           "k_range": [a.k_min, a.k_max, a.k_step], "nmf_rank": rank,
           "rank_at_boundary": bool(rank in (a.k_min, a.k_max - 1)),
           "rel_frobenius_error": round(_rel_frob, 6),
           "frac_sumsq_captured": round(1 - _rel_frob ** 2, 6),
           "mae": round(_mae, 8), "mae_zero_predictor": round(_mae_zero, 8),
           "beats_zero_predictor": bool(_mae < _mae_zero),
           "pct_zeros_input": round(float(100 * (X == 0).mean()), 4),
           "max_iter": a.max_iter, "seed": a.seed, "liana": li.__version__, "git_sha": _sha,
           "wall_min": round((time.time() - t0) / 60, 1)},
          open(os.path.join(OUT, "run_manifest.json"), "w"), indent=2)
log(f"DONE in {(time.time()-t0)/60:.1f} min -> {OUT}")
