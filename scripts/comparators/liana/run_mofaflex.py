#!/usr/bin/env python
"""MOFA-Flex on a finished LIANA+ ``inflow`` output.

Follows ``liana-py/docs/notebooks/inflow_mofaflex.ipynb`` cell by cell, starting from an
already-computed inflow AnnData (we do NOT re-run ``li.mt.inflow``).

Everything is parameterised; nothing about the GBM dataset is hardcoded except the
defaults you can override on the command line.

Example
-------
    python run_mofaflex.py \
        --lrdata results/comparators/liana/GBM/cellchatdb2_inflow/data/inflow_lrdata.h5ad \
        --outdir results/comparators/liana/GBM/mofaflex_inflow \
        --region-key grade --sample-key tma_id --drop-obs motif patch_id
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import colors  # noqa: E402


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lrdata", required=True, help="Path to the finished inflow lrdata .h5ad")
    p.add_argument("--outdir", required=True, help="Output directory (data/, plots/, models/ are created)")

    # tutorial cell 18/19 + 21
    p.add_argument("--nzf-threshold", type=float, default=0.01,
                   help="Keep features with var['nonzero_fraction'] > this (tutorial cell 19 uses 0.01)")
    p.add_argument("--nzf-mode", choices=["global", "reachability"], default="global",
                   help="'global' = the tutorial's criterion, nonzero_fraction > thr over ALL "
                        "cells. 'reachability' divides by P(>=1 neighbour of that feature's "
                        "sender type) first, so the threshold is scale-free in sender abundance. "
                        "DEVIATION from the tutorial; default stays 'global' so the "
                        "tutorial-faithful run remains reproducible.")
    p.add_argument("--bandwidth", type=float, default=13.1454,
                   help="gaussian sigma of the kernel the inflow scores were built with; used "
                        "only to rebuild the support radius for --nzf-mode reachability")
    p.add_argument("--cutoff", type=float, default=0.1,
                   help="kernel cutoff; support radius R = bandwidth * sqrt(-2 ln cutoff)")
    p.add_argument("--bandwidth-support-radius", type=float, default=None,
                   help="override R directly instead of deriving it from --bandwidth/--cutoff")
    p.add_argument("--xy-sep", default="^", help="separator in '<sender>^<lig>^<rec>' feature names")
    p.add_argument("--moran-fdr", type=float, default=0.05, help="SVI filter FDR cutoff (tutorial cell 21)")
    p.add_argument("--moran-i", type=float, default=0.01, help="SVI filter Moran's I cutoff (tutorial cell 21)")
    p.add_argument("--skip-svi", action="store_true", help="Skip the Moran SVI filter entirely")

    # tutorial cell 23
    p.add_argument("--min-features", type=int, default=25, help="lrdata_to_mudata min_features")
    p.add_argument("--min-cells", type=int, default=5, help="lrdata_to_mudata min_cells")

    # tutorial cell 25
    p.add_argument("--n-factors", type=int, default=20,
                   help="CEILING on the number of factors, not a selection -- MOFA-Flex fits "
                        "exactly this many and inactive ones are pruned afterwards by --r2-floor. "
                        "The default 20 was chosen because inflow_mofaflex.ipynb cell 25 uses "
                        "n_factors=20; it coincidentally equals ALARMIST's K=20 for this dataset, "
                        "but that is NOT the justification. Check the ceiling-binding warning "
                        "printed after the fit before quoting the active count as a result.")
    p.add_argument("--batch-size", type=int, default=2048)
    p.add_argument("--n-particles", type=int, default=1)
    p.add_argument("--lr", type=float, default=0.005)
    p.add_argument("--max-epochs", type=int, default=1000)
    p.add_argument("--patience", type=int, default=50)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--n-inducing", type=int, default=100, help="GaussianProcess n_inducing (mofaflex default 100)")
    p.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    p.add_argument("--densify", default="auto", choices=["auto", "yes", "no"],
                   help="densify the inflow matrix before training; 'auto' densifies on MPS, "
                        "which cannot hold sparse CSR tensors")

    # our-data adaptations
    p.add_argument("--spatial-key", default="spatial", help="obsm key with the coordinates (tutorial: X_spatial_coords)")
    p.add_argument("--celltype-key", default="cell_type")
    p.add_argument("--region-key", default="grade", help="tutorial uses major_brain_region; ours is grade")
    p.add_argument("--sample-key", default="tma_id", help="replicate unit for the group-level test")
    p.add_argument("--drop-obs", nargs="*", default=["motif", "patch_id"],
                   help="obs columns to drop at load (ALARMIST leakage)")
    p.add_argument("--drop-uns", nargs="*", default=["motif_colors", "moranI"],
                   help="uns keys to drop at load (leakage / stale gene-level moranI)")

    # reporting
    p.add_argument("--r2-floor", type=float, default=0.02, help="active-factor R2 floor (tutorial cell 38)")
    p.add_argument("--leiden-resolution", type=float, default=0.4)
    p.add_argument("--n-focus-factors", type=int, default=2)
    p.add_argument("--top-weights-n", type=int, default=10,
                   help="features per facet in top_weights.png (tutorial cell 43 uses 5)")
    p.add_argument("--circle-all-factors", action=argparse.BooleanOptionalAction, default=True,
                   help="emit circle_plot_Factor<N>.png for EVERY active factor. The tutorial "
                        "draws them only for its hand-picked focus factors.")
    p.add_argument("--circle-top-n", type=int, default=10,
                   help="interactions per circle plot, by |loading| (tutorial cell 59 uses 10)")
    p.add_argument("--lr-of-interest", nargs="*", default=["GRN^SORT1", "ANXA1^FPR1"],
                   help="'<ligand>^<receptor>' strings to trace through the loadings")
    p.add_argument("--tag", default="mofaflex_inflow", help="model file stem")
    p.add_argument("--determinism-check", action="store_true",
                   help="fit twice for --determinism-epochs and compare weights (MPS has no deterministic-algorithm guarantee)")
    p.add_argument("--determinism-epochs", type=int, default=20)
    p.add_argument("--refit", action="store_true", help="delete an existing model file and refit")
    return p.parse_args(argv)


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def peak_rss_gb() -> float:
    """macOS reports ru_maxrss in bytes, Linux in kB."""
    v = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return v / 1024**3 if sys.platform == "darwin" else v / 1024**2


def save_gg(gg, path, width=8, height=6, dpi=200):
    """plotnine ggplot -> png."""
    gg.save(str(path), width=width, height=height, dpi=dpi, verbose=False)
    print(f"  [plot] {path}")


def save_fig(fig, path, dpi=200):
    fig.savefig(str(path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] {path}")


def save_current(path, dpi=200):
    """scanpy/liana plotters that draw onto the current MATPLOTLIB figure.

    Refuses to write a blank canvas. Several li.pl.* entry points (dotplot, connectivity)
    return a PLOTNINE ggplot and never touch plt.gcf() -- saving the current figure there
    silently produces an all-white PNG. Use save_gg() for those instead.
    """
    fig = plt.gcf()
    has = any(ax.has_data() or ax.get_images() or ax.collections or ax.patches or ax.texts
              for ax in fig.axes)
    if not fig.axes or not has:
        plt.close(fig)
        print(f"  [plot] BLANK, not written: {path}  "
              f"(plotnine return-type trap? use save_gg)")
        return False
    fig.savefig(str(path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] {path}")
    return True


def save_gg(gg, path, width=18, height=6, dpi=200):
    """plotnine ggplot returned by li.pl.dotplot / li.pl.connectivity -- must use its own .save()."""
    gg.save(str(path), dpi=dpi, width=width, height=height, verbose=False)
    plt.close("all")
    print(f"  [plot] {path}")
    return True


def git_sha(repo: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


# --------------------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)
    t_start = time.time()

    outdir = Path(args.outdir).resolve()
    data_dir, plot_dir, model_dir = outdir / "data", outdir / "plots", outdir / "models"
    for d in (data_dir, plot_dir, model_dir):
        d.mkdir(parents=True, exist_ok=True)

    # ---- imports (torch env knobs must precede the torch import) ----------------------
    # Tutorial cell 6 sets CUBLAS_WORKSPACE_CONFIG=":4096:8"; that is a CUDA-only knob.
    # We set it anyway (harmless off-CUDA) and record the MPS situation in the manifest.
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    # MPS has no kernel for `aten::_pdist_forward`, which mofaflex's GP prior calls once at
    # init to set the kernel lengthscale from the inducing-point spread. Without this the fit
    # dies with NotImplementedError. Only that op falls back to CPU. Recorded as a deviation.
    os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

    import anndata as ad
    import torch

    # Tutorial: torch.use_deterministic_algorithms(True). On MPS several ops have no
    # deterministic implementation registered, which turns that into a hard error mid-fit,
    # so we use warn_only=True. Recorded as a deviation.
    torch.use_deterministic_algorithms(True, warn_only=True)

    import liana as li
    import mofaflex as mfl
    import scanpy as sc
    import squidpy as sq
    from liana.multi import lrdata_to_mudata

    from anndata._warnings import ImplicitModificationWarning
    warnings.filterwarnings("ignore", category=ImplicitModificationWarning)
    warnings.filterwarnings("ignore", message=".*pull_on_update.*")
    sc.set_figure_params(dpi=80, dpi_save=200, facecolor="white")

    if args.device != "auto":
        device = args.device
    elif torch.cuda.is_available():
        device = "cuda"
    elif torch.backends.mps.is_available():
        device = "mps"
    else:
        device = "cpu"
    print(f"Using device: {device}")

    # gpytorch/linear_operator run the variational-strategy Cholesky in float64 by default.
    # MPS has no float64, so on MPS we drop those two linalg dtypes to float32.
    # Recorded as a deviation (reduced precision in the GP inducing-point solve).
    linalg_float32 = False
    if device == "mps":
        from linear_operator import settings as _lo_settings
        _lo_settings._linalg_dtype_cholesky._set_value(torch.float32)
        _lo_settings._linalg_dtype_symeig._set_value(torch.float32)
        linalg_float32 = True

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    manifest: dict = {
        "script": str(Path(__file__).resolve()),
        "argv": sys.argv,
        "args": vars(args),
        "device": device,
        "linalg_dtype_float32": linalg_float32,
        "env": {k: os.environ.get(k) for k in ("CUBLAS_WORKSPACE_CONFIG", "PYTORCH_ENABLE_MPS_FALLBACK")},
        "versions": {
            "python": sys.version.split()[0],
            "liana": li.__version__,
            "mofaflex": mfl.__version__,
            "torch": torch.__version__,
            "scanpy": sc.__version__,
            "squidpy": sq.__version__,
            "anndata": ad.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "git_sha": git_sha(Path(__file__).resolve().parents[3]),
    }

    # ---------------------------------------------------------------- load (tutorial 9-14)
    print(f"\n[1] loading {args.lrdata}")
    lrdata = ad.read_h5ad(args.lrdata)
    manifest["input"] = {"path": str(Path(args.lrdata).resolve()), "shape_in": list(lrdata.shape)}

    dropped_obs = [c for c in args.drop_obs if c in lrdata.obs.columns]
    lrdata.obs = lrdata.obs.drop(columns=dropped_obs)
    dropped_uns = [k for k in args.drop_uns if k in lrdata.uns]
    for k in dropped_uns:
        del lrdata.uns[k]
    print(f"    dropped obs {dropped_obs}, uns {dropped_uns}")
    manifest["input"]["dropped_obs"] = dropped_obs
    manifest["input"]["dropped_uns"] = dropped_uns

    for k in args.drop_obs:
        assert k not in lrdata.obs.columns, f"leakage column {k} still present"

    # tutorial cell 10: MPS does not support float64. Our inflow matrix is also float64
    # (the tutorial's is float32 because its input h5ad was), so cast X too.
    lrdata.obsm[args.spatial_key] = np.asarray(lrdata.obsm[args.spatial_key]).astype("float32")
    manifest["input"]["X_dtype_in"] = str(lrdata.X.dtype)
    lrdata.X = lrdata.X.astype("float32")

    # ------------------------------------------------------ feature QC (tutorial 18-19)
    print("\n[2] feature QC on the inflow score")
    values = lrdata.X.data if hasattr(lrdata.X, "data") else np.asarray(lrdata.X).ravel()
    nzf = lrdata.var["nonzero_fraction"].values
    thr = args.nzf_threshold

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].hist(values, bins=100)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("inflow score (non-zero entries)")
    axes[0].set_ylabel("# cell x feature entries (log)")
    axes[0].set_title("Distribution of inflow values")
    pos = nzf[nzf > 0]
    axes[1].hist(nzf, bins=np.logspace(np.log10(pos.min()), 0, 60))
    axes[1].set_xscale("log")
    axes[1].set_xlabel("feature non-zero fraction")
    axes[1].set_ylabel("# features")
    axes[1].set_title("Per-feature sparsity")
    axes[1].axvline(thr, color="crimson", ls="--", lw=1, label=f"{thr * 100:.2f}% cutoff")
    axes[1].legend()
    fig.tight_layout()
    save_fig(fig, plot_dir / "qc_inflow_distributions.png")

    n_before_nzf = lrdata.n_vars

    if args.nzf_mode == "global":
        # tutorial cell 19, verbatim
        keep = lrdata.var["nonzero_fraction"].values > thr
        print(f"    nonzero_fraction > {thr}: {n_before_nzf} -> {int(keep.sum())} features")
    else:
        # ---- DEVIATION: reachability-normalised nonzero_fraction -------------------------
        # An inflow feature <sender>^<lig>^<rec> can only be non-zero for a cell that has a
        # cell of that SENDER inside the kernel support. So its nonzero_fraction is bounded
        # above by the sender's REACHABILITY, P(>=1 sender-s neighbour), which scales with the
        # sender's abundance and falls with its spatial clustering. A single GLOBAL threshold
        # therefore imposes a different effective stringency per sender: measured on this
        # dataset, surviving-feature count vs sender abundance is Spearman rho = 0.917, and
        # Lymphoid (2.54% reachable, max attainable nonzero_fraction 0.0043) cannot satisfy a
        # 1% global cut for arithmetic rather than biological reasons.
        #
        # Dividing by reachability asks the intended question -- "is this feature non-zero in
        # an appreciable share of the cells that COULD receive it" -- and is scale-free in
        # sender abundance.
        from sklearn.neighbors import radius_neighbors_graph

        R = args.bandwidth_support_radius or (args.bandwidth * np.sqrt(-2 * np.log(args.cutoff)))
        xy = np.asarray(lrdata.obsm[args.spatial_key], dtype=float)
        ct_arr = lrdata.obs[args.celltype_key].astype(str).values
        G = radius_neighbors_graph(xy, radius=R, mode="connectivity", include_self=True)
        senders = np.array([f.split(args.xy_sep)[0] for f in lrdata.var_names.astype(str)])

        reach = {s: float((G @ (ct_arr == s).astype(float) > 0).mean()) for s in set(senders)}
        reach_vec = np.array([reach[s] for s in senders])
        nzf_norm = lrdata.var["nonzero_fraction"].values / np.maximum(reach_vec, 1e-12)
        keep = nzf_norm > thr
        lrdata.var["reachability"] = reach_vec
        lrdata.var["nonzero_fraction_reach_norm"] = nzf_norm

        print(f"    reachability-normalised nonzero_fraction > {thr} "
              f"(support radius {R:.4f} um): {n_before_nzf} -> {int(keep.sum())} features")
        print(f"    {'sender':>16} {'reach':>7} {'kept':>6} {'/ total':>8}")
        for s in sorted(reach, key=lambda x: -reach[x]):
            m = senders == s
            print(f"    {s:>16} {reach[s]:7.3f} {int(keep[m].sum()):6d} {int(m.sum()):8d}")
        # Honest caveat, printed so it lands in the log: features from a low-reachability
        # sender are supported by proportionally fewer cells, so that view's factors are
        # estimated from a smaller effective sample. No absolute cell-count floor is applied
        # here (deliberately -- the change requested was the normalisation alone); read the
        # per-view R2 with the reachability column above in mind.
        low = [s for s in reach if reach[s] < 0.05]
        if low:
            print(f"    NOTE: reachability < 0.05 for {low} -- their features rest on a small "
                  f"cell base; check per-view R2 before interpreting those factors.")

    manifest["input"]["nzf_mode"] = args.nzf_mode
    lrdata = lrdata[:, keep].copy()

    # -------------------------------------------------------- SVI filter (tutorial 21)
    n_before_svi = lrdata.n_vars
    if args.skip_svi:
        svi_kept = n_before_svi
        print("    SVI filter SKIPPED by request")
    else:
        sq.gr.spatial_autocorr(lrdata, mode="moran", use_raw=False, seed=args.seed)
        moran = lrdata.uns["moranI"]
        svis = moran.index[(moran["pval_norm_fdr_bh"] <= args.moran_fdr) & (moran["I"] > args.moran_i)]
        lrdata = lrdata[:, svis].copy()
        svi_kept = lrdata.n_vars
        moran.to_csv(data_dir / "moranI_features.csv")
        print(f"    SVI filter: {n_before_svi} -> {svi_kept} features "
              f"({n_before_svi - svi_kept} dropped)")

    manifest["feature_qc"] = {
        "n_features_input": n_before_nzf,
        "nzf_threshold": thr,
        "n_after_nzf": n_before_svi,
        "n_after_svi": svi_kept,
        "svi_is_noop": bool(svi_kept == n_before_svi),
    }

    # --------------------------------------------------------- to MuData (tutorial 23)
    # MPS cannot hold sparse CSR tensors (`new_compressed_tensor` is unimplemented), and
    # mofaflex's collate_fn moves the batch to the device *before* densifying. Our inflow
    # matrix is sparse CSR; densify it here. Recorded as a deviation.
    densify = args.densify == "yes" or (args.densify == "auto" and device == "mps")
    from scipy import sparse as _sp
    if densify and _sp.issparse(lrdata.X):
        nbytes = lrdata.n_obs * lrdata.n_vars * 4
        print(f"    densifying X for the {device} backend ({nbytes / 1024**3:.2f} GB)")
        lrdata.X = np.asarray(lrdata.X.todense(), dtype="float32")
    manifest["densified"] = bool(densify)

    print("\n[3] lrdata -> MuData (one view per sender cell type)")
    obs_keys = [k for k in [args.celltype_key, args.region_key, args.sample_key] if k in lrdata.obs.columns]
    mdata = lrdata_to_mudata(
        lrdata, min_features=args.min_features, min_cells=args.min_cells, obs_keys=obs_keys, verbose=True
    )
    view_counts = {m: int(mdata.mod[m].n_vars) for m in mdata.mod}
    print(f"    views kept: {view_counts}")

    senders_in = pd.Series([v.split("^")[0] for v in lrdata.var_names]).value_counts().to_dict()
    dropped_views = sorted(set(senders_in) - set(view_counts))
    print(f"    views DROPPED (< min_features={args.min_features}): {dropped_views}")
    pd.DataFrame(
        {"sender": list(senders_in), "n_features_after_qc": [senders_in[s] for s in senders_in],
         "kept_as_view": [s in view_counts for s in senders_in]}
    ).sort_values("n_features_after_qc", ascending=False).to_csv(data_dir / "view_feature_counts.csv", index=False)

    manifest["mudata"] = {
        "n_views": len(view_counts), "view_feature_counts": view_counts,
        "dropped_views": dropped_views, "n_obs": int(mdata.n_obs),
        "features_per_sender_after_qc": {k: int(v) for k, v in senders_in.items()},
    }

    # make sure the GP prior can find the coordinates in every view
    for m in mdata.mod:
        mdata.mod[m].obsm[args.spatial_key] = lrdata[mdata.mod[m].obs_names].obsm[args.spatial_key].astype("float32")
    mdata.obsm[args.spatial_key] = lrdata[mdata.obs_names].obsm[args.spatial_key].astype("float32")

    # --------------------------------------------------------------- fit (tutorial 25)
    save_path = model_dir / f"{args.tag}.hdf5"
    if args.refit and save_path.exists():
        save_path.unlink()

    def build_model():
        return mfl.terms.MofaFlex(
            n_factors=args.n_factors,
            factor_prior=mfl.priors.GaussianProcess(
                covariates_mkey=args.spatial_key,
                independent_lengthscales=True,
                n_inducing=args.n_inducing,
            ),
            weight_prior=mfl.priors.Horseshoe(),
            nonnegative_weights=False,
        )

    fit_seconds = None
    if not save_path.exists():
        print(f"\n[4] fitting MOFA-Flex (n_factors={args.n_factors}, max_epochs={args.max_epochs}) -> {save_path}")
        model = build_model()
        t0 = time.time()
        model.fit(
            mdata,
            batch_size=args.batch_size,
            n_particles=args.n_particles,
            lr=args.lr,
            max_epochs=args.max_epochs,
            save_path=str(save_path),
            early_stopper_patience=args.patience,
            seed=args.seed,
            device=device,
        )
        fit_seconds = time.time() - t0
        print(f"    fit finished in {fit_seconds / 60:.1f} min")
        save_current(plot_dir / "data_overview.png")
    else:
        print(f"\n[4] reusing cached model at {save_path}")

    # optional determinism probe (tutorial cell 6 is CUDA-specific and does not apply here)
    determinism = None
    if args.determinism_check:
        print("\n[4b] determinism probe: two short fits with the same seed")
        ws = []
        probe_dir = model_dir / "determinism_probe"
        probe_dir.mkdir(exist_ok=True)
        for rep in range(2):
            m = build_model()
            # save_path=None makes mofaflex dump a timestamped .h5 into the CWD; give it a path.
            m.fit(mdata, batch_size=args.batch_size, n_particles=args.n_particles, lr=args.lr,
                  max_epochs=args.determinism_epochs,
                  save_path=str(probe_dir / f"{args.tag}_probe_rep{rep}.hdf5"),
                  early_stopper_patience=10**9,
                  seed=args.seed, device=device, plot_data_overview=False)
            ws.append(pd.concat(m.get_weights().values(), axis=1))
            plt.close("all")
        a, b = ws[0], ws[1].reindex(index=ws[0].index, columns=ws[0].columns)
        d = np.abs(a.values - b.values)
        determinism = {"epochs": args.determinism_epochs, "max_abs_weight_diff": float(np.nanmax(d)),
                       "bitwise_identical": bool(np.nanmax(d) == 0.0)}
        print(f"    {determinism}")
    manifest["determinism_probe"] = determinism

    # ------------------------------------------------------------- reload (tutorial 28)
    model = mfl.MOFAFLEX.load(str(save_path))

    # ---------------------------------------------------- diagnostics (tutorial 31-35)
    print("\n[5] diagnostics")
    save_gg(mfl.pl.factor_correlation(model, figsize=(6, 6)), plot_dir / "factor_correlation.png", 6, 6)
    save_gg(mfl.pl.variance_explained(model, figsize=(6, 6)), plot_dir / "variance_explained.png", 6, 6)
    save_gg(mfl.pl.variance_explained(model, group_by="view", figsize=(9, 6)),
            plot_dir / "variance_explained_by_view.png", 9, 6)

    r2_df = model.get_r2(type="term", term=None)
    r2_df.to_csv(data_dir / "r2_per_factor_view.csv", index=False)

    r2_by_factor = r2_df.groupby("component")["R2"]
    active = r2_by_factor.sum()[r2_by_factor.max() >= args.r2_floor].sort_values(ascending=False)
    factors = active.index.tolist()
    print(f"    active factors (R2 >= {args.r2_floor} in >=1 view): {len(factors)} / {args.n_factors}")
    # K IS A CEILING, NOT A SELECTION. MOFA-Flex does not choose the number of factors; it
    # fits exactly n_factors and we prune afterwards with the R2 floor. If nearly all of them
    # survive, the ceiling is binding and the "number of programs" is an artefact of the
    # argument, not a property of the data -- the same failure that produced the k_range
    # confound in run_nmf.py, which is why that runner warns on a boundary rank.
    _frac = len(factors) / max(args.n_factors, 1)
    manifest["n_factors_requested"] = int(args.n_factors)
    manifest["n_factors_active"] = int(len(factors))
    manifest["ceiling_binding"] = bool(_frac >= 0.9)
    if _frac >= 0.9:
        print(f"    *** WARNING: {len(factors)}/{args.n_factors} factors are active "
              f"({100*_frac:.0f}%). n_factors is a CEILING and it is BINDING -- the model would "
              f"very likely use more if given more, so this factor count must NOT be reported "
              f"as a discovered number. Re-fit with a larger --n-factors to find where the "
              f"active count saturates. ***")
    print(f"    {factors}")

    pd.DataFrame({"factor": active.index, "total_R2": active.values,
                  "max_view_R2": r2_by_factor.max().reindex(active.index).values}).to_csv(
        data_dir / "active_factors.csv", index=False)

    # variance explained per view (= per sender)
    view_r2 = r2_df.groupby("view")["R2"].sum().sort_values(ascending=False)
    view_r2.to_csv(data_dir / "r2_per_view.csv")
    print(f"    R2 per view (sender):\n{view_r2}")

    manifest["r2"] = {
        "n_active_factors": len(factors),
        "active_factors": factors,
        "total_r2_per_active_factor": {k: float(v) for k, v in active.items()},
        "r2_per_view": {k: float(v) for k, v in view_r2.items()},
    }

    if len(factors) == 0:
        print("\nNO ACTIVE FACTORS -- stopping before the downstream analysis.")
        manifest["wall_seconds"] = time.time() - t_start
        manifest["peak_rss_gb"] = peak_rss_gb()
        with open(outdir / "run_manifest.json", "w") as fh:
            json.dump(manifest, fh, indent=2, default=str)
        return

    # --------------------------------------------------- factor AnnData (tutorial 38)
    fac = model.get_factors()
    group_name = next(iter(fac))
    factor_adata = ad.AnnData(fac[group_name])[:, factors].copy()
    factor_adata.obsm["spatial"] = np.asarray(
        lrdata[factor_adata.obs_names].obsm[args.spatial_key]
    ).astype("float32")
    meta_keys = [k for k in [args.celltype_key, args.region_key, args.sample_key] if k in lrdata.obs.columns]
    shared = lrdata.obs_names.intersection(factor_adata.obs_names)
    for k in meta_keys:
        factor_adata.obs[k] = pd.Series(lrdata.obs.loc[shared, k].astype(str), index=shared).reindex(
            factor_adata.obs_names)
        factor_adata.obs[k] = factor_adata.obs[k].astype("category")

    factor_adata.to_df().to_csv(data_dir / "factor_scores.csv.gz", compression="gzip")

    # ----------------------------------------------------- top weights (tutorial 43)
    # NOTE on how to read this figure, because it is easy to misread
    # (mofaflex/pl/_plotting.py:1139-1169):
    #   * x is | Weight | -- the ABSOLUTE loading. The SIGN is carried by the glyph only,
    #     (+) for weight >= 0 and (-) for weight < 0, so a strongly negative feature sits
    #     just as far right as a strongly positive one.
    #   * within each facet the top n_features are taken by |weight| and sorted ascending,
    #     so the largest sits at the TOP of the panel.
    #   * facet_wrap uses scales="free", so EVERY PANEL HAS ITS OWN X-SCALE -- bar lengths
    #     are not comparable across factors.
    #   * these are RAW weights. They are dominated by high-prevalence adhesion pairs
    #     (NCAM1^NCAM1 occupies a top-5 slot in 17 of 20 factors, and more than one slot in
    #     11 of 20). ALARMIST divides its factor by each LRI's global prevalence before
    #     ranking (V* = V/(mean_LR+1)); there is no equivalent step here.
    _tw_h = max(8.0, 8.0 * args.top_weights_n / 5.0)   # 20 facets x n rows needs room
    save_gg(mfl.pl.top_weights(model, n_features=args.top_weights_n, factors=None,
                               figsize=(24, _tw_h)),
            plot_dir / "top_weights.png", 24, _tw_h)
    print(f"    top_weights.png: top {args.top_weights_n} features per factor")

    # ------------------------------------------- clustering / embedding (tutorial 45-52)
    print("\n[6] factor-space clustering")
    X = np.asarray(factor_adata.X, dtype="float64")
    z = (X - X.mean(0)) / X.std(0)
    factor_adata.obsm["X_clipped"] = np.clip(z, -5, 5)
    sc.pp.neighbors(factor_adata, use_rep="X_clipped", random_state=args.seed)
    sc.tl.umap(factor_adata, neighbors_key="neighbors", random_state=args.seed)
    sc.tl.leiden(factor_adata, resolution=args.leiden_resolution, flavor="igraph",
                 random_state=args.seed, n_iterations=2, directed=False)

    sc.pl.umap(factor_adata, color=["leiden", *meta_keys], size=8, wspace=0.8, show=False)
    save_current(plot_dir / "umap_leiden_annotations.png")

    sc.pl.embedding(factor_adata, basis="spatial", color=["leiden", *meta_keys], s=8, wspace=0.3, show=False)
    save_current(plot_dir / "spatial_leiden_annotations.png")

    focus_factors = factors[: args.n_focus_factors]
    print(f"    focus factors (top total R2): {focus_factors}")

    absmax = float(np.nanpercentile(np.abs(factor_adata[:, factors].X), 99))
    sc.pl.umap(factor_adata, color=focus_factors, cmap="RdBu_r", vcenter=0, vmin=-absmax, vmax=absmax,
               size=8, ncols=2, wspace=0.1, show=False)
    save_current(plot_dir / "umap_focus_factors.png")

    vals = np.hstack([np.asarray(factor_adata[:, f].X).flatten() for f in focus_factors])
    p_low, p_high = np.nanpercentile(vals, 1), np.nanpercentile(vals, 99)
    am = max(abs(p_low), abs(p_high))
    norm = colors.TwoSlopeNorm(vcenter=0, vmin=-am, vmax=am)
    sc.pl.embedding(factor_adata, basis="spatial", color=focus_factors, cmap="RdBu_r", norm=norm,
                    s=10, ncols=2, wspace=0.15, show=False)
    save_current(plot_dir / "spatial_focus_factors.png")

    # ---------------------------------------------- liana bridge (tutorial 55, 57, 59)
    print("\n[7] MOFA-Flex weights -> liana interaction table")
    keys = ["source", "ligand_complex", "receptor_complex"]
    loadings = li.ut.get_variable_loadings(
        loadings=model.get_weights(), variable_sep="^", var_names=keys
    )
    # NB: get_variable_loadings re-sorts rows by |Factor 1|. NEVER zip positionally against
    # anything else - always merge on `keys`.
    loadings.to_csv(data_dir / "mofaflex_loadings.csv", index=False)

    inflow_means = (
        lrdata.to_df()
        .groupby(lrdata.obs[args.celltype_key], observed=True)
        .mean()
        .T
        .reset_index(names="feature")
        .melt(id_vars="feature", var_name="target", value_name="inflow")
    )
    inflow_means[keys] = inflow_means["feature"].str.split("^", expand=True)
    inflow_means.to_csv(data_dir / "inflow_means_by_receiver.csv.gz", index=False, compression="gzip")

    def factor_interactions(factor, top_n=12):
        lo = loadings[[*keys, factor]].rename(columns={factor: "loading"})
        lo = lo.reindex(lo["loading"].abs().sort_values(ascending=False).index).head(top_n)
        return lo.merge(inflow_means[[*keys, "target", "inflow"]], on=keys, how="left")

    for f in factors:
        factor_interactions(f, top_n=25).to_csv(
            data_dir / f"factor_interactions_{f.replace(' ', '')}.csv", index=False)

    # dotplot (tutorial 57)
    top_feats = pd.concat(
        [loadings.reindex(loadings[f].abs().sort_values(ascending=False).index).head(8)[keys]
         for f in focus_factors]
    ).drop_duplicates()
    loadings_long = (
        loadings.merge(top_feats, on=keys)
        .melt(id_vars=keys, value_vars=focus_factors, var_name="target", value_name="loading")
    )
    loadings_long["abs_loading"] = loadings_long["loading"].abs()
    # li.pl.dotplot returns a PLOTNINE ggplot, so it must be saved with its own .save().
    # The first pass used save_current() and wrote a 4,377-byte all-white PNG -- the same
    # return-type trap that produced the blank connectivity.png in plot_liana_full.py.
    _gg = li.pl.dotplot(liana_res=loadings_long, colour="loading", size="abs_loading",
                        cmap="RdBu_r", figure_size=(18, 6), return_fig=True)
    save_gg(_gg, plot_dir / "dotplot_focus_factors.png", width=18, height=6)

    # circle plots (tutorial 59)
    # The tutorial draws these only for its two hand-picked focus factors. We emit one per
    # ACTIVE factor: with 19-20 active factors the focus pair is an arbitrary slice, and the
    # sender->receiver topology is the main thing these plots carry, so it should be available
    # for every program rather than only the two with the largest total R2.
    circle_targets = factors if args.circle_all_factors else focus_factors
    print(f"\n    circle plots for {len(circle_targets)} factor(s)"
          f"{' (ALL active)' if args.circle_all_factors else ' (focus only)'}")
    for f in circle_targets:
        try:
            lrdata.uns["liana_res"] = factor_interactions(f, top_n=args.circle_top_n)
            ax = li.pl.circle_plot(lrdata, groupby=args.celltype_key, source_key="source",
                                   target_key="target", score_key="inflow", pivot_mode="mean",
                                   figure_size=(6, 6))
            ax.set_title(f"{f}  (top {args.circle_top_n} interactions by |loading|)")
            save_current(plot_dir / f"circle_plot_{f.replace(' ', '')}.png")
        except Exception as e:
            plt.close("all")
            print(f"      circle_plot {f} FAILED: {type(e).__name__}: {e}")

    # ------------------------------------------------------ LR pairs of interest
    print("\n[8] tracing the ligand-receptor pairs of interest")
    loi_rows = []
    for pair in args.lr_of_interest:
        lig, rec = pair.split("^")
        sub = loadings[(loadings["ligand_complex"] == lig) & (loadings["receptor_complex"] == rec)]
        if sub.empty:
            print(f"    {pair}: ABSENT from the fitted feature set (removed by QC or view drop)")
            loi_rows.append({"ligand_complex": lig, "receptor_complex": rec, "source": None,
                             "status": "absent"})
            continue
        for _, r in sub.iterrows():
            rec_row = {"source": r["source"], "ligand_complex": lig, "receptor_complex": rec,
                       "status": "present"}
            for f in factors:
                rec_row[f] = float(r[f])
            fabs = {f: abs(float(r[f])) for f in factors}
            best = max(fabs, key=fabs.get)
            rec_row["top_factor"] = best
            rec_row["top_loading"] = float(r[best])
            # rank of this feature within the top factor, by |loading|
            col = loadings[best].abs().sort_values(ascending=False)
            rec_row["rank_in_top_factor"] = int(
                np.where(col.index == r.name)[0][0] + 1)
            rec_row["n_features"] = int(loadings.shape[0])
            loi_rows.append(rec_row)
    loi = pd.DataFrame(loi_rows)
    loi.to_csv(data_dir / "lr_of_interest_loadings.csv", index=False)
    print(loi[[c for c in ["source", "ligand_complex", "receptor_complex", "status",
                           "top_factor", "top_loading", "rank_in_top_factor"] if c in loi.columns]].to_string(index=False))

    # ---------------------------------------- grade association at the replicate unit
    print(f"\n[9] {args.region_key} association, aggregated to {args.sample_key}")
    grade_stats = None
    if args.sample_key in factor_adata.obs and args.region_key in factor_adata.obs:
        from scipy.stats import mannwhitneyu
        from statsmodels.stats.multitest import multipletests

        df = factor_adata.to_df()
        df[args.sample_key] = factor_adata.obs[args.sample_key].values
        df[args.region_key] = factor_adata.obs[args.region_key].values
        punch = df.groupby(args.sample_key, observed=True)[factors].mean()
        punch_grade = df.groupby(args.sample_key, observed=True)[args.region_key].agg(
            lambda s: s.value_counts().idxmax())
        punch.assign(**{args.region_key: punch_grade}).to_csv(data_dir / "factor_means_per_punch.csv")

        levels = sorted(punch_grade.unique())
        assert len(levels) == 2, f"expected 2 {args.region_key} levels, got {levels}"
        a_lv, b_lv = levels
        rows = []
        for f in factors:
            a = punch.loc[punch_grade == a_lv, f].values
            b = punch.loc[punch_grade == b_lv, f].values
            u, p = mannwhitneyu(a, b, alternative="two-sided")
            rows.append({"factor": f, f"n_{a_lv}": len(a), f"n_{b_lv}": len(b),
                         f"mean_{a_lv}": float(a.mean()), f"mean_{b_lv}": float(b.mean()),
                         "delta": float(b.mean() - a.mean()), "U": float(u), "p": float(p)})
        grade_stats = pd.DataFrame(rows)
        grade_stats["p_bh"] = multipletests(grade_stats["p"], method="fdr_bh")[1]
        grade_stats = grade_stats.sort_values("p")
        grade_stats.to_csv(data_dir / "factor_grade_punch_mannwhitney.csv", index=False)
        print(grade_stats.to_string(index=False))

        # bar plot, tutorial cell 63 analogue
        region_means = (
            factor_adata.to_df()[focus_factors]
            .groupby(factor_adata.obs[args.region_key], observed=True)
            .mean()
            .sort_values(focus_factors[0], ascending=False)
        )
        ax = region_means.plot(kind="bar", figsize=(6, 4))
        ax.set_ylabel("Mean factor score")
        ax.set_xlabel(args.region_key)
        ax.set_title(f"Focus factor activity by {args.region_key}")
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_current(plot_dir / "focus_factors_by_region.png")

        # per-punch strip plot for every active factor
        n = len(factors)
        ncol = min(5, n)
        nrow = int(np.ceil(n / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(3 * ncol, 3 * nrow), squeeze=False)
        for i, f in enumerate(factors):
            ax = axs[i // ncol][i % ncol]
            for j, lv in enumerate(levels):
                y = punch.loc[punch_grade == lv, f].values
                ax.scatter(np.full(len(y), j) + np.random.uniform(-0.08, 0.08, len(y)), y, s=22)
                ax.hlines(y.mean(), j - 0.2, j + 0.2, color="k")
            ax.set_xticks(range(len(levels)))
            ax.set_xticklabels(levels)
            pv = float(grade_stats.loc[grade_stats["factor"] == f, "p"].iloc[0])
            ax.set_title(f"{f}\np={pv:.3f}", fontsize=9)
        for i in range(n, nrow * ncol):
            axs[i // ncol][i % ncol].axis("off")
        fig.suptitle(f"Per-punch mean factor score by {args.region_key} (n={len(punch)} punches)", y=1.001)
        fig.tight_layout()
        save_fig(fig, plot_dir / "factor_by_punch_grade.png")

        manifest["grade_association"] = {
            "n_punches": int(len(punch)),
            "levels": {lv: int((punch_grade == lv).sum()) for lv in levels},
            "min_possible_two_sided_p": float(
                mannwhitneyu(np.arange(int((punch_grade == levels[0]).sum())),
                             np.arange(int((punch_grade == levels[0]).sum()),
                                       int(len(punch))), alternative="two-sided")[1]),
            "n_sig_bh_0.05": int((grade_stats["p_bh"] < 0.05).sum()),
        }

    # ------------------------------------------------------------------- manifest
    # A REPLOT (cached model reused, fit_seconds is None) must NOT clobber the fit provenance
    # written by the run that actually did the fitting. On 2026-08-04 a replot to regenerate two
    # blank figures overwrote fit_seconds / n_epochs / determinism_probe with nulls, leaving the
    # 70.5 min fit recoverable only from logs/. Carry those keys forward instead.
    wall = time.time() - t_start
    _prev = {}
    _mpath = outdir / "run_manifest.json"
    if _mpath.exists():
        try:
            _prev = json.load(open(_mpath))
        except Exception:
            _prev = {}
    FIT_KEYS = ("fit_seconds", "n_epochs", "determinism_probe", "peak_rss_gb", "wall_seconds")
    if fit_seconds is None and _prev.get("fit_seconds") is not None:
        for k in FIT_KEYS:
            if _prev.get(k) is not None:
                manifest[k] = _prev[k]
        manifest["last_replot"] = dict(
            wall_seconds=wall, peak_rss_gb=peak_rss_gb(),
            timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
            note="model reloaded from cache; figures/tables regenerated. Fit provenance above "
                 "belongs to the original fit and was carried forward, not remeasured.")
        print("  [manifest] replot: fit provenance carried forward from the previous manifest")
    else:
        manifest["fit_seconds"] = fit_seconds
        manifest["wall_seconds"] = wall
        manifest["peak_rss_gb"] = peak_rss_gb()
    manifest["model_path"] = str(save_path)
    manifest["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    with open(_mpath, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    print(f"\nDONE in {wall / 60:.1f} min, peak RSS {manifest['peak_rss_gb']:.1f} GB")
    print(f"manifest -> {outdir / 'run_manifest.json'}")


if __name__ == "__main__":
    main()
