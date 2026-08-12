#!/usr/bin/env python
"""LIANA+ MISTy (LR-MISTy) on a spatial slide, following the authors' misty.ipynb tutorial.

Tutorial reference (liana-py docs/notebooks/misty.ipynb, section "Ligand-Receptor MISTy"):

    misty = li.mt.lrMistyData(adata[:, hvg], bandwidth=200, set_diag=False,
                              cutoff=0.01, nz_threshold=0.1)
    misty(model=LinearModel, bypass_intra=True, verbose=True)

The *generic* section of the same notebook (a different model, different views) is the only
place `RandomForestModel` appears; it is NOT the LR-MISTy configuration.  This script therefore
runs **LinearModel as the primary** and exposes `--model rf` only as an explicitly-labelled
secondary / rate probe.

Everything is parameterised; nothing about the GBM slide is hardcoded.  Every argument default
below is either the tutorial's own value or the installed function's default -- the run manifest
records which, plus every deviation.

Verified against the installed liana 1.8.1 source before use:
  * `lrMistyData(adata, resource_name, resource, nz_threshold, use_raw, layer, spatial_key,
                 kernel, bandwidth, set_diag, cutoff, zoi, verbose)`
  * `MistyData.__call__(model, bypass_intra, predict_self, maskby, k_cv, alphas, seed,
                        inplace, verbose, **kwargs)`
  * `_make_view` for the intra view is called with a hardcoded `nz_threshold=0`
    (_misty_constructs.py ~L240), so every receptor in the resource survives as a target;
    only the *extra* (ligand) view is thinned by `--nz-threshold`.
  * `spatial_neighbors` is called by `lrMistyData` without `max_neighbours`, i.e. at its own
    default of 100 -- the connectivity graph is a 100-NN graph reweighted by the kernel.
  * the target loop in `MistyData.__call__` is a plain `for target in intra.var_names` with no
    shared state, so wall time is exactly linear in the number of targets.

Usage
-----
    python run_misty.py --h5ad data/xenium_mm_final_cell_id.h5ad \
        --cellchatdb data/LRdatabase/CellChatDBv2.0.human.csv \
        --outdir results/comparators/liana/GBM/misty --model linear
"""
from __future__ import annotations

import argparse
import json
import os
import resource as _resource
import subprocess
import sys
import time
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

import liana as li
from liana.method import MistyData, lrMistyData
from liana.method.sp import LinearModel, RandomForestModel

# --------------------------------------------------------------------------- args
p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--h5ad", required=True)
p.add_argument("--outdir", required=True)
p.add_argument("--cellchatdb", default=None,
               help="CSV with ligand/receptor columns. If omitted, LIANA's own resource is used.")
p.add_argument("--resource-name", default="consensus",
               help="LIANA resource name, only used when --cellchatdb is omitted.")
p.add_argument("--count-layer", default="counts",
               help="layer holding raw counts, or 'X' if .X already holds them")
p.add_argument("--cell-type-col", default="cell_type")
p.add_argument("--spatial-key", default="spatial")
p.add_argument("--group-col", default=None,
               help="obs column used only for reporting per-group cell counts (e.g. tma_id)")
p.add_argument("--subset-col", default=None, help="obs column to subset on")
p.add_argument("--subset-val", default=None, help="value of --subset-col to keep")

# --- tutorial LR-MISTy configuration -------------------------------------------
p.add_argument("--bandwidth", type=float, default=200.0, help="tutorial value (pkg default 100)")
p.add_argument("--cutoff", type=float, default=0.01, help="tutorial value (pkg default 0.1)")
p.add_argument("--nz-threshold", type=float, default=0.1, help="tutorial value = pkg default")
p.add_argument("--set-diag", type=int, default=0, help="tutorial value = pkg default (False)")
p.add_argument("--kernel", default="misty_rbf", help="pkg default, tutorial does not override")
p.add_argument("--zoi", type=float, default=0.0, help="pkg default")
p.add_argument("--model", choices=["linear", "rf"], default="linear",
               help="linear = the tutorial's LR-MISTy model; rf = the generic-section model")
p.add_argument("--bypass-intra", type=int, default=1, help="tutorial value for LR-MISTy")
p.add_argument("--k-cv", type=int, default=10, help="pkg default")
p.add_argument("--predict-self", type=int, default=0, help="pkg default")
p.add_argument("--seed", type=int, default=1337, help="pkg default (liana DefaultValues.seed)")
p.add_argument("--n-jobs", type=int, default=-1, help="only forwarded for --model rf")

# --- preprocessing -------------------------------------------------------------
p.add_argument("--qc", type=int, default=1,
               help="filter_cells(min_genes) + filter_genes(min_cells), the tutorial QC")
p.add_argument("--min-genes", type=int, default=10)
p.add_argument("--min-cells", type=int, default=3)
p.add_argument("--drop-underscore-genes", type=int, default=1,
               help="drop var_names containing '_' (LIANA's complex separator)")
p.add_argument("--hvg", type=int, default=0,
               help="tutorial pre-step; 0 = OFF. See DEVIATION 1 in the manifest.")
p.add_argument("--n-top-genes", type=int, default=1521, help="only used when --hvg 1")

# --- run control ---------------------------------------------------------------
p.add_argument("--targets", default="all",
               help="'all' (tutorial) or a comma-separated subset -- used for rate probes only")
p.add_argument("--report-genes", default="",
               help="comma-separated genes whose non-zero fraction to report explicitly")
p.add_argument("--report-subset-col", default=None,
               help="obs column for a secondary non-zero-fraction report (e.g. tma_id)")
p.add_argument("--report-subset-val", default=None)
p.add_argument("--focus-targets", default="",
               help="comma-separated targets to profile in detail in the manifest")
p.add_argument("--focus-predictors", default="",
               help="comma-separated predictors to look up for --focus-targets")
p.add_argument("--construct-only", type=int, default=0,
               help="build the views, record their shapes, and stop (no model fit)")
p.add_argument("--figures", type=int, default=1)
p.add_argument("--top-n-plot", type=int, default=30)
p.add_argument("--tag", default=None, help="subdirectory under --outdir; default = --model")
a = p.parse_args()

T0 = time.time()
LOG = []


def log(*m):
    s = f"[{time.strftime('%H:%M:%S')}] " + " ".join(str(x) for x in m)
    print(s, flush=True)
    LOG.append(s)


def rss_gb():  # macOS ru_maxrss is bytes
    return _resource.getrusage(_resource.RUSAGE_SELF).ru_maxrss / 1e9


TAG = a.tag or a.model
OUT = os.path.abspath(a.outdir)
RUN = os.path.join(OUT, TAG)
DATA = os.path.join(RUN, "data")
PLOTS = os.path.join(RUN, "plots")
for d in (DATA, PLOTS):
    os.makedirs(d, exist_ok=True)

DEV: list[dict] = []          # deviations from the tutorial, each with the number that forced it
MAN: dict = {"script": os.path.basename(__file__)}

# --------------------------------------------------------------------------- load
log("liana", li.__version__, "| python", sys.version.split()[0])
log("reading", a.h5ad)
adata = sc.read_h5ad(a.h5ad)
log(f"loaded {adata.shape[0]:,} x {adata.shape[1]:,}")
n0 = adata.shape

if a.count_layer != "X":
    if a.count_layer not in adata.layers:
        raise SystemExit(f"ERROR: layer '{a.count_layer}' absent; have {list(adata.layers)}")
    adata.X = adata.layers[a.count_layer]      # no .copy(): this layer is a DENSE int64
adata.layers.clear()                           # 100197 x 5119 int64 = 4.1 GB -- drop the alias
adata.raw = None
adata.uns.pop("log1p", None)
if not sp.issparse(adata.X):                   # sparsify immediately, memory not semantics
    adata.X = sp.csr_matrix(adata.X)
adata.X = adata.X.astype(np.float32)
# drop precomputed embeddings / graphs from the source file; they are irrelevant here and are
# duplicated by every AnnData .copy() downstream
for k in [k for k in adata.obsm if k != a.spatial_key]:
    del adata.obsm[k]
for k in list(adata.obsp):
    del adata.obsp[k]
for k in list(adata.varm):
    del adata.varm[k]
for k in list(adata.varp):
    del adata.varp[k]

if a.drop_underscore_genes:
    bad = [g for g in adata.var_names.astype(str) if "_" in g]
    if bad:
        adata = adata[:, [g for g in adata.var_names.astype(str) if "_" not in g]].copy()
        log(f"dropped {len(bad)} '_'-containing genes (LIANA complex separator) -> {adata.shape[1]}")
    MAN["n_underscore_genes_dropped"] = len(bad)

if a.subset_col:
    msk = adata.obs[a.subset_col].astype(str) == str(a.subset_val)
    adata = adata[msk].copy()
    log(f"subset {a.subset_col}=={a.subset_val} -> {adata.shape[0]:,} cells")

# ALARMIST outputs must never leak into a comparator
for c in ("motif", "patch_id"):
    if c in adata.obs:
        del adata.obs[c]

n_pre = adata.shape[0]
if a.qc:
    sc.pp.filter_cells(adata, min_genes=a.min_genes)
    sc.pp.filter_genes(adata, min_cells=a.min_cells)
    log(f"QC: filter_cells(min_genes={a.min_genes}) + filter_genes(min_cells={a.min_cells}) "
        f"-> dropped {n_pre - adata.shape[0]} cells, {adata.shape}")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
log(f"normalize_total(1e4) + log1p; X max {adata.X.max():.2f}")

if a.cell_type_col in adata.obs:
    adata.obs["cell_type"] = adata.obs[a.cell_type_col].astype(str).astype("category")
adata.obsm[a.spatial_key] = np.asarray(adata.obsm[a.spatial_key], dtype=float)
xy = adata.obsm[a.spatial_key]
log(f"spatial extent {np.ptp(xy[:, 0]):.0f} x {np.ptp(xy[:, 1]):.0f} (units of obsm)")

MAN.update({"h5ad": a.h5ad, "shape_loaded": list(n0), "shape_after_prep": list(adata.shape),
            "qc": bool(a.qc), "min_genes": a.min_genes, "min_cells": a.min_cells,
            "n_cells_dropped_by_qc": int(n_pre - adata.shape[0]),
            "subset": None if not a.subset_col else f"{a.subset_col}=={a.subset_val}"})
if a.group_col and a.group_col in adata.obs:
    MAN["group_sizes"] = adata.obs[a.group_col].astype(str).value_counts().sort_index().to_dict()
if "cell_type" in adata.obs:
    MAN["cell_type_counts"] = adata.obs["cell_type"].value_counts().to_dict()


# ------------------------------------------------------- non-zero fraction report
def nz_frac(ad, genes):
    """Non-zero fraction per gene, exactly as liana's `_get_props` computes it."""
    out = {}
    for g in genes:
        if g not in ad.var_names:
            out[g] = None
            continue
        col = ad[:, g].X
        out[g] = float(col.getnnz() / ad.n_obs) if hasattr(col, "getnnz") else \
            float((np.asarray(col).ravel() != 0).sum() / ad.n_obs)
    return out


report_genes = [g for g in a.report_genes.split(",") if g]
nz_rows = []
if report_genes:
    fr = nz_frac(adata, report_genes)
    log("non-zero fraction (whole object): " +
        ", ".join(f"{g}={('NA' if v is None else f'{100*v:.2f}%')}" for g, v in fr.items()))
    nz_rows += [{"scope": "run_object", "n_cells": int(adata.n_obs), "gene": g, "nz_frac": v}
                for g, v in fr.items()]
    MAN["nz_fraction_run_object"] = fr
    if a.report_subset_col and a.report_subset_col in adata.obs:
        m = adata.obs[a.report_subset_col].astype(str) == str(a.report_subset_val)
        sub = adata[m]
        fr2 = nz_frac(sub, report_genes)
        log(f"non-zero fraction ({a.report_subset_col}=={a.report_subset_val}, "
            f"n={sub.n_obs:,}): " +
            ", ".join(f"{g}={('NA' if v is None else f'{100*v:.2f}%')}" for g, v in fr2.items()))
        nz_rows += [{"scope": f"{a.report_subset_col}=={a.report_subset_val}",
                     "n_cells": int(sub.n_obs), "gene": g, "nz_frac": v} for g, v in fr2.items()]
        MAN[f"nz_fraction_{a.report_subset_col}_{a.report_subset_val}"] = fr2
        if "cell_type" in sub.obs:
            MAN[f"cell_type_counts_{a.report_subset_col}_{a.report_subset_val}"] = \
                sub.obs["cell_type"].value_counts().to_dict()
    pd.DataFrame(nz_rows).to_csv(os.path.join(DATA, "nz_fractions.csv"), index=False)

# --------------------------------------------------------------------- resource
if a.cellchatdb:
    db = pd.read_csv(a.cellchatdb)
    lr = db[["ligand", "receptor"]].astype(str).drop_duplicates().reset_index(drop=True)
    log(f"resource {os.path.basename(a.cellchatdb)}: {lr.shape[0]} unique ligand-receptor pairs")
    MAN.update({"resource": os.path.basename(a.cellchatdb), "resource_pairs": int(lr.shape[0])})
    res_kw = {"resource": lr}
else:
    log(f"resource: LIANA '{a.resource_name}'")
    MAN.update({"resource": f"liana:{a.resource_name}", "resource_pairs": None})
    res_kw = {"resource_name": a.resource_name}

# ---------------------------------------------------------------- HVG (deviation)
view_input = adata
if a.hvg:
    sc.pp.highly_variable_genes(adata, n_top_genes=a.n_top_genes)
    view_input = adata[:, adata.var.highly_variable]
    log(f"HVG pre-step ON: {view_input.shape[1]} genes")
else:
    DEV.append({
        "id": 1,
        "deviation": "tutorial's `adata[:, hvg]` pre-step NOT applied; the full panel is used",
        "tutorial_text": "the notebook frames HVG selection as 'for the sake of computational "
                         "speed' on a genome-wide Visium slide",
        "number_that_forced_it": "this is a 5,119-gene TARGETED panel, not a genome-wide assay, "
                                 "so HVG selection removes signal rather than noise. Measured "
                                 "with --hvg 1 --n-top-genes 1521 --construct-only 1: FULL SLIDE "
                                 "intra 100190x82 / extra 100190x13, against 100190x382 / "
                                 "100190x37 without it -- 78% of the receptor targets and 65% of "
                                 "the ligand predictors are discarded. On a single 5,363-cell "
                                 "punch the same pre-step collapses the extra view to ONE "
                                 "predictor (intra 5363x54, extra 5363x1), which makes the model "
                                 "vacuous. The speed argument does not apply here either: the "
                                 "full-panel LinearModel fit takes 2.7 min.",
    })
    log("DEVIATION 1: HVG pre-step deliberately OFF (would leave a 1-ligand extra view)")

# ----------------------------------------------------------------- build views
log(f"lrMistyData(bandwidth={a.bandwidth}, set_diag={bool(a.set_diag)}, cutoff={a.cutoff}, "
    f"nz_threshold={a.nz_threshold}, kernel={a.kernel})  [all tutorial values]")
t = time.perf_counter()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    misty = lrMistyData(view_input, nz_threshold=a.nz_threshold, bandwidth=a.bandwidth,
                        set_diag=bool(a.set_diag), cutoff=a.cutoff, kernel=a.kernel,
                        zoi=a.zoi, spatial_key=a.spatial_key, verbose=True, **res_kw)
construct_s = time.perf_counter() - t
conn = misty.mod["extra"].obsp[f"{a.spatial_key}_connectivities"]
deg = np.asarray((conn > 0).sum(axis=1)).ravel()
log(f"views built in {construct_s:.1f}s | intra {tuple(misty.mod['intra'].shape)} (receptors) | "
    f"extra {tuple(misty.mod['extra'].shape)} (ligands)")
log(f"connectivity: {conn.nnz:,} nnz | degree median {np.median(deg):.0f} max {deg.max()} "
    f"(spatial_neighbors max_neighbours default = 100, not overridden by lrMistyData)")
MAN.update({"construct_s": construct_s,
            "intra_shape": list(misty.mod["intra"].shape),
            "extra_shape": list(misty.mod["extra"].shape),
            "conn_nnz": int(conn.nnz),
            "conn_degree_median": float(np.median(deg)),
            "conn_degree_max": int(deg.max()),
            "conn_max_neighbours_cap": 100,
            "conn_cap_binds_frac": float((deg >= 100).mean())})

pd.DataFrame({"view": ["intra", "extra"],
              "role": ["targets (receptors)", "predictors (ligands, spatially weighted)"],
              "nz_threshold": [0.0, a.nz_threshold],
              "n_features": [misty.mod["intra"].n_vars, misty.mod["extra"].n_vars],
              "features": [";".join(misty.mod["intra"].var_names),
                           ";".join(misty.mod["extra"].var_names)]}
             ).to_csv(os.path.join(DATA, "view_features.csv"), index=False)

if a.construct_only:
    MAN.update({"construct_only": True, "hvg": bool(a.hvg), "n_top_genes": a.n_top_genes,
                "wall_min": (time.time() - T0) / 60.0, "peak_rss_gb": rss_gb(),
                "liana_version": li.__version__, "deviations": DEV})
    with open(os.path.join(RUN, "run_manifest.json"), "w") as fh:
        json.dump(MAN, fh, indent=2, default=str)
    with open(os.path.join(RUN, "run.log"), "w") as fh:
        fh.write("\n".join(LOG) + "\n")
    log("construct-only: stopping before the fit")
    print("MANIFEST " + json.dumps({"hvg": bool(a.hvg), "n_genes_in": int(view_input.shape[1]),
                                    "intra_shape": MAN["intra_shape"],
                                    "extra_shape": MAN["extra_shape"]}))
    sys.exit(0)

# optional target subsetting (rate probes only)
n_targets_full = misty.mod["intra"].n_vars
if a.targets != "all":
    tg = [t_ for t_ in a.targets.split(",") if t_]
    missing = [t_ for t_ in tg if t_ not in misty.mod["intra"].var_names]
    if missing:
        raise SystemExit(f"ERROR: targets absent from intra view: {missing}")
    intra = misty.mod["intra"][:, tg].copy()
    misty = MistyData(data={"intra": intra, "extra": misty.mod["extra"].copy()},
                      obs=intra.obs, spatial_key=a.spatial_key)
    log(f"TARGET SUBSET (rate probe): {len(tg)} of {n_targets_full} targets")
MAN.update({"n_targets_available": int(n_targets_full),
            "n_targets_run": int(misty.mod["intra"].n_vars),
            "targets_arg": a.targets})

# ------------------------------------------------------------------------- fit
Model = LinearModel if a.model == "linear" else RandomForestModel
kw = {} if a.model == "linear" else {"n_jobs": a.n_jobs}
log(f"misty(model={Model.__name__}, bypass_intra={bool(a.bypass_intra)}, k_cv={a.k_cv}, "
    f"seed={a.seed}{'' if not kw else ', ' + str(kw)}) on {misty.mod['intra'].n_vars} targets")
t = time.perf_counter()
misty(model=Model, bypass_intra=bool(a.bypass_intra), predict_self=bool(a.predict_self),
      k_cv=a.k_cv, seed=a.seed, verbose=True, **kw)
run_s = time.perf_counter() - t
log(f"fit done in {run_s/60:.2f} min ({run_s/misty.mod['intra'].n_vars:.4f} s/target)")

tm = misty.uns["target_metrics"].copy()
inter = misty.uns["interactions"].copy()
tm.to_csv(os.path.join(DATA, "target_metrics.csv"), index=False)
inter.to_csv(os.path.join(DATA, "interactions.csv"), index=False)
view_cols = [c for c in tm.columns if c in ("intra", "extra", "juxta", "para")]
tm[["target", *view_cols]].melt(id_vars="target", var_name="view", value_name="contribution") \
    .to_csv(os.path.join(DATA, "contributions.csv"), index=False)

MAN.update({"model": Model.__name__, "bypass_intra": bool(a.bypass_intra),
            "k_cv": a.k_cv, "seed": a.seed, "predict_self": bool(a.predict_self),
            "n_jobs": (a.n_jobs if a.model == "rf" else "n/a (LinearModel default -1)"),
            "run_s": run_s, "run_min": run_s / 60.0,
            "s_per_target": run_s / misty.mod["intra"].n_vars,
            "peak_rss_gb": rss_gb(),
            "target_metrics_shape": list(tm.shape), "interactions_shape": list(inter.shape),
            "target_metrics_cols": list(tm.columns), "interactions_cols": list(inter.columns)})

if bool(a.bypass_intra):
    DEV.append({
        "id": 0,
        "deviation": "none -- but a required reading note",
        "tutorial_text": "LR-MISTy is run with bypass_intra=True",
        "number_that_forced_it": "with bypass_intra=True liana sets intra_R2 := 0 "
                                 "(_Misty.py `intra_r2 = ... if not bypass_intra else 0`), so "
                                 "gain_R2 == multi_R2 == the cross-validated R2 of the ligand "
                                 "extra view ALONE. It is not a gain over an intrinsic model. "
                                 "Also, with a single view the meta-model coefficient is "
                                 "normalised to itself, so the 'extra' contribution is 1.0 for "
                                 "every target by construction.",
    })

# ------------------------------------------------------------------- summarise
g = tm["gain_R2"]
n_nan = int(g.isna().sum())
zero_exact = int((g == 0).sum())
summary = {
    "n_targets": int(tm.shape[0]),
    "gain_R2_nan": n_nan,
    "gain_R2_exactly_zero": zero_exact,
    "gain_R2_min": float(np.nanmin(g)) if n_nan < len(g) else None,
    "gain_R2_median": float(np.nanmedian(g)) if n_nan < len(g) else None,
    "gain_R2_mean": float(np.nanmean(g)) if n_nan < len(g) else None,
    "gain_R2_q75": float(np.nanquantile(g, 0.75)) if n_nan < len(g) else None,
    "gain_R2_q90": float(np.nanquantile(g, 0.90)) if n_nan < len(g) else None,
    "gain_R2_q99": float(np.nanquantile(g, 0.99)) if n_nan < len(g) else None,
    "gain_R2_max": float(np.nanmax(g)) if n_nan < len(g) else None,
    "n_gain_gt_0.001": int((g > 0.001).sum()),
    "n_gain_gt_0.01": int((g > 0.01).sum()),
    "n_gain_gt_0.05": int((g > 0.05).sum()),
    "n_gain_gt_0.1": int((g > 0.1).sum()),
}
ranked = tm.sort_values("gain_R2", ascending=False, na_position="last").reset_index(drop=True)
ranked["rank"] = np.arange(1, len(ranked) + 1)
ranked.to_csv(os.path.join(DATA, "target_metrics_ranked.csv"), index=False)
summary["top20"] = ranked.head(20)[["rank", "target", "gain_R2", "multi_R2"]] \
    .to_dict(orient="records")
MAN["target_metrics_summary"] = summary
log("gain_R2: n=%d  median=%.5f  max=%.5f  >0.01: %d  ==0: %d  NaN: %d" %
    (summary["n_targets"], summary["gain_R2_median"] or float("nan"),
     summary["gain_R2_max"] or float("nan"), summary["n_gain_gt_0.01"], zero_exact, n_nan))
log("top10 by gain_R2: " + ", ".join(f"{r['target']}={r['gain_R2']:.4f}"
                                     for r in summary["top20"][:10]))

# per-target predictor importances for the focus targets
focus_t = [x for x in a.focus_targets.split(",") if x]
focus_p = [x for x in a.focus_predictors.split(",") if x]
focus_out = {}
for tgt in focus_t:
    row = ranked[ranked["target"] == tgt]
    if row.empty:
        focus_out[tgt] = {"present": False,
                          "reason": "not a receptor in the resource-filtered intra view"}
        log(f"FOCUS {tgt}: absent from intra view")
        continue
    r = row.iloc[0]
    sub = inter[(inter["target"] == tgt) & (inter["view"] == "extra")] \
        .dropna(subset=["importances"]).copy()
    sub["abs_importance"] = sub["importances"].abs()
    sub = sub.sort_values("abs_importance", ascending=False).reset_index(drop=True)
    sub["rank"] = np.arange(1, len(sub) + 1)
    sub.to_csv(os.path.join(DATA, f"interactions_{tgt}.csv"), index=False)
    d = {"present": True, "gain_R2": float(r["gain_R2"]), "multi_R2": float(r["multi_R2"]),
         "rank_by_gain_R2": int(r["rank"]), "of_n_targets": int(len(ranked)),
         "n_predictors": int(len(sub)),
         "top5_predictors": sub.head(5)[["predictor", "importances"]].to_dict(orient="records")}
    for pr in focus_p:
        q = sub[sub["predictor"] == pr]
        d[f"predictor_{pr}"] = ({"present": False} if q.empty else
                                {"present": True,
                                 "importance_tvalue": float(q.iloc[0]["importances"]),
                                 "rank_by_abs_importance": int(q.iloc[0]["rank"]),
                                 "of_n_predictors": int(len(sub))})
    focus_out[tgt] = d
    log(f"FOCUS {tgt}: gain_R2={d['gain_R2']:.5f} rank {d['rank_by_gain_R2']}/{d['of_n_targets']}"
        f" | top predictors " + ", ".join(f"{x['predictor']}({x['importances']:.1f})"
                                          for x in d["top5_predictors"]))
    for pr in focus_p:
        v = d[f"predictor_{pr}"]
        if v["present"]:
            log(f"    predictor {pr}: t={v['importance_tvalue']:.2f}, "
                f"rank {v['rank_by_abs_importance']}/{v['of_n_predictors']}")
        else:
            log(f"    predictor {pr}: ABSENT from the extra view")
MAN["focus"] = focus_out

# ---------------------------------------------------------------------- figures
figs = []


def _p9save(pl, name):
    fp = os.path.join(PLOTS, name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pl.save(fp, dpi=200, verbose=False)
    figs.append(name)
    log("  figure:", name)


if a.figures and tm.shape[0] > 1:
    try:
        _p9save(li.pl.target_metrics(misty=misty, stat="gain_R2", top_n=a.top_n_plot,
                                     figure_size=(10, 4), return_fig=True),
                "target_metrics_gain_R2_top.png")
    except Exception as e:
        log(f"  FAILED target_metrics plot: {type(e).__name__}: {e}")
    try:
        _p9save(li.pl.target_metrics(misty=misty, stat="multi_R2", top_n=a.top_n_plot,
                                     figure_size=(10, 4), return_fig=True),
                "target_metrics_multi_R2_top.png")
    except Exception as e:
        log(f"  FAILED target_metrics multi_R2 plot: {type(e).__name__}: {e}")
    try:
        top_names = list(ranked.head(a.top_n_plot)["target"])
        _p9save(li.pl.contributions(target_metrics=tm[tm["target"].isin(top_names)],
                                    view_names=view_cols, figure_size=(10, 4), return_fig=True),
                "contributions_top.png")
    except Exception as e:
        log(f"  FAILED contributions plot: {type(e).__name__}: {e}")
    try:
        _p9save(li.pl.interactions(misty=misty, view="extra", top_n=200, figure_size=(10, 8),
                                   return_fig=True),
                "interactions_extra_top200.png")
    except Exception as e:
        log(f"  FAILED interactions plot: {type(e).__name__}: {e}")

    # honest distributional view -- the headline of this run
    try:
        fig, ax = plt.subplots(figsize=(5, 3.4), dpi=200)
        ax.hist(g.dropna(), bins=60, color="#4C72B0")
        ax.axvline(0.01, color="crimson", ls="--", lw=1,
                   label=f"0.01  (n={summary['n_gain_gt_0.01']}/{summary['n_targets']})")
        ax.set_xlabel("gain_R2 (= extra-view R$^2$; intra bypassed)")
        ax.set_ylabel("targets")
        ax.set_title(f"LR-MISTy {Model.__name__}, {summary['n_targets']} receptor targets")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(PLOTS, "gain_R2_hist.png"))
        plt.close(fig)
        figs.append("gain_R2_hist.png")
        log("  figure: gain_R2_hist.png")
    except Exception as e:
        log(f"  FAILED hist: {type(e).__name__}: {e}")

    for tgt, d in focus_out.items():
        if not d.get("present"):
            continue
        try:
            sub = pd.read_csv(os.path.join(DATA, f"interactions_{tgt}.csv")).head(15)[::-1]
            fig, ax = plt.subplots(figsize=(4.6, 4.2), dpi=200)
            cols = ["#C44E52" if v < 0 else "#4C72B0" for v in sub["importances"]]
            ax.barh(sub["predictor"], sub["importances"], color=cols)
            ax.set_xlabel("importance (OLS t-value)")
            ax.set_title(f"{tgt}  gain_R2={d['gain_R2']:.4f} "
                         f"(rank {d['rank_by_gain_R2']}/{d['of_n_targets']})", fontsize=9)
            fig.tight_layout()
            fig.savefig(os.path.join(PLOTS, f"importances_{tgt}.png"))
            plt.close(fig)
            figs.append(f"importances_{tgt}.png")
            log(f"  figure: importances_{tgt}.png")
        except Exception as e:
            log(f"  FAILED importances_{tgt}: {type(e).__name__}: {e}")

MAN["figures"] = figs
MAN["n_figures"] = len(figs)

# --------------------------------------------------------------------- manifest
try:
    MAN["git_sha"] = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                    capture_output=True, text=True,
                                    cwd=os.path.dirname(os.path.abspath(__file__))
                                    ).stdout.strip()
except Exception:
    MAN["git_sha"] = None
MAN.update({
    "liana_version": li.__version__,
    "python": sys.version.split()[0],
    "tutorial_reference": "liana-py docs/notebooks/misty.ipynb, 'Ligand-Receptor MISTy' section",
    "tutorial_call": "lrMistyData(adata[:, hvg], bandwidth=200, set_diag=False, cutoff=0.01, "
                     "nz_threshold=0.1); misty(bypass_intra=True, model=LinearModel, verbose=True)",
    "params_provenance": {
        "bandwidth": f"{a.bandwidth} -- TUTORIAL value (package default is 100). NOT the "
                     f"13.1454 um bandwidth used by this repo's bivariate/inflow LIANA runs.",
        "cutoff": f"{a.cutoff} -- tutorial value (package default 0.1)",
        "nz_threshold": f"{a.nz_threshold} -- tutorial value, equals the package default; "
                        f"applies to the EXTRA (ligand) view only, the intra view is built "
                        f"with a hardcoded nz_threshold=0",
        "set_diag": f"{bool(a.set_diag)} -- tutorial value = package default",
        "kernel": f"{a.kernel} -- package default; tutorial does not override it",
        "max_neighbours": "100 -- spatial_neighbors default; lrMistyData does not expose it",
        "k_cv": f"{a.k_cv} -- package default",
        "seed": f"{a.seed} -- package default (liana DefaultValues.seed)",
    },
    "deviations": DEV,
    "wall_min": (time.time() - T0) / 60.0,
    "peak_rss_gb": rss_gb(),
})
with open(os.path.join(RUN, "run_manifest.json"), "w") as fh:
    json.dump(MAN, fh, indent=2, default=str)
with open(os.path.join(RUN, "run.log"), "w") as fh:
    fh.write("\n".join(LOG) + "\n")
log(f"DONE  wall {MAN['wall_min']:.2f} min  peak RSS {MAN['peak_rss_gb']:.2f} GB -> {RUN}")
print("MANIFEST " + json.dumps({k: MAN[k] for k in
                                ("n_targets_run", "run_min", "s_per_target", "peak_rss_gb",
                                 "wall_min")}))
