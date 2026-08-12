#!/usr/bin/env python
"""LIANA+ spatially-informed bivariate LR scoring, following the authors' bivariate tutorial.

Tutorial: /Users/jiayifan/tansey_lab/liana-py/docs/notebooks/bivariate.ipynb
Call contract: scripts/comparators/liana/NOTES.md   Deviations: DEVIATIONS.md

Runs WHOLE-SLIDE (no core split) -- verified safe: the kernel support here is 28.2 um while the
minimum inter-core cell-cell distance is 222.9 um, giving zero cross-core pairs.

Bandwidth is set by equal-area correspondence to an s x s patch (s = 50 um):
    support radius = s / sqrt(pi)                    = 28.2095 um
    bandwidth      = support / sqrt(-2 ln cutoff)    = s / 3.804 = 13.1454 um

Usage: python run_liana.py --h5ad PATH --out-dir DIR [options]
"""
import argparse, json, os, subprocess, time
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import liana as li

p = argparse.ArgumentParser()
p.add_argument("--h5ad", required=True)
p.add_argument("--out-dir", required=True)
p.add_argument("--cell-type-col", default="cell_type")
p.add_argument("--count-layer", default="counts")
p.add_argument("--db", default="consensus",
               help="'consensus' = LIANA's own resource; else path to CellChatDBv2.0.human.csv")
p.add_argument("--bandwidth", type=float, default=13.1454,
               help="gaussian kernel bandwidth in COORDINATE UNITS. Set by equal-area correspondence to an s x s patch: support radius = s/sqrt(pi) = 28.2 um for s=50, and bandwidth = support / sqrt(-2 ln cutoff) = s/3.804 = 13.1454 um. (The earlier 18.75, derived from the bivariate tutorial's ~6-neighbour \"first ring\" rule, is VOID: first ring is a hexagonal-lattice concept with no referent on irregularly packed single cells.)")
p.add_argument("--cutoff", type=float, default=0.1)
p.add_argument("--kernel", default="gaussian")
p.add_argument("--set-diag", action=argparse.BooleanOptionalAction, default=True)
p.add_argument("--max-neighbours", type=int, default=100)
p.add_argument("--local-name", default="cosine")
p.add_argument("--global-name", default="morans")
p.add_argument("--n-perms", type=int, default=100)
p.add_argument("--nz-prop", type=float, default=0.02,
               help="tutorial uses 0.2 = 20%% of Visium SPOTS; the binomial spot->cell "
                    "equivalent is ~1.5%% of cells, so 0.02 (see DEVIATIONS.md)")
p.add_argument("--requested-lrs", default="GRN^SORT1,ANXA1^FPR1")
p.add_argument("--n-top-plots", type=int, default=6)
p.add_argument("--seed", type=int, default=0)
p.add_argument("--punch-col", default="tma_id",
               help="obs column with the TMA punch/core ID; written to cell_meta.csv "
                    "so the cross-punch interaction filter can run before NMF")
p.add_argument("--plots-only", action="store_true", default=False,
               help="skip scoring entirely: rebuild lrdata from <out-dir>/data/ and redraw the "
                    "figures. Nothing under data/ is written and run_manifest.json is NOT "
                    "touched; a separate plot_manifest.json records the replot.")
a = p.parse_args()
np.random.seed(a.seed)

OUT = a.out_dir
PLOTS, REQ, DATA = (os.path.join(OUT, "plots"), os.path.join(OUT, "plots", "requested"),
                    os.path.join(OUT, "data"))
for d in (PLOTS, REQ, DATA): os.makedirs(d, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
SAVED, BLANK = [], []


def _is_blank(fig):
    """Refuse to write a figure with no drawn content. plotnine-backed liana plotters leave the
    matplotlib canvas empty, which used to yield silent blank PNGs elsewhere in this suite."""
    if fig is None or not fig.axes:
        return True
    return not any(ax.has_data() or ax.get_images() or ax.collections or ax.patches or ax.texts
                   for ax in fig.axes)


def savefig(nm, sub=PLOTS):
    fig = plt.gcf() if plt.get_fignums() else None
    if _is_blank(fig):
        plt.close("all"); BLANK.append(nm); log(f"  BLANK, not written: {nm}"); return False
    fig.savefig(os.path.join(sub, nm + ".png"), dpi=200, bbox_inches="tight")
    plt.close("all"); SAVED.append(nm); return True

def try_plot(fn, nm, sub=PLOTS):
    try:
        if fn() is not False and savefig(nm, sub): log(f"  plot ok: {nm}")
    except Exception as e: plt.close("all"); log(f"  plot FAILED {nm}: {type(e).__name__}: {e}")

# ------------------------------------------------------------------ load
log("reading", a.h5ad)
adata = sc.read_h5ad(a.h5ad)
log(f"loaded {adata.shape[0]} x {adata.shape[1]}")

# Tutorial recipe: layers['counts'] = X.copy(); normalize_total(target_sum=1e4); log1p.
# Our .X is already log-normalised by an unknown recipe, so we start from genuine raw counts
# and reproduce the tutorial's two steps exactly.
if a.count_layer != "X":
    if a.count_layer not in adata.layers:
        raise SystemExit(f"ERROR: layer '{a.count_layer}' absent; have {list(adata.layers)}")
    adata.X = adata.layers[a.count_layer].copy()
adata.layers.clear(); adata.raw = None
adata.uns.pop("log1p", None)          # stale key -> spurious "already log-transformed" warning
# LIANA uses '_' as its heteromeric-subunit separator and warns about gene names containing it.
# On the Xenium 5K panel the offenders are the Intergenic_Region_* genomic CONTROL probes --
# not real genes, absent from every LR resource. Dropped BEFORE normalisation so library size
# reflects the real panel (same treatment as the stLearn run).
bad = [g for g in adata.var_names.astype(str) if "_" in g]
if bad:
    adata = adata[:, [g for g in adata.var_names.astype(str) if "_" not in g]].copy()
    log(f"dropped {len(bad)} control-probe genes containing '_' -> {adata.shape[1]} genes")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
log(f"normalize_total(1e4) + log1p (tutorial recipe); X max = {adata.X.max():.2f}")

if a.cell_type_col not in adata.obs:
    raise SystemExit(f"ERROR: '{a.cell_type_col}' not in obs")
adata.obs["cell_type"] = adata.obs[a.cell_type_col].astype(str).astype("category")
adata.obsm["spatial"] = np.asarray(adata.obsm["spatial"], dtype=float)
xy = adata.obsm["spatial"]
log(f"extent {np.ptp(xy[:,0]):.0f} x {np.ptp(xy[:,1]):.0f} um, whole slide (no core split)")

# ------------------------------------------------- step 1: spatial neighbors
R = a.bandwidth * np.sqrt(-2 * np.log(a.cutoff))
log(f"spatial_neighbors(bandwidth={a.bandwidth}, cutoff={a.cutoff}, kernel={a.kernel}) "
    f"-> effective radius {R:.1f} um")
li.ut.spatial_neighbors(adata, bandwidth=a.bandwidth, cutoff=a.cutoff, kernel=a.kernel,
                        set_diag=a.set_diag, max_neighbours=a.max_neighbours)
conn = adata.obsp["spatial_connectivities"]
deg = np.asarray((conn > 0).sum(axis=1)).ravel()
log(f"  connectivities: {conn.nnz:,} nonzeros | neighbours per cell median {np.median(deg):.0f}, "
    f"max {deg.max()} | max_neighbours cap binds for {100*(deg >= a.max_neighbours).mean():.4f}%")

# ------------------------------------------------------- step 2: LR resource
res_fingerprint = None
if a.db == "consensus":
    resource, res_label, res_kw = None, "LIANA consensus", dict(resource_name="consensus")
    res_n = int(li.rs.select_resource("consensus").shape[0])
    log(f"resource: {res_label} -> {res_n} rows")
else:
    db = pd.read_csv(a.db)
    # LIANA's resource is a DataFrame with ['ligand','receptor']; heteromeric subunits are
    # '_'-joined, which is exactly our CellChatDB v2 encoding -> direct, LOSSLESS handover.
    # LIANA names a pair '<ligand>^<receptor>'.
    resource = (db[["ligand", "receptor"]].astype(str)
                  .drop_duplicates().reset_index(drop=True))
    res_label, res_kw = os.path.basename(a.db), dict(resource=resource, resource_name=None)
    res_n = int(resource.shape[0])
    # Subunit ORDER discriminates the 2026-07-28 re-export from CellChatDBv2.0.human.old.csv, and
    # the two are NOT join-compatible (~1,120 of 3,218 LR keys differ). Record which one this run
    # consumed so a downstream join against ALARMIST results cannot silently mismatch.
    _key = (resource["ligand"] + "^" + resource["receptor"])
    res_fingerprint = {"path": os.path.abspath(a.db), "rows": int(db.shape[0]),
                       "unique_pairs": res_n,
                       "has_TGFBR1_TGFBR2": bool(_key.str.contains("TGFBR1_TGFBR2").any()),
                       "has_TGFBR2_TGFBR1": bool(_key.str.contains("TGFBR2_TGFBR1").any())}
    log(f"resource: {res_label} -> {res_n} unique LR pairs | fingerprint {res_fingerprint}")

# ----------------------------------------------- step 3: bivariate scoring
if a.plots_only:
    # Rebuild lrdata from the persisted matrices. Nothing under data/ is rewritten -- this path
    # exists so the figure suite can be extended without re-running the permutation test or
    # touching a result that other outputs (NMF, METHODS.md numbers) were derived from.
    _z = np.load(os.path.join(DATA, "local_scores.npz"), allow_pickle=True)
    lrdata = ad.AnnData(X=np.asarray(_z["values"], dtype=np.float32),
                        obs=pd.DataFrame(index=_z["cells"].astype(str)),
                        var=pd.DataFrame(index=_z["pairs"].astype(str)))
    lrdata.obsm["spatial"] = xy
    for _lay, _f in (("pvals", "local_pvals.npz"), ("cats", "local_categories.npz")):
        _p = os.path.join(DATA, _f)
        if os.path.exists(_p):
            lrdata.layers[_lay] = np.asarray(np.load(_p, allow_pickle=True)["values"])
    gv = pd.read_csv(os.path.join(DATA, "global_scores.csv"), index_col=0)
    assert list(gv.index.astype(str)) == list(lrdata.var_names), \
        "global_scores.csv and local_scores.npz disagree on the LR-pair order"
    log(f"  PLOTS-ONLY: reloaded {lrdata.shape[0]} cells x {lrdata.shape[1]} LR pairs "
        f"from {DATA} (layers: {list(lrdata.layers)})")
else:
    log(f"bivariate(local={a.local_name}, global={a.global_name}, n_perms={a.n_perms}, "
        f"nz_prop={a.nz_prop})")
    lrdata = li.mt.bivariate(adata, local_name=a.local_name, global_name=a.global_name,
                             n_perms=a.n_perms, mask_negatives=False, add_categories=True,
                             nz_prop=a.nz_prop, use_raw=False, verbose=True, **res_kw)
    log(f"  lrdata: {lrdata.shape[0]} cells x {lrdata.shape[1]} LR pairs")

    # -------------------------------------------------------- persist QUANT
    gv = lrdata.var.copy()
    gv.to_csv(os.path.join(DATA, "global_scores.csv"))
    np.savez_compressed(os.path.join(DATA, "local_scores.npz"),
                        values=np.asarray(lrdata.X.todense() if hasattr(lrdata.X, "todense")
                                          else lrdata.X, dtype=np.float32),
                        pairs=np.asarray(lrdata.var_names).astype(str),
                        cells=np.asarray(lrdata.obs_names).astype(str))
    if "pvals" in lrdata.layers:
        pv = lrdata.layers["pvals"]
        np.savez_compressed(os.path.join(DATA, "local_pvals.npz"),
                            values=np.asarray(pv.todense() if hasattr(pv, "todense") else pv,
                                              dtype=np.float32))
    if "cats" in lrdata.layers:
        ct_ = lrdata.layers["cats"]
        np.savez_compressed(os.path.join(DATA, "local_categories.npz"),
                            values=np.asarray(ct_.todense() if hasattr(ct_, "todense") else ct_))
    _meta = {"cell": adata.obs_names, "x": xy[:, 0], "y": xy[:, 1],
             "cell_type": adata.obs["cell_type"].astype(str)}
    if a.punch_col and a.punch_col in adata.obs:
        _meta["core"] = adata.obs[a.punch_col].astype(str).values
    pd.DataFrame(_meta).to_csv(os.path.join(DATA, "cell_meta.csv"), index=False)
    log(f"  saved global_scores.csv ({gv.shape}) + local score/pval matrices")

# ------------------------------------------------------------------- PLOTS
sort_key = a.global_name if a.global_name in gv.columns else "mean"
ranked = gv.sort_values(sort_key, ascending=False).index.tolist()
top = ranked[:a.n_top_plots]
requested = [s.strip() for s in a.requested_lrs.split(",") if s.strip()]
log(f"top by {sort_key}: {top}")

TAG = lambda lr: lr.replace("^", "-").replace("/", "-")


def _sp(color, layer=None, cmap="magma", **kw):
    sc.pl.spatial(lrdata, color=color, layer=layer, size=1.4, cmap=cmap, spot_size=20,
                  show=False, **kw)


def _genes_of(lr):
    """Constituent ligand + receptor genes. Heteromeric subunits are '_'-joined by LIANA."""
    lig, rec = lr.split("^", 1)
    return [g for g in dict.fromkeys(lig.split("_") + rec.split("_")) if g in adata.var_names]


def _gene_maps(lr):
    """bivariate.ipynb plots the constituent genes beside the LR score map."""
    genes = _genes_of(lr)
    if not genes:
        log(f"  {lr}: no constituent gene on the panel -- gene map skipped"); return False
    sc.pl.spatial(adata, color=genes, size=1.4, cmap="viridis", spot_size=20,
                  ncols=min(3, len(genes)), vmin=0, vmax="p99.5", show=False)


try_plot(lambda: _sp(top[:4], vmax=1), f"top_{sort_key}")
for lr in top:
    # score, p-value and local-category maps -- the tutorial draws all three, and all three
    # matrices are persisted under data/. Previously only the score map was drawn for the
    # method's own top hits (the p-value map existed only for the requested pair).
    try_plot(lambda lr=lr: _sp([lr]), f"local_{TAG(lr)}")
    if "pvals" in lrdata.layers:
        try_plot(lambda lr=lr: _sp([lr], layer="pvals", cmap="magma_r"), f"local_{TAG(lr)}_pvals")
    if "cats" in lrdata.layers:
        try_plot(lambda lr=lr: _sp([lr], layer="cats", cmap="coolwarm"), f"local_{TAG(lr)}_cats")
    try_plot(lambda lr=lr: _gene_maps(lr), f"genes_{TAG(lr)}")

req_rank = {}
for lr in requested:
    if lr in gv.index:
        rk = ranked.index(lr) + 1
        row = gv.loc[lr]
        req_rank[lr] = dict(rank=rk, of=len(ranked),
                            **{k: (float(row[k]) if k in row else None)
                               for k in (sort_key, f"{sort_key}_pvals", "mean", "std")})
        log(f"  requested {lr}: rank {rk}/{len(ranked)} by {sort_key}")
        try_plot(lambda lr=lr: _sp([lr]), f"rank{rk}_{TAG(lr)}", REQ)
        if "pvals" in lrdata.layers:
            try_plot(lambda lr=lr: _sp([lr], layer="pvals", cmap="magma_r"),
                     f"rank{rk}_{TAG(lr)}_pvals", REQ)
        if "cats" in lrdata.layers:
            try_plot(lambda lr=lr: _sp([lr], layer="cats", cmap="coolwarm"),
                     f"rank{rk}_{TAG(lr)}_cats", REQ)
        try_plot(lambda lr=lr: _gene_maps(lr), f"rank{rk}_{TAG(lr)}_genes", REQ)
    else:
        req_rank[lr] = None
        log(f"  REQUESTED {lr}: absent -- filtered by nz_prop={a.nz_prop} or not in resource")

# ---------------- calibration + graph figures the bivariate tutorial draws ----------------
def _bandwidth_query():
    """bivariate.ipynb cell 19. NOTE the x-axis is a hard query RADIUS, not the gaussian sigma,
    so the comparable point for a sigma of `a.bandwidth` is the support radius R."""
    import plotnine as p9
    plot, dfb = li.ut.query_bandwidth(coordinates=xy, start=5, end=120, interval_n=40)
    dfb.to_csv(os.path.join(DATA if not a.plots_only else PLOTS, "bandwidth_query.csv"),
               index=False)
    (plot + p9.geom_vline(xintercept=float(R), linetype="dashed", color="red")
          + p9.annotate("text", x=float(R), y=dfb.neighbours.max(),
                        label=f"support radius = {R:.1f} um  (sigma = {a.bandwidth})",
                        angle=90, va="bottom", ha="right", color="red")
     ).save(os.path.join(PLOTS, "bandwidth_query.png"), dpi=200, width=7, height=4.5, verbose=False)
    log("  plot ok: bandwidth_query"); SAVED.append("bandwidth_query")
try:
    _bandwidth_query()
except Exception as e:
    log(f"  plot FAILED bandwidth_query: {type(e).__name__}: {e}")

def _connectivity():
    """bivariate.ipynb cell 23. li.pl.connectivity returns a PLOTNINE ggplot -- it must be saved
    with its own .save(); plt.gcf() would be an empty canvas and would write a blank PNG."""
    idx = int(np.argmin(np.abs(xy - xy.mean(0)).sum(1)))
    gg = li.pl.connectivity(adata, idx=idx, size=0.6, figure_size=(6, 6),
                            spatial_key="spatial", return_fig=True)
    gg.save(os.path.join(PLOTS, f"connectivity_idx{idx}.png"), dpi=200, width=6, height=6,
            verbose=False)
    log(f"  plot ok: connectivity_idx{idx}"); SAVED.append(f"connectivity_idx{idx}")
try:
    _connectivity()
except Exception as e:
    log(f"  plot FAILED connectivity: {type(e).__name__}: {e}")

# ------------------------------------------------------------------ manifests
try:
    _sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True).stdout.strip() or None
except Exception:
    _sha = None
_common = {"h5ad": a.h5ad, "dataset": os.path.basename(a.out_dir.rstrip("/")),
           "tier": "default" if a.db == "consensus" else "cellchatdb2",
           "resource": res_label, "resource_n_pairs": res_n,
           "resource_fingerprint": res_fingerprint,
           "n_pairs_tested": int(lrdata.shape[1]), "bandwidth": a.bandwidth,
           "cutoff": a.cutoff, "kernel": a.kernel, "effective_radius_um": round(float(R), 1),
           "max_neighbours": a.max_neighbours, "set_diag": bool(a.set_diag),
           "cap_binds_frac": float((deg >= a.max_neighbours).mean()),
           "connectivity_nnz": int(conn.nnz), "median_neighbours": float(np.median(deg)),
           "max_neighbours_observed": int(deg.max()),
           "local_name": a.local_name, "global_name": a.global_name, "n_perms": a.n_perms,
           "nz_prop": a.nz_prop, "n_cells": int(adata.shape[0]), "split": "none (whole slide)",
           "requested": req_rank, "liana": li.__version__, "seed": a.seed, "git_sha": _sha,
           "figures_saved": len(SAVED), "blank_figures_suppressed": BLANK,
           "wall_min": round((time.time() - t0) / 60, 1)}
# In --plots-only mode run_manifest.json describes a run we did NOT redo: leave it alone.
_name = "plot_manifest.json" if a.plots_only else "run_manifest.json"
json.dump({"script": "run_liana.py", "mode": "plots-only" if a.plots_only else "full", **_common},
          open(os.path.join(OUT, _name), "w"), indent=2, default=str)
log(f"DONE in {(time.time()-t0)/60:.1f} min -> {OUT} "
    f"({len(SAVED)} figures, {len(BLANK)} blank suppressed) [{_name}]")
