#!/usr/bin/env python
"""LIANA+ INFLOW score — the single-cell-resolution branch of the authors' decision tree.

Tutorial: /Users/jiayifan/tansey_lab/liana-py/docs/notebooks/inflow_score.ipynb

Why this and not `li.mt.bivariate`: LIANA's own README decision tree routes
    Spatially-resolved -> Resolution? -> Single-cell -> Interaction scoring -> INFLOW SCORE
while Local Bivariate Metrics sits under the *Spot-based* branch. Xenium is single-cell.
Correspondingly `_inflow.py` defaults `nz_prop=0.001` against bivariate's 0.05 -- the authors
addressed single-cell sparsity by writing a different method, not by asking users to lower a
threshold.

    Inflow_{i,s,l,r} = ( sum_j W_ij * L_{j,l} * C_{j,s} / sum_j W_ij ) * R_{i,r}

so a FEATURE is (sender cell type s) x (ligand l, receptor r), named '<cell_type>^<lig^rec>'.
Every term is non-negative on log1p input, so scores should be NMF-admissible; verified at runtime.

Usage: python run_inflow.py --h5ad PATH --out-dir DIR [options]
"""
import argparse, json, os, time
import numpy as np, pandas as pd
import scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import liana as li

p = argparse.ArgumentParser()
p.add_argument("--h5ad", required=True)
p.add_argument("--out-dir", required=True)
p.add_argument("--cell-type-col", default="cell_type")
p.add_argument("--count-layer", default="counts")
p.add_argument("--db", default="consensus")
p.add_argument("--bandwidth", type=float, default=13.1454,
               help="gaussian kernel bandwidth in COORDINATE UNITS. Set by equal-area correspondence to an s x s patch: support radius = s/sqrt(pi) = 28.2 um for s=50, and bandwidth = support / sqrt(-2 ln cutoff) = s/3.804 = 13.1454 um. (The earlier 18.75, derived from the bivariate tutorial's ~6-neighbour \"first ring\" rule, is VOID: first ring is a hexagonal-lattice concept with no referent on irregularly packed single cells.)")
p.add_argument("--cutoff", type=float, default=0.1)
p.add_argument("--nz-prop", type=float, default=0.001, help="inflow's own default")
p.add_argument("--svg-filter", action=argparse.BooleanOptionalAction, default=True,
               help="tutorial's Moran's I gene filter (FDR<0.05 & I>0.01) before inflow")
p.add_argument("--svi-filter", action="store_true", default=False,
               help="tutorial's OPTIONAL Moran's I filter on interactions after inflow")
p.add_argument("--seed", type=int, default=0)
p.add_argument("--punch-col", default="tma_id",
               help="obs column with the TMA punch/core ID; written to cell_meta.csv")
a = p.parse_args()
np.random.seed(a.seed)

OUT = a.out_dir
PLOTS, DATA = os.path.join(OUT, "plots"), os.path.join(OUT, "data")
for d in (PLOTS, DATA): os.makedirs(d, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)

log("reading", a.h5ad)
adata = sc.read_h5ad(a.h5ad)
if a.count_layer != "X":
    adata.X = adata.layers[a.count_layer].copy()
adata.layers.clear(); adata.raw = None; adata.uns.pop("log1p", None)
bad = [g for g in adata.var_names.astype(str) if "_" in g]
if bad:
    adata = adata[:, [g for g in adata.var_names.astype(str) if "_" not in g]].copy()
    log(f"dropped {len(bad)} '_'-containing control probes -> {adata.shape[1]} genes")
adata.obs["cell_type"] = adata.obs[a.cell_type_col].astype(str).astype("category")
adata.obsm["spatial"] = np.asarray(adata.obsm["spatial"], float)

# tutorial QC. NOTE run_liana.py does NOT do this, which is why the bivariate branch has 100,197
# cells against inflow's 100,190 -- any cross-branch comparison must join on cell ID, not position.
n_cells_pre = adata.shape[0]
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_genes(adata, min_cells=3)
log(f"QC dropped {n_cells_pre - adata.shape[0]} cells (min_genes=10)")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
log(f"after QC + normalize_total(1e4) + log1p: {adata.shape[0]} cells x {adata.shape[1]} genes")

li.ut.spatial_neighbors(adata, bandwidth=a.bandwidth, cutoff=a.cutoff, kernel="gaussian",
                        set_diag=True)
log(f"spatial_neighbors(bandwidth={a.bandwidth}) -> {adata.obsp['spatial_connectivities'].nnz:,} nnz")

n_genes_pre = adata.shape[1]
if a.svg_filter:
    import squidpy as sq
    sq.gr.spatial_autocorr(adata, mode="moran", use_raw=False, show_progress_bar=False)
    mi = adata.uns["moranI"]
    svgs = mi.index[(mi["pval_norm_fdr_bh"] < 0.05) & (mi["I"] > 0.01)]
    log(f"SVG filter (FDR<0.05 & I>0.01): {n_genes_pre} -> {len(svgs)} genes")
    adata.uns["moranI"].to_csv(os.path.join(DATA, "gene_moranI.csv"))
    adata = adata[:, svgs].copy()

if a.db == "consensus":
    resource, res_label = li.rs.select_resource("consensus"), "LIANA consensus"
else:
    db = pd.read_csv(a.db)
    resource = db[["ligand", "receptor"]].astype(str).drop_duplicates().reset_index(drop=True)
    res_label = os.path.basename(a.db)
log(f"resource: {res_label} ({resource.shape[0]} pairs)")

log(f"li.mt.inflow(groupby='cell_type', nz_prop={a.nz_prop})")
lrdata = li.mt.inflow(adata, groupby="cell_type", resource=resource, nz_prop=a.nz_prop,
                      use_raw=False, verbose=True)
log(f"  lrdata: {lrdata.shape[0]} cells x {lrdata.shape[1]} FEATURES")

if a.svi_filter:
    import squidpy as sq
    sq.gr.spatial_autocorr(lrdata, mode="moran", use_raw=False, show_progress_bar=False)
    mi = lrdata.uns["moranI"]
    svis = mi.index[(mi["pval_norm_fdr_bh"] <= 0.05) & (mi["I"] > 0.01)]
    log(f"  SVI filter: {lrdata.shape[1]} -> {len(svis)} features")
    lrdata.uns["moranI"].to_csv(os.path.join(DATA, "interaction_moranI.csv"))
    lrdata = lrdata[:, svis].copy()

# ---------------------------------------------------- REQUIRED DIAGNOSTICS
X = lrdata.X.toarray() if hasattr(lrdata.X, "toarray") else np.asarray(lrdata.X)
n_neg = int((X < 0).sum()); mn = float(X.min())
log("")
log("=== NMF ADMISSIBILITY ===")
log(f"  min value = {mn:.6g} | negative entries = {n_neg:,} of {X.size:,}")
log(f"  NMF-admissible (non-negative): {'YES' if n_neg == 0 else 'NO -- STOP'}")
log("=== FEATURE SPACE ===")
log(f"  n_features = {lrdata.shape[1]}   (bivariate gave 131)")
log(f"  var columns: {list(lrdata.var.columns)}")
log(f"  example var_names: {list(lrdata.var_names[:4])}")
n_ct = adata.obs['cell_type'].nunique()
log(f"  cell types = {n_ct} -> implied unique interactions = {lrdata.shape[1] / n_ct:.1f}")
log(f"  sparsity: {100*(X == 0).mean():.2f}% zeros")

lrdata.var.to_csv(os.path.join(DATA, "feature_var.csv"))
np.savez_compressed(os.path.join(DATA, "inflow_scores.npz"),
                    values=X.astype(np.float32),
                    features=np.asarray(lrdata.var_names).astype(str),
                    cells=np.asarray(lrdata.obs_names).astype(str))
_meta = {"cell": lrdata.obs_names,
         "x": lrdata.obsm["spatial"][:, 0] if "spatial" in lrdata.obsm else np.nan,
         "y": lrdata.obsm["spatial"][:, 1] if "spatial" in lrdata.obsm else np.nan,
         "cell_type": adata.obs["cell_type"].astype(str).values}
if a.punch_col and a.punch_col in adata.obs:
    _meta["core"] = adata.obs[a.punch_col].astype(str).values
pd.DataFrame(_meta).to_csv(os.path.join(DATA, "cell_meta.csv"), index=False)
lrdata.write(os.path.join(DATA, "inflow_lrdata.h5ad"))

_feat = np.asarray(lrdata.var_names).astype(str)
_lr_of = np.array(["^".join(f.split("^")[1:]) for f in _feat])
_send_of = np.array([f.split("^")[0] for f in _feat])
_per_sender = pd.Series(_send_of).value_counts()
try:
    import subprocess
    _sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True).stdout.strip() or None
except Exception:
    _sha = None
json.dump({"script": "run_inflow.py", "dataset": os.path.basename(a.out_dir.rstrip("/")),
           "tier": "default" if a.db == "consensus" else "cellchatdb2",
           "h5ad": a.h5ad, "db": (None if a.db == "consensus" else os.path.abspath(a.db)),
           "resource": res_label, "resource_n_pairs": int(resource.shape[0]),
           "bandwidth": a.bandwidth,
           "cutoff": a.cutoff, "nz_prop": a.nz_prop, "svg_filter": a.svg_filter,
           "svi_filter": a.svi_filter, "set_diag": True, "n_cells": int(lrdata.shape[0]),
           "n_cells_before_qc": int(n_cells_pre), "n_features": int(lrdata.shape[1]),
           "n_cell_types": int(n_ct),
           # the sender x LR expansion is RAGGED, not a complete grid -- do not report it as
           # n_unique_lr * n_cell_types
           "n_unique_lr_pairs": int(len(set(_lr_of))),
           "features_per_sender_min": int(_per_sender.min()),
           "features_per_sender_max": int(_per_sender.max()),
           "expansion_vs_unique_lr": round(float(lrdata.shape[1] / len(set(_lr_of))), 3),
           "n_genes_after_svg": int(adata.shape[1]), "min_value": mn,
           "n_negative": n_neg, "n_entries": int(X.size),
           "nmf_admissible": bool(n_neg == 0),
           "pct_zeros": float(100*(X == 0).mean()),
           "liana": li.__version__, "seed": a.seed, "git_sha": _sha,
           "wall_min": round((time.time()-t0)/60, 1)},
          open(os.path.join(OUT, "run_manifest.json"), "w"), indent=2, default=str)
log(f"DONE in {(time.time()-t0)/60:.1f} min -> {OUT}")
