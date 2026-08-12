#!/usr/bin/env python
"""stLearn CCI on a spatial h5ad, following the authors' Xenium tutorial one-to-one.

Tutorial: stlearn.readthedocs.io/en/latest/tutorials/cell_cell_interaction_xenium.html
Call contract: scripts/comparators/stlearn/NOTES.md    Deviations: DEVIATIONS.md

Usage:
  python run_stlearn.py --h5ad PATH --out-dir DIR [options]

Everything is parameterised; nothing about a dataset is hardcoded.
"""
import argparse, json, os, platform, random, sys, time

p = argparse.ArgumentParser()
p.add_argument("--h5ad", required=True)
p.add_argument("--out-dir", required=True)
p.add_argument("--cell-type-col", default="cell_type")
p.add_argument("--count-layer", default="counts",
               help="layer holding RAW counts; 'X' to use .X directly")
p.add_argument("--lrs", default="default",
               help="'default' = stLearn's connectomeDB2020_lit; else a file of L_R lines")
p.add_argument("--n-row", type=int, default=125, help="grid bins along y (imagerow)")
p.add_argument("--n-col", type=int, default=125, help="grid bins along x (imagecol)")
p.add_argument("--distance", type=float, default=250.0, help="tutorial default, physical units")
p.add_argument("--n-pairs", type=int, default=10000, help="tutorial recommends ~10000")
p.add_argument("--n-perms", type=int, default=1000, help="tutorial recommends ~1000")
p.add_argument("--min-spots", type=int, default=20)
p.add_argument("--seed", type=int, default=0)
p.add_argument("--n-top-plots", type=int, default=6)
p.add_argument("--requested-lrs", default="GRN_SORT1,ANXA1_FPR1",
               help="always plotted regardless of rank, into plots/requested/")
p.add_argument("--n-cpus", type=int, default=None)
p.add_argument("--library-id", default="sample", help="key for the Visium-style uns['spatial']")
a = p.parse_args()

# Tutorial cell 1: seeds + macOS BLAS/numba guard, BEFORE numpy/scanpy import.
seed = a.seed
random.seed(seed); os.environ["PYTHONHASHSEED"] = str(seed)
if platform.system() == "Darwin":
    # Tutorial cell 1 pins these to 1 on macOS ("BLAS/numba deadlocks are common") and must run
    # before numpy/scanpy import. Default reproduces that exactly. If --n-cpus is given we raise
    # the thread caps to match -- otherwise NUMBA_NUM_THREADS=1 would silently keep the njit
    # kernels single-threaded no matter what n_cpus we pass to stLearn. Performance only.
    n_cpus = 1 if a.n_cpus is None else a.n_cpus
    for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[v] = str(n_cpus)
else:
    n_cpus = a.n_cpus

import numpy as np, pandas as pd, scanpy as sc, stlearn as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
np.random.seed(seed)
import warnings; warnings.filterwarnings("ignore")

t0 = time.time()
OUT = a.out_dir
PLOTS = os.path.join(OUT, "plots"); REQ = os.path.join(PLOTS, "requested")
DATA = os.path.join(OUT, "data")
for d in (PLOTS, REQ, DATA): os.makedirs(d, exist_ok=True)
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)

def save(fig_or_none, name, sub=PLOTS):
    """Tutorial plots draw onto the current figure; save + close whatever is active."""
    f = plt.gcf()
    f.savefig(os.path.join(sub, name + ".png"), dpi=200, bbox_inches="tight")
    plt.close("all")

# ---------------------------------------------------------------- load
log("reading", a.h5ad)
adata = sc.read_h5ad(a.h5ad)
log(f"loaded {adata.shape[0]} cells x {adata.shape[1]} genes")

if a.count_layer != "X":
    if a.count_layer not in adata.layers:
        sys.exit(f"ERROR: layer '{a.count_layer}' not found; have {list(adata.layers)}")
    adata.X = adata.layers[a.count_layer].copy()
    log(f"X <- layers['{a.count_layer}'] (raw counts; tutorial forbids log1p before run())")
adata.layers.clear()
if adata.raw is not None: adata.raw = None

# stLearn encodes an LR pair as the string "LIGAND_RECEPTOR", so '_' is its L/R separator and
# st.tl.cci.run() hard-errors on any gene name containing one ("Recommend to rename adata.var_names
# or remove these genes"). On the Xenium 5K panel the offenders are the Intergenic_Region_* genomic
# CONTROL probes -- not real genes, absent from every LR database, and normally excluded upstream
# (the tutorial's 313 genes are real panel only, 10x keeps controls in separate feature types).
# Dropped BEFORE normalize_total so library size reflects the real panel, matching the tutorial.
bad = [g for g in adata.var_names.astype(str) if "_" in g]
if bad:
    log(f"dropping {len(bad)} genes whose names contain '_' (stLearn's L_R separator): "
        f"{bad[:3]}{' ...' if len(bad) > 3 else ''}")
    adata = adata[:, [g for g in adata.var_names.astype(str) if "_" not in g]].copy()
    log(f"  -> {adata.shape[1]} genes remain")

if a.cell_type_col not in adata.obs:
    sys.exit(f"ERROR: '{a.cell_type_col}' not in obs; have {list(adata.obs.columns)}")
adata.obs["cell_type"] = adata.obs[a.cell_type_col].astype(str).astype("category")
keep = ~adata.obs["cell_type"].isin(["nan", "NaN", "None", ""])
if (~keep).sum(): log(f"dropping {int((~keep).sum())} cells with no cell type"); adata = adata[keep].copy()

# stLearn expects the Visium-style adata.uns['spatial'] that its OWN readers build; both
# st.convert_scanpy() and st.tl.cci.grid() index into it (grid() copies it verbatim, line 172).
# We load a plain h5ad, so we replicate exactly what st.read_xenium() constructs:
#   uns['spatial'][lib] = {images:{quality:arr}, use_quality:quality, scalefactors:{...}}
# read_xenium ITSELF creates a BLANK placeholder image when no image file is supplied, sized
# (1.1*max_coord)^2 -- here 18,299^2 RGBA = 1.34 GB, which would also bloat grid.h5ad. We use a
# 1x1 placeholder instead: every plot below passes show_image=False, and the array is only read
# on the _add_image() path. Per stLearn's own source comment the scalefactor scales the IMAGE to
# the spots, never the spots to the image, so it cannot move a data point. scale=1 matches the
# tutorial's read_xenium(scale=1) for Xenium micron coordinates, which makes convert_scanpy's
# `image_coor = obsm['spatial'] * scale` the identity.
LIB = a.library_id
adata.uns["spatial"] = {LIB: {
    "images": {"hires": np.zeros((1, 1, 4), dtype=np.uint8)},
    "use_quality": "hires",
    "scalefactors": {"tissue_hires_scalef": 1.0, "spot_diameter_fullres": 15},
}}
adata = st.convert_scanpy(adata)          # -> obs['imagecol','imagerow'], identity at scalef=1
log("uns['spatial'] built (read_xenium-equivalent); convert_scanpy -> imagecol/imagerow")
assert np.allclose(adata.obs["imagecol"].values, np.asarray(adata.obsm["spatial"])[:, 0]), \
    "convert_scanpy moved the coordinates -- scalefactor is not 1"
xs, ys = adata.obs["imagecol"].values, adata.obs["imagerow"].values
log(f"extent {xs.max()-xs.min():.0f} x {ys.max()-ys.min():.0f} (coord units)")

# ------------------------------------------------- tutorial steps 1,2,5
st.pp.filter_genes(adata, min_counts=10)
st.pp.filter_cells(adata, min_counts=10)
log(f"after QC (min_counts=10 both): {adata.shape[0]} cells x {adata.shape[1]} genes")
adata.raw = adata
st.pp.normalize_total(adata)          # library-size ONLY -- no log1p (see NOTES.md)
log("normalize_total done (no log1p, per tutorial)")

# ------------------------------------------------------- step 6: grid
log(f"gridding n_row={a.n_row} (y) n_col={a.n_col} (x)")
grid = st.tl.cci.grid(adata, n_row=a.n_row, n_col=a.n_col, use_label="cell_type", n_cpus=n_cpus)
gx, gy = grid.obs["imagecol"].values, grid.obs["imagerow"].values
spot_x = (xs.max()-xs.min())/a.n_col; spot_y = (ys.max()-ys.min())/a.n_row
log(f"grid: {grid.shape[0]} occupied spots, spot {spot_x:.1f} x {spot_y:.1f} units, "
    f"{adata.shape[0]/grid.shape[0]:.1f} cells/occupied-spot")

# -------------------------------------------------- step 7: LR database
if a.lrs == "default":
    lrs = st.tl.cci.load_lrs(["connectomeDB2020_lit"], species="human")
    db_name = "connectomeDB2020_lit"
else:
    lrs = np.array([l.strip() for l in open(a.lrs) if l.strip()])
    db_name = os.path.basename(a.lrs)
log(f"LR database '{db_name}': {len(lrs)} pairs")

# ------------------------------------------ step 8: LR permutation test
log(f"run(min_spots={a.min_spots}, distance={a.distance}, n_pairs={a.n_pairs})")
st.tl.cci.run(grid, lrs, min_spots=a.min_spots, distance=a.distance,
              n_pairs=a.n_pairs, n_cpus=n_cpus, random_state=seed)
lr_summary = grid.uns["lr_summary"]
log(f"lr_summary: {lr_summary.shape[0]} LR pairs survived filtering")
lr_summary.to_csv(os.path.join(DATA, "lr_summary.csv"))

# ------------------------------------------------- step 9: cell-type CCI
# --- GUARD: stLearn maps cell types to deconvolution columns by SUBSTRING and takes the first
# hit (het_helpers.get_data_for_counting: `[ct in col for col in cols]`). A label that is a
# substring of another (here 'mGAM' vs 'non-mGAM') can therefore silently bind to the WRONG
# column, corrupting every interaction reported for it. Verify the mapping instead of trusting
# alphabetical luck.
_cols = grid.uns["cell_type"].columns.values.astype(str)
_all = np.unique(grid.obs["cell_type"].astype(object).values.astype(str))
_mis = {ct: [c for c in _cols if ct in c][0] for ct in _all
        if [c for c in _cols if ct in c][0] != ct}
if _mis:
    sys.exit(f"ERROR: stLearn substring cell-type mapping is wrong for {_mis}. "
             f"Rename the offending labels so none is a substring of another.")
log(f"cell-type -> deconvolution column mapping verified for all {len(_all)} labels "
    f"(substring hazard: {[ct for ct in _all if sum(ct in c for c in _cols) > 1]})")

# --- WORKAROUND (version gap, see DEVIATIONS.md): run_cci does
# `adata.obs[label].values.astype(str)` and passes the result into an @njit kernel. Under
# pandas 3 / numpy 2 that yields a StringDtype/object array, which numba cannot type
# ("non-precise type array(pyobject, 1d, C)"). Object-dtype python strings make the same
# expression return a plain <U array, which numba accepts. Labels are byte-identical; only the
# dtype changes. run_cci never uses .cat, but the plotting functions do, so restore afterwards.
_cats = list(grid.obs["cell_type"].cat.categories)
_before = grid.obs["cell_type"].astype(str).tolist()
grid.obs["cell_type"] = grid.obs["cell_type"].astype(object)

log(f"run_cci(n_perms={a.n_perms}, spot_mixtures=True, cell_prop_cutoff=0.1)")
st.tl.cci.run_cci(grid, "cell_type", min_spots=2, spot_mixtures=True, cell_prop_cutoff=0.1,
                  sig_spots=True, n_perms=a.n_perms, random_state=seed, n_cpus=n_cpus)

grid.obs["cell_type"] = pd.Categorical(grid.obs["cell_type"], categories=_cats)
assert grid.obs["cell_type"].astype(str).tolist() == _before, "dtype workaround altered labels"

# ----------------------------------------------------- persist QUANT
lr_names = list(lr_summary.index)
for key in ["lr_scores", "p_vals", "p_adjs", "-log10(p_adjs)", "lr_sig_scores"]:
    if key in grid.obsm:
        fn = key.replace("(", "").replace(")", "").replace("-", "neg")
        pd.DataFrame(np.asarray(grid.obsm[key]), index=grid.obs_names, columns=lr_names) \
          .to_csv(os.path.join(DATA, f"spot_{fn}.csv.gz"), compression="gzip")
        log(f"  saved spot x LR matrix: {key} {np.asarray(grid.obsm[key]).shape}")
pd.DataFrame({"spot": grid.obs_names, "imagecol": gx, "imagerow": gy,
              "n_cells": grid.obs.get("n_cells", pd.Series(index=grid.obs_names, dtype=float)),
              "dominant_cell_type": grid.obs["cell_type"].astype(str)}) \
  .to_csv(os.path.join(DATA, "spot_meta.csv"), index=False)
if "cell_type" in grid.uns:
    pd.DataFrame(grid.uns["cell_type"], index=grid.obs_names) \
      .to_csv(os.path.join(DATA, "spot_celltype_proportions.csv.gz"), compression="gzip")
for key in [f"lr_cci_cell_type", f"per_lr_cci_cell_type"]:
    if key not in grid.uns: continue
    v = grid.uns[key]
    if isinstance(v, dict):
        os.makedirs(os.path.join(DATA, key), exist_ok=True)
        for lr, m in v.items():
            pd.DataFrame(m).to_csv(os.path.join(DATA, key, f"{lr}.csv"))
        log(f"  saved {len(v)} per-LR cell-type interaction matrices")
    else:
        pd.DataFrame(v).to_csv(os.path.join(DATA, f"{key}.csv")); log(f"  saved {key}")
grid.write(os.path.join(DATA, "grid.h5ad"))

# ------------------------------------------------------------- PLOTS
top = lr_names[:a.n_top_plots]
requested = [s.strip() for s in a.requested_lrs.split(",") if s.strip()]
log(f"top LRs: {top}")

def try_plot(fn, name, sub=PLOTS):
    try:
        fn(); save(None, name, sub); log(f"  plot ok: {name}")
    except Exception as e:
        plt.close("all"); log(f"  plot FAILED {name}: {type(e).__name__}: {e}")

try_plot(lambda: st.pl.cluster_plot(grid, use_label="cell_type", size=10, show_image=False,
                                    show_plot=False), "grid_celltypes")
try_plot(lambda: st.pl.cluster_plot(adata, use_label="cell_type", size=1, show_image=False,
                                    show_plot=False), "cell_celltypes")
try_plot(lambda: st.pl.lr_summary(grid, n_top=min(50, len(lr_names)), figsize=(10, 3)), "lr_summary_top50")
try_plot(lambda: st.pl.lr_summary(grid, n_top=min(500, len(lr_names))), "lr_summary_all")
try_plot(lambda: st.pl.lr_diagnostics(grid, figsize=(10, 2.5)), "lr_diagnostics")
try_plot(lambda: st.pl.lr_n_spots(grid, n_top=min(50, len(lr_names)), figsize=(11, 3)), "lr_n_spots")
try_plot(lambda: st.pl.cci_check(grid, "cell_type", figsize=(16, 5)), "cci_check")

STATS = ["lr_scores", "-log10(p_adjs)", "lr_sig_scores"]
def lr_panel(lr, sub):
    fig, axes = plt.subplots(ncols=len(STATS), figsize=(20, 7))
    for i, s in enumerate(STATS):
        st.pl.lr_result_plot(grid, use_result=s, use_lr=lr, show_color_bar=False,
                             ax=axes[i], show_image=False)
        axes[i].set_aspect("equal"); axes[i].set_title(f"{lr}  {s}")
    return fig

for lr in top:
    try_plot(lambda lr=lr: lr_panel(lr, PLOTS), f"lr_result_{lr}")
for lr in requested:                       # ALWAYS, regardless of rank
    if lr in lr_names:
        rank = lr_names.index(lr) + 1
        try_plot(lambda lr=lr: lr_panel(lr, REQ), f"rank{rank}_{lr}", REQ)
    else:
        log(f"  REQUESTED LR '{lr}' absent from lr_summary -> not testable/significant here")

try_plot(lambda: st.pl.ccinet_plot(grid, "cell_type", min_counts=30), "ccinet_all")
for lr in top[:3]:
    try_plot(lambda lr=lr: st.pl.ccinet_plot(grid, "cell_type", lr, min_counts=2, figsize=(10, 7.5)),
             f"ccinet_{lr}")
for lr in requested:
    if lr in lr_names:
        try_plot(lambda lr=lr: st.pl.ccinet_plot(grid, "cell_type", lr, min_counts=2, figsize=(10, 7.5)),
                 f"ccinet_{lr}", REQ)
try_plot(lambda: st.pl.lr_chord_plot(grid, "cell_type"), "chord_all")
for lr in top[:3]:
    try_plot(lambda lr=lr: st.pl.lr_chord_plot(grid, "cell_type", lr), f"chord_{lr}")
for lr in requested:
    if lr in lr_names:
        try_plot(lambda lr=lr: st.pl.lr_chord_plot(grid, "cell_type", lr), f"chord_{lr}", REQ)
try_plot(lambda: st.pl.cci_map(grid, "cell_type"), "cci_map")
if len(top):        # lr_cci_map raises UnboundLocalError on an empty lrs list
    try_plot(lambda: st.pl.lr_cci_map(grid, "cell_type", lrs=list(top), min_total=10), "lr_cci_map")

# --------------------------------------------------------- manifest
json.dump({
    "script": "run_stlearn.py", "h5ad": a.h5ad, "out_dir": OUT,
    "lr_database": db_name, "n_lr_pairs_input": int(len(lrs)),
    "n_lr_pairs_tested": int(lr_summary.shape[0]),
    "cell_type_col": a.cell_type_col, "count_layer": a.count_layer,
    "n_cells_loaded": int(adata.raw.shape[0]) if adata.raw is not None else None,
    "n_cells_after_qc": int(adata.shape[0]), "n_genes_after_qc": int(adata.shape[1]),
    "n_occupied_spots": int(grid.shape[0]),
    "cells_per_occupied_spot": round(adata.shape[0] / grid.shape[0], 2),
    "n_row": a.n_row, "n_col": a.n_col,
    "spot_size": [round(float(spot_x), 1), round(float(spot_y), 1)],
    "distance": a.distance, "n_pairs": a.n_pairs, "n_perms": a.n_perms,
    "min_spots": a.min_spots, "seed": seed,
    "stlearn": st.__version__, "scanpy": sc.__version__, "python": platform.python_version(),
    "wall_min": round((time.time() - t0) / 60, 1),
}, open(os.path.join(OUT, "run_manifest.json"), "w"), indent=2)
log(f"DONE in {(time.time()-t0)/60:.1f} min -> {OUT}")
