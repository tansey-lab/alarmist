#!/usr/bin/env python
"""The stLearn Xenium CCI tutorial's figure set, reproduced call-for-call.

Tutorial: stlearn.readthedocs.io/en/latest/tutorials/cell_cell_interaction_xenium.html
(local mirror: /Users/jiayifan/tansey_lab/stLearn/stlearn.readthedocs.io/en/latest/tutorials/)

SCOPE RULE: this script emits EXACTLY the plots the Xenium tutorial emits, with the tutorial's
own argument values and its own multi-panel figure compositions -- nothing more. Verified
against the tutorial source: the Xenium vignette calls only

    cluster_plot, feat_plot, gene_plot, lr_summary, lr_result_plot, cci_check,
    ccinet_plot, lr_chord_plot        (+ the tl call st.tl.cci.adj_pvals)

It NEVER calls lr_diagnostics, lr_n_spots, cci_map, lr_cci_map, lr_plot or lr_go -- those come
from the *generic Visium* CCI vignette -- and neither vignette calls het_plot, grid_plot or
deconvolution_plot. Those live in the pre-existing plots/ and plots_full/ trees and are
deliberately not reproduced here. See METHODS.md for the provenance table.

Nothing is recomputed: run()/run_cci() results are read back from data/grid.h5ad. The
cell-level `adata` is rebuilt (QC + normalize_total, ~30 s, fully deterministic) because the
tutorial's CELL 40 / 42 / 44 figures are grid-vs-single-cell side-by-side panels and the
single-cell panel cannot be drawn from grid.h5ad alone.

Usage:
  python plot_stlearn_tutorial.py --h5ad PATH --run-dir DIR --out-dir DIR
"""
import argparse, json, os, platform, random, sys, time

p = argparse.ArgumentParser()
p.add_argument("--h5ad", required=True, help="the SAME source h5ad run_stlearn.py was given")
p.add_argument("--run-dir", required=True, help="run directory containing data/grid.h5ad")
p.add_argument("--out-dir", required=True)
p.add_argument("--cell-type-col", default="cell_type")
p.add_argument("--count-layer", default="counts")
p.add_argument("--library-id", default="sample")
p.add_argument("--seed", type=int, default=0)
p.add_argument("--requested-lrs", default="GRN_SORT1,ANXA1_FPR1",
               help="standing project request; plotted into out-dir/requested/ so they are "
                    "never confused with the method's own ranking")
p.add_argument("--n-cpus", type=int, default=None)
a = p.parse_args()

# Mirror run_stlearn.py's macOS thread guard so preprocessing is bit-identical.
seed = a.seed
random.seed(seed); os.environ["PYTHONHASHSEED"] = str(seed)
if platform.system() == "Darwin":
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
REQ = os.path.join(OUT, "requested")
for d in (OUT, REQ): os.makedirs(d, exist_ok=True)
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
LABEL = "cell_type"          # our stand-in for the tutorial's 'leiden' (DEVIATIONS.md row 1)

def save(name, sub=OUT):
    plt.gcf().savefig(os.path.join(sub, name + ".png"), dpi=200, bbox_inches="tight")
    plt.close("all")

def guard(fn, name, sub=OUT):
    try:
        fn(); save(name, sub); log(f"  ok: {name}")
        return True
    except Exception as e:
        plt.close("all"); log(f"  FAILED {name}: {type(e).__name__}: {e}")
        return False

# ============================================================== rebuild cell-level adata
# Byte-for-byte the same preprocessing run_stlearn.py performed (its lines 70-131).
log("rebuilding cell-level adata from", a.h5ad)
adata = sc.read_h5ad(a.h5ad)
if a.count_layer != "X":
    if a.count_layer not in adata.layers:
        sys.exit(f"ERROR: layer '{a.count_layer}' not found; have {list(adata.layers)}")
    adata.X = adata.layers[a.count_layer].copy()
adata.layers.clear()
if adata.raw is not None: adata.raw = None

bad = [g for g in adata.var_names.astype(str) if "_" in g]
if bad:
    adata = adata[:, [g for g in adata.var_names.astype(str) if "_" not in g]].copy()
    log(f"  dropped {len(bad)} genes containing '_' (stLearn's L/R separator) -> {adata.shape[1]}")

adata.obs[LABEL] = adata.obs[a.cell_type_col].astype(str).astype("category")
keep = ~adata.obs[LABEL].isin(["nan", "NaN", "None", ""])
if (~keep).sum(): adata = adata[keep].copy()

adata.uns["spatial"] = {a.library_id: {
    "images": {"hires": np.zeros((1, 1, 4), dtype=np.uint8)},
    "use_quality": "hires",
    "scalefactors": {"tissue_hires_scalef": 1.0, "spot_diameter_fullres": 15},
}}
adata = st.convert_scanpy(adata)
st.pp.filter_genes(adata, min_counts=10)      # tutorial CELL 12
st.pp.filter_cells(adata, min_counts=10)      # tutorial CELL 12
adata.raw = adata                             # tutorial CELL 18
st.pp.normalize_total(adata)                  # tutorial CELL 27 (no log1p -- see NOTES.md)
log(f"  adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

# ============================================================== load the fitted grid
log("loading", os.path.join(a.run_dir, "data", "grid.h5ad"))
grid = sc.read_h5ad(os.path.join(a.run_dir, "data", "grid.h5ad"))
log(f"  grid: {grid.shape[0]} spots x {grid.shape[1]} genes")

man_p = os.path.join(a.run_dir, "run_manifest.json")
if os.path.exists(man_p):
    man = json.load(open(man_p))
    assert adata.shape[0] == man["n_cells_after_qc"], (
        f"rebuilt adata has {adata.shape[0]} cells but the run recorded "
        f"{man['n_cells_after_qc']} -- preprocessing drifted, refusing to plot")
    assert adata.shape[1] == man["n_genes_after_qc"], "gene count drifted from the run manifest"
    assert grid.shape[0] == man["n_occupied_spots"], "grid spot count drifted from the run manifest"
    log(f"  verified against run_manifest.json: {adata.shape[0]} cells, {adata.shape[1]} genes, "
        f"{grid.shape[0]} spots")

# The categorical must carry the SAME categories in both objects, or list_clusters= and the
# colour mapping silently disagree between the grid and single-cell panels of one figure.
cats = list(grid.obs[LABEL].cat.categories)
adata.obs[LABEL] = pd.Categorical(adata.obs[LABEL].astype(str), categories=cats)
assert not adata.obs[LABEL].isna().any(), "a cell type present in adata is missing from the grid"
log(f"  {len(cats)} cell types, categories aligned across grid and adata")

# ---------------------------------------------------------------- figure canvas sizing
# The ONLY tutorial argument we do not copy literally. The tutorial's figsize=(20, 5) / (20, 8)
# are sized for its 1.37:1 breast section; combined with the tutorial's own set_aspect('equal')
# they squeeze our 2.20:1 TMA into ~10% of the canvas and push the legend over the tissue.
# figsize changes no datum, no statistic and no selection -- it is pure canvas -- so we keep
# every other argument verbatim and give each panel a box matching the tissue's aspect.
# Recorded in DEVIATIONS.md.
XY = np.asarray(grid.obsm["spatial"], float)
ASPECT = (XY[:, 1].max() - XY[:, 1].min()) / (XY[:, 0].max() - XY[:, 0].min())
PANEL_H = 10.0
PANEL_W = PANEL_H / ASPECT
def figsize(n_panels, legend_room=1.0):
    return (n_panels * PANEL_W * legend_room, PANEL_H)
log(f"  tissue aspect {ASPECT:.2f} -> panel {PANEL_W:.1f} x {PANEL_H:.1f} "
    f"(tutorial's literal (20,5) would letterbox it into ~10% of the canvas)")

# ============================================ CELL 23 -- single-cell labels, standalone figure
guard(lambda: st.pl.cluster_plot(adata, use_label=LABEL, size=1, show_image=False,
                                 bbox_to_anchor=(1.2, 1), figsize=figsize(1, 1.9)),
      "cell23_cell_labels")

# ================================= CELL 40 -- 2-panel: grid dominant spots | single-cell labels
def cell40():
    fig, axes = plt.subplots(ncols=2, figsize=figsize(2, 1.9))
    st.pl.cluster_plot(grid, use_label=LABEL, size=10, ax=axes[0], show_plot=False,
                       show_image=False, bbox_to_anchor=(1.1, 1))
    st.pl.cluster_plot(adata, use_label=LABEL, ax=axes[1], show_plot=False,
                       show_image=False, bbox_to_anchor=(1.1, 1))
    axes[0].set_title(f"Grid {LABEL} dominant spots"); axes[0].set_aspect("equal")
    axes[1].set_title(f"Cell {LABEL} labels");         axes[1].set_aspect("equal")
guard(cell40, "cell40_grid_vs_cell_labels")

# ==================== CELL 42 -- 3-panel per cell type: proportions | grid subset | cell subset
# The tutorial loops over 3 clusters chosen by POSITION (groups[6], groups[10], groups[11]) --
# the idc/dcis/stroma names are local variables, not a selection rule. There is no positional
# equivalent for our 9 annotated types, so we iterate the complete set: same call, same
# arguments, complete iteration. Argument values (vmax=1, show_color_bar=False) are the
# tutorial's -- vmax=1 is what makes the proportion panels comparable across cell types.
def cell42(ct):
    fig, axes = plt.subplots(ncols=3, figsize=figsize(3, 1.6))
    grid.obs["group"] = grid.uns[LABEL][ct].values
    st.pl.feat_plot(grid, feature="group", ax=axes[0], show_plot=False, vmax=1,
                    show_color_bar=False, show_image=False)
    st.pl.cluster_plot(grid, use_label=LABEL, list_clusters=[ct], ax=axes[1],
                       show_plot=False, show_image=False)
    st.pl.cluster_plot(adata, use_label=LABEL, list_clusters=[ct], ax=axes[2],
                       show_plot=False, show_image=False)
    axes[0].set_title(f"Grid {ct} proportions (max = 1)"); axes[0].set_aspect("equal")
    axes[1].set_title(f"Grid {ct} max spots");             axes[1].set_aspect("equal")
    axes[2].set_title(f"Individual cell {ct}");            axes[2].set_aspect("equal")
for ct in cats:
    guard(lambda ct=ct: cell42(ct), f"cell42_composition_{ct.replace('/', '-')}")

# ===================================== CELL 44 -- 2-panel gene expression: grid | single cell
# The tutorial plots one gene, CXCL12, a prior biological pick for its breast section. We have
# no tutorial-given gene, so we use the ligand and receptor of the top-ranked LR pair -- the
# tutorial's own `lr_summary.index[0]` idiom, applied to gene identity. The tutorial's
# vmax=80 on the single-cell panel is a display clip tuned to its counts and is NOT copied;
# recorded as a deviation.
summ = grid.uns["lr_summary"]
def cell44(gene):
    fig, axes = plt.subplots(ncols=2, figsize=figsize(2, 1.15))
    st.pl.gene_plot(grid, gene_symbols=gene, ax=axes[0], show_color_bar=False,
                    show_plot=False, show_image=False)
    st.pl.gene_plot(adata, gene_symbols=gene, ax=axes[1], show_color_bar=False,
                    show_plot=False, show_image=False)
    axes[0].set_aspect("equal"); axes[1].set_aspect("equal")
    axes[0].set_title(f"Grid {gene} expression")
    axes[1].set_title(f"Cell {gene} expression")

top_lr_pre = str(summ.index.values[0])
for gene in dict.fromkeys(top_lr_pre.split("_")):
    if gene in grid.var_names and gene in adata.var_names:
        guard(lambda g=gene: cell44(g), f"cell44_gene_{gene}")
    else:
        log(f"  gene {gene} absent from the panel -- skipped")

# ================================================================ CELL 58 -- LR ranking bars
guard(lambda: st.pl.lr_summary(grid, n_top=500), "cell58_lr_summary_n500")
guard(lambda: st.pl.lr_summary(grid, n_top=50, figsize=(10, 3)), "cell58_lr_summary_n50")

# ============================================================ CELL 62 -- adjust significance
# The tutorial's own comment is "### Can adjust significance thresholds." and it passes exactly
# adj_pvals' defaults, which st.tl.cci.run() has already applied (adj_pvals docstring,
# analysis.py:382: "Default settings of this function are already run in st.tl.cci.run").
# We make the call for literal fidelity and PROVE it changed no statistic. It does re-sort
# lr_summary with an unstable argsort, which permutes LRs that are tied on n_spots_sig; we
# assert the top of the ranking (everything we plot) is untouched.
lrs_before = list(summ.index)
mats_before = {k: np.asarray(grid.obsm[k]).copy()
               for k in ["lr_scores", "p_vals", "p_adjs", "-log10(p_adjs)", "lr_sig_scores"]}
st.tl.cci.adj_pvals(grid, correct_axis="spot", pval_adj_cutoff=0.05, adj_method="fdr_bh")
summ = grid.uns["lr_summary"]
lrs_after = list(summ.index)
pos = {lr: j for j, lr in enumerate(lrs_after)}
idx = np.array([pos[lr] for lr in lrs_before])
for k, before in mats_before.items():
    assert np.allclose(before, np.asarray(grid.obsm[k])[:, idx], equal_nan=True), \
        f"adj_pvals changed {k} -- it is NOT the no-op the docstring claims"
moved = sum(1 for x, y in zip(lrs_before, lrs_after) if x != y)
assert lrs_before[:3] == lrs_after[:3], "adj_pvals reordered the top-3 LRs -- ranking is not stable"
log(f"  adj_pvals: verified no-op on all 5 obsm matrices; {moved}/{len(lrs_after)} tied LRs "
    f"reordered; top-3 unchanged {lrs_after[:3]}")

# ========================================== CELL 67 -- per-LR spatial maps, 3 statistics
best_lr = str(summ.index.values[0])          # tutorial: "Just choosing one of the top"
STATS = ["lr_scores", "-log10(p_adjs)", "lr_sig_scores"]
def lr_result(lr):
    fig, axes = plt.subplots(ncols=len(STATS), figsize=figsize(len(STATS), 1.15))
    for i, stat in enumerate(STATS):
        st.pl.lr_result_plot(grid, use_result=stat, use_lr=lr, show_color_bar=False,
                             ax=axes[i], show_image=False)
        axes[i].set_aspect("equal"); axes[i].set_title(f"{lr} {stat}")
guard(lambda: lr_result(best_lr), f"cell67_lr_result_{best_lr}")

# ================================================================== CELL 79 -- CCI diagnostic
guard(lambda: st.pl.cci_check(grid, LABEL, figsize=(16, 5)), "cell79_cci_check")

# ====================================== CELL 81 -- cell-type interaction networks (shared pos)
# return_pos=True / pos=pos_1 is the tutorial's mechanism for giving every per-LR network the
# SAME node layout, which is what makes them visually comparable. run_stlearn.py omitted it.
pos_1 = None
try:
    pos_1 = st.pl.ccinet_plot(grid, LABEL, return_pos=True, min_counts=30)
    save("cell81_ccinet_all"); log("  ok: cell81_ccinet_all (node layout captured)")
except Exception as e:
    plt.close("all"); log(f"  FAILED cell81_ccinet_all: {type(e).__name__}: {e}")

lrs_top3 = [str(x) for x in summ.index.values[0:3]]     # tutorial: index.values[0:3]
for lr in lrs_top3[0:2]:                                # tutorial: for best_lr in lrs[0:2]
    guard(lambda lr=lr: st.pl.ccinet_plot(grid, LABEL, lr, min_counts=2, figsize=(10, 7.5),
                                          pos=pos_1),
          f"cell81_ccinet_{lr}")

# ============================================================== CELL 83 -- chord diagrams
guard(lambda: st.pl.lr_chord_plot(grid, LABEL), "cell83_chord_all")
for lr in lrs_top3[0:2]:
    guard(lambda lr=lr: st.pl.lr_chord_plot(grid, LABEL, lr), f"cell83_chord_{lr}")

# ===================================================== standing request -> requested/ only
# Not a tutorial step. The comparator-benchmark rule requires GRN->SORT1 and ANXA1->FPR1 to be
# plotted whatever their rank, in a SEPARATE directory so they cannot be mistaken for stLearn's
# own ranking. Same calls, same arguments as the tutorial cells above.
requested = [s.strip() for s in a.requested_lrs.split(",") if s.strip()]
lr_names = list(summ.index)
for lr in requested:
    if lr not in lr_names:
        log(f"  REQUESTED '{lr}' absent from lr_summary -> not testable on this panel")
        continue
    rank = lr_names.index(lr) + 1
    guard(lambda lr=lr: lr_result(lr), f"cell67_lr_result_rank{rank}_{lr}", REQ)
    guard(lambda lr=lr: st.pl.ccinet_plot(grid, LABEL, lr, min_counts=2, figsize=(10, 7.5),
                                          pos=pos_1), f"cell81_ccinet_rank{rank}_{lr}", REQ)
    guard(lambda lr=lr: st.pl.lr_chord_plot(grid, LABEL, lr), f"cell83_chord_rank{rank}_{lr}", REQ)

# ================================================================================ manifest
n_main = len([f for f in os.listdir(OUT) if f.endswith(".png")])
n_req = len([f for f in os.listdir(REQ) if f.endswith(".png")])
json.dump({
    "script": "plot_stlearn_tutorial.py",
    "scope": "the stLearn Xenium CCI tutorial's figure set, call-for-call; nothing else",
    "tutorial": "stlearn.readthedocs.io/en/latest/tutorials/cell_cell_interaction_xenium.html",
    "tutorial_cells_reproduced": ["23", "40", "42", "44", "58", "62", "67", "79", "81", "83"],
    "not_reproduced_generic_visium_vignette_only":
        ["lr_diagnostics", "lr_n_spots", "cci_map", "lr_cci_map", "lr_plot", "lr_go"],
    "not_reproduced_in_neither_vignette": ["het_plot", "grid_plot", "deconvolution_plot"],
    "h5ad": a.h5ad, "run_dir": a.run_dir, "out_dir": OUT,
    "n_cells": int(adata.shape[0]), "n_genes": int(adata.shape[1]),
    "n_spots": int(grid.shape[0]), "cell_types": cats,
    "top_lr": best_lr, "top3_lrs": lrs_top3,
    "adj_pvals_verified_noop": True, "tied_lrs_reordered_by_adj_pvals": int(moved),
    "gene_plot_genes": top_lr_pre.split("_"),
    "requested_lrs": requested,
    "n_png_main": n_main, "n_png_requested": n_req,
    "seed": seed, "stlearn": st.__version__, "scanpy": sc.__version__,
    "python": platform.python_version(),
    "wall_min": round((time.time() - t0) / 60, 2),
}, open(os.path.join(OUT, "plot_manifest.json"), "w"), indent=2)
log(f"DONE in {(time.time()-t0)/60:.2f} min -> {OUT} ({n_main} PNGs + {n_req} in requested/)")
