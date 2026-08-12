#!/usr/bin/env python
"""Completeness pass: every applicable stLearn plot not produced on the first run.

Replots from the saved grid.h5ad (no recomputation). First run covered cluster_plot,
lr_summary, lr_result_plot, cci_check, ccinet_plot, lr_chord_plot, cci_map, lr_cci_map,
lr_diagnostics, lr_n_spots. This adds:

  lr_plot            -- per-LR ligand/receptor expression + significance, the plot that
                        actually says WHICH genes in WHICH spots (the key omission)
  het_plot           -- cell-type heterogeneity / diversity per spot
  feat_plot          -- per-spot cell-type proportion (from spot_mixtures)
  gene_plot          -- ligand and receptor gene expression
  deconvolution_plot -- per-spot cell-type composition pies
  grid_plot          -- gridding sanity check
  lr_go              -- GO enrichment of LR results (needs a GO source; guarded)

Excluded as a different analysis, not this workflow: trajectory / pseudotime / subcluster /
interactive-bokeh / mask / qc.

Usage: python plot_stlearn_full.py --run-dir DIR --out-dir DIR
"""
import argparse, json, os, time, warnings
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import stlearn as st
warnings.filterwarnings("ignore")

p = argparse.ArgumentParser()
p.add_argument("--run-dir", required=True)
p.add_argument("--out-dir", required=True)
p.add_argument("--n-top", type=int, default=4)
p.add_argument("--requested-lrs", default="GRN_SORT1,ANXA1_FPR1")
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)

def save(nm):
    plt.gcf().savefig(os.path.join(a.out_dir, nm + ".png"), dpi=180, bbox_inches="tight")
    plt.close("all")
def guard(fn, nm):
    try: fn(); save(nm); log(f"  ok: {nm}")
    except Exception as e: plt.close("all"); log(f"  FAILED {nm}: {type(e).__name__}: {e}")

log("loading grid.h5ad")
grid = sc.read_h5ad(os.path.join(a.run_dir, "data", "grid.h5ad"))
log(f"  {grid.shape[0]} spots x {grid.shape[1]} genes")
xy = np.asarray(grid.obsm["spatial"], float)
ASPECT = (xy[:, 1].max() - xy[:, 1].min()) / (xy[:, 0].max() - xy[:, 0].min())
FS = (4.6, 4.6 * ASPECT)          # tall TMA -- respect the aspect ratio
log(f"  tissue aspect {ASPECT:.2f} -> figsize {FS[0]:.1f} x {FS[1]:.1f}")

summ = pd.read_csv(os.path.join(a.run_dir, "data", "lr_summary.csv"), index_col=0)
top = list(summ.index[:a.n_top])
requested = [s.strip() for s in a.requested_lrs.split(",") if s.strip()]
log(f"  top LRs: {top} | requested: {requested}")

# ---- lr_plot: the omission that matters -- ligand & receptor expression WITH significance
for lr in top + [r for r in requested if r in summ.index]:
    guard(lambda lr=lr: st.pl.lr_plot(grid, lr, outer_size_prop=1, outer_mode="binary",
                                      pt_scale=6, use_label=None, show_image=False,
                                      sig_spots=True, figsize=FS),
          f"lr_plot_binary_{lr}")
    guard(lambda lr=lr: st.pl.lr_plot(grid, lr, inner_size_prop=0.08, middle_size_prop=0.3,
                                      outer_size_prop=0.5, outer_mode="continuous", pt_scale=40,
                                      use_label=None, show_image=False, sig_spots=True,
                                      figsize=FS),
          f"lr_plot_continuous_{lr}")

# ---- cell-type context
# st.tl.cci.het.count() writes obs['het'] (cell-type diversity per spot); het_plot reads it.
guard(lambda: (st.tl.cci.het.count(grid, use_label="cell_type", use_het="het"),
               st.pl.het_plot(grid, use_label="cell_type", figsize=FS)), "het_plot_celltype")
guard(lambda: st.pl.cluster_plot(grid, use_label="cell_type", size=6, show_image=False,
                                 show_plot=False, figsize=FS), "cluster_plot_celltype")
guard(lambda: st.pl.grid_plot(grid, use_label="cell_type", n_row=40, n_col=18, size=4,
                              figsize=FS), "grid_plot")
# deconvolution_plot reads uns['deconvolution']; grid() stored the same frame under the label.
if "cell_type" in grid.uns and "deconvolution" not in grid.uns:
    grid.uns["deconvolution"] = grid.uns["cell_type"]
for ct in list(grid.obs["cell_type"].cat.categories)[:4]:
    # tutorial pattern: copy the per-spot proportion column into obs before feat_plot
    grid.obs[f"prop_{ct}"] = np.asarray(grid.uns["cell_type"][ct].values, dtype=float)
    guard(lambda ct=ct: st.pl.deconvolution_plot(grid, use_label="cell_type", celltype=ct,
                                                 figsize=FS), f"deconvolution_{ct}")
    guard(lambda ct=ct: st.pl.feat_plot(grid, feature=f"prop_{ct}", figsize=FS,
                                        show_color_bar=True), f"feat_plot_prop_{ct}")

# ---- ligand / receptor gene expression for the requested + top pairs
genes = []
for lr in top[:2] + [r for r in requested if r in summ.index]:
    genes += [g for g in lr.split("_") if g in grid.var_names]
for g in dict.fromkeys(genes):
    guard(lambda g=g: st.pl.gene_plot(grid, gene_symbols=g, show_image=False, figsize=FS),
          f"gene_plot_{g}")

# ---- GO enrichment of LR results (needs a GO source; may require network)
# st.tl.cci.run_lr_go(adata, r_path=...) REQUIRES a working R + clusterProfiler install.
# Our comp-stlearn env is Python-only, so this is genuinely unavailable, not skipped by choice.
log("  lr_go: SKIPPED -- st.tl.cci.run_lr_go() requires an R install (r_path=) with "
    "clusterProfiler; comp-stlearn is a Python-only env. Recorded in DEVIATIONS.md.")

grid.write(os.path.join(a.run_dir, "data", "grid.h5ad"))   # persist het/prop columns added here
log("  re-saved grid.h5ad with het + proportion columns")
n = len([f for f in os.listdir(a.out_dir) if f.endswith(".png")])
json.dump({"script": "plot_stlearn_full.py", "n_png": n, "aspect": round(float(ASPECT), 2),
           "top_lrs": top, "wall_min": round((time.time() - t0) / 60, 1)},
          open(os.path.join(a.out_dir, "plot_manifest.json"), "w"), indent=2)
log(f"DONE in {(time.time()-t0)/60:.1f} min -> {a.out_dir} ({n} PNGs)")
