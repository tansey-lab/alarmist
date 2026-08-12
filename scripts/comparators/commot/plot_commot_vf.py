#!/usr/bin/env python
"""Redraw COMMOT signalling vector fields with the TUTORIAL's argument values.

Why this exists: `run_commot.py`'s first pass called `ct.pl.plot_cell_communication` with the
function's own defaults (`plot_method='cell'`, `scale=1.0`, `ndsize=1`, no `normalize_v`). Those
render, but at 3k-26k cells per core the arrows are invisible and the colour scale is swamped by
near-zero values. The authors' notebooks use a different argument set entirely:

    ct.pl.plot_cell_communication(adata, database_name=..., lr_pair=..., plot_method='grid',
        background_legend=True, scale=0.00003, ndsize=8, grid_density=0.4, summary='sender',
        background='summary', clustering='leiden', cmap='Reds',
        normalize_v=True, normalize_v_quantile=0.995)
    (Basic_usage.ipynb; visium-mouse_brain.ipynb is the same call with pathway_name='PSAP')

This script reproduces that call. It reads each core's `adata_commot.h5ad` -- so it costs no OT
re-run, which is the entire point of having persisted `obsp`.

`scale` is the one value that cannot be copied verbatim. `_utils/_plotting.py:318-320` passes it
to `quiver(..., scale_units='x')`, so arrow length is `|v| / scale` in DATA units. The tutorial's
0.00003 is tuned to Visium full-resolution pixels (~9,000 units across); our cores are microns
(~1,000-3,000 across). Copying it would draw arrows orders of magnitude too long. We rescale to
preserve the tutorial's arrow-length-as-a-fraction-of-field: scale = 0.00003 * (9000 / x_extent).
Same class of unit trap as `dis_thr` -- see NOTES.md.

Usage: python plot_commot_vf.py --run-dir DIR [--cores 1,13] [--out-sub plots] [--scale-ref 9000]
"""
import argparse, glob, json, os, time
import numpy as np
import anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import commot as ct

p = argparse.ArgumentParser()
p.add_argument("--run-dir", required=True)
p.add_argument("--cores", default=None, help="comma-separated; default = every core dir found")
p.add_argument("--out-sub", default="plots")
p.add_argument("--db-name", default="cellchat")
p.add_argument("--clustering", default="cell_type")
p.add_argument("--requested-lrs", default="GRN_SORT1,ANXA1_FPR1")
p.add_argument("--n-top-pathways", type=int, default=5)
p.add_argument("--scale-ref", type=float, default=9000.0,
               help="x-extent of the tutorial's Visium field, in ITS units. The tutorial's "
                    "scale=0.00003 is defined against this; we hold arrow length as a constant "
                    "fraction of the field. 9000 ~= 65 hex rows x 137 units (NOTES.md).")
p.add_argument("--tutorial-scale", type=float, default=0.00003)
p.add_argument("--dry-run", action="store_true")
a = p.parse_args()

DBN, CLUST = a.db_name, a.clustering
requested = [s.strip() for s in a.requested_lrs.split(",") if s.strip()]
cores = (a.cores.split(",") if a.cores else
         sorted(os.path.basename(os.path.dirname(f))
                for f in glob.glob(os.path.join(a.run_dir, "*", "adata_commot.h5ad"))))
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
log(f"cores: {cores}")

# tutorial argument values, verbatim except `scale` (see module docstring)
TUT = dict(plot_method="grid", background_legend=True, ndsize=8, grid_density=0.4,
           normalize_v=True, normalize_v_quantile=0.995)

made, failed = 0, []
for c in cores:
    apath = os.path.join(a.run_dir, c, "adata_commot.h5ad")
    if not os.path.exists(apath):
        failed.append((c, "no adata_commot.h5ad")); continue
    A = ad.read_h5ad(apath)
    PLOTS = os.path.join(a.run_dir, c, a.out_sub); os.makedirs(PLOTS, exist_ok=True)
    pts = np.asarray(A.obsm["spatial"], float)
    extent = float(pts[:, 0].max() - pts[:, 0].min())
    scale = a.tutorial_scale * (a.scale_ref / extent)
    log(f"--- core {c}: {A.shape[0]} cells, x-extent {extent:.0f} um -> scale {scale:.3e}")
    if a.dry_run: continue

    # the keys we actually have cluster/vector results for: top-N pathways + requested pairs
    rcv = A.obsm[f"commot-{DBN}-sum-receiver"]
    info = A.uns[f"commot-{DBN}-info"]["df_ligrec"]
    pair_tags = {f"{r[0]}-{r[1]}" for r in info.itertuples(index=False)}
    pw_cols = [x for x in rcv.columns if x.startswith("r-") and x[2:] not in pair_tags
               and x[2:] != "total-total"]
    top_pw = [x[2:] for x in rcv[pw_cols].sum(axis=0).sort_values(ascending=False).index
              ][:a.n_top_pathways]
    lr_of = {f"{r[0]}_{r[1]}": (r[0], r[1]) for r in info.itertuples(index=False)}

    targets = [(dict(pathway_name=pw), f"pathway_{pw}") for pw in top_pw]
    targets += [(dict(lr_pair=lr_of[lr]), f"requested_{lr}") for lr in requested if lr in lr_of]

    for kw, label in targets:
        for smry in ("sender", "receiver"):
            # tutorial variant 1: background='summary', cmap='Reds' -- signal magnitude
            # tutorial variant 2: background='image' + clustering + cmap='Alphabet'. We have no
            # H&E, so we use background='cluster', the function's own no-image equivalent, to
            # get the cell-type-coloured panel the tutorial gets from the image.
            for bg, cmap, sfx in (("summary", "Reds", "sig"), ("cluster", "Alphabet", "ct")):
                nm = f"native_vf_{smry}_{label}_{sfx}.png"
                try:
                    ct.tl.communication_direction(A, database_name=DBN, k=5, **kw)
                    ct.pl.plot_cell_communication(
                        A, database_name=DBN, summary=smry, background=bg, cmap=cmap,
                        clustering=CLUST, scale=scale,
                        filename=os.path.join(PLOTS, nm), **TUT, **kw)
                    plt.close("all"); made += 1
                except Exception as e:
                    plt.close("all"); failed.append((c, nm, f"{type(e).__name__}: {e}"))
                    log(f"  FAILED {nm}: {type(e).__name__}: {e}")
    del A

log(f"DONE: {made} figures, {len(failed)} failures")
for f in failed: log("  FAIL", f)
json.dump({"script": "plot_commot_vf.py", "run_dir": a.run_dir, "cores": cores,
           "tutorial_args": TUT, "tutorial_scale": a.tutorial_scale, "scale_ref": a.scale_ref,
           "n_made": made, "failures": failed},
          open(os.path.join(a.run_dir, "vf_replot_manifest.json"), "w"), indent=2)
