#!/usr/bin/env python
"""Persist the run_cci() quantitative outputs that run_stlearn.py left inside grid.h5ad.

NOT an analysis step and NOT a tutorial step -- the stLearn tutorial is a notebook and
persists nothing. This exists purely to satisfy the comparator-benchmark output contract
("Persist quantitative outputs, not just plot images ... score matrices, p-values, and
per-cell/per-spot assignments must land on disk in a re-readable form"). Nothing is
recomputed: every value written here was already produced by the original run and is read
back verbatim from data/grid.h5ad.

Two gaps it closes:

  1. `run_stlearn.py:206` loops only over ["lr_cci_<label>", "per_lr_cci_<label>"], so the
     stage-2 SIGNIFICANCE (`per_lr_cci_pvals_<label>`) and the uncorrected counts
     (`per_lr_cci_raw_<label>`, `lr_cci_raw_<label>`) never reached disk. The per-LR CSVs on
     disk therefore carried interaction counts with no p-values.
  2. `run_stlearn.py:156` writes lr_summary.csv BEFORE run_cci(), so the on-disk copy has 3
     columns while the object ends with 6 -- the three `*_<label>` columns run_cci() appends
     (analysis.py:770-772) are the stage-2 ranking and were missing.

grid.h5ad is opened READ-ONLY and never rewritten.

Usage: python export_stlearn_quant.py --run-dir DIR [--label cell_type]
"""
import argparse, json, os, time, warnings

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

p = argparse.ArgumentParser()
p.add_argument("--run-dir", required=True, help="run directory containing data/grid.h5ad")
p.add_argument("--label", default="cell_type", help="the use_label passed to run_cci()")
a = p.parse_args()

t0 = time.time()
DATA = os.path.join(a.run_dir, "data")
GRID = os.path.join(DATA, "grid.h5ad")
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)

log("reading", GRID, "(read-only)")
grid = sc.read_h5ad(GRID)
written = {}

# ---- 1. lr_summary with the run_cci() columns included -------------------------------
summ = grid.uns["lr_summary"]
dest = os.path.join(DATA, "lr_summary.csv")
prev = pd.read_csv(dest, index_col=0) if os.path.exists(dest) else None
summ.to_csv(dest)
written["lr_summary.csv"] = list(summ.columns)
log(f"  lr_summary.csv: {summ.shape[0]} LRs x {summ.shape[1]} cols {list(summ.columns)}"
    + (f"  (was {prev.shape[1]} cols: {list(prev.columns)})" if prev is not None else ""))
if prev is not None:
    shared = [c for c in prev.columns if c in summ.columns]
    same = all(np.array_equal(prev[c].values, summ.loc[prev.index, c].values) for c in shared)
    assert same, "existing lr_summary.csv disagrees with grid.h5ad on a shared column"
    log(f"  verified: all {len(shared)} pre-existing columns unchanged (pure column addition)")

# ---- 2. cell-type x cell-type matrices, pooled over LRs -------------------------------
for key in [f"lr_cci_{a.label}", f"lr_cci_raw_{a.label}"]:
    if key not in grid.uns:
        log(f"  {key}: ABSENT in grid.h5ad -- skipped")
        continue
    fn = os.path.join(DATA, f"{key}.csv")
    df = pd.DataFrame(grid.uns[key])
    df.to_csv(fn)
    written[f"{key}.csv"] = list(df.shape)
    log(f"  {key}.csv: {df.shape[0]} x {df.shape[1]}")

# ---- 3. per-LR cell-type x cell-type matrices ----------------------------------------
for key in [f"per_lr_cci_{a.label}", f"per_lr_cci_pvals_{a.label}", f"per_lr_cci_raw_{a.label}"]:
    if key not in grid.uns:
        log(f"  {key}: ABSENT in grid.h5ad -- skipped")
        continue
    d = grid.uns[key]
    sub = os.path.join(DATA, key)
    os.makedirs(sub, exist_ok=True)
    for lr, m in d.items():
        pd.DataFrame(m).to_csv(os.path.join(sub, f"{lr}.csv"))
    written[key] = len(d)
    log(f"  {key}/: {len(d)} per-LR matrices")

# ---- 4. manifest ----------------------------------------------------------------------
json.dump(
    {
        "script": "export_stlearn_quant.py",
        "run_dir": a.run_dir,
        "label": a.label,
        "source": "data/grid.h5ad (read-only; nothing recomputed)",
        "written": written,
        "wall_min": round((time.time() - t0) / 60, 2),
    },
    open(os.path.join(DATA, "export_manifest.json"), "w"),
    indent=2,
)
log(f"DONE in {(time.time()-t0)/60:.2f} min -> {DATA}")
