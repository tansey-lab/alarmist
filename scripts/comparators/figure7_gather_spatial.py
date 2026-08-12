#!/usr/bin/env python
"""Extract every method's spatial map of the two motif-1 arms, for ALL 13 TMA punches.

Ids run 1-14 with 7 absent. Each method scores a DIFFERENT quantity on a DIFFERENT
support, so maps are comparable only in shape, never in value; each is stored raw
here and percentile-scaled at plot time.

  CytoSignal  cell   LRscore, diffusion slot           whole-slide run, split by tma_id
  stLearn     spot   co-expression, 51.3 um grid       whole-slide run, split by tma_id
  SpatialDM   cell   local bivariate Moran z           per-core run; ONLY for globally
                                                       selected pairs -- GRN in 7 cores,
                                                       ANXA1 in 6, both in only 1/8/13
  COMMOT      cell   received transport mass           per-core run
  NICHES      cell   neighbourhood mechanism score     per-core run, ALRA tier
  LIANA+      cell   local bivariate cosine            whole-slide run, split by core
  CellChat    --     NO SPATIAL OUTPUT EXISTS
  ALARMIST    cell   motif-1 loading (one map, not two)

Whole-slide methods are cropped to a punch, never recomputed on it.

Env: /Users/jiayifan/anaconda3/envs/bptf/bin/python
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path("/Users/jiayifan/tansey_lab/alarmist")
CMP = ROOT / "results/comparators"
OUT = ROOT / "results/comparators/_benchmark/figure7"
OUT.mkdir(parents=True, exist_ok=True)

ARMS = ("GRN_SORT1", "ANXA1_FPR1")
store: dict[str, np.ndarray] = {}


def put(core: str, method: str, arm: str, x, y, v) -> None:
    store[f"{core}|{method}|{arm}|x"] = np.asarray(x, np.float32)
    store[f"{core}|{method}|{arm}|y"] = np.asarray(y, np.float32)
    store[f"{core}|{method}|{arm}|v"] = np.asarray(v, np.float32)


# ---------------------------------------------------------------- reference cells
with h5py.File(ROOT / "data/xenium_mm_final_cell_id.h5ad", "r") as f:
    def cat(name):
        g = f["obs"][name]
        if isinstance(g, h5py.Group):
            cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                               for c in g["categories"][:]])
            return cats[g["codes"][:]]
        return np.asarray([a.decode() if isinstance(a, bytes) else str(a) for a in g[:]])

    idx_key = f["obs"].attrs.get("_index", "_index")
    barcodes = np.asarray([b.decode() if isinstance(b, bytes) else str(b)
                           for b in f["obs"][idx_key][:]])
    tma = cat("tma_id")
    grade = cat("grade")
    gx = f["obs"]["centroid_x"][:].astype(float)
    gy = f["obs"]["centroid_y"][:].astype(float)

CORES = sorted(np.unique(tma), key=int)
core_grade = {c: grade[tma == c][0] for c in CORES}
print("cores:", ", ".join(f"{c}({core_grade[c][:2]},{int((tma == c).sum()):,})" for c in CORES))

h5_pos = {b: i for i, b in enumerate(barcodes)}

# ---------------------------------------------------------------- ALARMIST
U = np.load(ROOT / "results/GBM/single_cell/cell_loadings.npy")
for c in CORES:
    m = tma == c
    put(c, "ALARMIST", "motif1", gx[m], gy[m], U[m, 1])

# ---------------------------------------------------------------- CytoSignal
# 1.87 GB MatrixMarket, cells x interactions (89,035 x 919), 1-indexed. Parsed ONCE.
q = CMP / "cytosignal/GBM/cellchatdb2/run_full/quant"
cs_cells = pd.read_csv(q / "score_diffusion_Raw_smooth.cells.tsv", header=None)[0].astype(str).to_numpy()
COLS = {423: 0, 410: 1}                       # 423 GRN-SORT1, 410 ANXA1-FPR1
cs_v = np.zeros((len(cs_cells), 2), np.float32)
print("parsing CytoSignal mtx (63.5M entries)...")
with open(q / "score_diffusion_Raw_smooth.mtx") as fh:
    fh.readline(); fh.readline()
    for line in fh:
        i, j, val = line.split()
        k = COLS.get(int(j))
        if k is not None:
            cs_v[int(i) - 1, k] = float(val)
cs_gi = np.array([h5_pos.get(b, -1) for b in cs_cells])
cs_ok = cs_gi >= 0
cs_tma = np.full(len(cs_cells), "", object)
cs_tma[cs_ok] = tma[cs_gi[cs_ok]]
for c in CORES:
    m = cs_tma == c
    if m.sum():
        for arm, k in zip(ARMS, (0, 1)):
            put(c, "CytoSignal", arm, gx[cs_gi[m]], gy[cs_gi[m]], cs_v[m, k])

# ---------------------------------------------------------------- stLearn
# Spots carry no core label; assign each to the punch of its nearest real cell.
sm = pd.read_csv(CMP / "stlearn/GBM/cellchatdb2/data/spot_meta.csv")
sc = pd.read_csv(CMP / "stlearn/GBM/cellchatdb2/data/spot_lr_scores.csv.gz", index_col=0)
_, nn = cKDTree(np.c_[gx, gy]).query(np.c_[sm["imagecol"], sm["imagerow"]])
sm_tma = tma[nn]
for c in CORES:
    m = (sm_tma == c)
    if m.sum():
        for arm in ARMS:
            put(c, "stLearn", arm, sm.loc[m, "imagecol"], sm.loc[m, "imagerow"],
                sc[arm].to_numpy()[m])

# ---------------------------------------------------------------- LIANA+
lz = np.load(CMP / "liana/GBM/cellchatdb2/data/local_scores.npz", allow_pickle=True)
lm = pd.read_csv(CMP / "liana/GBM/cellchatdb2/data/cell_meta.csv")
lg = pd.read_csv(CMP / "liana/GBM/cellchatdb2/data/global_scores.csv")
LCOL = {a: int(lg.index[lg["interaction"] == k][0])
        for a, k in zip(ARMS, ("GRN^SORT1", "ANXA1^FPR1"))}
for c in CORES:
    m = (lm["core"].astype(str) == c).to_numpy()
    if m.sum():
        for arm in ARMS:
            put(c, "LIANA+", arm, lm.loc[m, "x"], lm.loc[m, "y"], lz["values"][m, LCOL[arm]])

# ---------------------------------------------------------------- per-core methods
missing_spatialdm: dict[str, list[str]] = {}
for c in CORES:
    # SpatialDM -- local statistic exists only for globally selected pairs
    d = CMP / f"spatialdm/GBM/cellchatdb2/{c}/data"
    z = np.load(d / "local_z.npz", allow_pickle=True)
    meta = pd.read_csv(d / "cell_meta.csv")
    pairs = list(pd.read_csv(d / "local_n_spots.csv")["interaction_name"])
    for arm in ARMS:
        if arm in pairs:
            put(c, "SpatialDM", arm, meta["x"], meta["y"], z["values"][pairs.index(arm)])
        else:
            missing_spatialdm.setdefault(c, []).append(arm)

    # COMMOT
    d = CMP / f"commot/GBM/cellchatdb2/{c}/data"
    rcv = pd.read_csv(d / "sum_receiver.csv.gz", index_col=0)
    meta = pd.read_csv(d / "cell_meta.csv")
    for arm, col in zip(ARMS, ("r-GRN-SORT1", "r-ANXA1-FPR1")):
        put(c, "COMMOT", arm, meta["x"], meta["y"],
            rcv[col].to_numpy() if col in rcv.columns else np.zeros(len(meta)))

    # NICHES (ALRA tier, NeighborhoodToCell)
    nq = CMP / f"niches/GBM/cellchatdb2/alra/core{c}/quant"
    nf = pd.read_csv(nq / "NeighborhoodToCell_features.tsv", header=None)[0].astype(str).to_numpy()
    nm = pd.read_csv(nq / "NeighborhoodToCell_metadata.csv")
    ROWS = {int(np.where(nf == "GRN—SORT1")[0][0]) + 1: 0,
            int(np.where(nf == "ANXA1—FPR1")[0][0]) + 1: 1}
    nv = np.zeros((len(nm), 2), np.float32)
    with open(nq / "NeighborhoodToCell_scores.mtx") as fh:
        fh.readline(); fh.readline()
        for line in fh:
            i, j, val = line.split()
            k = ROWS.get(int(i))
            if k is not None:
                nv[int(j) - 1, k] = float(val)
    for arm, k in zip(ARMS, (0, 1)):
        put(c, "NICHES", arm, nm["x"], nm["y"], nv[:, k])

    print(f"  core {c:>2} done")

np.savez_compressed(OUT / "spatial_maps_all_cores.npz", **store)
(OUT / "provenance_all_cores.json").write_text(json.dumps({
    "cores": CORES,
    "grade": core_grade,
    "n_cells": {c: int((tma == c).sum()) for c in CORES},
    "spatialdm_pairs_not_globally_selected": missing_spatialdm,
    "cellchat": "no per-cell / per-spot / per-core output exists; no map is drawable",
    "stlearn_core_assignment": "each spot inherits the tma_id of its nearest cell (cKDTree)",
    "whole_slide_methods_cropped_not_refit": ["CytoSignal", "stLearn", "LIANA+"],
    "warning": "values are NOT comparable across methods; percentile-scale at plot time",
}, indent=2))
print(f"\nwrote {OUT}/spatial_maps_all_cores.npz  ({len(store)//3} maps)")
