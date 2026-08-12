#!/usr/bin/env python
"""Choose the inflow kernel bandwidth from LIANA-internal evidence, per the tutorial's own recipe.

Why this script exists
----------------------
The bandwidth in force (13.1454 um, support 28.21 um) was set by EQUAL-AREA correspondence to
ALARMIST's 50 um patch. `comparator-benchmark/SKILL.md:45-46` forbids exactly that ("do not match
it to ALARMIST's patch size"), so the value has no LIANA-internal justification and is an open
contract violation. `inflow_mofaflex.ipynb` states the criteria the authors actually want:

  * biologically, the bandwidth "should reflect the typical range of molecular signaling"
  * technically, too wide blurs fine spatial pattern, too narrow misses signalling gradients

and demonstrates `li.ut.query_bandwidth()` for exploring the trade-off. This script supplies the
three measurements needed to make that choice on THIS tissue, so the bandwidth can be re-derived
without reference to ALARMIST.

A UNIT TRAP THAT THE TUTORIAL DOES NOT FLAG
-------------------------------------------
`li.ut.query_bandwidth` counts neighbours inside a HARD radius -- it calls
`BallTree.query_radius(..., r=max_distance)` and its returned column is literally spelled
`bandwith`. But `li.ut.spatial_neighbors(bandwidth=b)` treats `b` as the gaussian sigma and
truncates at `cutoff`, giving an effective reach of

    R = b * sqrt(-2 ln cutoff) = 2.145966 * b     (at cutoff = 0.1)

So a value read off the query_bandwidth curve is an R, and must be DIVIDED by 2.146 before being
passed as `bandwidth=`. Reading the curve and passing the number straight through inflates the
neighbourhood area by 2.146^2 = 4.6x. Everything below is reported in BOTH units.

Run:
    /Users/jiayifan/anaconda3/envs/comp-liana/bin/python \\
        scripts/comparators/liana/choose_bandwidth.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree, radius_neighbors_graph

ROOT = Path("/Users/jiayifan/tansey_lab/alarmist")

p = argparse.ArgumentParser()
p.add_argument("--inflow-dir", default=str(ROOT / "results/comparators/liana/GBM/cellchatdb2_inflow"))
p.add_argument("--lric-dir", default=str(ROOT / "results/comparators/liana/GBM/lric_percore"))
p.add_argument("--out-dir", default=str(ROOT / "results/comparators/liana/GBM/bandwidth_choice"))
p.add_argument("--cutoff", type=float, default=0.1)
p.add_argument("--start", type=float, default=5.0)
p.add_argument("--end", type=float, default=120.0,
               help="tutorial uses 35; extended so the plateau and the 57.9 um inflow-tutorial "
                    "scale are both on the curve")
p.add_argument("--interval-n", type=int, default=60)
p.add_argument("--current-bandwidth", type=float, default=13.1454)
a = p.parse_args()

K = float(np.sqrt(-2 * np.log(a.cutoff)))          # sigma -> support-radius multiplier
OUT = Path(a.out_dir); (OUT / "data").mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(Path(a.inflow_dir) / "data" / "cell_meta.csv")
xy = meta[["x", "y"]].values.astype(float)
ct = meta["cell_type"].astype(str).values
core = meta["core"].astype(str).values if "core" in meta.columns else None
n = len(xy)
print(f"{n:,} cells | sigma->R multiplier k = {K:.6f} | current sigma {a.current_bandwidth} "
      f"-> R {a.current_bandwidth*K:.2f} um\n")

# ---------------------------------------------------------------- 1. neighbours vs radius
tree = BallTree(xy, metric="euclidean")
radii = np.linspace(a.start, a.end, a.interval_n)
counts = [float(np.mean(tree.query_radius(xy, r=r, count_only=True)) - 1) for r in radii]
curve = pd.DataFrame({"radius_um": radii, "mean_neighbours": counts,
                      "implied_sigma_um": radii / K})
curve.to_csv(OUT / "data" / "neighbours_vs_radius.csv", index=False)

# ---------------------------------------------------------- 2. cell size: the juxtacrine anchor
d1 = tree.query(xy, k=2)[0][:, 1]          # distance to the single nearest other cell
nn = dict(median=float(np.median(d1)), q25=float(np.percentile(d1, 25)),
          q75=float(np.percentile(d1, 75)), q95=float(np.percentile(d1, 95)))
print("=== cell spacing (nearest-neighbour distance) ===")
print(f"  median {nn['median']:.2f} um  IQR {nn['q25']:.2f}-{nn['q75']:.2f}  p95 {nn['q95']:.2f}")
print(f"  -> a strictly JUXTACRINE reach (touching cells only) is about {nn['median']:.0f}-"
      f"{nn['q75']:.0f} um, i.e. sigma ~ {nn['median']/K:.1f}-{nn['q75']/K:.1f} um")

print("\n=== mean neighbours at candidate scales ===")
print(f"{'R (um)':>9} {'sigma (um)':>11} {'mean nb':>9}   note")
NOTES = {28.21: "IN FORCE (equal-area to ALARMIST 50 um patch) -- contract violation",
         57.94: "inflow tutorial's own scale (bandwidth 27 um, MERFISH mouse brain)"}
for r in [10, 15, 20, 25, 28.21, 35, 40, 50, 57.94, 70, 90, 120]:
    c = float(np.mean(tree.query_radius(xy, r=r, count_only=True)) - 1)
    key = min(NOTES, key=lambda k: abs(k - r))
    note = NOTES[key] if abs(key - r) < 0.5 else ""
    print(f"{r:9.2f} {r/K:11.2f} {c:9.1f}   {note}")

# ------------------------------------------- 3. reachability of each sender vs radius
print("\n=== reachability: fraction of cells with >=1 neighbour of each sender type ===")
print("    (this is the ceiling on that sender's nonzero_fraction, i.e. what drives view attrition)")
types = sorted(set(ct))
rows = []
for r in [20, 28.21, 40, 57.94, 80]:
    G = radius_neighbors_graph(xy, radius=r, mode="connectivity", include_self=True)
    row = {"radius_um": r}
    for s in types:
        row[s] = float((G @ (ct == s).astype(float) > 0).mean())
    rows.append(row)
reach = pd.DataFrame(rows).set_index("radius_um")
reach.to_csv(OUT / "data" / "reachability_vs_radius.csv")
order = sorted(types, key=lambda s: -(ct == s).mean())
print(reach[order].round(3).to_string())

# ------------------------------------------- 4. empirical signalling range, from LRIC
lric = Path(a.lric_dir) / "combined" / "aggregate_per_bin.csv"
emp = None
if lric.exists():
    g = pd.read_csv(lric)
    print(f"\n=== empirical signalling range from LRIC g(r) ({lric.name}) ===")
    print(f"  columns: {list(g.columns)[:8]}")
    emp = g.head(20).to_dict("records")
else:
    print(f"\n=== LRIC aggregate not found at {lric} — skipping the empirical-range anchor ===")

# ------------------------------------------- 5. the tutorial's own figure, done correctly
# inflow_mofaflex.ipynb:
#     plot, df = li.ut.query_bandwidth(coordinates=..., start=5, end=35, interval_n=40)
#     plot + p9.scale_y_continuous(breaks=range(int(df.neighbours.min()),
#                                               int(df.neighbours.max())+1))
# Reproduced verbatim, then repeated over a wider range. The vertical guides are drawn at the
# SUPPORT RADIUS R, not at sigma -- query_bandwidth's x-axis is a hard query radius, so marking
# a sigma on it compares two different quantities and understates the kernel by 2.146x. The
# existing results/.../cellchatdb2_inflow/plots/global/bandwidth_query.png has that error.
import liana as li  # noqa: E402
import plotnine as p9  # noqa: E402

FIG = OUT / "figures"; FIG.mkdir(exist_ok=True)
GUIDES = [(a.current_bandwidth * K, f"in force: sigma {a.current_bandwidth} -> R {a.current_bandwidth*K:.1f}", "#c0392b"),
          (27.0 * K, f"inflow tutorial: sigma 27 -> R {27.0*K:.1f}", "#1f4e79")]

for tag, (lo, hi, ni) in {"tutorial_5_35": (5, 35, 40),
                          "extended_5_120": (a.start, a.end, a.interval_n)}.items():
    try:
        plot, df = li.ut.query_bandwidth(coordinates=xy, start=lo, end=hi, interval_n=ni)
        df.to_csv(OUT / "data" / f"query_bandwidth_{tag}.csv", index=False)
        # NOTE: the returned column is spelled 'bandwith' (missing d) in liana 1.8.1.
        g = plot + p9.scale_y_continuous(
            breaks=range(int(df.neighbours.min()), int(df.neighbours.max()) + 1,
                         max(1, int((df.neighbours.max() - df.neighbours.min()) // 12))))
        span = df.neighbours.max() - df.neighbours.min()
        for i, (x, lab, col) in enumerate(GUIDES):
            if lo <= x <= hi:
                # place the label INSIDE the panel and stagger the two guides, otherwise the
                # text is clipped at the top edge (the first version read "in fo...")
                g = (g + p9.geom_vline(xintercept=x, linetype="dashed", colour=col)
                       + p9.annotate("text", x=x - (hi - lo) * 0.012,
                                     y=df.neighbours.min() + span * (0.42 - 0.22 * i),
                                     label=lab, angle=90, ha="center", va="bottom",
                                     size=7.5, colour=col))
        g = g + p9.labs(
            x="query radius (um)   —   a HARD radius, not the gaussian sigma",
            # query_bandwidth.py:78 -> avg_nn = np.ceil(np.median(num_neighbors)); ... avg_nn - 1
            # so this is ceil(MEDIAN) - 1, despite the source calling it 'avg_nn'.
            y="ceil(median neighbours) − 1\n(as liana computes it — not the mean)",
            title=f"li.ut.query_bandwidth, {lo}–{hi} µm   (n={n:,} cells)")
        g = g + p9.theme(plot_title=p9.element_text(size=10),
                         axis_title=p9.element_text(size=9),
                         plot_margin=0.035)
        g.save(FIG / f"query_bandwidth_{tag}.png", dpi=200, width=8.2, height=5.0, verbose=False)
        print(f"  [plot] {FIG / f'query_bandwidth_{tag}.png'}")
    except Exception as e:
        print(f"  query_bandwidth {tag} FAILED: {type(e).__name__}: {e}")

json.dump({"k_sigma_to_R": K, "n_cells": n, "nn_distance_um": nn,
           "current_sigma": a.current_bandwidth, "current_R": a.current_bandwidth * K,
           "inflow_tutorial_sigma": 27.0, "inflow_tutorial_R": 27.0 * K,
           "reachability": reach.to_dict()},
          open(OUT / "bandwidth_choice.json", "w"), indent=2, default=str)
print(f"\nwrote {OUT}/bandwidth_choice.json")
