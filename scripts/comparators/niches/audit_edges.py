#!/usr/bin/env python
"""Audit the NICHES spatial edgelist: self-edges, effective neighbourhood radius, cross-core.

NICHES defines neighbours by a *mutual k-nearest-neighbour* graph with k = 4 and no distance
cutoff, so unlike every other comparator its neighbourhood has no radius parameter to quote.
This recovers the implied radius empirically from the edge list that NICHES actually built
(the CellToCellSpatial column names are '<sending cell>—<receiving cell>'), so NICHES can be
placed in METHODS.md's cross-method neighbourhood table in microns.

Also quantifies the self-edge behaviour: ComputeEdgelist.R:46-47 takes order(dist)[1:(k+1)],
which includes the cell itself at distance 0, and the mutual-NN symmetrisation keeps it. So
every cell is its own neighbour and NeighborhoodToCell mixes autocrine self-signal into the
niche average.

Env: bptf.  Usage:
    python audit_edges.py --root results/comparators/niches/GBM
"""
import argparse
import os

import numpy as np
import pandas as pd

EMDASH = "—"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", default="results/comparators/niches/GBM")
    p.add_argument("--tier", default="noimpute", help="edgelist is imputation-independent")
    a = p.parse_args()

    base = os.path.join(a.root, "cellchatdb2")
    cores = pd.read_csv(os.path.join(a.root, "input", "cores.csv"))

    rows, all_d = [], []
    for _, c in cores.iterrows():
        core = c["core"]
        meta = pd.read_csv(os.path.join(a.root, "input", core, "meta.csv"))
        xy = dict(zip(meta.cell_id.astype(str), zip(meta.x, meta.y)))

        colf = os.path.join(base, a.tier, core, "quant", "CellToCellSpatial_columns.tsv")
        if not os.path.exists(colf):
            continue
        cols = [l.rstrip("\n") for l in open(colf) if l.strip()]

        n_self, d = 0, []
        for col in cols:
            s, r = col.split(EMDASH)
            if s == r:
                n_self += 1
                continue
            (x1, y1), (x2, y2) = xy[s], xy[r]
            d.append(np.hypot(x1 - x2, y1 - y2))
        d = np.asarray(d)
        all_d.append(d)
        rows.append(dict(
            core=core, n_cells=int(c.n_cells), grade=c["grade"], n_edges=len(cols),
            n_self_edges=n_self,
            pct_self=round(100 * n_self / len(cols), 2),
            n_real_edges=len(d),
            real_edges_per_cell=round(len(d) / c.n_cells, 2),
            min_um=round(float(d.min()), 2), median_um=round(float(np.median(d)), 2),
            p95_um=round(float(np.percentile(d, 95)), 2), max_um=round(float(d.max()), 2),
        ))
        print(f"{core}: {len(cols)} edges ({n_self} self, {len(d)} real), "
              f"median {np.median(d):.1f} µm, p95 {np.percentile(d, 95):.1f} µm, max {d.max():.1f} µm")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(base, "summary_edge_audit.csv"), index=False)

    d = np.concatenate(all_d)
    print("\n=== pooled over all cores ===")
    print(f"real (non-self) edges      : {len(d):,}")
    print(f"self-edges                 : {int(df.n_self_edges.sum()):,} "
          f"({100 * df.n_self_edges.sum() / df.n_edges.sum():.1f}% of all edges)")
    print(f"median edge length         : {np.median(d):.1f} µm")
    print(f"95th percentile            : {np.percentile(d, 95):.1f} µm")
    print(f"99th percentile            : {np.percentile(d, 99):.1f} µm")
    print(f"max                        : {d.max():.1f} µm")
    print("\ncross-core edges           : 0 by construction "
          "(one NICHES run per core; no edge can span cores)")


if __name__ == "__main__":
    main()
