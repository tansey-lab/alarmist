#!/usr/bin/env python
"""Per-niche cell-type composition as absolute cell counts, not proportions.

`vessel_niche_cluster.py` writes `neighborhood_cluster_cell_types`, where every bar is
normalised to 100 % and a niche of 1 vessel looks as tall as one of 104. This redraws the same
data with bar height = the number of ring cells in that niche, so niche size is visible.

The per-vessel runs store proportions in `neighborhood_counts.npy`, so counts are recovered as
proportion x n_ring from `niche_labels.csv` (exact: max deviation from integer 1.4e-14, and the
total matches sum(n_ring), 20,279 cells over the 300 vessels).

Colours come from the run's own palette (celltype_palette), so the percentage figure is
re-rendered alongside the new one rather than left on plain tab20.

The same is done for the motif k-means vessel clusterings (`cluster_celltype_composition` ->
`cluster_celltype_counts`). Those runs cluster the SAME 300 vessels on the same 30 um rings, so
the per-vessel counts are taken from the niche run and joined on vessel_id; the vessel sets are
asserted identical. Their figures live in the sibling `*_composition` folder.

Outputs per run folder: neighborhood_cluster_cell_counts.{png,pdf,svg},
    neighborhood_cluster_cell_types.{png,pdf,svg} (re-rendered), niche_cell_type_counts.csv;
    for k-means runs: cluster_celltype_counts.{png,pdf,svg},
    cluster_celltype_composition.{png,pdf,svg} (re-rendered), cluster_celltype_counts.csv

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/plot_niche_cell_counts.py
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from celltype_palette import PRESENTATION, celltype_colors

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
# celltype_palette.PRESENTATION is exactly what was asked for here: Tumor_epi in the deep
# purple #54278f, plus the 3-cycle SMC pink / Plasma green / Langhans_cell grey. Reused rather
# than copied, so these figures cannot drift from the H&E and composition ones.
SWAP = PRESENTATION
CLUSTER_COL = "Neighborhood cluster"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.bbox": "tight",
})


def per_group(counts, labels, k):
    return np.vstack([counts[labels == i].sum(0) for i in range(k)])


def stacked(counts, cell_types, colours, ylabel, title, path_stem, xlabel=CLUSTER_COL,
            n_vessels=None):
    fig, ax = plt.subplots(figsize=(max(4, .55 * counts.shape[0] + 2), 4))
    x = np.arange(counts.shape[0])
    bottom = np.zeros(counts.shape[0])
    for i, ct in enumerate(cell_types):
        ax.bar(x, counts[:, i], bottom=bottom, label=ct, color=colours[ct])
        bottom += counts[:, i]
    if n_vessels is not None:                 # vessels per cluster, above its bar
        top = bottom.max()
        for xi, n in zip(x, n_vessels):
            ax.text(xi, bottom[xi] + .015 * top, f"{n} vessels", ha="center", va="bottom",
                    fontsize=6.5, color="#444444")
        ax.set_ylim(0, top * 1.10)
    ax.set_xticks(x, [str(i) for i in x])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(title="Cell type", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{path_stem}.{ext}", dpi=400 if ext == "png" else None)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", default=f"{BASE}/vessel_niche_ring30_pervessel*",
                    help="glob of niche run folders")
    ap.add_argument("--kmeans-runs", default=f"{BASE}/vessel_kmeans_ring30_*",
                    help="glob of motif k-means run folders (the *_composition siblings "
                         "receive the figures)")
    a = ap.parse_args()
    per_vessel = None        # counts table shared with the k-means runs
    for run in sorted(glob.glob(a.runs)):
        L = pd.read_csv(f"{run}/niche_labels.csv")
        P = np.load(f"{run}/neighborhood_counts.npy")
        cols = list(pd.read_csv(f"{run}/niche_cell_type_proportions.csv", index_col=0).columns)
        assert P.shape[1] == len(cols), (P.shape, len(cols))
        counts = P * L.n_ring.to_numpy()[:, None]            # proportions -> cells
        assert np.abs(counts - np.round(counts)).max() < 1e-6
        counts = np.round(counts)
        lab = L[CLUSTER_COL].to_numpy()
        k = lab.max() + 1
        per_niche = np.vstack([counts[lab == i].sum(0) for i in range(k)])
        colours = celltype_colors(cols, override=SWAP)

        stacked(per_niche, cols, colours, "Cells in the 30 µm rings",
                f"{CLUSTER_COL} cell counts", f"{run}/neighborhood_cluster_cell_counts",
                n_vessels=[int((lab == i).sum()) for i in range(k)])
        stacked(per_niche / per_niche.sum(1, keepdims=True), cols, colours,
                "Cell type proportions", f"{CLUSTER_COL} cell type proportions",
                f"{run}/neighborhood_cluster_cell_types")
        pd.DataFrame(per_niche.astype(int), columns=cols,
                     index=[f"niche{i}" for i in range(k)]).assign(
            n_vessels=[int((lab == i).sum()) for i in range(k)],
            n_cells=per_niche.sum(1).astype(int)).to_csv(f"{run}/niche_cell_type_counts.csv")
        print(f"  {os.path.basename(run)}: {k} niches, {int(per_niche.sum()):,} ring cells "
              f"(largest {int(per_niche.sum(1).max()):,}, smallest {int(per_niche.sum(1).min()):,})")
        if per_vessel is None:
            per_vessel = pd.DataFrame(counts, columns=cols, index=L.vessel_id)

    # ---- the motif k-means runs, same vessels and rings ---------------------------------
    cols = list(per_vessel.columns)
    colours = celltype_colors(cols, override=SWAP)
    for run in sorted(glob.glob(a.kmeans_runs)):
        if run.endswith("_composition") or not os.path.exists(f"{run}/vessels.csv"):
            continue
        V = pd.read_csv(f"{run}/vessels.csv")
        if "cluster" not in V.columns:
            print(f"  {os.path.basename(run)}: no cluster column, skipped")
            continue
        assert set(V.vessel_id) == set(per_vessel.index), os.path.basename(run)
        out = f"{run}_composition"
        os.makedirs(out, exist_ok=True)
        C = per_vessel.loc[V.vessel_id].to_numpy()
        lab = V.cluster.to_numpy()
        k = int(lab.max()) + 1
        tot = per_group(C, lab, k)
        stacked(tot, cols, colours, "Cells in the 30 µm rings", "Cluster cell counts",
                f"{out}/cluster_celltype_counts", xlabel="cluster",
                n_vessels=[int((lab == i).sum()) for i in range(k)])
        stacked(tot / tot.sum(1, keepdims=True), cols, colours, "Cell type proportions",
                "Cluster cell type proportions", f"{out}/cluster_celltype_composition",
                xlabel="cluster")
        pd.DataFrame(tot.astype(int), columns=cols,
                     index=[f"cluster{i}" for i in range(k)]).assign(
            n_vessels=[int((lab == i).sum()) for i in range(k)],
            n_cells=tot.sum(1).astype(int)).to_csv(f"{out}/cluster_celltype_counts.csv")
        print(f"  {os.path.basename(run)}: {k} clusters, largest {int(tot.sum(1).max()):,} "
              f"cells, smallest {int(tot.sum(1).min()):,}")


if __name__ == "__main__":
    main()
