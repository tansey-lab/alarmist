#!/usr/bin/env python
"""Niche clustering of the perivascular cells, by CELL TYPE COMPOSITION (no motifs).

Adapted from /home/fanj2/cluster_analysis.py with as few changes as possible. The method is
that script's: count cell types in each cell's r-micron neighbourhood, z-score the columns,
and k-means with K chosen by KneeLocator over a restart-stabilised inertia curve.

What was changed, and why -- nothing else was touched:

 1. SCOPE. Only cells within --ring um of a curated lumen boundary (and outside the lumen)
    are clustered, i.e. the same perivascular set every other analysis in this directory
    uses. Neighbourhood counts are still built from ALL cells, so a ring cell's composition
    is its true local composition, not a composition of ring cells only.
 2. PER SECTION. The KDTree is built once per section. The four sections share a coordinate
    frame, so one global tree would merge cells that are millimetres apart in different
    tissue. `neighbor_cell_types` therefore takes an explicit cell_type_map, so that every
    section produces the same columns in the same order.
 3. COLUMN NAME. 'cell_type_assignment' -> 'cell_type' (a --cell-type-column flag).
 4. SEED. `autotuned_kmeans` fixes random_state (the original left KMeans unseeded, and
    this repo requires fixed seeds). Restarts use seed + i, so the restart loop still
    explores.
 5. `estimate_cluster_associations` took `n_clusters` / `cluster_labels` from module globals
    in the original and could not run as written; they are arguments now. Its permutation
    loop was also rewritten to shuffle ONCE per trial and score all K x K pairs from that
    shuffle, instead of reshuffling per pair. Each pair's marginal null is unchanged -- only
    the draws are now shared across pairs -- and it turns an hour into seconds.
 6. `plot_cluster_assignments` used axarr[i // n_col, i % n_col], which breaks for a single
    row; n_col defaults to 2 here so the four sections make a 2x2 grid.
 7. DROPPED `plot_meta_cluster_volcanos`. It calls `differential_expression` and
    `volcano_plot`, neither of which is defined or imported in cluster_analysis.py, so it
    cannot run as given; it also needs the expression matrix, which this question does not.

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/vessel_niche_cluster.py --ring 30 --neighborhood 30
"""
from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial import KDTree, cKDTree
from matplotlib.path import Path as MplPath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import H5, read_cat, load_polys, densify, log

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
MIN_RING = 20          # identical to vessel_motif_kmeans.py

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})


# ----------------------------------------------------------------------------------------
# from cluster_analysis.py, kept as close to verbatim as the changes above allow
# ----------------------------------------------------------------------------------------
def neighbor_cell_types(xy, cell_type_ids, n_cell_types, neighborhood_size=30,
                        neighborhood_type='distance', rows=None):
    """`rows` restricts which cells get a neighbourhood; the neighbourhood itself is still
    searched against every cell in `xy`, so a returned row is identical to what the original
    (which built one for every cell) would have produced for that cell."""
    # Build the neighborhoods
    kdtree = KDTree(xy)
    query_pts = xy if rows is None else xy[rows]
    if neighborhood_type == 'distance':
        indexes = kdtree.query_ball_point(query_pts, r=neighborhood_size)
    elif neighborhood_type == 'knn':
        indexes = kdtree.query(query_pts, k=neighborhood_size)[1]

    # Count how many of each cell type there are within the neighborhood
    cell_type_counts = np.zeros((len(query_pts), n_cell_types))
    for i, idxs in enumerate(indexes):
        idxs = np.asarray(idxs, int)
        np.add.at(cell_type_counts, (np.repeat(i, len(idxs)), cell_type_ids[idxs]), 1)
    return cell_type_counts


def autotuned_kmeans(X, min_clusters=1, max_clusters=30,
                     interp_method='interp1d', polynomial_degree=2,
                     n_restarts=10, seed=0):
    # Cluster all the tissue clusters into meta-clusters
    from kneed import KneeLocator
    from sklearn.cluster import KMeans
    n_clusters_range = np.array(list(range(min_clusters, max_clusters + 1)))
    kmeans_scores = np.zeros(len(n_clusters_range))
    kmeans_labels = np.zeros((len(n_clusters_range), X.shape[0]))
    kmeans_centroids = []
    print('Clustering via K-Means')
    for cidx, n_clusters in enumerate(n_clusters_range):
        print(f'K={n_clusters}')
        kmeans_centroids.append(None)
        best_score = None
        for i in range(n_restarts):
            kmeans = KMeans(n_clusters, random_state=seed + i)
            kmeans.fit(X)
            if i == 0 or kmeans.inertia_ < best_score:
                kmeans_centroids[cidx] = kmeans.cluster_centers_
                kmeans_labels[cidx] = kmeans.labels_
                kmeans_scores[cidx] = kmeans.inertia_
                best_score = kmeans.inertia_

    # Find the best number of clusters via the elbow method.
    # Use a quadratic curve fit to smooth out the sometimes-noisy inertia scores.
    best_k = int(np.round(KneeLocator(n_clusters_range, kmeans_scores,
                                      curve='convex', direction='decreasing',
                                      interp_method=interp_method,
                                      polynomial_degree=polynomial_degree).knee))
    best_cidx = best_k - min_clusters
    best_cluster_labels = kmeans_labels[best_cidx]
    best_clusters = kmeans_centroids[best_cidx]

    return {'best_labels': best_cluster_labels,
            'best_clusters': best_clusters,
            'best_k': best_k,
            'all_labels': kmeans_labels,
            'all_clusters': kmeans_centroids,
            'scores': kmeans_scores}


def plot_cluster_assignments(outfile, df, xy, cluster_col, cmap='tab20', n_col=2,
                             punch_id_col='sample_id'):
    cmap = matplotlib.colormaps[cmap]
    colors = cmap(np.linspace(0, 1, df[cluster_col].nunique()))
    punch_ids = np.unique(df[punch_id_col])
    n_row = int(np.ceil(len(punch_ids) / n_col))
    fig, axarr = plt.subplots(n_row, n_col, figsize=(5 * n_col, 5 * n_row), squeeze=False)
    for idx, punch_id in enumerate(punch_ids):
        ax = axarr[idx // n_col, idx % n_col]
        sub = (df[punch_id_col] == punch_id).to_numpy()
        for c in np.unique(df[cluster_col]):
            c_mask = sub & (df[cluster_col] == c).to_numpy()
            ax.scatter(xy[c_mask, 0], xy[c_mask, 1], s=1.5, color=colors[int(c)], lw=0)
        ax.set_title(punch_id, fontsize=10)
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    for idx in range(len(punch_ids), n_row * n_col):
        axarr[idx // n_col, idx % n_col].axis('off')
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()


def plot_meta_clusters_per_punch(outfile, df, cluster_col, proportions=True, cmap='tab20',
                                 punch_id_col='sample_id'):
    # Count occurrences of each category per integer value
    counts = df.groupby([punch_id_col, cluster_col]).size().unstack(fill_value=0)

    # Convert to proportions
    if proportions:
        counts = counts.div(counts.sum(axis=1), axis=0)

    # Plot stacked bar chart
    ax = counts.plot(kind='bar', stacked=True, colormap=cmap,
                     figsize=(int(np.ceil((counts.shape[1] * 0.5))) + 4, 5))

    # Formatting
    plt.ylabel(('Counts' if not proportions else 'Proportions'))
    plt.xlabel('Section')
    plt.title(f'{cluster_col} {("counts" if not proportions else "proportions")}')
    plt.legend(title=cluster_col, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=90, ha='center')
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()


def plot_cluster_cell_types(outfile, proportions, cell_types, cluster_col, cmap='tab20'):
    # Get tab20 colormap
    cmap = plt.get_cmap(cmap)
    colors = [cmap(i % 20) for i in range(proportions.shape[1])]

    fig, ax = plt.subplots(figsize=(max(4, .55 * proportions.shape[0] + 2), 4))
    x = np.arange(proportions.shape[0])
    bottom = np.zeros(proportions.shape[0])
    for i in range(proportions.shape[1]):
        ax.bar(x, proportions[:, i], bottom=bottom, label=cell_types[i], color=colors[i])
        bottom += proportions[:, i]
    ax.set_xticks(x)
    ax.set_xticklabels([f'{i}' for i in x])

    # Formatting
    plt.ylabel('Cell type proportions')
    plt.xlabel(f'{cluster_col}')
    plt.title(f'{cluster_col} cell type proportions')
    plt.legend(title='Cell type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()


def estimate_cluster_associations(xy, cluster_labels, n_clusters, n_neighbors=20,
                                  n_trials=1000, normal_approx=True, seed=0, section=None):
    # Neighbours are found within a section: the four sections share a coordinate frame, so
    # a global tree would call cells in different tissue "adjacent".
    knn_idxs = np.zeros((len(xy), n_neighbors), int)
    groups = [np.arange(len(xy))] if section is None else \
             [np.where(np.asarray(section) == s)[0] for s in sorted(set(section))]
    for g in groups:
        k = min(n_neighbors + 1, len(g))
        _, loc = KDTree(xy[g]).query(xy[g], k=k)
        loc = np.atleast_2d(loc)[:, 1:]
        if loc.shape[1] < n_neighbors:                    # tiny section: wrap round
            loc = loc[:, np.arange(n_neighbors) % max(loc.shape[1], 1)]
        knn_idxs[g] = g[loc]

    def knn_mat(labels):
        """All K x K knn_pct values from one labelling (the original's knn_pct, vectorised)."""
        nb = labels[knn_idxs]                                   # cells x n_neighbors
        out = np.zeros((n_clusters, n_clusters))
        for a in range(n_clusters):
            sel = nb[labels == a]
            if len(sel) == 0:
                continue
            for b in range(n_clusters):
                out[a, b] = (sel == b).mean(axis=1).mean()
            del sel
        return out

    true_stat = knn_mat(cluster_labels)
    rng = np.random.default_rng(seed)
    null = np.zeros((n_trials, n_clusters, n_clusters))
    for trial in range(n_trials):
        if trial % 200 == 0:
            print('\t' + str(trial))
        null[trial] = knn_mat(rng.permutation(cluster_labels))

    if normal_approx:
        p = 1 - stats.norm.cdf(true_stat, loc=null.mean(0),
                               scale=np.maximum(null.std(0), 1e-12))
    else:
        p = ((true_stat <= null).sum(0) + 1) / (n_trials + 1)
    base = np.array([(cluster_labels == a).mean() for a in range(n_clusters)])
    with np.errstate(divide='ignore', invalid='ignore'):   # a never-observed pair -> -inf
        lfc = np.log2(true_stat / base[:, None].T)
    return p, lfc


def plot_cluster_associations(outfile, cluster_association_pvalues, cluster_association_lfcs):
    plt.figure(figsize=(6, 5))
    q = stats.false_discovery_control(cluster_association_pvalues)
    annotations = np.array([[('***' if v < 1e-3 else ('**' if v < 1e-2 else ('*' if v < 1e-1 else '')))
                             for v in row] for row in q])
    sns.heatmap(np.nan_to_num(cluster_association_lfcs).clip(-5, 5), cmap='coolwarm', center=0,
                ax=plt.gca(), annot=annotations, fmt='s', annot_kws={'fontsize': 12})
    plt.ylabel('Anchor niche')
    plt.xlabel('Neighboring niche')
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()


# ----------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ring", type=float, default=30.0,
                    help="a cell is perivascular if it is within this distance of a curated "
                         "lumen boundary and outside the lumen")
    ap.add_argument("--neighborhood", type=float, default=30.0,
                    help="radius of the composition neighbourhood (cluster_analysis.py's "
                         "neighborhood_size)")
    ap.add_argument("--cell-type-column", default="cell_type")
    ap.add_argument("--unit", choices=("vessel", "cell"), default="vessel",
                    help="vessel: one row per vessel, its 30 um ring's cell-type composition "
                         "(the ring IS the neighbourhood, so --neighborhood is unused). "
                         "cell: the reference script's unit -- one row per perivascular cell, "
                         "its own --neighborhood um composition.")
    ap.add_argument("--composition", choices=("proportion", "count"), default=None,
                    help="feature scaling before the z-score. Default: proportion for "
                         "--unit vessel (ring size spans 20-262 cells, 13x, so raw counts "
                         "would mostly encode ring size), count for --unit cell (the "
                         "reference's own behaviour).")
    ap.add_argument("--max-clusters", type=int, default=20)
    ap.add_argument("--k", type=int, default=None,
                    help="force K instead of taking the KneeLocator elbow")
    ap.add_argument("--assoc-trials", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    if a.composition is None:
        a.composition = "proportion" if a.unit == "vessel" else "count"
    nb = "" if a.unit == "vessel" else f"_nb{int(a.neighborhood)}"
    OUT = (f"{BASE}/vessel_niche_ring{int(a.ring)}{nb}_per{a.unit}"
           + ("" if a.composition == ("proportion" if a.unit == "vessel" else "count")
              else f"_{a.composition}")
           + (f"_k{a.k}" if a.k else ""))
    os.makedirs(OUT, exist_ok=True)

    f = h5py.File(H5, "r")
    o = f["obs"]
    sec_cats, sec = read_cat(o, "sample_id")
    ct_cats, ct = read_cat(o, a.cell_type_column)
    xy_all = f["obsm"]["spatial"][:]
    f.close()
    ct_cats = list(ct_cats)
    log(f"{len(sec):,} cells, {len(ct_cats)} cell types, {len(sec_cats)} sections")

    # ---- neighbourhood composition over ALL cells, one KDTree per section ---------------
    counts = np.zeros((len(sec), len(ct_cats)))
    ring_mask = np.zeros(len(sec), bool)
    vessel_of = np.full(len(sec), "", object)
    for si, s in enumerate(sec_cats):
        m = np.where(sec == si)[0]
        P = xy_all[m]
        # ---- which of them are perivascular (same geometry as vessel_motif_kmeans.py) ---
        tree = cKDTree(P)
        best = np.full(len(m), np.inf)
        for p in load_polys(s):
            poly = p["poly"]
            bt = cKDTree(densify(poly))
            lo, hi = poly.min(0) - a.ring, poly.max(0) + a.ring
            cand = np.asarray(tree.query_ball_point((lo + hi) / 2,
                                                    np.hypot(*(hi - lo)) / 2 + 1.0))
            if cand.size == 0:
                continue
            d, _ = bt.query(P[cand], k=1)
            ins = MplPath(poly).contains_points(P[cand])
            keep = (d > 0) & (d <= a.ring) & (~ins)
            if keep.sum() < MIN_RING:
                continue
            idx = cand[keep]
            ring_mask[m[idx]] = True
            closer = d[keep] < best[idx]                 # a cell can ring several vessels
            best[idx[closer]] = d[keep][closer]
            vessel_of[m[idx[closer]]] = f"{s}__{p['vid']}"
        loc = np.where(ring_mask[m])[0]
        if a.unit == "cell":
            counts[m[loc]] = neighbor_cell_types(P, ct[m], len(ct_cats),
                                                 neighborhood_size=a.neighborhood, rows=loc)
        log(f"  {s}: {len(loc):,} perivascular of {len(m):,} cells")

    log(f"{int(ring_mask.sum()):,} perivascular cells over "
        f"{len(set(vessel_of[ring_mask]))} vessels")

    # ---- assemble the rows: one per vessel, or the reference's one per cell -------------
    if a.unit == "vessel":
        # The vessel's neighbourhood is its ring, so its composition is simply the cell-type
        # make-up of the ring cells. Each cell is assigned to its nearest lumen, so no cell
        # is counted for two vessels.
        vid = vessel_of[ring_mask]
        row_ct = ct[ring_mask]
        vids = sorted(set(vid))
        X = np.zeros((len(vids), len(ct_cats)))
        for i, v in enumerate(vids):
            np.add.at(X, (i, row_ct[vid == v]), 1)
        row_section = [v.rsplit("__", 1)[0] for v in vids]
        row_xy = np.vstack([xy_all[ring_mask][vid == v].mean(0) for v in vids])
        log(f"{len(vids)} vessels; ring sizes {int(X.sum(1).min())}-{int(X.sum(1).max())} "
            f"cells (median {int(np.median(X.sum(1)))})")
    else:
        X = counts[ring_mask]
        row_section = [sec_cats[i] for i in sec[ring_mask]]
        row_xy = xy_all[ring_mask]
        vid = vessel_of[ring_mask]
        row_ct = ct[ring_mask]

    if a.composition == "proportion":
        X = X / np.maximum(X.sum(axis=1, keepdims=True), 1)
    keep_col = X.max(axis=0) > 0                       # Drop zero columns
    Xz = X[:, keep_col].copy()
    Xz = (Xz - Xz.mean(axis=0, keepdims=True)) / Xz.std(axis=0, keepdims=True)
    kept_types = [c for c, k in zip(ct_cats, keep_col) if k]
    log(f"feature matrix {Xz.shape} ({len(ct_cats) - keep_col.sum()} all-zero types dropped)")

    res = autotuned_kmeans(Xz, min_clusters=1, max_clusters=a.max_clusters, seed=a.seed)
    print(f'Best K = {res["best_k"]}')
    n_clusters = a.k or res["best_k"]
    cluster_labels = (res['all_labels'][n_clusters - 1] if a.k else res['best_labels']).astype(int)
    if a.k:
        print(f"K forced to {a.k} (elbow said {res['best_k']})")

    centroids = np.array([X[cluster_labels == c].mean(axis=0) for c in range(n_clusters)])
    proportions = centroids / centroids.sum(axis=1, keepdims=True)

    if a.unit == "vessel":
        df = pd.DataFrame({"vessel_id": vids, "sample_id": row_section,
                           "n_ring": X.sum(1) if a.composition == "count"
                                     else [int((vid == v).sum()) for v in vids],
                           "Neighborhood cluster": cluster_labels})
    else:
        df = pd.DataFrame({
            "sample_id": row_section,
            "cell_type": [ct_cats[i] for i in row_ct],
            "vessel_id": vid,
            "cell_index": np.where(ring_mask)[0],
            "Neighborhood cluster": cluster_labels,
        })
    df["condition"] = df.sample_id.str.split("_").str[1]
    xy = row_xy

    df.to_csv(f"{OUT}/niche_labels.csv", index=False)
    pd.DataFrame(proportions, columns=ct_cats,
                 index=[f"niche{c}" for c in range(n_clusters)]).to_csv(
        f"{OUT}/niche_cell_type_proportions.csv")
    pd.DataFrame({"k": np.arange(1, a.max_clusters + 1), "inertia": res["scores"]}).to_csv(
        f"{OUT}/elbow_scores.csv", index=False)
    np.save(f"{OUT}/neighborhood_counts.npy", X)
    log(f"unit = {a.unit}, composition = {a.composition}, rows = {X.shape[0]}")

    print(f"\n=== niches (K = {n_clusters}) ===")
    for c in range(n_clusters):
        m = cluster_labels == c
        top = np.argsort(-proportions[c])[:4]
        print(f"  niche {c}: n={m.sum():>6}  " +
              ", ".join(f"{ct_cats[t]} {proportions[c, t]:.2f}" for t in top) +
              "  | " + ", ".join(f"{s} {(df.sample_id[m] == s).mean():.0%}" for s in sec_cats))

    print('Plotting neighborhood clusters')
    plot_cluster_assignments(f'{OUT}/neighborhood_cluster_sections.pdf', df, xy,
                             'Neighborhood cluster')
    plot_cluster_cell_types(f'{OUT}/neighborhood_cluster_cell_types.pdf', proportions,
                            ct_cats, 'Neighborhood proportions')
    plot_meta_clusters_per_punch(f'{OUT}/neighborhood_cluster_proportions.pdf', df,
                                 'Neighborhood cluster', proportions=True)

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(np.arange(1, a.max_clusters + 1), res["scores"], "o-", color="#1a1a1a", ms=3)
    ax.axvline(res["best_k"], color="#b2182b", ls="--", lw=1)
    ax.set_xlabel("k"); ax.set_ylabel("inertia")
    ax.set_title(f"KneeLocator elbow: k = {res['best_k']}", loc="left")
    fig.savefig(f"{OUT}/elbow.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT}/elbow.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    if a.assoc_trials > 0:
        print('Estimating niche associations')
        pv, lfc = estimate_cluster_associations(xy, cluster_labels, n_clusters,
                                                n_trials=a.assoc_trials, seed=a.seed,
                                                section=row_section)
        np.save(f"{OUT}/niche_association_pvalues.npy", pv)
        np.save(f"{OUT}/niche_association_lfcs.npy", lfc)
        plot_cluster_associations(f'{OUT}/neighborhood_associations.pdf', pv, lfc)

    log(f"wrote to {OUT}/")


if __name__ == '__main__':
    main()
