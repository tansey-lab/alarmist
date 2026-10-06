#!/usr/bin/env python
"""Pick K for the vessel niches by per-cluster bootstrap Jaccard, not by a global index.

The rule, in the order it is applied:

  1. scan K from KMIN to the elbow's K (default 2..9);
  2. at each K, take SUB (default 100) subsamples of FRAC (default 0.8) of the vessels,
     WITHOUT replacement, recluster each one, and score every reference cluster by its
     best Jaccard against the subsample's clusters;
  3. choose the LARGEST K at which every cluster's mean Jaccard clears THRESH (0.75) --
     not the most stable K, which is almost always K=2 and throws away resolution;
  4. inspect the small clusters: low Jaccard means noise (evidence for a smaller K); high
     Jaccard but tiny means check whether it is one section only, or low-count rings (a QC
     artefact) before believing it is a rare but real population;
  5. clustree: how each cluster at K splits into K+1, to see whether a new small cluster is
     genuinely new structure or just an outlier being shaved off a large cluster's edge.

Features are the per-vessel 30 um ring compositions of vessel_niche_ring30_pervessel
(proportions over 19 cell types, z-scored, clipped at --clip because unclipped z lets one
16-sigma peribronchial vessel own a cluster).

The same subsamples also give the consensus matrix at every K (how often each pair of vessels
lands together), so the Jaccard numbers and the consensus panels are one experiment, not two.

Outputs (figures_niche_robustness/): niche_jaccard_by_k, niche_clustree,
    niche_consensus_all_k, niche_consensus_k<K> for every K,
    niche_jaccard_per_cluster.csv, niche_k_choice.csv

Usage:
    /home/fanj2/.conda/envs/spatial/bin/python \
        /home/fanj2/alarmist/scripts/research/luad_revision/niche_k_jaccard.py
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import log

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
SRC = f"{BASE}/vessel_niche_ring30_pervessel"
OUT = f"{BASE}/figures_niche_robustness"
MM = 1 / 25.4

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def save(fig, stem):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{OUT}/{stem}.{ext}", dpi=400 if ext == "png" else None)
    plt.close(fig)
    print(f"  {stem}.{{png,pdf,svg}}")


def zscore(P, clip):
    Z = (P - P.mean(0)) / P.std(0)
    return np.clip(Z, -clip, clip)


def fit(Z, k, seed=0):
    return KMeans(k, random_state=seed, n_init=10).fit_predict(Z)


def contingency(a, b, ka, kb):
    return np.bincount(a * kb + b, minlength=ka * kb).reshape(ka, kb)


def jaccard(Z, k, ref, rng, nboot, frac):
    """Per-cluster bootstrap Jaccard and the consensus matrix, from the same subsamples.

    Returns (mean Jaccard, sd Jaccard, consensus) where consensus[i, j] is the fraction of the
    subsamples containing both vessels in which they landed in the same cluster.
    """
    n, m = len(Z), int(round(frac * len(Z)))
    J = np.full((nboot, k), np.nan)
    co, seen = np.zeros((n, n)), np.zeros((n, n))
    for b in range(nboot):
        idx = rng.choice(n, m, replace=False)
        lb = fit(Z[idx], k, seed=b)
        C = contingency(ref[idx], lb, k, k).astype(float)
        rows, cols = C.sum(1, keepdims=True), C.sum(0, keepdims=True)
        with np.errstate(invalid="ignore"):
            Jm = C / (rows + cols - C)          # Jaccard of every (ref, boot) pair
        present = rows[:, 0] > 0
        J[b, present] = np.nanmax(Jm[present], axis=1)
        co[np.ix_(idx, idx)] += (lb[:, None] == lb[None, :])
        seen[np.ix_(idx, idx)] += 1
    return np.nanmean(J, 0), np.nanstd(J, 0), co / np.maximum(seen, 1)


def draw_consensus(ax, M, ref, k, title):
    order = np.argsort(ref, kind="stable")
    im = ax.imshow(M[np.ix_(order, order)], cmap="Blues", vmin=0, vmax=1,
                   interpolation="nearest")
    for e in np.cumsum(np.bincount(ref, minlength=k))[:-1]:
        ax.axhline(e - .5, color="#d55e00", lw=.5)
        ax.axvline(e - .5, color="#d55e00", lw=.5)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_linewidth(.4); sp.set_color("#999999")
    ax.set_title(title, loc="left", pad=3, fontsize=7)
    return im


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=9, help="the elbow's K; scan stops here")
    ap.add_argument("--sub", type=int, default=100, help="subsamples per K")
    ap.add_argument("--frac", type=float, default=0.8)
    ap.add_argument("--thresh", type=float, default=0.75)
    ap.add_argument("--clip", type=float, default=5.0)
    ap.add_argument("--small", type=int, default=10, help="flag clusters with <= this many vessels")
    a = ap.parse_args()

    P = np.load(f"{SRC}/neighborhood_counts.npy")
    L = pd.read_csv(f"{SRC}/niche_labels.csv")
    types = pd.read_csv(f"{SRC}/niche_cell_type_proportions.csv", index_col=0).columns.to_numpy()
    Z = zscore(P, a.clip)
    rng = np.random.default_rng(0)
    KS = list(range(a.kmin, a.kmax + 1))
    log(f"{len(Z)} vessels x {P.shape[1]} cell types; K {a.kmin}..{a.kmax}; "
        f"{a.sub} subsamples at {a.frac:.0%}, threshold {a.thresh}")

    rows, labels, cons = [], {}, {}
    for k in KS:
        ref = fit(Z, k)
        labels[k] = ref
        mu, sd, M = jaccard(Z, k, ref, rng, a.sub, a.frac)
        cons[k] = M
        for i in range(k):
            m = ref == i
            sec = L.sample_id[m].value_counts()
            comp = P[m].mean(0)
            top = np.argsort(comp)[::-1][:3]
            rows.append(dict(
                k=k, cluster=i, n_vessels=int(m.sum()), jaccard=mu[i], jaccard_sd=sd[i],
                n_sections=int(sec.size), top_section=sec.index[0],
                top_section_frac=sec.iloc[0] / m.sum(),
                median_n_ring=float(np.median(L.n_ring[m])),
                top_cell_types="; ".join(f"{types[t]} {comp[t]:.2f}" for t in top)))
        log(f"  K={k}: min Jaccard {mu.min():.2f} (cluster {int(mu.argmin())}, "
            f"n={int((ref == mu.argmin()).sum())}), smallest cluster "
            f"{int(np.bincount(ref).min())} vessels")
    R = pd.DataFrame(rows)

    # ---- the choice ---------------------------------------------------------------------
    per_k = R.groupby("k").jaccard.min()
    ok = per_k[per_k >= a.thresh]
    chosen = int(ok.index.max()) if len(ok) else int(per_k.idxmax())
    if len(ok):
        log(f"largest K with every cluster >= {a.thresh}: K = {chosen}")
    else:
        log(f"no K clears {a.thresh}; falling back to the best, K = {chosen}")

    # ---- figure 1: per-cluster Jaccard against K ----------------------------------------
    fig, ax = plt.subplots(figsize=(78 * MM, 50 * MM))
    for k in KS:
        s = R[R.k == k]
        ax.scatter(np.full(len(s), k), s.jaccard, s=6 + 42 * s.n_vessels / R.n_vessels.max(),
                   facecolor="none", edgecolor="#8aa8c2", lw=.7, zorder=3)
        tiny = s[s.n_vessels <= a.small]
        ax.scatter(np.full(len(tiny), k), tiny.jaccard, s=6 + 42 * tiny.n_vessels / R.n_vessels.max(),
                   facecolor="none", edgecolor="#d55e00", lw=.9, zorder=4)
    ax.plot(per_k.index, per_k.values, color="#0072b2", lw=1.3, marker="o", ms=3.2, zorder=5,
            label="worst cluster at that K")
    ax.axhline(a.thresh, color="#999999", lw=.7, ls="--")
    ax.axvline(chosen, color="#cccccc", lw=.7, zorder=0)
    ax.text(KS[0] - .35, a.thresh + .012, f"{a.thresh:g}", fontsize=6, color="#666666", va="bottom")
    ax.text(chosen, 1.04, f"K = {chosen}", ha="center", fontsize=6.5, color="#0072b2")
    ax.set(xlabel="K", ylabel="per-cluster Jaccard",
           ylim=(0, 1.1), xticks=KS, yticks=np.arange(0, 1.01, .25))
    ax.set_title(f"{a.sub} subsamples of {a.frac:.0%} of the vessels, without replacement",
                 loc="left", pad=4, fontsize=6.5)
    ax.scatter([], [], facecolor="none", edgecolor="#8aa8c2", lw=.7, s=20, label="one cluster")
    ax.scatter([], [], facecolor="none", edgecolor="#d55e00", lw=.9, s=12,
               label=f"cluster ≤ {a.small} vessels")
    ax.legend(frameon=False, loc="lower left", handletextpad=.4, borderpad=0)
    save(fig, "niche_jaccard_by_k")

    # ---- figure 2: clustree -------------------------------------------------------------
    jmap = {(r.k, r.cluster): r.jaccard for r in R.itertuples()}
    pos = {}
    for n, k in enumerate(KS):
        ref = labels[k]
        if n == 0:
            o = np.argsort(np.bincount(ref, minlength=k))[::-1]    # biggest cluster first
        else:                                                      # barycentre of the parents
            k0 = KS[n - 1]
            C = contingency(labels[k0], ref, k0, k).astype(float)
            w = C / np.maximum(C.sum(0, keepdims=True), 1)
            bary = np.array([sum(w[i, j] * pos[(k0, i)] for i in range(k0)) for j in range(k)])
            o = np.argsort(bary, kind="stable")
        for x, c in enumerate(o):
            pos[(k, c)] = x - (k - 1) / 2
    fig, ax = plt.subplots(figsize=(80 * MM, 14 * len(KS) * MM))
    cmap = plt.get_cmap("viridis")
    for k0, k1 in zip(KS[:-1], KS[1:]):
        C = contingency(labels[k0], labels[k1], k0, k1)
        for i in range(k0):
            for j in range(k1):
                if not C[i, j]:
                    continue
                f = C[i, j] / C[i].sum()
                ax.plot([pos[(k0, i)], pos[(k1, j)]], [-k0, -k1], color="#777777",
                        lw=.3 + 2.2 * f, alpha=.15 + .55 * f, zorder=1,
                        solid_capstyle="round")
    nmax = R.n_vessels.max()
    for k in KS:
        for c in range(k):
            n = int((labels[k] == c).sum())
            j = jmap[(k, c)]
            ax.scatter(pos[(k, c)], -k, s=22 + 220 * n / nmax, color=cmap(j), zorder=3,
                       edgecolor="#333333" if j < 0.75 else "none", lw=.6)
            ax.text(pos[(k, c)], -k, str(n), ha="center", va="center", fontsize=5,
                    color="white" if j < .6 else "#222222", zorder=4)
    ax.set_yticks([-k for k in KS]); ax.set_yticklabels([f"K = {k}" for k in KS])
    ax.set_xticks([]); ax.spines["left"].set_visible(False); ax.spines["bottom"].set_visible(False)
    ax.tick_params(left=False)
    ax.set_xlim(-max(KS) / 2 - .6, max(KS) / 2 + .6)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=ax, fraction=.03, pad=.02)
    cb.set_label("per-cluster Jaccard"); cb.outline.set_linewidth(.4)
    ax.set_title("node = cluster (number = vessels), edge width = share of the K cluster "
                 "going to K+1", loc="left", pad=6, fontsize=6.5)
    save(fig, "niche_clustree")

    # ---- figure 3: consensus matrices, every K ------------------------------------------
    ncol = 4
    nrow = int(np.ceil(len(KS) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 38 * MM, nrow * 40 * MM))
    for ax in np.ravel(axes)[len(KS):]:
        ax.set_axis_off()
    for ax, k in zip(np.ravel(axes), KS):
        im = draw_consensus(ax, cons[k], labels[k], k,
                            f"K = {k}" + ("  ←" if k == chosen else ""))
    cb = fig.colorbar(im, ax=axes, fraction=.02, pad=.015)
    cb.set_label("co-clustering frequency"); cb.outline.set_linewidth(.4)
    save(fig, "niche_consensus_all_k")

    for k in KS:
        fig, ax = plt.subplots(figsize=(58 * MM, 54 * MM))
        im = draw_consensus(ax, cons[k], labels[k], k, f"consensus, K = {k}")
        cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.03)
        cb.set_label("co-clustering frequency"); cb.outline.set_linewidth(.4)
        save(fig, f"niche_consensus_k{k}")

    # ---- tables --------------------------------------------------------------------------
    os.makedirs(OUT, exist_ok=True)
    R.round(3).to_csv(f"{OUT}/niche_jaccard_per_cluster.csv", index=False)
    choice = per_k.rename("min_jaccard").reset_index()
    choice["all_clusters_stable"] = choice.min_jaccard >= a.thresh
    choice["smallest_cluster"] = [int(np.bincount(labels[k]).min()) for k in choice.k]
    choice["consensus_within"] = [
        float(np.mean([cons[k][np.ix_(labels[k] == i, labels[k] == i)].mean() for i in range(k)]))
        for k in choice.k]
    choice["consensus_between"] = [
        float(cons[k][labels[k][:, None] != labels[k][None, :]].mean()) for k in choice.k]
    choice["chosen"] = choice.k == chosen
    choice.round(3).to_csv(f"{OUT}/niche_k_choice.csv", index=False)
    print(choice.round(3).to_string(index=False))
    print()
    print("clusters at the chosen K and every small cluster in the scan:")
    show = R[(R.k == chosen) | (R.n_vessels <= a.small)]
    print(show.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
