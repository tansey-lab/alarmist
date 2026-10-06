#!/usr/bin/env python
"""Is the chosen K for the vessel niches defensible? Stability, consensus and perturbations.

The inertia curve this clustering was originally cut with is nearly linear, so KneeLocator's
K=9 is an artefact of its sensitivity parameter, and from K=4 one vessel (P17_AIS__6, a
peribronchial outlier) takes a cluster of its own. This measures what can actually be claimed:

  stability    mean adjusted Rand index between the reference labels and labels refitted on
               BOOT bootstrap resamples, per K -- how reproducible the partition is
  silhouette   per K, on the same features -- how separated the clusters are at all
  consensus    at the chosen K, how often each pair of vessels lands together over the
               resamples (the 300 x 300 matrix, ordered by cluster)
  perturbation ARI of the chosen K's labels against refits under: other seeds, three feature
               definitions (raw z, z clipped at 5, rare types pooled), leaving each section
               out, and 80 % subsamples

Features are the per-vessel ring compositions of vessel_niche_ring30_pervessel (proportions,
z-scored); --clip 5 is the default because unclipped z lets one 16-sigma vessel own a cluster.

Outputs (figures_niche_robustness/): niche_k_stability, niche_k_silhouette,
    niche_consensus_k<K>, niche_robustness.csv

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/niche_k_robustness.py
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
from sklearn.metrics import adjusted_rand_score, silhouette_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import log

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
SRC = f"{BASE}/vessel_niche_ring30_pervessel"
OUT = f"{BASE}/figures_niche_robustness"
MM = 1 / 25.4
BOOT = 200
KS = range(2, 13)

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


def features(P, share, mode, clip):
    if mode == "pooled":
        rare = share < 0.005
        X = np.column_stack([P[:, ~rare], P[:, rare].sum(1)])
    else:
        X = P
    Z = (X - X.mean(0)) / X.std(0)
    return np.clip(Z, -clip, clip) if mode in ("clipped", "pooled") else Z


def fit(Z, k, seed=0):
    return KMeans(k, random_state=seed, n_init=10).fit_predict(Z)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--k", type=int, default=3, help="the K being defended")
    ap.add_argument("--clip", type=float, default=5.0)
    ap.add_argument("--boot", type=int, default=BOOT)
    a = ap.parse_args()

    P = np.load(f"{SRC}/neighborhood_counts.npy")
    L = pd.read_csv(f"{SRC}/niche_labels.csv")
    share = (P * L.n_ring.to_numpy()[:, None]).sum(0)
    share = share / share.sum()
    Z = features(P, share, "clipped", a.clip)
    rng = np.random.default_rng(0)

    # ---- stability and silhouette per K ------------------------------------------------
    rows, consensus = [], None
    for k in KS:
        ref = fit(Z, k)
        aris, co, seen = [], np.zeros((len(Z), len(Z))), np.zeros((len(Z), len(Z)))
        for b in range(a.boot):
            idx = np.unique(rng.choice(len(Z), len(Z), replace=True))
            lb = fit(Z[idx], k, seed=b + 1)
            full = KMeans(k, random_state=b + 1, n_init=10).fit(Z[idx]).predict(Z)
            aris.append(adjusted_rand_score(ref, full))
            same = (lb[:, None] == lb[None, :]).astype(float)
            co[np.ix_(idx, idx)] += same
            seen[np.ix_(idx, idx)] += 1
        rows.append(dict(k=k, silhouette=silhouette_score(Z, ref),
                         stability_ari=float(np.mean(aris)),
                         stability_sd=float(np.std(aris)),
                         smallest_cluster=int(np.bincount(ref).min())))
        if k == a.k:
            consensus = co / np.maximum(seen, 1)
            ref_k = ref
        log(f"  K={k}: silhouette {rows[-1]['silhouette']:.3f}, stability "
            f"{rows[-1]['stability_ari']:.2f}, smallest cluster {rows[-1]['smallest_cluster']}")
    R = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(70 * MM, 48 * MM))
    ax.errorbar(R.k, R.stability_ari, yerr=R.stability_sd, color="#0072b2", lw=1.2,
                marker="o", ms=3, capsize=2, elinewidth=.6)
    ax.axvline(a.k, color="#999999", lw=.6)
    ax.set(xlabel="K", ylabel=f"stability (ARI, {a.boot} resamples)", ylim=(0, 1), xticks=list(KS))
    save(fig, "niche_k_stability")

    fig, ax = plt.subplots(figsize=(70 * MM, 48 * MM))
    ax.plot(R.k, R.silhouette, color="#d55e00", lw=1.2, marker="o", ms=3)
    ax.axhline(0.25, color="#999999", lw=.6)
    ax.axvline(a.k, color="#999999", lw=.6)
    ax.set(xlabel="K", ylabel="silhouette", ylim=(0, .35), xticks=list(KS))
    ax.text(KS[-1], .255, "0.25: weak structure below", ha="right", va="bottom", fontsize=6,
            color="#666666")
    save(fig, "niche_k_silhouette")

    # ---- consensus matrix at the chosen K ----------------------------------------------
    order = np.argsort(ref_k, kind="stable")
    fig, ax = plt.subplots(figsize=(62 * MM, 58 * MM))
    im = ax.imshow(consensus[np.ix_(order, order)], cmap="Blues", vmin=0, vmax=1)
    edges = np.cumsum(np.bincount(ref_k))[:-1]
    for e in edges:
        ax.axhline(e - .5, color="#d55e00", lw=.7)
        ax.axvline(e - .5, color="#d55e00", lw=.7)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"consensus, K = {a.k}", loc="left", pad=4)
    cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.03)
    cb.set_label("co-clustering frequency")
    cb.outline.set_linewidth(.4)
    save(fig, f"niche_consensus_k{a.k}")
    within = np.mean([consensus[np.ix_(ref_k == i, ref_k == i)].mean() for i in range(a.k)])
    between = consensus[(ref_k[:, None] != ref_k[None, :])].mean()
    log(f"consensus at K={a.k}: within-cluster {within:.2f}, between-cluster {between:.2f}")

    # ---- perturbations ------------------------------------------------------------------
    pert = []
    for s in range(1, 6):
        pert.append(("seed " + str(s), adjusted_rand_score(ref_k, fit(Z, a.k, seed=s))))
    for mode in ("raw", "pooled"):
        Zi = features(P, share, mode, a.clip)
        pert.append((f"features: {mode}", adjusted_rand_score(ref_k, fit(Zi, a.k))))
    for s in sorted(L.sample_id.unique()):
        m = (L.sample_id != s).to_numpy()
        lb = KMeans(a.k, random_state=0, n_init=10).fit(Z[m]).predict(Z)
        pert.append((f"drop {s}", adjusted_rand_score(ref_k, lb)))
    sub = []
    for b in range(50):
        idx = rng.choice(len(Z), int(.8 * len(Z)), replace=False)
        lb = KMeans(a.k, random_state=b, n_init=10).fit(Z[idx]).predict(Z)
        sub.append(adjusted_rand_score(ref_k, lb))
    pert.append(("80 % subsample (mean of 50)", float(np.mean(sub))))
    PT = pd.DataFrame(pert, columns=["perturbation", "ARI_vs_reference"])
    os.makedirs(OUT, exist_ok=True)
    R.round(4).to_csv(f"{OUT}/niche_robustness.csv", index=False)
    PT.round(3).to_csv(f"{OUT}/niche_perturbations_k{a.k}.csv", index=False)
    print(R.round(3).to_string(index=False))
    print(PT.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
