#!/usr/bin/env python
"""Do the cell-type niches line up with the four vasculature motifs? One heatmap per K.

For every K in the Jaccard scan, the niche labels are refit exactly as in `niche_k_jaccard.py`
(ring proportions, z-scored, not clipped, k-means seed 0) and each cluster is correlated
with each motif's per-vessel ring positive fraction:

    r(cluster i, motif k) = Pearson r between the 0/1 membership of cluster i and f_k(v)

over the 300 curated vessels -- the point-biserial correlation, i.e. how much being in that
niche raises that motif's share of the ring. The best-matching motif of each cluster is boxed,
and the per-K headline is the mean of those best matches: how well, at that K, the niches can
be said to correspond to the motifs at all.

Motifs: the four picked vasculature motifs, labelled by name, not index -- SMC vascular
stabilization (m1), Vascular homeostasis (m23), Healthy alveolar (m24), Tumor vasculature (m10)
-- in the usual figure order and with the usual colours (`motif_state_palette.ANCHOR` / `NAME`). Fractions come from
vessel_kmeans_ring30_m0-1-6-10-23-24_frac/features.npy, joined to the niches by vessel_id.

Colour: Reds, 0 to the largest |r| in the scan, shared by every panel. Negative correlations
are left white -- the number is still printed in the cell.

Outputs (figures_niche_robustness/): niche_motif_corr_k<K> for every K,
    niche_motif_corr_all_k, niche_motif_best_match (K x motif, each cell the r of whichever
    niche matches that motif best at that K; _slim is the same without the cluster sizes),
    niche_motif_corr.csv, niche_motif_best_match.csv

Usage:
    /home/fanj2/.conda/envs/spatial/bin/python \
        /home/fanj2/alarmist/scripts/research/luad_revision/plot_niche_motif_correlation.py
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
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import log
from motif_state_palette import ANCHOR, NAME, SHORT
from niche_k_jaccard import OUT, SRC, fit, save, zscore

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
FRAC = f"{BASE}/vessel_kmeans_ring30_m0-1-6-10-23-24_frac"
FEAT_MOTIFS = [0, 1, 6, 10, 23, 24]          # column order of features.npy
MOTIFS = [1, 23, 24, 10]                     # what is drawn, in the usual order
LABEL = {m: NAME[m].replace(" motif", "") for m in MOTIFS}
MM = 1 / 25.4

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def corr(lab, k, Fm):
    """Pearson r of every (cluster, motif) pair; clusters are 0/1 indicators."""
    R = np.empty((k, Fm.shape[1]))
    for i in range(k):
        x = (lab == i).astype(float)
        R[i] = [np.corrcoef(x, Fm[:, j])[0, 1] for j in range(Fm.shape[1])]
    return R


def heat(ax, R, n, vmax, title, ylab=True, xlab=True):
    im = ax.imshow(R, cmap="Reds", vmin=0, vmax=vmax, aspect="auto")
    best = R.argmax(1)
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            ax.text(j, i, f"{R[i, j]:.2f}", ha="center", va="center", fontsize=5.8,
                    color="white" if R[i, j] > .62 * vmax else "#333333")
        ax.add_patch(Rectangle((best[i] - .5, i - .5), 1, 1, fill=False,
                               edgecolor="#111111", lw=1.1, zorder=4))
    ax.set_xticks(range(len(MOTIFS)))
    if xlab:
        ax.set_xticklabels([LABEL[m] for m in MOTIFS], fontsize=6, rotation=38,
                           ha="right", rotation_mode="anchor")
        for t, m in zip(ax.get_xticklabels(), MOTIFS):
            t.set_color(ANCHOR[m]); t.set_fontweight("bold")
    else:
        ax.set_xticklabels([])
    ax.set_yticks(range(R.shape[0]))
    ax.set_yticklabels([f"niche {i}  (n={n[i]})" for i in range(R.shape[0])] if ylab
                       else [f"{i}" for i in range(R.shape[0])])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_linewidth(.4); sp.set_color("#999999")
    ax.set_title(title, loc="left", pad=4, fontsize=7.5)
    return im


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=9)
    ap.add_argument("--clip", type=float, default=float("inf"),
                    help="clip z at +/- this; default is no clipping, as in niche_k_jaccard")
    a = ap.parse_args()
    KS = list(range(a.kmin, a.kmax + 1))

    P = np.load(f"{SRC}/neighborhood_counts.npy")
    L = pd.read_csv(f"{SRC}/niche_labels.csv")
    V = pd.read_csv(f"{FRAC}/vessels.csv")
    F = pd.DataFrame(np.load(f"{FRAC}/features.npy"),
                     columns=[f"m{k}" for k in FEAT_MOTIFS])
    F["vessel_id"] = V.vessel_id.to_numpy()
    F = F.set_index("vessel_id").loc[L.vessel_id].reset_index()     # join by id, not position
    Fm = F[[f"m{m}" for m in MOTIFS]].to_numpy()
    Z = zscore(P, a.clip)
    log(f"{len(Z)} vessels; motifs {MOTIFS}; K {a.kmin}..{a.kmax}")

    Rs, ns, rows = {}, {}, []
    for k in KS:
        lab = fit(Z, k)
        Rs[k], ns[k] = corr(lab, k, Fm), np.bincount(lab, minlength=k)
        for i in range(k):
            best = int(Rs[k][i].argmax())
            rows.append(dict(k=k, niche=i, n_vessels=int(ns[k][i]),
                             best_motif=NAME[MOTIFS[best]],
                             best_r=Rs[k][i][best],
                             **{NAME[m]: Rs[k][i][j] for j, m in enumerate(MOTIFS)}))
        log(f"  K={k}: mean best-match r {Rs[k].max(1).mean():.2f}, "
            f"best motifs {[SHORT[MOTIFS[j]] for j in Rs[k].argmax(1)]}")
    R = pd.DataFrame(rows)
    vmax = float(np.ceil(max(np.abs(Rs[k]).max() for k in KS) * 20) / 20)

    for k in KS:
        fig, ax = plt.subplots(figsize=(52 * MM, (14 + 7 * k) * MM))
        im = heat(ax, Rs[k], ns[k], vmax,
                  f"K = {k}   mean best match r = {Rs[k].max(1).mean():.2f}")
        cb = fig.colorbar(im, ax=ax, fraction=.055, pad=.04)
        cb.set_label("Pearson r"); cb.outline.set_linewidth(.4)
        save(fig, f"niche_motif_corr_k{k}")

    ncol = 4
    nrow = int(np.ceil(len(KS) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 48 * MM, nrow * 56 * MM))
    for ax in np.ravel(axes)[len(KS):]:
        ax.set_axis_off()
    for n, (ax, k) in enumerate(zip(np.ravel(axes), KS)):
        im = heat(ax, Rs[k], ns[k], vmax,
                  f"K = {k}   mean best r = {Rs[k].max(1).mean():.2f}", ylab=False,
                  xlab=n >= len(KS) - ncol)
        ax.set_ylabel("niche")
    fig.subplots_adjust(hspace=.22)
    cb = fig.colorbar(im, ax=axes, fraction=.02, pad=.015)
    cb.set_label("Pearson r"); cb.outline.set_linewidth(.4)
    fig.suptitle("niche x motif correlation; box = that niche's best-matching motif",
                 x=.005, ha="left", fontsize=7.5)
    save(fig, "niche_motif_corr_all_k")

    # ---- the best match per motif, every K in one heatmap --------------------------------
    B = np.array([[Rs[k][:, j].max() for j in range(len(MOTIFS))] for k in KS])
    who = np.array([[int(Rs[k][:, j].argmax()) for j in range(len(MOTIFS))] for k in KS])
    for stem, with_n in (("niche_motif_best_match", True),
                         ("niche_motif_best_match_slim", False)):
        fig, ax = plt.subplots(figsize=(62 * MM, (16 + 8 * len(KS)) * MM))
        im = ax.imshow(B, cmap="Reds", vmin=0, vmax=vmax, aspect="auto")
        for r in range(len(KS)):
            for c in range(len(MOTIFS)):
                dark = B[r, c] > .62 * vmax
                sub = f"niche {who[r, c]}" + (f", n={ns[KS[r]][who[r, c]]}" if with_n else "")
                ax.text(c, r - .12, f"{B[r, c]:.2f}", ha="center", va="center", fontsize=6.4,
                        color="white" if dark else "#333333")
                ax.text(c, r + .22, sub, ha="center", va="center",
                        fontsize=4.6 if with_n else 5.2,
                        color="#f0f0f0" if dark else "#777777")
        ax.set_xticks(range(len(MOTIFS)))
        ax.set_xticklabels([LABEL[m] for m in MOTIFS], fontsize=6, rotation=38, ha="right",
                           rotation_mode="anchor")
        for t, m in zip(ax.get_xticklabels(), MOTIFS):
            t.set_color(ANCHOR[m]); t.set_fontweight("bold")
        ax.set_yticks(range(len(KS)))
        ax.set_yticklabels([f"K = {k}" for k in KS])
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_linewidth(.4); sp.set_color("#999999")
        ax.set_title("best-matching niche per motif", loc="left", pad=4, fontsize=7.5)
        cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.04)
        cb.set_label("Pearson r of that niche"); cb.outline.set_linewidth(.4)
        save(fig, stem)
    BT = pd.DataFrame([dict(k=k, motif=NAME[m], best_niche=int(who[i, j]),
                            n_vessels=int(ns[k][who[i, j]]), best_r=B[i, j])
                       for i, k in enumerate(KS) for j, m in enumerate(MOTIFS)])
    BT.round(3).to_csv(f"{OUT}/niche_motif_best_match.csv", index=False)

    os.makedirs(OUT, exist_ok=True)
    R.round(3).to_csv(f"{OUT}/niche_motif_corr.csv", index=False)
    print(R.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
