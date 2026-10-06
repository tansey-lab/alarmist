#!/usr/bin/env python
"""Per motif: positive fraction against vessel size, and against distance to the vessel.

Two figures per motif, 25 motifs each.

  size_vs_positive/m<k>      one point per vessel (n=300): x = lumen area (log, it spans
                             147-333,400 um2), y = fraction of that vessel's 30 um ring
                             cells positive for the motif. Colour AND marker encode section,
                             so identity is never colour-alone. Spearman is reported overall
                             and per section, because section is a known confound.

  distance_vs_positive/m<k>  one point per CELL: x = distance to the nearest curated lumen
                             boundary, y = the GMM state (0/1). Smoothed with a binomial GLM
                             on a natural-cubic-spline basis (df=6), fitted per section.
                             Cells are binned at 2 um and the GLM is fitted to the
                             (positive, total) counts per bin -- binomial sufficiency makes
                             that identical to fitting every cell, and much faster.

The same fit is also run over NON-VASCULAR cells only (Vas_Endo, Lym_Endo, SMC, Pericyte
dropped) and reported in distance_vs_positive.csv, but is no longer drawn. It still matters
when reading the figures: 83 % of cells 0-10 um from a lumen ARE vessel-wall cells, against a
16.5 % tissue baseline, so a curve that rises toward the lumen may only be tracking cell-type
composition. The CSV's delta_nonvasc column is what separates the two cases.

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/plot_motif_size_distance.py
"""
from __future__ import annotations

import argparse
import textwrap
import os
import sys

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from patsy import dmatrix, build_design_matrices
from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_state_palette import ANCHOR as STATE_COL, NAME as STATE_NAME
from motif_vessel_alignment import (H5, NM, SECTIONS, VASCULAR, read_cat, load_polys,
                                    densify, signed_distance, log)
from lumen_geometry import perimeter

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
OUT = f"{BASE}/figures_motif_profiles"
MM = 1 / 25.4
SINGLE = 89 * MM
R_RING, MIN_RING = 30.0, 20
ANN = {0: "T cell", 1: "SMC", 6: "macrophage", 10: "tumour-endo",
       23: "endothelial", 24: "alveolar"}
NAME = {1: "SMC vascular stabilization motif", 10: "Tumor vasculature motif",     # user-given names,
        23: "Vascular homeostasis motif", 24: "Healthy alveolar motif"}  # 2026-09-18
SECCOL = {"P17_AIS": "#0072b2", "P17_LUAD": "#e69f00",
          "P21_AIS": "#009e73", "P21_LUAD": "#cc79a7"}      # Okabe-Ito, validated
SECMK = {"P17_AIS": "o", "P17_LUAD": "s", "P21_AIS": "^", "P21_LUAD": "D"}
DBIN, DMAX, SPL_DF = 2.0, 300.0, 6
XMAX_DIST = 100.0     # displayed x range of the per-motif distance figures; fit uses DMAX
# Vessel size metric. Area shrinks when a lumen is flattened during sectioning (veins collapse
# more than arteries); the traced perimeter does not, so P/pi -- the diameter of the circle with
# the same circumference -- is the squish-invariant size. --size switches between them.
SZ = "area"
SIZE_LAB = r"lumen area ($\mu$m$^2$, log scale)"
MIN_BIN = 30            # a 2 um bin thinner than this in a section is dropped from the fit

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
    "axes.linewidth": .5, "xtick.major.width": .5, "ytick.major.width": .5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def save(fig, sub, stem):
    d = f"{OUT}/{sub}"
    os.makedirs(d, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{d}/{stem}.{ext}", dpi=400 if ext == "png" else None)
    plt.close(fig)


def spline_fit(d, y, keep):
    """P(positive | distance) by binomial GLM on a natural-cubic-spline basis.

    Cells are binned at DBIN um first; a binomial GLM on the per-bin (successes, trials) is
    the same likelihood as one row per cell, since the counts are sufficient.
    Returns (grid, p, lo, hi) with the GLM's pointwise 95 % CI for the mean, or None when
    the section has too little data. That CI is the SPLINE FIT's uncertainty only. It is
    narrow here because each section contributes 70-110 k cells, and it does NOT account for
    spatial autocorrelation between neighbouring cells, so it understates the real
    uncertainty. The spread BETWEEN the four section curves is the larger and more honest
    one, and is drawn separately as a band on the mean.
    """
    d, y = d[keep], y[keep]
    m = (d > 0) & (d <= DMAX)
    if m.sum() < 500:
        return None
    b = np.floor(d[m] / DBIN).astype(int)
    nb = int(DMAX / DBIN)
    tot = np.bincount(b, minlength=nb)[:nb].astype(float)
    pos = np.bincount(b, weights=y[m].astype(float), minlength=nb)[:nb]
    ok = tot >= MIN_BIN
    if ok.sum() < SPL_DF + 2:
        return None
    x = (np.arange(nb) + .5) * DBIN
    X = dmatrix(f"cr(x, df={SPL_DF})", {"x": x[ok]}, return_type="dataframe")
    try:
        r = sm.GLM(np.column_stack([pos[ok], tot[ok] - pos[ok]]), X,
                   family=sm.families.Binomial()).fit()
    except Exception:
        return None
    g = np.linspace(x[ok].min(), x[ok].max(), 300)
    Xg = build_design_matrices([X.design_info], {"x": g})[0]
    try:
        sf = r.get_prediction(Xg).summary_frame(alpha=0.05)
        lo = np.asarray(sf["mean_ci_lower"]); hi = np.asarray(sf["mean_ci_upper"])
    except Exception:
        lo = hi = None
    return g, np.asarray(r.predict(Xg)), lo, hi


VASC_MOTIFS = [0, 1, 6, 10, 23, 24]
MCOL = ["#0072b2", "#e69f00", "#009e73", "#cc79a7", "#56b4e9", "#d55e00"]  # validated PASS
JIT, RUG_Y = 0.045, -0.14     # y-jitter for the 0/1 cloud; height of the x-axis rug
NPOINTS = 20000               # cells drawn; the FIT always uses all of them
THRESH = 0.5                  # ring positive fraction above which a vessel is "positive"
POSC, NEGC = "#b2182b", "#2166ac"   # validated PASS, --pairs all
KDE_MIN = 15                  # below this a group gets a rug only; a KDE would invent shape


def pooled_distance_figure(D, S, k, col, n_points=NPOINTS, seed=0):
    """All four sections pooled into one curve, over the RAW per-cell observations.

    Every cell is a point at y = its GMM state (0 or 1), y-jittered so the two bands have
    thickness, plus a rug of the same points' x-projections on the axis. Binned means are
    deliberately NOT drawn: they lie on the fitted curve by construction and show nothing
    the curve does not.

    The fit uses every cell in range (341,289). The point cloud is a seeded random subsample
    of n_points, because 341 k marks render as two solid bars and as a ~50 MB vector file;
    the scatter and rug are rasterized so axis text stays vector. Subsampling affects only
    what is drawn.

    Pooling weights sections by their cell count -- P17_LUAD contributes 110 k of the
    341 k -- so the curve leans on it. The per-section figures are what show the
    disagreement.
    """
    r = spline_fit(D, S[k], np.ones(len(D), bool))
    if r is None:
        return None
    g, p, lo, hi = r

    idx = np.where((D > 0) & (D <= DMAX))[0]
    rng = np.random.default_rng(seed)
    sub = idx if len(idx) <= n_points else rng.choice(idx, n_points, replace=False)
    xs = D[sub]
    ys = S[k][sub].astype(float) + rng.uniform(-JIT, JIT, len(sub))

    fig, ax = plt.subplots(figsize=(SINGLE, 62 * MM))
    ax.scatter(xs, ys, s=1.1, color="#3f3f3f", alpha=.20, lw=0, rasterized=True, zorder=1)
    # at 20 k marks over 300 um any opaque rug is a solid bar; alpha this low is what makes
    # it read as a density instead
    ax.scatter(xs, np.full(len(xs), RUG_Y), s=11, marker="|", color="#1a1a1a", alpha=.025,
               lw=.3, rasterized=True, zorder=1)
    if lo is not None:
        ax.fill_between(g, lo, hi, color=col, alpha=.35, lw=0, zorder=3)
    ax.plot(g, p, color=col, lw=2.0, zorder=4)
    ax.set_xlim(0, DMAX); ax.set_ylim(RUG_Y - .05, 1.12)
    ax.set_yticks([0, .5, 1])
    ax.set_xlabel(r"distance from the lumen boundary ($\mu$m)")
    ax.set_ylabel("P(motif positive)")
    ax.set_title(f"m{k}" + (f" {ANN[k]}" if k in ANN else ""), loc="left", pad=5)
    save(fig, "distance_vs_positive_pooled", f"m{k}")
    return dict(motif=f"m{k}", annotation=ANN.get(k, ""), n_cells=int(len(idx)),
                n_drawn=int(len(sub)),
                p_at_10um=float(np.interp(10, g, p)),
                p_at_150um=float(np.interp(150, g, p)),
                delta=float(np.interp(10, g, p) - np.interp(150, g, p)))


def _pfmt(pv):
    """P as a journal would print it."""
    if not np.isfinite(pv):
        return ""
    if pv < 1e-3:
        e = int(np.floor(np.log10(pv)))
        return rf"$P$ = {pv / 10 ** e:.0f}$\times$10$^{{{e}}}$"
    return rf"$P$ = {pv:.3f}"


def size_kde(x, y, k, lab):
    """Density of log lumen area for vessels called positive / negative at THRESH.

    Same x axis as the size scatter. A group smaller than KDE_MIN gets a rug only: a
    Gaussian KDE over a handful of points draws a shape the data does not support, and one
    motif has such a group (m23: 8 negative vessels of 298).
    """
    pos = y >= THRESH
    n1, n2 = int(pos.sum()), int((~pos).sum())
    # Mann-Whitney U on log area. It is rank-based, so U and P are IDENTICAL on raw area
    # (verified over all 25 motifs); the log only sets the axis. AUC = U / (n1 n2) is
    # P(a positive vessel is larger than a negative one), the two-sample counterpart of the
    # Spearman rho reported on the scatter.
    if n1 >= 2 and n2 >= 2:
        u = mannwhitneyu(x[pos], x[~pos])
        U, auc, pv = float(u.statistic), float(u.statistic / (n1 * n2)), float(u.pvalue)
    else:
        U = auc = pv = float("nan")
    grid = np.linspace(x.min() - .15, x.max() + .15, 400)
    fig, ax = plt.subplots(figsize=(SINGLE, 56 * MM))
    peak, out = 0.0, {}
    for name, m, c in (("positive", pos, POSC), ("negative", ~pos, NEGC)):
        n = int(m.sum())
        out[f"n_{name}"] = n
        out[f"median_{SZ}_{name}"] = float(10 ** np.median(x[m])) if n else float("nan")
        if n >= KDE_MIN and np.ptp(x[m]) > 0:
            d = gaussian_kde(x[m])(grid)
            peak = max(peak, float(d.max()))
            ax.plot(grid, d, color=c, lw=1.6, label=f"{name} (n={n})", zorder=3)
        elif n:
            ax.plot([], [], color=c, lw=1.6, label=f"{name} (n={n})")
    peak = peak or 1.0
    for j, (name, m, c) in enumerate((("positive", pos, POSC), ("negative", ~pos, NEGC))):
        ax.plot(x[m], np.full(int(m.sum()), -(.045 + .05 * j) * peak), "|", color=c,
                ms=4, mew=.6, zorder=2)
    ax.set_ylim(-.13 * peak, peak * 1.10)
    tk = np.arange(np.floor(x.min()), np.ceil(x.max()) + 1)
    ax.set_xticks(tk)
    ax.set_xticklabels([("%g" % (10 ** t)) if t < 3 else f"$10^{{{int(t)}}}$" for t in tk])
    ax.set_xlim(grid[0], grid[-1])
    ax.set_xlabel(SIZE_LAB)
    ax.set_ylabel("density")
    ax.set_title(lab, loc="left", pad=5)
    if np.isfinite(auc):
        ax.set_title(f"AUC {auc:.2f}   " + _pfmt(pv), loc="right", fontsize=6.5,
                     color="#444444", pad=5)
    ax.legend(frameon=False, loc="best", handlelength=1.4)
    save(fig, "size_kde", f"m{k}")
    out["mannwhitney_u"], out["auc"], out["p"] = U, auc, pv
    return dict(motif=f"m{k}", annotation=ANN.get(k, ""), **out)


def density_by_motif(V, RF, motifs, thresh, out):
    """One column per motif: the violin cut in half and laid on its side.

    x is log lumen area, y the density of that motif's positive vessels (all sections
    pooled). Every column shares the same x range and the same y range, so the panels are
    directly comparable. Vessels are drawn as ticks under each curve and the median as a
    black bar. A group below KDE_MIN gets ticks only -- no KDE is drawn.
    """
    x = np.log10(V['size'].to_numpy())
    grid = np.linspace(x.min() - .15, x.max() + .15, 400)
    dens, vals, peak = {}, {}, 0.0
    for k in motifs:
        v = x[RF[:, k] >= thresh]
        vals[k] = v
        if len(v) >= KDE_MIN and np.ptp(v) > 0:
            dens[k] = gaussian_kde(v)(grid)
            peak = max(peak, float(dens[k].max()))
    peak = peak or 1.0
    fig, axes = plt.subplots(1, len(motifs), figsize=(len(motifs) * 62 * MM, 30 * MM),
                             sharex=True, sharey=True, squeeze=False)
    rows = []
    for j, k in enumerate(motifs):
        ax = axes[0, j]
        c = STATE_COL.get(k, MCOL[j % len(MCOL)])     # same anchor as the UpSet / maps
        v, d = vals[k], dens.get(k)
        if d is not None:
            ax.fill_between(grid, 0, d, color=c, alpha=.45, lw=0, zorder=2)
            ax.plot(grid, d, color=c, lw=1.1, zorder=3)
        if len(v):
            ax.plot(v, np.full(len(v), -.09 * peak), "|", color=c, ms=3.6, mew=.6, alpha=.8,
                    zorder=3)
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            for q, lw, h in ((q1, .8, .10), (med, 1.4, .16), (q3, .8, .10)):
                ax.plot([q] * 2, [-.03 * peak, h * peak], color="#1a1a1a", lw=lw, zorder=4)
        ax.set_ylim(-.15 * peak, peak * 1.12)
        ax.set_xlim(grid[0], grid[-1])
        tk = np.arange(np.floor(x.min()), np.ceil(x.max()) + 1)
        ax.set_xticks(tk)
        ax.set_xticklabels([("%g" % (10 ** t)) if t < 3 else f"$10^{{{int(t)}}}$" for t in tk])
        ax.set_title(f"m{k}   n = {len(v)}", loc="left", fontsize=6.5, pad=3)
        ax.set_xlabel(SIZE_LAB, fontsize=6.5)
        if j == 0:
            ax.set_ylabel("density")
        rows.append(dict(motif=f"m{k}", annotation=ANN.get(k, ""), name=NAME.get(k, ""),
                         n_positive=len(v),
                         **{f"q1_{SZ}": float(10 ** np.percentile(v, 25)) if len(v) else np.nan,
                            f"median_{SZ}": float(10 ** np.median(v)) if len(v) else np.nan,
                            f"q3_{SZ}": float(10 ** np.percentile(v, 75)) if len(v) else np.nan},
                         kde_drawn=d is not None))
    fig.subplots_adjust(bottom=.34, wspace=.20)
    save(fig, "size_kde", f"density_thresh{thresh:g}")
    T = pd.DataFrame(rows)
    T.round(4).to_csv(f"{out}/size_kde_density.csv", index=False)
    print(T.round(2).to_string(index=False))
    return T


def combined_size_kde(V, RF, motifs, thresh, out):
    """Lumen-area density of the POSITIVE vessels only, one curve per motif, one axes.

    Same x axis and KDE as size_kde(), but the negative group is not drawn and the cut is
    --kde-thresh (0.9) rather than THRESH (0.5). Curve colours are the VASC_MOTIFS colours
    of the pooled distance figure, so a motif keeps its colour across the two figures.
    """
    x = np.log10(V['size'].to_numpy())
    grid = np.linspace(x.min() - .15, x.max() + .15, 400)
    fig, ax = plt.subplots(figsize=(SINGLE, 56 * MM))
    peak, rows = 0.0, []
    for j, k in enumerate(motifs):
        c = MCOL[VASC_MOTIFS.index(k)] if k in VASC_MOTIFS else MCOL[j % len(MCOL)]
        m = RF[:, k] >= thresh
        n = int(m.sum())
        lab = f"m{k}" + (f" {ANN[k]}" if k in ANN else "") + f" (n={n})"
        if n >= KDE_MIN and np.ptp(x[m]) > 0:
            d = gaussian_kde(x[m])(grid)
            peak = max(peak, float(d.max()))
            ax.plot(grid, d, color=c, lw=1.6, label=lab, zorder=3)
        elif n:
            ax.plot([], [], color=c, lw=1.6, label=lab)
        rows.append(dict(motif=f"m{k}", annotation=ANN.get(k, ""), threshold=thresh,
                         n_positive=n, n_vessels=int(len(x)),
                         **{f"median_{SZ}_positive": float(10 ** np.median(x[m])) if n
                            else float("nan")}))
    peak = peak or 1.0
    for j, k in enumerate(motifs):
        c = MCOL[VASC_MOTIFS.index(k)] if k in VASC_MOTIFS else MCOL[j % len(MCOL)]
        m = RF[:, k] >= thresh
        ax.plot(x[m], np.full(int(m.sum()), -(.045 + .045 * j) * peak), "|", color=c,
                ms=4, mew=.6, alpha=.7, zorder=2)
    ax.set_ylim(-(.06 + .045 * len(motifs)) * peak, peak * 1.10)
    ax.set_yticks([t for t in ax.get_yticks() if t >= 0])   # the rugs sit below zero
    tk = np.arange(np.floor(x.min()), np.ceil(x.max()) + 1)
    ax.set_xticks(tk)
    ax.set_xticklabels([("%g" % (10 ** t)) if t < 3 else f"$10^{{{int(t)}}}$" for t in tk])
    ax.set_xlim(grid[0], grid[-1])
    ax.set_xlabel(SIZE_LAB)
    ax.set_ylabel("density")
    ax.set_title(f"positive fraction $\geq$ {thresh:g}", loc="left", pad=5)
    ax.legend(frameon=False, loc="best", handlelength=1.4)
    save(fig, "size_kde", f"combined_thresh{thresh:g}")
    # ---- the same groups as violins ------------------------------------------------
    groups = [x[RF[:, k] >= thresh] for k in motifs]
    cols = [MCOL[VASC_MOTIFS.index(k)] if k in VASC_MOTIFS else MCOL[j % len(MCOL)]
            for j, k in enumerate(motifs)]
    pos = np.arange(len(motifs)) + 1

    def draw_violins(groups, cols, labels, slot=26, height=56):
        fig, ax = plt.subplots(figsize=(max(SINGLE, slot * MM * len(groups)), height * MM))
        pos = np.arange(len(groups)) + 1
        parts = ax.violinplot([g for g in groups], positions=pos, widths=.8,
                              showextrema=False, showmedians=False)
        for b, c in zip(parts["bodies"], cols):
            b.set_facecolor(c); b.set_edgecolor(c); b.set_alpha(.35); b.set_lw(.8)
        rng = np.random.RandomState(0)      # jitter only; nothing statistical rides on it
        for j, (g, c) in enumerate(zip(groups, cols)):
            ax.plot(pos[j] + rng.uniform(-.13, .13, len(g)), g, "o", ms=1.8, mfc=c,
                    mec="none", alpha=.45, zorder=2)
            q1, med, q3 = np.percentile(g, [25, 50, 75])
            ax.plot([pos[j], pos[j]], [q1, q3], color="#1a1a1a", lw=1.1, zorder=3,
                    solid_capstyle="butt")
            ax.plot(pos[j], med, "_", color="#1a1a1a", ms=9, mew=1.4, zorder=4)
        ax.set_xticks(pos)
        ax.set_xticklabels(labels, fontsize=6.5, linespacing=1.3)
        tk = np.arange(np.floor(x.min()), np.ceil(x.max()) + 1)
        ax.set_yticks(tk)
        ax.set_yticklabels([("%g" % (10 ** t)) if t < 3 else f"$10^{{{int(t)}}}$"
                            for t in tk])
        ax.set_ylabel(SIZE_LAB)
        ax.set_title(f"positive fraction $\geq$ {thresh:g}", loc="left", pad=5)
        return fig, ax

    fig, ax = draw_violins(groups, cols,
                           [f"m{k}" + (f"\n{ANN[k]}" if k in ANN else "") + f"\n(n={len(g)})"
                            for k, g in zip(motifs, groups)])
    save(fig, "size_kde", f"combined_thresh{thresh:g}_violin")

    # the _sig figure: violins sorted by median lumen area (largest left), user names
    srt = sorted(range(len(motifs)), key=lambda j: -np.median(groups[j]))
    groups, cols, motifs = [groups[j] for j in srt], [cols[j] for j in srt], \
        [motifs[j] for j in srt]
    pos = np.arange(len(motifs)) + 1

    # ---- same figure, Holm-corrected pairwise Mann-Whitney brackets -----------------
    # Mann-Whitney on log area = Mann-Whitney on raw area (rank based); Holm over the
    # len(motifs) choose 2 pairs, the same multipletests used in core/glm.py.
    pairs, pv = [], []
    for i in range(len(motifs)):
        for j in range(i + 1, len(motifs)):
            n1, n2 = len(groups[i]), len(groups[j])
            u = mannwhitneyu(groups[i], groups[j])
            pairs.append(dict(motif_a=f"m{motifs[i]}", motif_b=f"m{motifs[j]}",
                              name_a=NAME.get(motifs[i], ""), name_b=NAME.get(motifs[j], ""),
                              i=i, j=j, n_a=n1, n_b=n2,
                              mannwhitney_u=float(u.statistic),
                              auc=float(u.statistic / (n1 * n2)),
                              **{f"median_{SZ}_a": float(10 ** np.median(groups[i])),
                                 f"median_{SZ}_b": float(10 ** np.median(groups[j]))},
                              p=float(u.pvalue)))
            pv.append(float(u.pvalue))
    padj = multipletests(pv, method="holm")[1]
    for d, q in zip(pairs, padj):
        d["p_holm"] = float(q)

    def stars(q):
        return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < .05 else "ns"

    fig, ax = draw_violins(groups, cols,
                           [textwrap.fill(NAME.get(k, f"m{k}"), 18) + f"\n(m{k}, n={len(g)})"
                            for k, g in zip(motifs, groups)], slot=32, height=74)
    sig = sorted(pairs,                  # every tested pair; non-significant ones read "ns"
                 key=lambda d: (d["j"] - d["i"], d["i"]))
    top = max(g.max() for g in groups)
    step = .14 * (top - min(g.min() for g in groups))   # one bracket + its stars per level
    for lvl, d in enumerate(sig):
        y = top + step * (lvl + 1)
        a_, b_ = pos[d["i"]], pos[d["j"]]
        ax.plot([a_, a_, b_, b_], [y - step * .22, y, y, y - step * .22],
                color="#1a1a1a", lw=.7, clip_on=False, zorder=5)
        ax.text((a_ + b_) / 2, y + step * .05, stars(d["p_holm"]), ha="center",
                va="bottom", fontsize=7.5, color="#1a1a1a", zorder=5)
    if sig:
        ax.set_ylim(top=top + step * (len(sig) + .9))
        tk2 = [t for t in np.arange(np.floor(x.min()), np.ceil(x.max()) + 1) if t <= top]
        ax.set_yticks(tk2)                       # no labelled tick inside the bracket zone
        ax.set_yticklabels([("%g" % (10 ** t)) if t < 3 else f"$10^{{{int(t)}}}$"
                            for t in tk2])
    save(fig, "size_kde", f"combined_thresh{thresh:g}_violin_sig")
    PW = pd.DataFrame(pairs).drop(columns=["i", "j"])
    PW.round(5).to_csv(f"{out}/size_kde_combined_pairwise.csv", index=False)
    print(PW.round(4).to_string(index=False))

    density_by_motif(V, RF, motifs, thresh, out)

    T = pd.DataFrame(rows)
    T.round(4).to_csv(f"{out}/size_kde_combined.csv", index=False)
    print(T.round(3).to_string(index=False))
    return T


def main():
    global OUT, SZ, SIZE_LAB
    ap = argparse.ArgumentParser()
    ap.add_argument("--motifs", default="all")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--kde-combined", action="store_true",
                    help="only the positive-vessel lumen-area densities of --kde-motifs, "
                         "all in one axes, at --kde-thresh")
    ap.add_argument("--kde-thresh", type=float, default=0.9,
                    help="ring positive fraction for --kde-combined")
    ap.add_argument("--kde-motifs", default="1,10,23,24")
    ap.add_argument("--size", choices=("area", "diameter"), default="area",
                    help="vessel size metric on the x axis: traced lumen area, or the "
                         "squish-invariant perimeter diameter P/pi")
    ap.add_argument("--pooled", action="store_true",
                    help="only the six vasculature motifs, all sections pooled into one "
                         "curve, with the binned observations drawn as dots")
    a = ap.parse_args()
    OUT = a.out
    if a.size == "diameter":
        SZ, SIZE_LAB = "diam", r"lumen diameter P/$\pi$ ($\mu$m, log scale)"
        OUT = f"{OUT}_diam" if a.out == ap.get_default("out") else OUT
    use = list(range(NM)) if a.motifs == "all" else [int(x) for x in a.motifs.split(",")]

    f = h5py.File(H5, "r")
    o = f["obs"]
    sec_cats, sec = read_cat(o, "sample_id")
    ct_cats, ct = read_cat(o, "cell_type")
    xy = f["obsm"]["spatial"][:]
    S = np.vstack([(read_cat(o, f"motif_{k}_state")[1] == 1) for k in range(NM)])
    f.close()
    ct_cats = list(ct_cats)
    nonvasc = ~np.isin(ct, [ct_cats.index(t) for t in VASCULAR if t in ct_cats])
    log(f"{len(sec):,} cells; {nonvasc.sum():,} non-vascular")

    # ---- per-vessel ring fractions and lumen areas (figure 1) -----------------------
    rows, RF = [], []
    for s in SECTIONS:
        m = np.where(sec == sec_cats.index(s))[0]
        P = xy[m]
        tree = cKDTree(P)
        for p in load_polys(s):
            poly = p["poly"]
            bt = cKDTree(densify(poly))
            lo, hi = poly.min(0) - R_RING, poly.max(0) + R_RING
            cand = np.asarray(tree.query_ball_point((lo + hi) / 2,
                                                    np.hypot(*(hi - lo)) / 2 + 1.0))
            if cand.size == 0:
                continue
            dd, _ = bt.query(P[cand], k=1)
            ins = MplPath(poly).contains_points(P[cand])
            ring = cand[(dd > 0) & (dd <= R_RING) & (~ins)]
            if ring.size < MIN_RING or not p["area"]:
                continue
            rows.append(dict(section=s, area=float(p["area"]), n=int(ring.size),
                             d_perim=perimeter(poly) / np.pi))
            RF.append(S[:, m[ring]].mean(1))
    V = pd.DataFrame(rows)
    V["size"] = V.area if SZ == "area" else V.d_perim
    RF = np.vstack(RF)
    log(f"{len(V)} vessels with a ring and a lumen; size metric = "
        f"{'lumen area' if SZ == 'area' else 'perimeter diameter P/pi'}")

    if a.kde_combined:
        combined_size_kde(V, RF, [int(t) for t in a.kde_motifs.split(",")],
                          a.kde_thresh, OUT)
        log(f"wrote to {OUT}/size_kde/")
        return

    # ---- distance field, once per section (figure 2) --------------------------------
    D = np.full(len(sec), np.nan)
    for s in SECTIONS:
        m = np.where(sec == sec_cats.index(s))[0]
        D[m] = signed_distance(xy[m], load_polys(s))[0]
        log(f"  {s}: {(D[m] > 0).sum():,} cells outside a lumen, "
            f"{((D[m] > 0) & (D[m] <= DMAX)).sum():,} within {DMAX:.0f} um")

    if a.pooled:
        rows2 = [pooled_distance_figure(D, S, k, MCOL[i])
                 for i, k in enumerate(VASC_MOTIFS)]
        rows2 = [r for r in rows2 if r]
        pd.DataFrame(rows2).round(4).to_csv(
            f"{OUT}/distance_vs_positive_pooled.csv", index=False)
        print(pd.DataFrame(rows2).round(3).to_string(index=False))
        log(f"wrote to {OUT}/distance_vs_positive_pooled/")
        return

    size_stats, dist_stats, kde_stats = [], [], []
    for k in use:
        lab = f"m{k}" + (f" {ANN[k]}" if k in ANN else "")

        # ---------------------------------------------------------------- figure 1
        y = RF[:, k]
        x = np.log10(V['size'].to_numpy())
        rho, pv = spearmanr(x, y)
        per = {s: spearmanr(x[(V.section == s).to_numpy()],
                            y[(V.section == s).to_numpy()])[0] for s in SECTIONS}
        fig, ax = plt.subplots(figsize=(SINGLE, 62 * MM))
        for s in SECTIONS:
            m = (V.section == s).to_numpy()
            ax.scatter(x[m], y[m], s=9, marker=SECMK[s], color=SECCOL[s], lw=.3,
                       edgecolor="w", label=s, zorder=3)
        lo = sm.nonparametric.lowess(y, x, frac=.6, return_sorted=True)
        ax.plot(lo[:, 0], lo[:, 1], color="#1a1a1a", lw=1.4, zorder=4)
        tk = np.arange(np.floor(x.min()), np.ceil(x.max()) + 1)
        ax.set_xticks(tk)
        ax.set_xticklabels([("%g" % (10 ** t)) if t < 3 else f"$10^{{{int(t)}}}$" for t in tk])
        ax.set_ylim(-.03, 1.03)
        ax.set_xlabel(SIZE_LAB)
        ax.set_ylabel(f"fraction of ring cells positive")
        ax.set_title(lab, loc="left", pad=5)
        fig.legend(*ax.get_legend_handles_labels(), frameon=False, ncol=4,
                   handletextpad=.2, columnspacing=1.0, loc="upper center",
                   bbox_to_anchor=(.5, -(5 * MM) / fig.get_figheight()),
                   bbox_transform=fig.transFigure)
        save(fig, "size_vs_positive", f"m{k}")
        kde_stats.append(size_kde(x, y, k, lab))
        size_stats.append(dict(motif=f"m{k}", annotation=ANN.get(k, ""), rho=rho, p=pv,
                               **{f"rho_{s}": per[s] for s in SECTIONS}))

        # ---------------------------------------------------------------- figure 2
        fig, ax = plt.subplots(figsize=(SINGLE, 62 * MM))
        curves, curves_nv = [], []
        for s in SECTIONS:
            m = sec == sec_cats.index(s)
            r = spline_fit(D, S[k], m)
            if r is None:
                continue
            if r[2] is not None:
                ax.fill_between(r[0], r[2], r[3], color=SECCOL[s], alpha=.22, lw=0)
            ax.plot(r[0], r[1], color=SECCOL[s], lw=.9, alpha=.9, label=s)
            curves.append(r)
            rnv = spline_fit(D, S[k], m & nonvasc)
            if rnv is not None:
                curves_nv.append(rnv)
        if curves:
            g = curves[0][0]
            M = np.vstack([np.interp(g, c[0], c[1]) for c in curves])
            mean, sd = M.mean(0), M.std(0)
            # the band on the mean is the BETWEEN-SECTION spread, a different and larger
            # quantity than the per-section fit CI shaded above
            ax.fill_between(g, np.clip(mean - sd, 0, 1), np.clip(mean + sd, 0, 1),
                            color="#1a1a1a", alpha=.13, lw=0, zorder=4)
            ax.plot(g, mean, color="#1a1a1a", lw=1.8, zorder=5,
                    label=r"mean $\pm$ 1 SD")
            if curves_nv:
                # the non-vascular fit is still computed and still reported in
                # distance_vs_positive.csv -- it is only no longer drawn on the figure
                mnv = np.mean([np.interp(g, c[0], c[1]) for c in curves_nv], 0)
                dist_stats.append(dict(motif=f"m{k}", annotation=ANN.get(k, ""),
                                       p_at_10um=float(np.interp(10, g, mean)),
                                       p_at_150um=float(np.interp(150, g, mean)),
                                       delta=float(np.interp(10, g, mean)
                                                   - np.interp(150, g, mean)),
                                       p_at_10um_nonvasc=float(np.interp(10, g, mnv)),
                                       delta_nonvasc=float(np.interp(10, g, mnv)
                                                           - np.interp(150, g, mnv))))
        ax.set_xlim(0, XMAX_DIST); ax.set_ylim(-.03, 1.03)
        ax.set_xlabel(r"distance from the lumen boundary ($\mu$m)")
        ax.set_ylabel("P(motif positive)")
        ax.set_title(lab, loc="left", pad=5)
        ax.legend(frameon=False, loc="best", handlelength=1.6, ncol=2)
        save(fig, "distance_vs_positive", f"m{k}")
        print(f"  m{k:<3} size rho {rho:+.2f} | P(pos) 10um {dist_stats[-1]['p_at_10um']:.2f}"
              f" -> 150um {dist_stats[-1]['p_at_150um']:.2f}"
              f" (non-vasc delta {dist_stats[-1]['delta_nonvasc']:+.2f})"
              if dist_stats and dist_stats[-1]["motif"] == f"m{k}" else f"  m{k:<3} size rho {rho:+.2f}")

    os.makedirs(OUT, exist_ok=True)
    pd.DataFrame(size_stats).round(4).to_csv(f"{OUT}/size_vs_positive.csv", index=False)
    pd.DataFrame(dist_stats).round(4).to_csv(f"{OUT}/distance_vs_positive.csv", index=False)
    K = pd.DataFrame(kde_stats)
    K.round(4).to_csv(f"{OUT}/size_kde.csv", index=False)
    print(f"\n=== lumen {'area' if SZ == 'area' else 'diameter'} of positive vs negative "
          f"vessels (threshold {THRESH}) ===")
    print(K[["motif", "annotation", "n_positive", "n_negative", f"median_{SZ}_positive",
             f"median_{SZ}_negative", "mannwhitney_u", "auc",
             "p"]].round(3).to_string(index=False))
    log(f"wrote to {OUT}/")


if __name__ == "__main__":
    main()
