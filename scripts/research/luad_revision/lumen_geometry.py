#!/usr/bin/env python
"""Vessel size metrics that survive sectioning: perimeter-derived diameter, not area.

A lumen is flattened when the block is cut and mounted. What that does to each metric:

  area A          falls -- a collapsed vein loses lumen area, an artery much less
  perimeter P     conserved: flattening an ellipse does not change the wall's circumference
  D_area  = 2 sqrt(A / pi)       the diameter of a circle of the same area -- shrinks with A
  D_perim = P / pi               the diameter of a circle of the same circumference -- stable
  circularity = 4 pi A / P^2     1.0 = circle, -> 0 as the lumen flattens; the squish measure

So D_perim is the size metric to prefer, and D_area / D_perim = sqrt(circularity) is exactly how
much that vessel was squashed. NOTE that D_area is a monotone function of A, so on a log axis it
is the area plot with a rescaled axis and *identical* Spearman statistics -- only D_perim can
change a ranking, and only the comparison against it tells you whether the choice matters.

Perimeter is computed on the curated polygon as traced. Lumens whose polygon failed to parse
fall back to the vessel's bounding box (`src = 1`): a rectangle, whose perimeter is not the
vessel's, so they are written out but excluded from the summary.

Outputs: lumen_geometry.csv (one row per curated lumen) and, in figures_vessel_size/,
    lumen_size_metrics -- D_perim vs D_area, the circularity distribution, and whether
    squishing tracks the mural-cell fraction (arteries squash less than veins).

Usage:
    /home/fanj2/.conda/envs/spatial/bin/python \
        /home/fanj2/alarmist/scripts/research/luad_revision/lumen_geometry.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import SECTIONS, load_polys, log

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
OUT = f"{BASE}/figures_vessel_size"
FRAC = f"{BASE}/vessel_kmeans_ring30_m0-1-6-10-23-24_frac"
MM = 1 / 25.4
SECCOL = {"P17_AIS": "#0072b2", "P17_LUAD": "#e69f00",
          "P21_AIS": "#009e73", "P21_LUAD": "#cc79a7"}
SECMRK = {"P17_AIS": "o", "P17_LUAD": "s", "P21_AIS": "^", "P21_LUAD": "D"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
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


def shoelace(P):
    x, y = P[:, 0], P[:, 1]
    return .5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def perimeter(P):
    Q = np.vstack([P, P[:1]])
    return float(np.hypot(*np.diff(Q, axis=0).T).sum())


def feret(P, n=180):
    """Max and min caliper width over `n` orientations (the min is the squashed axis)."""
    th = np.linspace(0, np.pi, n, endpoint=False)
    proj = P @ np.vstack([np.cos(th), np.sin(th)])        # (vertices, angles)
    w = proj.max(0) - proj.min(0)
    return float(w.max()), float(w.min())


def geometry():
    rows = []
    for sec in SECTIONS:
        for p in load_polys(sec):
            P, A, Per = p["poly"], shoelace(p["poly"]), perimeter(p["poly"])
            fmax, fmin = feret(P)
            rows.append(dict(
                vessel_id=f"{sec}__{p['vid']}", section=sec, src=p["src"],
                flagged=p["flagged"], area_curated=p["area"], area=A, perimeter=Per,
                d_area=2 * np.sqrt(A / np.pi), d_perim=Per / np.pi,
                feret_max=fmax, feret_min=fmin,
                circularity=4 * np.pi * A / Per ** 2, aspect=fmax / max(fmin, 1e-9)))
    return pd.DataFrame(rows)


def main():
    G = geometry()
    G.to_csv(f"{BASE}/lumen_geometry.csv", index=False)
    log(f"{len(G)} lumens; {int((G.src == 1).sum())} fell back to the bounding box "
        f"(excluded below)")
    g = G[G.src == 0].copy()
    ok = g.area_curated.notna()
    log(f"shoelace vs curated area: max relative difference "
        f"{np.nanmax(np.abs(g.area[ok] - g.area_curated[ok]) / g.area_curated[ok]):.2e}")
    log(f"circularity: median {g.circularity.median():.2f}, "
        f"10th pct {g.circularity.quantile(.1):.2f}, "
        f"{(g.circularity < .5).mean():.0%} below 0.5")
    rs = spearmanr(g.area, g.d_perim).statistic
    rp = spearmanr(g.d_area, g.d_perim).statistic
    log(f"Spearman area vs d_perim {rs:.3f} (d_area vs d_perim is the same, {rp:.3f}: "
        f"d_area is a monotone function of area)")
    rk = pd.DataFrame({"a": g.area.rank(), "d": g.d_perim.rank()})
    log(f"rank shift between the two: median {np.median(np.abs(rk.a - rk.d)):.0f} places, "
        f"max {np.max(np.abs(rk.a - rk.d)):.0f} of {len(g)}")

    V = pd.read_csv(f"{FRAC}/vessels.csv").set_index("vessel_id")
    g = g.join(V[["frac_SMC", "n_ring"]], on="vessel_id")

    fig, axes = plt.subplots(1, 4, figsize=(198 * MM, 46 * MM))
    ax = axes[0]
    for s in SECTIONS:
        m = g.section == s
        ax.scatter(g.d_area[m], g.d_perim[m], s=7, marker=SECMRK[s], facecolor="none",
                   edgecolor=SECCOL[s], lw=.6, label=s)
    lim = [g[["d_area", "d_perim"]].to_numpy().min() * .9,
           g[["d_area", "d_perim"]].to_numpy().max() * 1.1]
    ax.plot(lim, lim, color="#999999", lw=.6, ls="--", zorder=0)
    ax.set(xscale="log", yscale="log", xlim=lim, ylim=lim,
           xlabel=r"area-equivalent diameter 2$\sqrt{A/\pi}$ ($\mu$m)",
           ylabel=r"perimeter diameter P/$\pi$ ($\mu$m)")
    ax.legend(frameon=False, loc="upper left", handletextpad=.2, borderpad=0, labelspacing=.2)
    ax.set_title(f"Spearman {rp:.3f}", loc="left", pad=3, fontsize=6.5)

    ax = axes[1]
    ax.hist(g.circularity, bins=np.linspace(0, 1, 31), color="#4e79a7", lw=0)
    ax.axvline(g.circularity.median(), color="#d55e00", lw=.9)
    ax.text(g.circularity.median() - .02, ax.get_ylim()[1] * .95,
            f"median {g.circularity.median():.2f}", ha="right", va="top", fontsize=6,
            color="#d55e00")
    ax.set(xlabel=r"circularity 4$\pi$A/P$^2$  (1 = circle)", ylabel="lumens", xlim=(0, 1))

    ax = axes[2]
    for s in SECTIONS:
        m = g.section == s
        ax.scatter(g.d_perim[m], g.circularity[m], s=7, marker=SECMRK[s], facecolor="none",
                   edgecolor=SECCOL[s], lw=.6)
    r = spearmanr(g.d_perim, g.circularity)
    ax.set(xscale="log", xlabel=r"lumen diameter P/$\pi$ ($\mu$m)",
           ylabel=r"circularity 4$\pi$A/P$^2$", ylim=(0, 1))
    ax.set_title(f"Spearman {r.statistic:.2f}", loc="left", pad=3, fontsize=6.5)

    ax = axes[3]
    for s in SECTIONS:
        m = g.section == s
        ax.scatter(g.frac_SMC[m], g.circularity[m], s=7, marker=SECMRK[s], facecolor="none",
                   edgecolor=SECCOL[s], lw=.6)
    r = spearmanr(g.frac_SMC, g.circularity, nan_policy="omit")
    g["band"] = pd.qcut(np.log10(g.area), 4, labels=False)
    inside = [spearmanr(d.frac_SMC, d.circularity, nan_policy="omit").statistic
              for _, d in g.groupby("band")]
    log("circularity vs SMC fraction: overall rho %+.2f (p=%.1g); within area quartiles %s"
        % (r.statistic, r.pvalue, ", ".join(f"{v:+.2f}" for v in inside)))
    ax.set(xlabel="SMC fraction of the 30 $\\mu$m ring (arterial proxy)",
           ylabel=r"circularity 4$\pi$A/P$^2$", ylim=(0, 1))
    ax.set_title(f"Spearman {r.statistic:+.2f} overall,\n"
                 f"{', '.join(f'{v:+.2f}' for v in inside)} within size quartiles",
                 loc="left", pad=3, fontsize=6.5)
    fig.subplots_adjust(wspace=.48)
    save(fig, "lumen_size_metrics")
    print(G.describe().round(2).to_string())


if __name__ == "__main__":
    main()
