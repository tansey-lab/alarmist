#!/usr/bin/env python
"""Figure 7 — one spatial panel per TMA punch, plus a cross-punch summary.

Per punch: columns are methods, rows are the two arms of ALARMIST motif 1, and
ALARMIST occupies a single cell spanning both rows because for it the two arms are
not two objects.

rho under each column is the Spearman correlation between that method's OWN two maps
in that method's OWN units on that punch -- a within-method statistic, so it carries
none of the cross-method scale problems the maps themselves have.

Colour is the percentile within each map. Values are NOT comparable across methods
(LRscore / co-expression / Moran z / transport mass / mechanism score / cosine).

Writes  results/comparators/_benchmark/figure7/punches/figure7_core<N>.{png,pdf,svg}
        results/comparators/_benchmark/figure7/figure7_rho_summary.{png,pdf,svg}

Env: /Users/jiayifan/anaconda3/envs/bptf/bin/python
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, spearmanr

from figure6_plot import ACCENT, BLUE, GREY, INK, U, panel_letter, save_all_formats, set_style

ROOT = Path("/Users/jiayifan/tansey_lab/alarmist")
DATA = ROOT / "results/comparators/_benchmark/figure7"
FIG = ROOT / "results/comparators/_benchmark/figure7"
PUNCH = FIG / "punches"
PUNCH.mkdir(parents=True, exist_ok=True)

COLS = ["COMMOT", "NICHES", "stLearn", "CytoSignal", "CellChat", "LIANA+", "SpatialDM"]
ARMS = [("GRN_SORT1", "GRN→SORT1"), ("ANXA1_FPR1", "ANXA1→FPR1")]
CMAP = "magma"
MS = {"stLearn": 3.2}       # stLearn is a 51.3 um grid; everything else is single cells


def pct(v: np.ndarray) -> np.ndarray:
    """Rank-scale to [0, 1].

    method='min', not the default 'average': COMMOT leaves most cells at exactly zero
    and averaging their ranks would paint that block mid-scale instead of at the floor.
    """
    return (rankdata(v, method="min") - 1) / max(len(v) - 1, 1)


def blank(ax, xlim, ylim, msg: str, title: str | None) -> None:
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.add_patch(plt.Rectangle((0.06, 0.06), 0.88, 0.88, transform=ax.transAxes,
                               facecolor="#F4F4F4", edgecolor=GREY, lw=0.6, ls=(0, (2, 2))))
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center",
            fontsize=4.6, color=GREY, style="italic")
    if title:
        ax.set_title(title, fontsize=5.5, pad=3)


def draw(ax, x, y, v, s, xlim, ylim) -> None:
    ok = ~np.isnan(v)          # SpatialDM leaves exactly one NaN cell in 5 core/arm combos
    x, y, v = x[ok], y[ok], v[ok]
    o = np.argsort(v)
    ax.scatter(x[o], y[o], c=pct(v)[o], s=s, cmap=CMAP, vmin=0, vmax=1,
               linewidths=0, rasterized=True)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def one_punch(z, core: str, grade: str, n_cells: int) -> dict[str, float]:
    ax_ = z[f"{core}|ALARMIST|motif1|x"]; ay_ = z[f"{core}|ALARMIST|motif1|y"]
    pad = 0.02 * max(np.ptp(ax_), np.ptp(ay_))
    cx, cy = (ax_.min() + ax_.max()) / 2, (ay_.min() + ay_.max()) / 2
    half = max(np.ptp(ax_), np.ptp(ay_)) / 2 + pad        # square box, so every punch is
    xlim = (cx - half, cx + half)                      # drawn at the same scale
    ylim = (cy - half, cy + half)

    fig = plt.figure(figsize=(5 * U, 2.35 * U))
    gs = fig.add_gridspec(2, len(COLS) + 1, wspace=0.10, hspace=0.10,
                          width_ratios=[1] * len(COLS) + [1.30])

    rho: dict[str, float] = {}
    flat: dict[str, bool] = {}
    for j, m in enumerate(COLS):
        if m == "CellChat":
            for i in range(2):
                blank(fig.add_subplot(gs[i, j]), xlim, ylim,
                      "no spatial\noutput\nexists", m if i == 0 else None)
            continue

        have = [f"{core}|{m}|{a}|v" in z for a, _ in ARMS]
        if all(have):
            v1 = z[f"{core}|{m}|{ARMS[0][0]}|v"]; v2 = z[f"{core}|{m}|{ARMS[1][0]}|v"]
            ok = ~(np.isnan(v1) | np.isnan(v2))
            v1, v2 = v1[ok], v2[ok]
            if len(v1) > 2 and v1.std() > 0 and v2.std() > 0:
                rho[m] = float(spearmanr(v1, v2).statistic)
            else:
                flat[m] = True   # a map with no variance -- almost always all-zero

        for i, (key, _) in enumerate(ARMS):
            axm = fig.add_subplot(gs[i, j])
            k = f"{core}|{m}|{key}|v"
            if k in z:
                draw(axm, z[f"{core}|{m}|{key}|x"], z[f"{core}|{m}|{key}|y"], z[k],
                     MS.get(m, 0.55), xlim, ylim)
                if i == 0:
                    axm.set_title(m, fontsize=5.5, pad=3)
            else:
                blank(axm, xlim, ylim, "not globally\nselected in\nthis punch",
                      m if i == 0 else None)
            if i == 1:
                lab = (f"$\\rho$ = {rho[m]:+.2f}" if m in rho
                       else ("all zero" if flat.get(m) else "not selected"))
                axm.text(0.5, -0.06, lab, transform=axm.transAxes, ha="center", va="top",
                         fontsize=5, color=INK if m in rho else GREY)

    axa = fig.add_subplot(gs[:, len(COLS)])
    draw(axa, ax_, ay_, z[f"{core}|ALARMIST|motif1|v"], 1.4, xlim, ylim)
    axa.set_title("ALARMIST", fontsize=5.5, pad=3, color=ACCENT, fontweight="bold")
    axa.text(0.5, -0.045, "motif 1 — a single object\nspanning 712 LRIs",
             transform=axa.transAxes, ha="center", va="top", fontsize=5, color=ACCENT)
    for sp in axa.spines.values():
        sp.set_visible(True); sp.set_color(ACCENT); sp.set_linewidth(0.9)

    bx = xlim[0] + 0.04 * (xlim[1] - xlim[0])
    by = ylim[0] + 0.055 * (ylim[1] - ylim[0])
    axa.plot([bx, bx + 500], [by, by], color=INK, lw=1.3, solid_capstyle="butt", zorder=6)
    axa.text(bx + 250, by + 0.015 * (ylim[1] - ylim[0]), "500 µm", ha="center", va="bottom",
             fontsize=4.8, zorder=6)

    for i, (_, lab) in enumerate(ARMS):
        fig.text(0.055, 0.72 - 0.42 * i, lab, rotation=90, ha="right", va="center", fontsize=5.5)
    fig.text(0.055, 0.98, f"punch {core}  ·  {grade} grade  ·  {n_cells:,} cells",
             ha="left", va="top", fontsize=6, fontweight="bold")

    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=fig.axes, fraction=0.008, pad=0.008, aspect=34, shrink=0.72)
    cb.set_ticks([0, 1]); cb.set_ticklabels(["low", "high"])
    cb.set_label("percentile within each map", fontsize=5, labelpad=1)
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, labelsize=5)

    save_all_formats(fig, PUNCH / f"figure7_core{core}")
    plt.close(fig)
    return {"rho": rho, "flat": flat}


def summary(rho_all: dict[str, dict[str, float]], cores, grade) -> None:
    """Heatmap of the within-method rho across every punch."""
    methods = [m for m in COLS if m != "CellChat"]
    M = np.full((len(methods), len(cores)), np.nan)
    Z = np.zeros((len(methods), len(cores)), bool)
    for j, c in enumerate(cores):
        for i, m in enumerate(methods):
            if m in rho_all[c]["rho"]:
                M[i, j] = rho_all[c]["rho"][m]
            elif rho_all[c]["flat"].get(m):
                Z[i, j] = True

    fig, ax = plt.subplots(figsize=(3.4 * U, 1.5 * U))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(len(cores)))
    ax.set_xticklabels([f"{c}" for c in cores])
    ax.set_yticks(range(len(methods)))
    ax.set_yticklabels(methods)
    ax.set_xlabel("TMA punch")
    for j, c in enumerate(cores):
        ax.text(j, -0.72, grade[c][:2], ha="center", va="center", fontsize=4.6,
                color=ACCENT if grade[c] == "high" else BLUE)
    for i in range(len(methods)):
        for j in range(len(cores)):
            if np.isnan(M[i, j]):
                ax.text(j, i, "0" if Z[i, j] else "ns", ha="center", va="center",
                        fontsize=4.3, color=GREY, style="italic")
            else:
                ax.text(j, i, f"{M[i, j]:.2f}".lstrip("0").replace("-0", "-"),
                        ha="center", va="center", fontsize=4.3,
                        color="white" if abs(M[i, j]) > 0.38 else INK)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.015)
    cb.set_label("Spearman $\\rho$ between the two arms,\nin each method's own units",
                 fontsize=5, labelpad=2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5, length=1.5)
    ax.text(0.0, 1.20, "ns = pair not globally selected in that punch    0 = map has no "
            "variance (all zero)", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=4.6, color=GREY)
    save_all_formats(fig, FIG / "figure7_rho_summary")
    plt.close(fig)

    med = {m: float(np.nanmedian(M[i])) for i, m in enumerate(methods)}
    (DATA / "within_method_rho_all_cores.json").write_text(json.dumps(
        {"per_core": rho_all, "median_across_cores": med}, indent=2))
    print("\nmedian rho across punches:", {k: round(v, 3) for k, v in med.items()})


def main() -> None:
    set_style()
    z = np.load(DATA / "spatial_maps_all_cores.npz")
    prov = json.loads((DATA / "provenance_all_cores.json").read_text())
    cores, grade, ncell = prov["cores"], prov["grade"], prov["n_cells"]

    rho_all = {c: one_punch(z, c, grade[c], ncell[c]) for c in cores}
    summary(rho_all, cores, grade)


if __name__ == "__main__":
    main()
