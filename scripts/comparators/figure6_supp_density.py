#!/usr/bin/env python
"""Supplementary figure — grade and cellularity are collinear in the 13-core GBM TMA.

The grade panel of Figure 6 rests on a 7-vs-6 punch comparison. High-grade punches
hold 3.4x more cells than low-grade ones, and an ALARMIST loading is a projection
over a 50 um neighbourhood, so loading tracks cellularity almost perfectly. This
figure states that limit rather than hiding it.

It does NOT show that the grade result is an artifact: rho(grade, cellularity) = 0.78,
so at n = 13 the two cannot be separated. Residualising on density removes the grade
signal for every motif, which is what collinearity does regardless of the truth.

Reads panel_f_*.csv from figure6_gather.py. Computes nothing new.

Env: /Users/jiayifan/anaconda3/envs/bptf/bin/python
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from figure6_plot import ACCENT, BLUE, GREY, INK, LIGHT, U, panel_letter, save_all_formats, set_style

ROOT = Path("/Users/jiayifan/tansey_lab/alarmist")
DATA = ROOT / "results/comparators/_benchmark/figure6"
FIG = ROOT / "results/comparators/_benchmark/figure6"

HI, LO = ACCENT, BLUE
MOTIF = 1  # the mGAM <-> MES-like loop


def main() -> None:
    set_style()
    cores = pd.read_csv(DATA / "panel_f_cores.csv")
    dens = pd.read_csv(DATA / "panel_f_density.csv")
    meta = json.loads((DATA / "panel_f_meta.json").read_text())

    hi = (cores["grade"] == "high").to_numpy()
    n_cells = cores["n_cells"].to_numpy(float)
    load = cores[f"mean_loading_m{MOTIF}"].to_numpy(float)

    fig = plt.figure(figsize=(5 * U, 1.65 * U))
    gs = fig.add_gridspec(1, 4, width_ratios=[0.85, 1.15, 1.15, 1.0], wspace=0.62)

    # ---- a: the design imbalance
    ax = fig.add_subplot(gs[0, 0])
    for j, (mask, lab, col) in enumerate(((hi, "high", HI), (~hi, "low", LO))):
        xj = np.full(mask.sum(), j) + np.linspace(-0.13, 0.13, mask.sum())
        ax.scatter(xj, n_cells[mask], s=16, facecolor=col, edgecolor="none", zorder=3)
        ax.plot([j - 0.26, j + 0.26], [n_cells[mask].mean()] * 2, color=INK, lw=1.2, zorder=4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"high\n(n = {meta['n_high']})", f"low\n(n = {meta['n_low']})"])
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("Cells per punch (n)")
    ax.set_yscale("log")
    ax.set_ylim(500, 40000)
    ax.set_yticks([1000, 10000])
    ax.set_yticklabels(["1,000", "10,000"])
    ax.text(0.5, 0.98, f"{meta['mean_cells_high'] / meta['mean_cells_low']:.1f}$\\times$",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.5)
    ax.set_xlabel("Tumour grade")
    panel_letter(ax, "a", dx=-0.45)

    # ---- b: loading tracks cellularity
    axb = fig.add_subplot(gs[0, 1])
    for mask, col in ((hi, HI), (~hi, LO)):
        axb.scatter(n_cells[mask], load[mask], s=18, facecolor=col, edgecolor="none", zorder=3)
    b = np.polyfit(np.log10(n_cells), load, 1)
    xs = np.linspace(np.log10(n_cells.min()), np.log10(n_cells.max()), 50)
    axb.plot(10 ** xs, np.polyval(b, xs), color=GREY, lw=0.9, ls=(0, (3, 2)), zorder=2)
    r = dens.loc[dens["motif"] == MOTIF, "rho_density"].iloc[0]
    p = dens.loc[dens["motif"] == MOTIF, "p_density"].iloc[0]
    axb.text(0.03, 0.97, f"Spearman $\\rho$ = {r:.2f}\n$p$ = {p:.1e}", transform=axb.transAxes,
             ha="left", va="top", fontsize=5)
    axb.set_xscale("log")
    axb.set_xlim(500, 40000)
    axb.set_xticks([1000, 10000])
    axb.set_xticklabels(["1,000", "10,000"])
    axb.set_xlabel("Cells per punch (n)")
    axb.set_ylabel(f"Mean motif-{MOTIF} loading")
    panel_letter(axb, "b", dx=-0.35)

    # ---- c: every motif does this, not just motif 1
    axc = fig.add_subplot(gs[0, 2])
    order = np.argsort(dens["rho_density"].to_numpy())
    yy = np.arange(len(order))
    cols = [ACCENT if dens["motif"].iloc[i] == MOTIF else GREY for i in order]
    axc.barh(yy, dens["rho_density"].to_numpy()[order], color=cols, height=0.7, lw=0)
    lab_i = int(np.where(dens["motif"].to_numpy()[order] == MOTIF)[0][0])
    axc.text(dens["rho_density"].to_numpy()[order][lab_i] + 0.03, lab_i, f"motif {MOTIF}",
             va="center", ha="left", fontsize=5, color=ACCENT)
    axc.set_yticks([])
    axc.set_xlim(0, 1.28)
    axc.set_xticks([0, 0.5, 1.0])
    axc.set_xlabel("Spearman $\\rho$ (loading vs cellularity)")
    axc.set_ylabel(f"ALARMIST motifs (n = {len(dens)})", labelpad=2)
    axc.spines["left"].set_visible(False)
    panel_letter(axc, "c", dx=-0.22)

    # ---- d: grade p before and after removing the cellularity trend
    axd = fig.add_subplot(gs[0, 3])
    raw = -np.log10(dens["p_grade_raw"].to_numpy())
    adj = -np.log10(dens["p_grade_density_adjusted"].to_numpy())
    for i in range(len(raw)):
        col = ACCENT if dens["motif"].iloc[i] == MOTIF else GREY
        axd.plot([0, 1], [raw[i], adj[i]], color=col, lw=1.1 if col == ACCENT else 0.6,
                 alpha=1.0 if col == ACCENT else 0.45,
                 zorder=3 if col == ACCENT else 1)
        axd.scatter([0, 1], [raw[i], adj[i]], s=9 if col == ACCENT else 5, color=col,
                    zorder=4 if col == ACCENT else 2,
                    alpha=1.0 if col == ACCENT else 0.45)
    axd.axhline(-np.log10(0.05), color=INK, lw=0.8, ls=(0, (3, 2)), zorder=0)
    axd.text(1.48, -np.log10(0.05), "$p$ = 0.05", fontsize=5, va="center", ha="left")
    axd.set_xticks([0, 1])
    axd.set_xticklabels(["as\nreported", "cellularity\nremoved"])
    axd.set_xlim(-0.35, 1.42)
    axd.set_ylabel("Grade association, $-$log$_{10}$ $p$")
    axd.spines["bottom"].set_bounds(0, 1)
    n_raw = int((dens["p_grade_raw"] < 0.05).sum())
    n_adj = int((dens["p_grade_density_adjusted"] < 0.05).sum())
    axd.text(0.0, axd.get_ylim()[1], f"{n_raw}/20", ha="center", va="bottom", fontsize=5)
    axd.text(1.0, axd.get_ylim()[1], f"{n_adj}/20", ha="center", va="bottom", fontsize=5)
    panel_letter(axd, "d", dx=-0.42)

    handles = [Line2D([], [], marker="o", ls="none", color=HI, ms=4, label="high grade"),
               Line2D([], [], marker="o", ls="none", color=LO, ms=4, label="low grade")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.055, 1.0), ncol=2,
               handletextpad=0.3, columnspacing=1.4)

    save_all_formats(fig, FIG / "figure6_supp_density_collinearity")
    plt.close(fig)


if __name__ == "__main__":
    main()
