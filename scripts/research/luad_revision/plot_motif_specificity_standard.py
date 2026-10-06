#!/usr/bin/env python
"""Visualise the selection rule of find_motif_specific_vessels.py.

Rule (per target motif k, 300 vessels, 30 um ring, GMM state): eligible = f_k >= --min-frac;
among eligible vessels, order by the highest positive fraction of the other motifs
(m23 excluded unless k = 23) ascending, ties -> higher f_k, more ring cells, vessel id;
the first --top are the picks. The picks come from find_motif_specific_vessels.select()
and are asserted identical to example_vessels.csv.

Outputs (figures_motif_specific_vessels/):
    m<k>_standard           all vessels: x = f_k, y = max other; eligible band shaded,
                            picks numbered
    standard_max_other_counts   eligible vessels with max other <= t, per motif
    standard_picks          picked vessels x six motifs (positive fraction), target outlined

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/plot_motif_specificity_standard.py
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
from find_motif_specific_vessels import SRC, OUT, MM, USE, ANN, UBIQ, select  # noqa: E402

SECCOL = {"P17_AIS": "#0072b2", "P17_LUAD": "#e69f00",
          "P21_AIS": "#009e73", "P21_LUAD": "#cc79a7"}
SECMK = {"P17_AIS": "o", "P17_LUAD": "s", "P21_AIS": "^", "P21_LUAD": "D"}
MCOL = ["#0072b2", "#e69f00", "#009e73", "#cc79a7", "#56b4e9", "#d55e00"]  # validated PASS
PICKCOL = "#b2182b"
OFFS = [(5, 5, "left", "bottom"), (-5, 5, "right", "bottom"),
        (0, -6, "center", "top")]   # pick-number label positions, 1 / 2 / 3

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def save(fig, out, stem):
    os.makedirs(out, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{out}/{stem}.{ext}", dpi=400 if ext == "png" else None)
    plt.close(fig)
    print(f"  {stem}.{{png,pdf,svg}}")


def mlab(k):
    return f"m{k} {ANN[k]}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--min-frac", type=float, default=0.5)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    V = pd.read_csv(f"{SRC}/vessels.csv")
    F = pd.DataFrame(np.load(f"{SRC}/features.npy"), columns=[f"m{k}" for k in USE])
    sec = V.section.to_numpy()

    sel = {k: select(V, F, k, a.min_frac, a.top) for k in USE}
    ref = pd.read_csv(f"{a.out}/example_vessels.csv")
    for k in USE:
        want = ref.loc[ref.motif == f"m{k}", "vessel_id"].tolist()
        got = V.vessel_id.iloc[sel[k][0]].tolist()
        assert got == want, (k, got, want)
    print("picks identical to example_vessels.csv")

    for k in USE:
        picks, mo, elig = sel[k]
        x = F[f"m{k}"].to_numpy()
        fig, ax = plt.subplots(figsize=(62 * MM, 58 * MM))
        ax.axvspan(a.min_frac, 1.03, color="#bdbdbd", alpha=.35, lw=0, zorder=0)
        ax.plot([0, 1], [0, 1], color="#bbbbbb", lw=.5, zorder=1)
        for s in SECCOL:
            m = sec == s
            ax.scatter(x[m], mo[m], s=10, marker=SECMK[s], color=SECCOL[s], lw=.25,
                       edgecolor="white", alpha=.85, zorder=2)
        for n, v in enumerate(picks, 1):
            ax.scatter(x[v], mo[v], s=48, facecolor="none", edgecolor=PICKCOL, lw=.9,
                       zorder=3)
            dx, dy, ha, va = OFFS[(n - 1) % len(OFFS)]
            ax.annotate(str(n), (x[v], mo[v]), xytext=(dx, dy), textcoords="offset points",
                        ha=ha, va=va, fontsize=6.5, color=PICKCOL, zorder=4)
        ax.set(xlim=(-.02, 1.03), ylim=(-.02, 1.03), xticks=[0, .5, 1], yticks=[0, .5, 1],
               xlabel=f"m{k} positive fraction",
               ylabel="max positive fraction, other motifs"
                      + ("" if k == UBIQ else " (excl. m23)"))
        ax.set_aspect("equal")
        ax.set_title(mlab(k), loc="left", pad=5)
        ax.set_title(f"eligible n = {int(elig.sum())}", loc="right", fontsize=6.5,
                     color="#444444", pad=5)
        h = [plt.Line2D([], [], ls="", marker=SECMK[s], mfc=SECCOL[s], mec="white", mew=.3,
                        ms=4, label=s) for s in SECCOL]
        h.append(plt.Line2D([], [], ls="", marker="o", mfc="none", mec=PICKCOL, mew=.9,
                            ms=6, label="picked"))
        ax.legend(handles=h, loc="upper center", bbox_to_anchor=(.5, -.2), ncol=5,
                  frameon=False, handletextpad=.2, columnspacing=.8)
        save(fig, a.out, f"m{k}_standard")

    ts = np.linspace(0, 1, 501)
    fig, ax = plt.subplots(figsize=(70 * MM, 52 * MM))
    for i, k in enumerate(USE):
        _, mo, elig = sel[k]
        ax.plot(ts, [(elig & (mo <= t)).sum() for t in ts], color=MCOL[i], lw=1,
                drawstyle="steps-post", label=mlab(k))
    ax.set(xlim=(0, 1), yscale="log", ylim=(.8, 320),
           xlabel="max positive fraction, other motifs",
           ylabel=f"eligible vessels (target ≥ {a.min_frac:g})")
    ax.set_yticks([1, 10, 100, 300], ["1", "10", "100", "300"])
    ax.legend(loc="center left", bbox_to_anchor=(1.01, .5), frameon=False, handlelength=1.2)
    save(fig, a.out, "standard_max_other_counts")

    order = [(k, v) for k in USE for v in sel[k][0]]
    ylab = [f"m{k} #{n}  {V.vessel_id.iloc[v]}" for k in USE
            for n, v in enumerate(sel[k][0], 1)]
    H = np.array([F.iloc[v].to_numpy() for _, v in order])
    fig, ax = plt.subplots(figsize=(52 * MM, (len(order) * 3.4 + 14) * MM))
    im = ax.imshow(H, cmap="Reds", vmin=0, vmax=1, aspect="auto")
    for r, (k, _) in enumerate(order):
        ax.add_patch(Rectangle((USE.index(k) - .5, r - .5), 1, 1, fill=False, ec="black",
                               lw=1))
    starts = np.cumsum([len(sel[k][0]) for k in USE])[:-1]
    for b in starts:
        ax.axhline(b - .5, color="white", lw=1.5)
    ax.set_xticks(range(len(USE)), [mlab(k) for k in USE], rotation=45, ha="right",
                  rotation_mode="anchor")
    ax.set_yticks(range(len(order)), ylab)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=.08, pad=.03, shrink=.35, anchor=(0, 1))
    cb.set_label("positive fraction")
    cb.outline.set_linewidth(.4)
    cb.ax.tick_params(labelsize=6)
    save(fig, a.out, "standard_picks")


if __name__ == "__main__":
    main()
