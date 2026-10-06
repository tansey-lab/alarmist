#!/usr/bin/env python
"""UpSet plot of the four vessel motifs: which cells / vessels carry which combination.

Two units, one figure each:
  cells    every cell of the four sections; a cell belongs to a motif's set when its
           obs motif_<k>_state is "positive" (the GMM ON/OFF call)
  vessels  the 300 curated vessels; a vessel belongs to a motif's set when at least
           --vessel-thresh of its 30 um ring cells are positive (0.5, the cut used by the
           size_kde figures)

Layout is the standard UpSet: intersection sizes on top, set sizes on the left, membership
matrix below. Intersections are exclusive -- a cell counts in exactly one column, the column
of the full combination it carries -- so the columns sum to the number of cells in at least
one set. The empty combination is reported in the right-hand title, not drawn.

Outputs (figures_motif_upset/): motif_upset_cells, motif_upset_vessels,
    motif_upset_counts.csv (every one of the 16 combinations, both units, plus the 0.9
    vessel threshold)

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/plot_motif_upset.py
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
import textwrap
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import H5, SECTIONS, read_cat, log
from motif_state_palette import (ANCHOR as COL, MOTIFS, NAME, SHORT,  # one shared map
                                 state_code, state_label, state_palette)

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
OUT = f"{BASE}/figures_motif_upset"
SRC = f"{BASE}/vessel_kmeans_ring30_m0-1-6-10-23-24_frac"
RING_MOTIFS = [0, 1, 6, 10, 23, 24]        # column order of features.npy
MM = 1 / 25.4
BARCOL, OFFCOL = "#4d4d4d", "#d9d9d9"
PAL = state_palette()          # the 16 OKLab state colours the section maps use

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
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


def combo_counts(M):
    """M: n x len(MOTIFS) boolean. Returns counts of all 2**k exclusive combinations."""
    return np.bincount(state_code(M), minlength=1 << M.shape[1])


def upset(counts, total, unit, out, stem, note=""):
    """standard UpSet: intersection bars on top, set sizes left, membership matrix below"""
    k = len(MOTIFS)
    combos = [c for c in range(1, 1 << k) if counts[c] > 0]
    combos.sort(key=lambda c: -counts[c])
    set_size = {j: sum(counts[c] for c in range(1 << k) if c >> j & 1) for j in range(k)}
    order = sorted(range(k), key=lambda j: -set_size[j])          # biggest set on top
    fmt = (lambda v: f"{v/1000:,.0f}k") if total > 20_000 else (lambda v: f"{v:,.0f}")

    fig = plt.figure(figsize=(max(118, 7 * len(combos) + 62) * MM, 118 * MM))
    gs = GridSpec(2, 2, width_ratios=[.75, 2.3], height_ratios=[2.1, 1], hspace=.06,
                  wspace=.36, figure=fig)      # the gap holds the motif names
    ax_b = fig.add_subplot(gs[0, 1])
    ax_m = fig.add_subplot(gs[1, 1], sharex=ax_b)
    ax_s = fig.add_subplot(gs[1, 0], sharey=ax_m)
    fig.add_subplot(gs[0, 0]).axis("off")

    x = np.arange(len(combos))
    ax_b.bar(x, [counts[c] for c in combos], width=.68,
             color=[PAL[c] for c in combos], edgecolor="#4d4d4d", lw=.3)
    ax_b.set_ylabel(f"{unit} in combination")
    ax_b.set_yticks(ax_b.get_yticks(), [fmt(v) for v in ax_b.get_yticks()])
    ax_b.set_ylim(0, max(counts[c] for c in combos) * 1.08)
    ax_b.tick_params(labelbottom=False, bottom=False)
    ax_b.spines["bottom"].set_visible(False)

    for r, j in enumerate(order):
        y = k - 1 - r
        ax_m.axhspan(y - .5, y + .5, color="#f2f2f2" if r % 2 else "white", lw=0, zorder=0)
        for i, c in enumerate(combos):
            on = bool(c >> j & 1)
            ax_m.plot(i, y, "o", ms=5.2, color=COL[MOTIFS[j]] if on else OFFCOL, zorder=3)
    for i, c in enumerate(combos):
        ys = [k - 1 - r for r, j in enumerate(order) if c >> j & 1]
        if len(ys) > 1:
            ax_m.plot([i, i], [min(ys), max(ys)], color=BARCOL, lw=1.4, zorder=2)
    ax_m.set_xlim(-.6, len(combos) - .4)
    ax_m.set_ylim(-.5, k - .5)
    ax_m.set_yticks(range(k), [textwrap.fill(NAME[MOTIFS[order[k - 1 - i]]], 18)
                               for i in range(k)])
    ax_m.tick_params(axis="y", pad=2)
    ax_m.set_xticks([])
    ax_m.tick_params(length=0)
    for sp in ax_m.spines.values():
        sp.set_visible(False)

    sizes = [set_size[order[k - 1 - i]] for i in range(k)]
    ax_s.barh(range(k), sizes, height=.5,
              color=[COL[MOTIFS[order[k - 1 - i]]] for i in range(k)])
    ax_s.invert_xaxis()
    ax_s.set_xlabel(f"{unit} positive")
    ax_s.set_xticks(ax_s.get_xticks(), [fmt(v) for v in ax_s.get_xticks()])
    ax_s.tick_params(labelleft=False, left=False)
    ax_s.spines["left"].set_visible(False)
    ax_s.set_ylim(-.5, k - .5)

    ax_b.set_title(note, loc="left", pad=5, fontsize=8)
    ax_b.set_title(f"none of the four: {counts[0]:,} of {total:,}", loc="right", pad=5,
                   fontsize=6.5, color="#444444")
    save(fig, out, stem)
    return combos


def load_cells(with_sections=False):
    """n_cells x len(MOTIFS) boolean of the GMM motif states, four sections pooled"""
    f = h5py.File(H5, "r")
    o = f["obs"]
    sec_cats, sec = read_cat(o, "sample_id")
    keep = np.isin(sec, [list(sec_cats).index(s) for s in SECTIONS])
    S = []
    for k in MOTIFS:
        cats, cc = read_cat(o, f"motif_{k}_state")
        assert list(cats) == ["negative", "positive"], (k, list(cats))
        S.append((cc == 1)[keep])
    f.close()
    M = np.column_stack(S)
    if not with_sections:
        return M
    return M, np.array(sec_cats)[sec[keep]]


def load_vessels(thresh=0.5):
    """n_vessels x len(MOTIFS) boolean: ring positive fraction >= thresh"""
    F = pd.DataFrame(np.load(f"{SRC}/features.npy"),
                     columns=[f"m{k}" for k in RING_MOTIFS])
    return np.column_stack([(F[f"m{k}"] >= thresh).to_numpy() for k in MOTIFS])


def upset_by_section(M, sections, out, stem, note):
    """same columns, one bar row per tissue section, each bar a % of that section's cells"""
    k = len(MOTIFS)
    counts = np.bincount(state_code(M), minlength=1 << k)
    combos = sorted([c for c in range(1, 1 << k) if counts[c] > 0], key=lambda c: -counts[c])
    per = {s: np.bincount(state_code(M[sections == s]), minlength=1 << k) for s in SECTIONS}
    tot = {s: int((sections == s).sum()) for s in SECTIONS}
    frac = {s: 100 * per[s] / tot[s] for s in SECTIONS}
    top = max(frac[s][c] for s in SECTIONS for c in combos) * 1.10
    set_order = sorted(range(k), key=lambda j: -sum(counts[c] for c in range(1 << k)
                                                    if c >> j & 1))

    fig = plt.figure(figsize=(max(104, 6.4 * len(combos) + 44) * MM, 138 * MM))
    gs = GridSpec(len(SECTIONS) + 1, 1, height_ratios=[1] * len(SECTIONS) + [1.6],
                  hspace=.30, figure=fig)
    axes = []
    for i, s in enumerate(SECTIONS):
        ax = fig.add_subplot(gs[i], sharex=axes[0] if axes else None)
        axes.append(ax)
        ax.bar(np.arange(len(combos)), [frac[s][c] for c in combos], width=.62,
               color=[PAL[c] for c in combos], edgecolor="#4d4d4d", lw=.3)
        ax.set_ylim(0, top)
        ax.text(.005, .95, s, transform=ax.transAxes, ha="left", va="top", fontsize=7)
        ax.tick_params(labelbottom=False, bottom=False, labelsize=6)
        ax.spines["bottom"].set_visible(False)
        ax.text(.995, .92, f"none of the four {frac[s][0]:.0f}%   ·   n = {tot[s]:,}",
                transform=ax.transAxes, ha="right", va="top", fontsize=6, color="#444444")
        if i == 0:
            ax.set_title(note, loc="left", pad=4, fontsize=8)
    fig.supylabel("% of the section's cells", fontsize=7, x=.04)

    ax_m = fig.add_subplot(gs[len(SECTIONS)], sharex=axes[0])
    for r, j in enumerate(set_order):
        y = k - 1 - r
        ax_m.axhspan(y - .5, y + .5, color="#f2f2f2" if r % 2 else "white", lw=0, zorder=0)
        for i, c in enumerate(combos):
            on = bool(c >> j & 1)
            ax_m.plot(i, y, "o", ms=5, color=COL[MOTIFS[j]] if on else OFFCOL, zorder=3)
    for i, c in enumerate(combos):
        ys = [k - 1 - r for r, j in enumerate(set_order) if c >> j & 1]
        if len(ys) > 1:
            ax_m.plot([i, i], [min(ys), max(ys)], color=BARCOL, lw=1.3, zorder=2)
    ax_m.set_xlim(-.6, len(combos) - .4)
    ax_m.set_ylim(-.5, k - .5)
    ax_m.set_yticks(range(k), [SHORT[MOTIFS[set_order[k - 1 - i]]] for i in range(k)],
                    fontsize=6.5)
    ax_m.set_xticks([])
    ax_m.tick_params(length=0, pad=2)
    for sp in ax_m.spines.values():
        sp.set_visible(False)
    save(fig, out, stem)
    return pd.DataFrame([dict(section=s, state=bin(c)[2:].zfill(k),
                              combination=state_label(c, short=False),
                              n_cells=int(per[s][c]), pct_of_section=float(frac[s][c]))
                         for s in SECTIONS for c in range(1 << k)])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vessel-thresh", type=float, default=0.5)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    f = h5py.File(H5, "r")
    o = f["obs"]
    sec_cats, sec = read_cat(o, "sample_id")
    keep = np.isin(sec, [list(sec_cats).index(s) for s in SECTIONS])
    S = []
    for k in MOTIFS:
        cats, cc = read_cat(o, f"motif_{k}_state")
        assert list(cats) == ["negative", "positive"], (k, list(cats))
        S.append((cc == 1)[keep])
    f.close()
    Mc = np.column_stack(S)
    cc_counts = combo_counts(Mc)
    log(f"{len(Mc):,} cells; positive per motif: "
        + ", ".join(f"m{k} {Mc[:, j].sum():,}" for j, k in enumerate(MOTIFS)))

    V = pd.read_csv(f"{SRC}/vessels.csv")
    F = pd.DataFrame(np.load(f"{SRC}/features.npy"),
                     columns=[f"m{k}" for k in RING_MOTIFS])
    Mv = np.column_stack([(F[f"m{k}"] >= a.vessel_thresh).to_numpy() for k in MOTIFS])
    vc_counts = combo_counts(Mv)
    Mv9 = np.column_stack([(F[f"m{k}"] >= .9).to_numpy() for k in MOTIFS])
    v9_counts = combo_counts(Mv9)
    log(f"{len(Mv)} vessels; ring fraction >= {a.vessel_thresh}: "
        + ", ".join(f"m{k} {Mv[:, j].sum()}" for j, k in enumerate(MOTIFS)))

    upset(cc_counts, len(Mc), "cells", a.out, "motif_upset_cells",
          note="cells, motif state")
    Ms, secs = load_cells(with_sections=True)
    upset_by_section(Ms, secs, a.out, "motif_upset_cells_by_section",
                     "cells, motif state").round(4).to_csv(
        f"{a.out}/motif_upset_by_section.csv", index=False)
    upset(vc_counts, len(Mv), "vessels", a.out, "motif_upset_vessels",
          note=f"vessels, ring positive fraction $\\geq$ {a.vessel_thresh:g}")

    rows = []
    for c in range(1 << len(MOTIFS)):
        members = [NAME[MOTIFS[j]] for j in range(len(MOTIFS)) if c >> j & 1]
        rows.append(dict(combination=" + ".join(members) if members else "(none)",
                         n_motifs=len(members), n_cells=int(cc_counts[c]),
                         frac_cells=float(cc_counts[c] / len(Mc)),
                         n_vessels_thr0_5=int(vc_counts[c]),
                         n_vessels_thr0_9=int(v9_counts[c])))
    T = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    os.makedirs(a.out, exist_ok=True)
    T.round(5).to_csv(f"{a.out}/motif_upset_counts.csv", index=False)
    print(T.to_string(index=False))
    log(f"wrote to {a.out}/")


if __name__ == "__main__":
    main()
