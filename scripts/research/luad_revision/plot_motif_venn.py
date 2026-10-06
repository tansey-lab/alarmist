#!/usr/bin/env python
"""Four-set Venn of the vessel motifs, in the same 16-state colours as the UpSet and the maps.

Four sets need ellipses, not circles: no arrangement of four circles produces all 15
intersections. The four congruent ellipses here are the classic Venn layout, and every region
is filled with its state colour from motif_state_palette, so a region's colour matches the
same combination in <section>_motif_states and motif_upset_*. The background is the
all-negative colour.

**The areas are not proportional to the counts** -- for four sets that is impossible; the
Venn shows which combinations exist and the printed count gives the size. The UpSet plot is
the one to read sizes off.

Regions are filled by rasterising the four ellipses on a grid and colouring each pixel by its
state, and each count is placed at the region's pole of inaccessibility (the pixel furthest
from that region's edge, via a distance transform), so labels never fall outside their region.

Outputs (figures_motif_upset/): motif_venn_cells, motif_venn_vessels

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/plot_motif_venn.py
"""
from __future__ import annotations

import argparse
import os
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Ellipse
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import log
from motif_state_palette import ANCHOR, MOTIFS, NAME, SHORT, state_code, state_palette
from plot_motif_upset import OUT, load_cells, load_vessels

MM = 1 / 25.4
# classic four-ellipse Venn: same ellipse, two mirrored pairs; (cx, cy, width, height, angle)
# searched over width / height / angle / offsets for a layout where all 15 regions exist and
# the smallest is as large as possible (6,137 px of a 700x700 grid; the common pyvenn layout
# gives 3,399), so every count has room inside its own region
ELLIPSES = [(0.365, 0.47, 0.76, 0.40, 145.0),
            (0.455, 0.55, 0.76, 0.40, 145.0),
            (0.545, 0.55, 0.76, 0.40, 35.0),
            (0.635, 0.47, 0.76, 0.40, 35.0)]
GRID = 1400

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 8, "legend.fontsize": 6.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def save(fig, out, stem):
    os.makedirs(out, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{out}/{stem}.{ext}", dpi=400 if ext == "png" else None)
    plt.close(fig)
    print(f"  {stem}.{{png,pdf,svg}}")


def ellipse_masks(n=GRID):
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x)
    out = []
    for cx, cy, w, h, ang in ELLIPSES:
        t = np.deg2rad(ang)
        dx, dy = X - cx, Y - cy
        u = dx * np.cos(t) + dy * np.sin(t)
        v = -dx * np.sin(t) + dy * np.cos(t)
        out.append((u / (w / 2)) ** 2 + (v / (h / 2)) ** 2 <= 1)
    return np.array(out)


def fmt(v, total):
    if total > 20_000 and v >= 10_000:
        return f"{v / 1000:.0f}k"
    return f"{v:,}"


def venn(counts, total, unit, out, stem, note):
    pal = state_palette()
    M = ellipse_masks()
    code = np.zeros(M.shape[1:], np.int64)
    for j in range(len(MOTIFS)):
        code += M[j].astype(np.int64) << j
    img = np.zeros(code.shape + (3,))
    for c in range(1 << len(MOTIFS)):
        img[code == c] = to_rgb(pal[c])
    missing = [c for c in range(1, 1 << len(MOTIFS)) if not (code == c).any()]
    assert not missing, f"layout misses regions: {missing}"

    fig, ax = plt.subplots(figsize=(108 * MM, 96 * MM))
    ax.imshow(img, extent=(0, 1, 0, 1), origin="lower", interpolation="nearest", zorder=1)
    for (cx, cy, w, h, ang), k in zip(ELLIPSES, MOTIFS):
        ax.add_patch(Ellipse((cx, cy), w, h, angle=ang, fc="none", ec=ANCHOR[k], lw=1.3,
                             zorder=3))
    for c in range(1, 1 << len(MOTIFS)):
        m = code == c
        d = distance_transform_edt(np.pad(m, 1))[1:-1, 1:-1]
        iy, ix = np.unravel_index(np.argmax(d), d.shape)
        L = 0.299 * img[iy, ix, 0] + 0.587 * img[iy, ix, 1] + 0.114 * img[iy, ix, 2]
        ax.text(ix / (GRID - 1), iy / (GRID - 1), fmt(counts[c], total), ha="center",
                va="center", fontsize=6.2, zorder=4,
                color="white" if L < 0.55 else "#1a1a1a")
    ax.text(.02, .03, f"none of the four\n{fmt(counts[0], total):>}", fontsize=6.2,
            color="#444444", ha="left", va="bottom", zorder=4)
    # set names, placed outside the four ellipse tips
    for (cx, cy, w, h, ang), k, (tx, ty, ha, va) in zip(
            ELLIPSES, MOTIFS, [(.02, .80, "left", "bottom"), (.16, .97, "left", "bottom"),
                               (.84, .97, "right", "bottom"), (.98, .80, "right", "bottom")]):
        ax.text(tx, ty, textwrap.fill(SHORT[k], 16), color=ANCHOR[k], fontsize=6.5,
                ha=ha, va=va, zorder=4, linespacing=1.25)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.06)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(note, loc="left", pad=4)
    ax.set_title(f"n = {total:,}", loc="right", pad=4, fontsize=6.5, color="#444444")
    save(fig, out, stem)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vessel-thresh", type=float, default=0.5)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    Mc = load_cells()
    cc = np.bincount(state_code(Mc), minlength=1 << len(MOTIFS))
    log(f"{len(Mc):,} cells")
    venn(cc, len(Mc), "cells", a.out, "motif_venn_cells", "cells, motif state")

    Mv = load_vessels(a.vessel_thresh)
    vc = np.bincount(state_code(Mv), minlength=1 << len(MOTIFS))
    log(f"{len(Mv)} vessels")
    venn(vc, len(Mv), "vessels", a.out, "motif_venn_vessels",
         f"vessels, ring positive fraction $\\geq$ {a.vessel_thresh:g}")
    log(f"wrote to {a.out}/")


if __name__ == "__main__":
    main()
