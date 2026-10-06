#!/usr/bin/env python
"""Whole tissue sections coloured by the 16-state combination of the four vessel motifs.

Every cell is drawn as its Xenium polygon and coloured by which of the four motifs it is
positive for -- one of 16 states -- using the shared OKLab palette of motif_state_palette.py
(singles keep their established motif colour, multi-positives are the OKLab mean darkened by
LSTEP per extra motif, all-negative is very light grey). All four sections share the palette,
so colours are comparable between them; the same {state: colour} map is what the UpSet plot
uses. Cells are drawn in order of how many motifs they carry, so multi-positive cells are not
hidden under single-positive ones. Curated lumens keep the bright cyan outline of the
per-motif maps.

Outputs (figures_section_motif_morphology/): <section>_motif_states,
    motif_state_palette.csv (the 16 colours), motif_state_counts_by_section.csv

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/plot_section_motif_states.py
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
from matplotlib.collections import PolyCollection
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import H5, SECTIONS, read_cat, load_polys, log
from motif_state_palette import (MOTIFS, NAME, draw_state_legend, palette_separation,
                                 state_code, state_label, state_palette)
from plot_section_motif_morphology import LUMEN_COL, OUT, WIDTH_MM, section_polygons

MM = 1 / 25.4

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 9, "legend.fontsize": 6.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def save(fig, out, stem, dpi):
    os.makedirs(out, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{out}/{stem}.{ext}", dpi=dpi)
    plt.close(fig)
    print(f"  {stem}.{{png,pdf,svg}}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sections", default="all")
    ap.add_argument("--lstep", type=float, default=None,
                    help="OKLab lightness removed per extra motif; 0 = the plain mean")
    ap.add_argument("--chroma", choices=("average", "restore"), default="average")
    ap.add_argument("--dpi", type=int, default=400)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    secs = SECTIONS if a.sections == "all" else a.sections.split(",")
    pal = (state_palette(chroma=a.chroma) if a.lstep is None
           else state_palette(chroma=a.chroma, lstep=a.lstep))
    sep, pair, _ = palette_separation(pal)
    log(f"palette: closest two states {sep:.3f} apart in OKLab "
        f"[{state_label(pair[0])}] vs [{state_label(pair[1])}]")

    f = h5py.File(H5, "r")
    o = f["obs"]
    sec_cats, sec = read_cat(o, "sample_id")
    cid_cats, cid_codes = read_cat(o, "cell_id")
    S = []
    for k in MOTIFS:
        cats, cc = read_cat(o, f"motif_{k}_state")
        assert list(cats) == ["negative", "positive"], (k, list(cats))
        S.append(cc == 1)
    f.close()
    sec_cats = list(sec_cats)
    cid_cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cid_cats])
    code_all = state_code(np.column_stack(S))

    rows = []
    for s in secs:
        polys, ids = section_polygons(s)
        m = np.flatnonzero(sec == sec_cats.index(s))
        take = pd.Series(m, index=cid_cats[cid_codes[m]]).reindex(ids)
        have = take.notna().to_numpy()
        polys = [p for p, h in zip(polys, have) if h]
        code = code_all[take.to_numpy(dtype="float64")[have].astype(int)]
        cnt = {c: int((code == c).sum()) for c in range(1 << len(MOTIFS))}
        log(f"{s}: {len(polys):,} cells drawn ({(~have).sum():,} not in the adata)")

        order = np.argsort([bin(c).count("1") for c in code], kind="stable")  #多 motif 在上
        xy_all = np.vstack([p.min(0) for p in polys] + [p.max(0) for p in polys])
        lo, hi = xy_all.min(0), xy_all.max(0)
        span = hi - lo
        fig = plt.figure(figsize=(WIDTH_MM * MM,
                                  WIDTH_MM * MM * span[1] / span[0] + 34 * MM))
        gs = GridSpec(2, 1, height_ratios=[span[1] / span[0] * WIDTH_MM, 30], hspace=.06,
                      figure=fig)
        ax = fig.add_subplot(gs[0])
        ax.add_collection(PolyCollection([polys[i] for i in order],
                                         facecolors=[pal[c] for c in code[order]],
                                         linewidths=0, antialiased=False, rasterized=True,
                                         zorder=2))
        for q in {f"{s}__{p['vid']}": p["poly"] for p in load_polys(s)}.values():
            ax.add_patch(Polygon(q, closed=True, fc="none", ec=LUMEN_COL, lw=.7,
                                 zorder=4))
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])   # y up, not image order
        ax.set_aspect("equal"); ax.axis("off")
        # inside the axes: below the data limits it is clipped by the gridspec, and just
        # outside it collides with the legend's count labels
        bx, by = lo[0] + span[0] * .02, lo[1] + span[1] * .03
        ax.plot([bx, bx + 1000], [by, by], color="#111111", lw=1.8, zorder=6)
        ax.text(bx, by + span[1] * .012, "1 mm", fontsize=6.5, va="bottom", zorder=6,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.75, pad=.6))
        ax.set_title(s, loc="left", pad=4)
        ax.set_title(f"{cnt[0]:,} of {len(code):,} cells carry none of the four", loc="right",
                     fontsize=7, color="#444444", pad=4)
        draw_state_legend(fig.add_subplot(gs[1]), pal, counts=cnt)
        save(fig, a.out, f"{s}_motif_states", a.dpi)
        for c in range(1 << len(MOTIFS)):
            rows.append(dict(section=s, state=bin(c)[2:].zfill(len(MOTIFS)),
                             combination=state_label(c, short=False), colour=pal[c],
                             n_motifs=bin(c).count("1"), n_cells=cnt[c],
                             frac_cells=cnt[c] / len(code)))
        del polys
    T = pd.DataFrame(rows)
    T.round(5).to_csv(f"{a.out}/motif_state_counts_by_section.csv", index=False)
    pd.DataFrame([dict(state=bin(c)[2:].zfill(len(MOTIFS)),
                       combination=state_label(c, short=False), colour=pal[c])
                  for c in sorted(pal)]).to_csv(f"{a.out}/motif_state_palette.csv",
                                                index=False)
    log(f"wrote to {a.out}/")


if __name__ == "__main__":
    main()
