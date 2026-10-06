#!/usr/bin/env python
"""Whole tissue section, every cell as its Xenium polygon, coloured by motif ON/OFF state.

One figure per section x motif. Cells positive for the motif (obs motif_<k>_state) are filled
solid, negative cells pale grey; the curated vessel lumens are outlined in cyan on top. No
H&E. Geometry comes from {XEN}/{section}_Xenium/cell_boundaries.parquet, which shares the
micron frame of obsm/spatial (verified by cell_id join, centroid difference 0.0).

The cell layer is rasterized inside the PDF/SVG -- a vector layer of ~600,000 polygons would
be unusable -- while the lumens and all text stay vector.

Outputs (figures_section_motif_morphology/): <section>_m<k>, section_motif_cell_counts.csv

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/plot_section_motif_morphology.py --motifs 1,10,23,24
"""
from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import H5, SECTIONS, read_cat, load_polys, log
from find_motif_specific_vessels import XEN, POSCOL, NEGCOL, ANN

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
OUT = f"{BASE}/figures_section_motif_morphology"
MM = 1 / 25.4
NAME = {1: "SMC vascular stabilization motif", 10: "Tumor vasculature motif",
        23: "Vascular homeostasis motif", 24: "Healthy alveolar motif"}
LUMEN_COL = "#00e0ff"
WIDTH_MM = 150

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 9, "legend.fontsize": 7,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def save(fig, out, stem, dpi):
    os.makedirs(out, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{out}/{stem}.{ext}", dpi=dpi)
    plt.close(fig)
    print(f"  {stem}.{{png,pdf,svg}}")


def section_polygons(section):
    """cell polygons of a whole section, in file order, with their cell_ids"""
    t = pq.read_table(f"{XEN}/{section}_Xenium/cell_boundaries.parquet",
                      columns=["cell_id", "vertex_x", "vertex_y"])
    codes, uniq = pd.factorize(t.column("cell_id").to_pandas(), sort=False)
    v = np.column_stack([t.column("vertex_x").to_numpy(), t.column("vertex_y").to_numpy()])
    starts = np.flatnonzero(np.r_[True, codes[1:] != codes[:-1]])
    return np.split(v, starts[1:]), np.asarray(uniq)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--motifs", default="1,10,23,24")
    ap.add_argument("--sections", default="all")
    ap.add_argument("--dpi", type=int, default=400)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    use = [int(t) for t in a.motifs.split(",")]
    secs = SECTIONS if a.sections == "all" else a.sections.split(",")

    f = h5py.File(H5, "r")
    o = f["obs"]
    sec_cats, sec = read_cat(o, "sample_id")
    cid_cats, cid_codes = read_cat(o, "cell_id")
    POS = {}
    for k in use:
        cats, cc = read_cat(o, f"motif_{k}_state")
        assert list(cats) == ["negative", "positive"], (k, list(cats))
        POS[k] = cc == 1
    f.close()
    sec_cats = list(sec_cats)
    cid_cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cid_cats])

    rows = []
    for s in secs:
        polys, ids = section_polygons(s)
        m = np.flatnonzero(sec == sec_cats.index(s))
        row_of = pd.Series(m, index=cid_cats[cid_codes[m]])
        take = row_of.reindex(ids)                       # NaN = cell not in the adata
        have = take.notna().to_numpy()
        rowi = take.to_numpy(dtype="float64")
        polys = [p for p, h in zip(polys, have) if h]
        rowi = rowi[have].astype(int)
        log(f"{s}: {len(ids):,} segmented cells, {have.sum():,} in the adata "
            f"({(~have).sum():,} filtered out and not drawn)")
        lum = {f"{s}__{p['vid']}": p["poly"] for p in load_polys(s)}
        xy_all = np.vstack([p.min(0) for p in polys] + [p.max(0) for p in polys])
        lo, hi = xy_all.min(0), xy_all.max(0)
        span = hi - lo
        for k in use:
            pos = POS[k][rowi]
            fig, ax = plt.subplots(figsize=(WIDTH_MM * MM,
                                            WIDTH_MM * MM * span[1] / span[0]))
            ax.add_collection(PolyCollection(
                polys, facecolors=np.where(pos, POSCOL, NEGCOL), linewidths=0,
                antialiased=False, rasterized=True, zorder=2))
            for q in lum.values():
                ax.add_patch(Polygon(q, closed=True, fc="none", ec=LUMEN_COL, lw=.7,
                                     zorder=4))
            ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])   # y up, not image order
            ax.set_aspect("equal"); ax.axis("off")
            x0, y0 = lo[0] + span[0] * .02, lo[1] + span[1] * .03
            ax.plot([x0, x0 + 1000], [y0, y0], color="#111111", lw=1.8, zorder=5)
            ax.text(x0, y0 + span[1] * .012, "1 mm", fontsize=6.5, va="bottom", zorder=5,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=.75, pad=.6))
            ax.set_title(f"{s}   {NAME.get(k, '')} (m{k})", loc="left", pad=4)
            ax.set_title(f"{pos.sum():,} / {len(pos):,} cells positive", loc="right",
                         fontsize=7, color="#444444", pad=4)
            ax.legend(handles=[Patch(facecolor=POSCOL, label="motif positive"),
                               Patch(facecolor=NEGCOL, label="negative"),
                               Line2D([], [], color=LUMEN_COL, lw=1.2,
                                      label=f"curated lumen (n={len(lum)})")],
                      loc="upper center", bbox_to_anchor=(.5, -.01), ncol=3, frameon=False)
            save(fig, a.out, f"{s}_m{k}", a.dpi)
            rows.append(dict(section=s, motif=f"m{k}", annotation=ANN.get(k, ""),
                             name=NAME.get(k, ""), n_cells_drawn=len(pos),
                             n_positive=int(pos.sum()),
                             frac_positive=float(pos.mean()), n_lumens=len(lum)))
        del polys
    T = pd.DataFrame(rows)
    T.round(4).to_csv(f"{a.out}/section_motif_cell_counts.csv", index=False)
    print(T.round(3).to_string(index=False))
    log(f"wrote to {a.out}/")


if __name__ == "__main__":
    main()
