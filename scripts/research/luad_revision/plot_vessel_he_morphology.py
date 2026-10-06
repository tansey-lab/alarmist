#!/usr/bin/env python
"""Example vessels on H&E, with translucent cell areas coloured by cell type.

One figure per motif, columns = that motif's example vessels (the picks of
find_motif_specific_vessels.py, read from example_vessels.csv). Background: the section's
aligned H&E (<section>_Xenium/HE/aligned_HE_rgb.ome.tif, 0.2125 um/pixel, already registered
to the Xenium frame -- the same frame as obsm/spatial and the boundary parquets). Overlay:
the Xenium cell polygon of every cell, filled translucent in its cell-type colour. The
palette is he_overlay_palette.HE_CELLTYPE_COLORS, built in OKLCH by make_he_palette.py so that
no cell colour falls in the H&E stain's own measured hue band. Nuclei are not drawn. The
curated lumen is cyan.

Outputs (figures_motif_specific_vessels/): m<k>_specific_examples_he

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/plot_vessel_he_morphology.py
"""
from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np
import pandas as pd
import tifffile
import zarr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import H5, read_cat, load_polys, log
from he_overlay_palette import HE_CELLTYPE_COLORS
from find_motif_specific_vessels import (OUT, MM, USE, ANN, XEN, FOV_MIN, FOV_MAX, CONTEXT_UM,
                                         LEG_MIN_FRAC, load_shapes, save)

HE_PX = 0.2125             # um per pixel of aligned_HE_rgb.ome.tif (OME PhysicalSizeX)
CELL_ALPHA = 0.78          # translucent, but the cell type should read at a glance
LUMEN_ALPHA = 0.6          # the curated lumen outline is translucent too
CT_COL = HE_CELLTYPE_COLORS     # the one dict; edit he_overlay_palette.py to tweak
TARGET_PX = 2000           # pyramid level chosen so a crop is about this wide
LUMEN_COL = "#00e0ff"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 8, "legend.fontsize": 6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def he_crop(section, lo, hi):
    """RGB crop of the aligned H&E covering [lo, hi] in microns, plus its true extent"""
    path = f"{XEN}/{section}_Xenium/HE/aligned_HE_rgb.ome.tif"
    with tifffile.TiffFile(path) as t:
        store = t.series[0].aszarr()
        g = zarr.open(store, mode="r")
        want = (hi[0] - lo[0]) / HE_PX
        lvl = 0
        while str(lvl + 1) in list(g.array_keys()) and want / 2 ** (lvl + 1) >= TARGET_PX:
            lvl += 1
        z = g[str(lvl)]
        px = HE_PX * 2 ** lvl
        x0, y0 = int(np.floor(lo[0] / px)), int(np.floor(lo[1] / px))
        x1, y1 = int(np.ceil(hi[0] / px)), int(np.ceil(hi[1] / px))
        x0, y0 = max(x0, 0), max(y0, 0)
        x1, y1 = min(x1, z.shape[1]), min(y1, z.shape[0])
        img = np.asarray(z[y0:y1, x0:x1, :3])
        store.close()
    return img, (x0 * px, x1 * px, y1 * px, y0 * px), lvl


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--motifs", default=",".join(str(k) for k in USE))
    ap.add_argument("--context", type=float, default=CONTEXT_UM)
    ap.add_argument("--vessels", default=None,
                    help="draw these vessel ids instead of the motif picks, e.g. "
                         "P17_AIS__6; the figure is named he_vessels_<ids>")
    ap.add_argument("--fov", type=float, default=None,
                    help="fixed field of view in um (default: lumen extent + 2 x context)")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    use = [int(t) for t in a.motifs.split(",")]

    if a.vessels:
        ids = [v.strip() for v in a.vessels.split(",")]
        area = {f"{s}__{q['vid']}": q["area"]
                for s in {v.rsplit("__", 1)[0] for v in ids} for q in load_polys(s)}
        E = pd.DataFrame([dict(motif="custom", vessel_id=v, section=v.rsplit("__", 1)[0],
                               lumen_area=area.get(v, float("nan"))) for v in ids])
        use = ["custom"]
    else:
        E = pd.read_csv(f"{a.out}/example_vessels.csv")
    f = h5py.File(H5, "r")
    o = f["obs"]
    sec_cats, sec = read_cat(o, "sample_id")
    ct_cats, ct = read_cat(o, "cell_type")
    cid_cats, cid_codes = read_cat(o, "cell_id")
    xy = f["obsm"]["spatial"][:]
    f.close()
    sec_cats, ct_cats = list(sec_cats), list(ct_cats)
    cid_cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cid_cats])
    CT = CT_COL
    assert not set(ct_cats) - set(CT), sorted(set(ct_cats) - set(CT))

    def ctcol(t):
        return CT[ct_cats[t]]

    POLY = {s: {f"{s}__{p['vid']}": p["poly"] for p in load_polys(s)} for s in set(E.section)}
    use = use if not a.vessels else ["custom"]

    for k in use:
        rows = E[E.motif == (k if isinstance(k, str) else f"m{k}")]
        if rows.empty:
            continue
        cols = []
        for _, r in rows.iterrows():
            poly = POLY[r.section][r.vessel_id]
            fov = (a.fov if a.fov else
                   float(np.clip(np.ptp(poly, 0).max() + 2 * a.context, FOV_MIN, FOV_MAX)))
            c = poly.mean(0)
            lo, hi = c - fov / 2, c + fov / 2
            m = np.flatnonzero((sec == sec_cats.index(r.section)) &
                               (xy[:, 0] > lo[0]) & (xy[:, 0] < hi[0]) &
                               (xy[:, 1] > lo[1]) & (xy[:, 1] < hi[1]))
            ids = cid_cats[cid_codes[m]]
            img, ext, lvl = he_crop(r.section, lo, hi)
            cols.append(dict(r=r, poly=poly, lo=lo, hi=hi, fov=fov, m=m, ids=ids,
                             shapes=load_shapes(r.section, lo, hi, set(ids.tolist()),
                                                kinds=("cell",)),
                             img=img, ext=ext))
            log(f"m{k} {r.vessel_id}: {len(m)} cells, {fov:.0f} um field, H&E level {lvl}, "
                f"{img.shape[1]}x{img.shape[0]} px")

        drawn = np.concatenate([c["m"] for c in cols])
        cnt = np.bincount(ct[drawn], minlength=len(ct_cats))
        share = cnt / max(cnt.sum(), 1)
        keep = [t for t in np.argsort(-cnt) if share[t] >= LEG_MIN_FRAC]

        fig, axes = plt.subplots(1, len(cols), figsize=(len(cols) * 55 * MM, 58 * MM),
                                 squeeze=False)
        for j, C in enumerate(cols):
            ax = axes[0, j]
            ax.imshow(C["img"], extent=C["ext"], interpolation="bilinear", zorder=1)
            P, FC = [], []
            for i, cid in enumerate(C["ids"]):
                for q in C["shapes"]["cell"].get(cid, ()):
                    P.append(q)
                    FC.append(ctcol(ct[C["m"][i]]))
            if P:
                ax.add_collection(PolyCollection(P, facecolors=FC, alpha=CELL_ALPHA,
                                                 linewidths=0, zorder=3))
            for vid, q in POLY[C["r"].section].items():
                if (q[:, 0].max() < C["lo"][0] or q[:, 0].min() > C["hi"][0] or
                        q[:, 1].max() < C["lo"][1] or q[:, 1].min() > C["hi"][1]):
                    continue
                subj = vid == C["r"].vessel_id
                ax.add_patch(Polygon(q, closed=True, fc="none", zorder=5,
                                     ec=LUMEN_COL if subj else "#9e9e9e",
                                     alpha=LUMEN_ALPHA, lw=1.8 if subj else .8))
            ax.set_xlim(C["lo"][0], C["hi"][0])
            ax.set_ylim(C["hi"][1], C["lo"][1])          # image row order: y increases down
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            L = 100.0 if C["fov"] >= 300 else 50.0
            x0 = C["lo"][0] + C["fov"] * .06
            y0 = C["hi"][1] - C["fov"] * .06
            ax.plot([x0, x0 + L], [y0, y0], color="#111111", lw=1.8, zorder=6)
            ax.text(x0, y0 - C["fov"] * .02, f"{L:.0f} " + r"$\mu$m", fontsize=6, zorder=6,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=.7, pad=.6))
            ax.set_title(f"{C['r'].vessel_id}\n{C['r'].lumen_area:,.0f} " + r"$\mu$m$^2$"
                         + f"   ·   {C['fov']:.0f} " + r"$\mu$m field", fontsize=6, pad=3,
                         linespacing=1.3)
        hs = [Patch(facecolor=ctcol(t), alpha=CELL_ALPHA, edgecolor="none",
                    label=f"{ct_cats[t]} {share[t]:.0%}") for t in keep]
        rest = [t for t in range(len(ct_cats)) if t not in keep and cnt[t]]
        if rest:
            hs.append(Line2D([], [], ls="", marker="", label=f"+{len(rest)} types <2% "
                             f"({share[rest].sum():.0%} of cells)"))
        hs += [Line2D([], [], color=LUMEN_COL, lw=1.8, alpha=LUMEN_ALPHA,
                      label="curated lumen")]
        fig.subplots_adjust(bottom=.08, top=.88)
        fig.legend(handles=hs, frameon=False, ncol=6, loc="upper center", fontsize=6,
                   handletextpad=.4, columnspacing=1.0, bbox_to_anchor=(.5, .04),
                   bbox_transform=fig.transFigure)
        if isinstance(k, str):      # the column titles already name the vessels
            save(fig, "he_vessels_" + "_".join(rows.vessel_id))
        else:
            fig.suptitle(f"m{k} {ANN.get(k, '')}", fontsize=8, x=.02, ha="left")
            save(fig, f"m{k}_specific_examples_he")
    log(f"wrote to {a.out}/")


if __name__ == "__main__":
    main()
