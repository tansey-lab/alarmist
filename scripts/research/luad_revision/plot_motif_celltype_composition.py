#!/usr/bin/env python
"""Cell-type composition of the cells positive for a motif (GMM ON/OFF state).

For each motif in --motifs: the composition of the cells whose motif_<k>_state is
"positive", pooled over the four sections. --cells ring restricts to cells 0 < d <= 30 um
outside a curated lumen; --cells all (default) uses every cell in the four sections.

Colours are tab20 indexed by the cell-type category order, the same mapping as
vessel_niche_cluster.plot_cluster_cell_types and the vessel crop figures, so a cell type
keeps its colour across all of them.

Outputs (figures_motif_celltype/):
    m<k>_celltype_pie          one pie per motif, wedges in category order
    motif_celltype_stacked     the same compositions as stacked proportions, one figure
    motif_celltype_composition.csv   proportions and counts, plus the all-cell background

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/plot_motif_celltype_composition.py
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
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import H5, SECTIONS, read_cat, load_polys, signed_distance, log
from he_overlay_palette import HE_CELLTYPE_COLORS

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
OUT = f"{BASE}/figures_motif_celltype"
MM = 1 / 25.4
ANN = {0: "T cell", 1: "SMC", 6: "macrophage", 10: "tumour-endo",
       23: "endothelial", 24: "alveolar"}
from motif_state_palette import NAME      # one source for the motif display names
VASC_MOTIFS = [0, 1, 6, 10, 23, 24]
MCOL = ["#0072b2", "#e69f00", "#009e73", "#cc79a7", "#56b4e9", "#d55e00"]  # validated PASS
BAR_CT = ["Tumor_epi", "Vas_Endo", "SMC", "Alveolar_epi"]   # the four asked for
R_RING = 30.0
PIE_LABEL_MIN = 0.05       # wedges below this share are not labelled on the pie

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
    return f"m{k}" + (f" {ANN[k]}" if k in ANN else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--motifs", default="1,23,24,10")   # the order used by the density plot
    ap.add_argument("--cells", choices=("all", "ring"), default="all",
                    help="all cells, or only cells within 30 um of a curated lumen")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    use = [int(t) for t in a.motifs.split(",")]
    out = a.out + ("_ring" if a.cells == "ring" else "")

    f = h5py.File(H5, "r")
    o = f["obs"]
    sec_cats, sec = read_cat(o, "sample_id")
    ct_cats, ct = read_cat(o, "cell_type")
    xy = f["obsm"]["spatial"][:] if a.cells == "ring" else None
    POS = {}
    for k in use:
        cats, cc = read_cat(o, f"motif_{k}_state")
        assert list(cats) == ["negative", "positive"], (k, list(cats))
        POS[k] = cc == 1
    f.close()
    ct_cats = list(ct_cats)
    keep = np.isin(sec, [sec_cats.index(s) for s in SECTIONS])
    if a.cells == "ring":
        D = np.full(len(sec), np.nan)
        for s in SECTIONS:
            m = np.where(sec == sec_cats.index(s))[0]
            D[m] = signed_distance(xy[m], load_polys(s))[0]
        keep &= (D > 0) & (D <= R_RING)
    log(f"{keep.sum():,} cells in scope ({a.cells}); {len(ct_cats)} cell types")

    CT = HE_CELLTYPE_COLORS      # same dict as the H&E overlays (he_overlay_palette.py)
    assert not set(ct_cats) - set(CT), sorted(set(ct_cats) - set(CT))
    col = [CT[c] for c in ct_cats]

    # ---- composition per motif --------------------------------------------------
    P, rows = {}, []
    bg = np.bincount(ct[keep], minlength=len(ct_cats)).astype(float)
    rows.append(dict(group="all cells", n_cells=int(keep.sum()),
                     **{ct_cats[i]: bg[i] / bg.sum() for i in range(len(ct_cats))}))
    for k in use:
        m = keep & POS[k]
        c = np.bincount(ct[m], minlength=len(ct_cats)).astype(float)
        P[k] = c / max(c.sum(), 1)
        rows.append(dict(group=f"m{k}", annotation=ANN.get(k, ""), n_cells=int(m.sum()),
                         **{ct_cats[i]: P[k][i] for i in range(len(ct_cats))}))
        top = np.argsort(-P[k])[:3]
        log(f"m{k} {ANN.get(k, ''):<12} {int(m.sum()):>8,} positive cells | "
            + ", ".join(f"{ct_cats[t]} {P[k][t]:.0%}" for t in top))

    # ---- pies -------------------------------------------------------------------
    for k in use:
        n = int((keep & POS[k]).sum())
        fig, ax = plt.subplots(figsize=(52 * MM, 52 * MM))
        # no cell-type names on the wedges: the colours are the stacked-bar legend's
        ax.pie(P[k], colors=col, startangle=90, counterclock=False,
               autopct=lambda v: f"{v:.0f}%" if v >= PIE_LABEL_MIN * 100 else "",
               pctdistance=.72, radius=1.0,
               textprops=dict(fontsize=6), wedgeprops=dict(lw=.3, edgecolor="white"))
        ax.set_title(f"{mlab(k)}   n = {n:,}", loc="left", pad=6)
        save(fig, out, f"m{k}_celltype_pie")

    # ---- stacked bar, format of vessel_niche_cluster.plot_cluster_cell_types -----
    prop = np.vstack([P[k] for k in use])
    fig, ax = plt.subplots(figsize=(max(4, .55 * len(use) + 2), 4))
    x = np.arange(len(use))
    bottom = np.zeros(len(use))
    for i in range(len(ct_cats)):
        ax.bar(x, prop[:, i], bottom=bottom, label=ct_cats[i], color=col[i], width=.75)
        bottom += prop[:, i]
    ax.set_xticks(x)
    ax.set_xticklabels([mlab(k) for k in use], rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Cell type proportions")
    ax.set_xlabel("motif")
    ax.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=6.5,
              title_fontsize=7, frameon=False)
    save(fig, out, "motif_celltype_stacked")

    os.makedirs(out, exist_ok=True)
    T = pd.DataFrame(rows)
    T.round(5).to_csv(f"{out}/motif_celltype_composition.csv", index=False)

    # ---- per section: % of each cell type positive for each motif ------------------
    ct_idx = {t: ct_cats.index(t) for t in BAR_CT if t in ct_cats}
    assert len(ct_idx) == len(BAR_CT), [t for t in BAR_CT if t not in ct_cats]
    bar = []
    for s in SECTIONS:
        ks = keep & (sec == sec_cats.index(s))
        for t, ti in ct_idx.items():
            base = ks & (ct == ti)
            for k in use:
                n = int(base.sum())
                bar.append(dict(section=s, cell_type=t, motif=f"m{k}",
                                name=NAME.get(k, f"m{k}"), n_cells=n,
                                n_positive=int((base & POS[k]).sum()),
                                frac_positive=float((POS[k][base]).mean()) if n else np.nan,
                                share_of_positive_cells=float(
                                    (base & POS[k]).sum() / max((ks & POS[k]).sum(), 1))))
    B = pd.DataFrame(bar)
    B.round(5).to_csv(f"{out}/celltype_positive_fraction.csv", index=False)
    # one column per motif; rows are the four cell types, pooled over sections.
    # No axes: the thin black frame is the 0-100 % scale, the fill is the value.
    tot = B.groupby(["cell_type", "motif"], observed=True)[["n_cells", "n_positive"]].sum()
    tot["frac_positive"] = tot.n_positive / tot.n_cells
    fig, axes = plt.subplots(1, len(use), figsize=(len(use) * 62 * MM, 30 * MM),
                             sharex=True, sharey=True, squeeze=False)
    for j, k in enumerate(use):
        ax = axes[0, j]
        for i, t in enumerate(BAR_CT):
            v = 100 * float(tot.loc[(t, f"m{k}"), "frac_positive"])
            ax.add_patch(Rectangle((0, i - .34), 100, .68, facecolor="none",
                                   edgecolor="#1a1a1a", lw=.5, zorder=3))
            ax.add_patch(Rectangle((0, i - .34), v, .68, facecolor=col[ct_cats.index(t)],
                                   edgecolor="none", zorder=2))
            inside = v > 72                      # long bars get the label inside the fill
            ax.text(v + (-1.8 if inside else 1.8), i, f"{v:.2f}%", va="center",
                    ha="right" if inside else "left", fontsize=5.8, color="#1a1a1a",
                    zorder=4)
        ax.set_xlim(-2, 102)
        ax.set_ylim(len(BAR_CT) - .5, -.5)
        ax.set_xticks([])
        ax.tick_params(length=0)
        if j == 0:                      # sharey: only the first column carries the labels
            ax.set_yticks(range(len(BAR_CT)), BAR_CT)
        else:
            ax.tick_params(labelleft=False)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(textwrap.fill(NAME.get(k, f"m{k}"), 18) + f"\n(m{k})", fontsize=6.5,
                     pad=4, linespacing=1.3)
    fig.subplots_adjust(wspace=.12)
    save(fig, out, "celltype_positive_fraction")
    log("% of each cell type positive (pooled; per-section numbers in the csv):\n"
        + (100 * tot.frac_positive).round(1).unstack().to_string())
    log(f"wrote to {out}/")


if __name__ == "__main__":
    main()
