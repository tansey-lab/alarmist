#!/usr/bin/env python
"""The cleanest single-motif vessels: high for one motif, low for the others.

Input: per-vessel 30 um ring positive fractions f_j(v) (GMM state) for the six vessel motifs.

Selection rule, per target motif k (raw fractions, no ranks):

    others(k)  = the other five motifs, minus m23 unless k = 23
    eligible   f_k(v) >= --min-frac (0.5, the same cut as the size_kde positive vessel)
    order      max_{j in others(k)} f_j(v) ascending;
               ties -> higher f_k, then more ring cells, then vessel id
    picked     the first --top eligible vessels

m23 (endothelial) has ring fraction >= 0.5 in 290 / 300 vessels (median 1.00), so as a
competitor it would dominate every other motif's ordering; it is left out of others(k) and
shown separately. Counts of eligible vessels with max-other <= 0.10 / 0.25 / 0.50 and
max-other < f_k are written to motif_specific_vessels.csv. The lowest max-other among
eligible vessels is 0.38 for m6 and 0.43 for m24: their picks are the least contaminated,
not clean.

Superseded (2026-09-17): a within-section percentile-rank margin. Rank ties at 0 and 1 and
per-section baselines let it pick vessels whose target fraction is below a competitor's
(e.g. m10 0.23 vs m1 0.95); outputs kept in superseded_rank_margin/.

Figures: per motif, rows = the six motifs (GMM ON/OFF state -- positive solid, negative
grey; loading not drawn) plus a cell-type row, columns = the picked vessels.
fov = clip(lumen extent + 2 x CONTEXT_UM, FOV_MIN, FOV_MAX); panels differ in scale, so every
column carries its own scale bar and its field size in the title.

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/find_motif_specific_vessels.py
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
import pyarrow.compute as pc
import pyarrow.dataset as pds
from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_vessel_alignment import H5, SECTIONS, read_cat, load_polys, log
from celltype_palette import celltype_colors

BASE = "/data1/tanseyw/projects/fanj2/alarmist_luad"
SRC = f"{BASE}/vessel_kmeans_ring30_m0-1-6-10-23-24_frac"
OUT = f"{BASE}/figures_motif_specific_vessels"
MM = 1 / 25.4
USE = [0, 1, 6, 10, 23, 24]
ANN = {0: "T cell", 1: "SMC", 6: "macrophage", 10: "tumour-endo",
       23: "endothelial", 24: "alveolar"}
UBIQ = 23          # positive in ~every vessel; excluded from the competitor set
MAX_OTHER_CUTS = (0.10, 0.25, 0.50)
POSCOL, NEGCOL = "#3b0f70", "#e2e2e2"   # GMM positive / negative
FOV_MIN, FOV_MAX, CONTEXT_UM = 350.0, 1200.0, 150.0
LEG_MIN_FRAC = 0.02        # a cell type below this share of the drawn cells joins "other"
XEN = "/data1/tanseyw/projects/spatial_data/linghua/j.ccell.2025.10.004"
CELL_A, NUC_A = 0.45, 0.95   # cell body / nucleus fill alpha in morphology mode

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def save(fig, stem):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{OUT}/{stem}.{ext}", dpi=400 if ext == "png" else None)
    plt.close(fig)
    print(f"  {stem}.{{png,pdf,svg}}")


def others(k):
    """competitor motifs for target k: the other five, minus m23 unless k is m23"""
    return [j for j in USE if j != k and (j != UBIQ or k == UBIQ)]


def select(V, F, k, min_frac=0.5, top=3):
    """eligible = f_k >= min_frac; order by max-other asc, then f_k desc, n_ring desc, id.

    Returns (picked row indices, max-other per vessel, eligible mask)."""
    fk = F[f"m{k}"].to_numpy()
    mo = F[[f"m{j}" for j in others(k)]].max(1).to_numpy()
    elig = fk >= min_frac
    idx = np.flatnonzero(elig)
    order = np.lexsort((V.vessel_id.to_numpy()[idx], -V.n_ring.to_numpy()[idx],
                        -fk[idx], mo[idx]))
    return idx[order][:top], mo, elig


def load_shapes(section, lo, hi, want, pad=30.0, kinds=("cell", "nucleus")):
    """Xenium cell and nucleus polygons inside the crop, keyed by cell_id.

    The boundary parquets are in the same micron frame as obsm/spatial (verified: centroid
    difference 0.0 for all 1,676,162 cells), so no transform is applied. A cell may carry
    several nucleus rings (multinucleate), one per label_id.
    """
    out = {}
    for kind in kinds:
        d = pds.dataset(f"{XEN}/{section}_Xenium/{kind}_boundaries.parquet", format="parquet")
        flt = ((pc.field("vertex_x") >= lo[0] - pad) & (pc.field("vertex_x") <= hi[0] + pad) &
               (pc.field("vertex_y") >= lo[1] - pad) & (pc.field("vertex_y") <= hi[1] + pad))
        t = d.to_table(columns=["cell_id", "vertex_x", "vertex_y", "label_id"],
                       filter=flt).to_pandas()
        t = t[t.cell_id.isin(want)]
        polys = {}
        for (cid, _), v in t.groupby(["cell_id", "label_id"], sort=False):
            polys.setdefault(cid, []).append(v[["vertex_x", "vertex_y"]].to_numpy())
        out[kind] = polys
    return out


def draw_morphology(ax, ids, shapes, colors):
    """cell bodies then nuclei, each filled with that cell's colour"""
    drawn = 0
    for kind, alpha, lw in (("cell", CELL_A, .12), ("nucleus", NUC_A, 0.0)):
        P, C = [], []
        for i, cid in enumerate(ids):
            for poly in shapes[kind].get(cid, ()):
                P.append(poly)
                C.append(colors[i])
        if kind == "cell":
            drawn = len({c for c in ids if shapes["cell"].get(c)})
        if P:
            ax.add_collection(PolyCollection(P, facecolors=C, alpha=alpha, linewidths=lw,
                                             edgecolors="white", zorder=2))
    return drawn


def main():
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=3, help="example vessels drawn per motif")
    ap.add_argument("--fov-min", type=float, default=FOV_MIN)
    ap.add_argument("--fov-max", type=float, default=FOV_MAX)
    ap.add_argument("--context", type=float, default=CONTEXT_UM,
                    help="um of tissue kept beyond the lumen on every side")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="eligible: target ring positive fraction >= this")
    ap.add_argument("--render", choices=("morphology", "points"), default="morphology",
                    help="cell/nucleus polygons from the Xenium boundary parquets, or dots")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    OUT = a.out

    V = pd.read_csv(f"{SRC}/vessels.csv")
    F = pd.DataFrame(np.load(f"{SRC}/features.npy"), columns=[f"m{k}" for k in USE])

    rows, picks = [], {}
    for k in USE:
        sel, mo, elig = select(V, F, k, a.min_frac, a.top)
        picks[k] = sel
        fk = F[f"m{k}"].to_numpy()
        rows.append(dict(
            motif=f"m{k}", annotation=ANN[k], n_eligible=int(elig.sum()),
            **{f"n_max_other_le_{c:.2f}": int((elig & (mo <= c)).sum()) for c in MAX_OTHER_CUTS},
            n_target_top=int((elig & (mo < fk)).sum()),
            best_vessel=V.vessel_id.iloc[sel[0]], best_frac=float(fk[sel[0]]),
            best_max_other=float(mo[sel[0]])))
        log(f"m{k} {ANN[k]:<12} eligible {rows[-1]['n_eligible']:>3} | max-other <=.25: "
            f"{rows[-1]['n_max_other_le_0.25']:>3} | target top: {rows[-1]['n_target_top']:>3}"
            f" | best {rows[-1]['best_vessel']} ({fk[sel[0]]:.2f} vs {mo[sel[0]]:.2f})")
    T = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    T.round(4).to_csv(f"{OUT}/motif_specific_vessels.csv", index=False)
    per = []
    for k in USE:
        _, mo, _ = select(V, F, k, a.min_frac, a.top)
        for n, v in enumerate(picks[k], 1):
            per.append(dict(motif=f"m{k}", pick=n, vessel_id=V.vessel_id.iloc[v],
                            section=V.section.iloc[v], target=float(F.iloc[v][f"m{k}"]),
                            max_other=float(mo[v]),
                            n_ring=int(V.n_ring.iloc[v]),
                            lumen_area=float(V.lumen_area.iloc[v]),
                            **{f"m{j}": float(F.iloc[v][f"m{j}"]) for j in USE}))
    pd.DataFrame(per).round(4).to_csv(f"{OUT}/example_vessels.csv", index=False)

    # ------------------------------------------------------------------ the crops
    f = h5py.File(H5, "r")
    o = f["obs"]
    sec_cats, sec = read_cat(o, "sample_id")
    ct_cats, ct = read_cat(o, "cell_type")
    cid_cats, cid_codes = read_cat(o, "cell_id")
    cid_cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cid_cats])
    xy = f["obsm"]["spatial"][:]
    # the GMM ON/OFF call straight from obs -- no loading, no percentile ramp
    POS = {}
    for k in USE:
        cats, cc = read_cat(o, f"motif_{k}_state")
        assert list(cats) == ["negative", "positive"], (k, list(cats))
        POS[k] = cc == 1
    f.close()
    ct_cats = list(ct_cats)
    CT = celltype_colors(ct_cats)      # the run's own palette, not plain tab20
    POLY = {s: {f"{s}__{p['vid']}": p["poly"] for p in load_polys(s)} for s in SECTIONS}

    def nice_bar(fov):
        """scale-bar length: the largest round number under a third of the field"""
        for L in (500.0, 200.0, 100.0, 50.0, 20.0):
            if L <= fov / 3:
                return L
        return 20.0

    rws = USE + [None]
    for i, k in enumerate(USE):
        sel = picks[k]
        # geometry first, so the cell-type legend can be built from what is actually drawn
        cols = []
        for v in sel:
            vid = V.vessel_id.iloc[v]
            sname = vid.rsplit("__", 1)[0]
            poly = POLY[sname][vid]
            ext = float(np.max(poly.max(0) - poly.min(0)))
            fov = float(np.clip(ext + 2 * a.context, a.fov_min, a.fov_max))
            cx, cy = poly.mean(0)
            lo = np.array([cx - fov / 2, cy - fov / 2]); hi = lo + fov
            m = np.where((sec == sec_cats.index(sname)) &
                         (xy[:, 0] > lo[0]) & (xy[:, 0] < hi[0]) &
                         (xy[:, 1] > lo[1]) & (xy[:, 1] < hi[1]))[0]
            cols.append(dict(v=v, vid=vid, s=sname, lo=lo, hi=hi, fov=fov, m=m, ext=ext))

        drawn = np.concatenate([c["m"] for c in cols])
        cnt = np.bincount(ct[drawn], minlength=len(ct_cats))
        share = cnt / max(cnt.sum(), 1)
        keep = [t for t in np.argsort(-cnt) if share[t] >= LEG_MIN_FRAC]

        fig, axes = plt.subplots(len(rws), len(cols),
                                 figsize=(len(cols) * 46 * MM, len(rws) * 43 * MM),
                                 squeeze=False)
        for col, C in enumerate(cols):
            m, lo, hi, fov = C["m"], C["lo"], C["hi"], C["fov"]
            dot = float(np.clip(2.6 * (350.0 / fov) ** .5, 1.0, 3.2))
            if a.render == "morphology":
                ids = cid_cats[cid_codes[m]]
                shapes = load_shapes(C["s"], lo, hi, set(ids.tolist()))
                got = len({c for c in ids if shapes["cell"].get(c)})
                if got < len(ids):
                    log(f"  {C['vid']}: {len(ids) - got} of {len(ids)} cells have no cell "
                        f"boundary in the parquet")
            for r, j in enumerate(rws):
                ax = axes[r, col]
                colours = ([CT[ct_cats[t]] for t in ct[m]] if j is None else
                           [POSCOL if v else NEGCOL for v in POS[j][m]])
                if a.render == "morphology":
                    draw_morphology(ax, ids, shapes, colours)
                elif j is None:
                    ax.scatter(xy[m, 0], xy[m, 1], s=dot, c=colours, lw=0)
                else:
                    pos = POS[j][m]
                    ax.scatter(xy[m][~pos, 0], xy[m][~pos, 1], s=dot, color=NEGCOL, lw=0)
                    ax.scatter(xy[m][pos, 0], xy[m][pos, 1], s=dot, color=POSCOL, lw=0,
                               zorder=3)
                for v2, q2 in POLY[C["s"]].items():
                    if (q2[:, 0].max() < lo[0] or q2[:, 0].min() > hi[0] or
                            q2[:, 1].max() < lo[1] or q2[:, 1].min() > hi[1]):
                        continue
                    subj = v2 == C["vid"]
                    ax.add_patch(Polygon(q2, closed=True, fc="none",
                                         ec="#00e0ff" if subj else "#9e9e9e",
                                         lw=1.8 if subj else .7, zorder=5))
                ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])
                ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
                if col == 0:
                    ax.set_ylabel("cell type" if j is None else f"m{j} {ANN[j]}",
                                  fontsize=6.5,
                                  color="#b2182b" if j == k else "#1a1a1a")
            # every column has its own scale, so every column gets its own bar
            ax = axes[-1, col]
            L = nice_bar(fov)
            x0, y0 = lo[0] + fov * .06, lo[1] + fov * .06
            ax.plot([x0, x0 + L], [y0, y0], color="#111", lw=1.6, zorder=6)
            ax.text(x0, y0 + fov * .03, f"{L:.0f} " + r"$\mu$m", fontsize=6, zorder=6)
            fr = F.iloc[C["v"]]
            oth = max(fr[f"m{j}"] for j in others(k))
            tail = "" if k == UBIQ else f" | m23 {fr[f'm{UBIQ}']:.2f}"
            axes[0, col].set_title(
                f"{C['vid']}\n"
                f"m{k} {fr[f'm{k}']:.2f} | max other {oth:.2f}{tail}\n"
                f"{V.lumen_area.iloc[C['v']]:,.0f} " + r"$\mu$m$^2$"
                + f"   ·   {fov:.0f} " + r"$\mu$m field",
                fontsize=5.8, pad=3, linespacing=1.35)
        from matplotlib.lines import Line2D
        hs = [Line2D([], [], ls="", marker="o", ms=3.4, color=CT[ct_cats[t]],
                     label=f"{ct_cats[t]} {share[t]:.0%}") for t in keep]
        rest = [t for t in range(len(ct_cats)) if t not in keep and cnt[t]]
        if rest:                       # listed with no swatch: a grey one would collide
            hs.append(Line2D([], [], ls="", marker="", label=f"+{len(rest)} types <2% "
                             f"({share[rest].sum():.0%} of cells)"))
        hs += [Line2D([], [], ls="", marker="o", ms=3.4, color=POSCOL, label="motif positive"),
               Line2D([], [], ls="", marker="o", ms=3.4, color=NEGCOL, label="negative")]
        fig.legend(handles=hs, frameon=False, ncol=6, loc="upper center", fontsize=5.8,
                   handletextpad=.2, columnspacing=.9,
                   bbox_to_anchor=(.5, -(4 * MM) / fig.get_figheight()),
                   bbox_transform=fig.transFigure)
        fig.suptitle(f"m{k} {ANN[k]}", fontsize=8, x=.02, ha="left")
        save(fig, f"m{k}_specific_examples")

    print("\n" + T.to_string(index=False))
    log(f"wrote to {OUT}/")


if __name__ == "__main__":
    main()
