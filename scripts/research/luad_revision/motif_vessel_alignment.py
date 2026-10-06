#!/usr/bin/env python
"""Which ALARMIST motif's spatial structure aligns with the curated blood-vessel annotation?

Model-free geometry against a model-independent annotation. For every cell we compute the
signed distance to the nearest curated H&E lumen boundary (negative = inside the lumen),
then ask, per section and per motif:

  * enrichment   P(motif positive | 0 < d <= 25 um) / P(motif positive | whole section)
  * AUROC        how well the continuous loading ranks "perivascular (d <= 25)" against
                 "far (d > 100)", over cells in those two strata only
  * specificity  what fraction of that motif's positive cells are perivascular
  * profile      mean loading and positive fraction in distance bands, so a sharp
                 perivascular peak can be told from a diffuse gradient
  * per-vessel   fraction of the 309 curated lumens whose 0-25 um ring is enriched

The two obvious confounds are handled explicitly rather than argued away:

  1. Endothelium is trivially at vessels. Every statistic is recomputed on NON-VASCULAR
     cells only (excluding Vas_Endo, Lym_Endo, SMC, Pericyte). A motif that only aligns
     because it is an endothelial motif loses its signal there.
  2. The annotation is sparse (71-83 lumens per section) and biased to large vessels, so
     "far from a curated lumen" is not "avascular". Enrichment is therefore reported
     against the section background, and also against a distance-matched null: random
     points drawn inside the tissue with the same distance-to-lumen distribution is not
     possible, so instead we report the tissue-layer composition of each band, which is
     what would otherwise drive a spurious result.

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/motif_vessel_alignment.py \
        --out /data1/tanseyw/projects/fanj2/alarmist_luad/motif_vessel_alignment.json
"""
from __future__ import annotations

import argparse
import json
import time

import h5py
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from scipy.spatial import cKDTree

H5 = "/data1/tanseyw/projects/fanj2/alarmist_luad/adata_with_region.h5ad"
VDIR = "/data1/tanseyw/projects/tosh/vessel_curator/figures/vessel_curation"
SECTIONS = ["P17_AIS", "P17_LUAD", "P21_AIS", "P21_LUAD"]
NM = 25
VASCULAR = ["Vas_Endo", "Lym_Endo", "SMC", "Pericyte"]

# distance bands in um from the lumen BOUNDARY; the first is inside the lumen
BANDS = [(-1e9, 0.0), (0.0, 10.0), (10.0, 25.0), (25.0, 50.0),
         (50.0, 100.0), (100.0, 200.0), (200.0, 1e9)]
BAND_NAMES = ["inside", "0-10", "10-25", "25-50", "50-100", "100-200", ">200"]
NEAR, FAR = 25.0, 100.0          # the two strata the AUROC contrasts (CLI overrides)
STEP = 1.5                       # um between densified boundary points


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def read_cat(obs, key):
    g = obs[key]
    cats = [x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]]
    return cats, g["codes"][:]


def parse_poly(v):
    """poly_world_json is sometimes None and sometimes the literal string 'null'."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    try:
        q = json.loads(v)
    except Exception:
        return None
    if not q or len(q) < 3:
        return None
    P = np.asarray(q, float)
    if P.ndim != 2 or P.shape[1] != 2 or not np.isfinite(P).all():
        return None
    if len(P) > 3 and np.allclose(P[0], P[-1]):
        P = P[:-1]
    return P if len(P) >= 3 else None


def load_polys(sec):
    """Curated lumens in world microns. Joined on the vessel_id suffix, never row order."""
    L = pd.read_parquet(f"{VDIR}/{sec}_lumens.parquet")
    V = pd.read_parquet(f"{VDIR}/{sec}_vessels.parquet")
    vidx = L["vessel_id"].str.rsplit("__", n=1).str[-1].astype(int).to_numpy()
    assert np.allclose(V.loc[vidx, "cx_um"].to_numpy(), L["cx_um"].to_numpy(), atol=1e-6)
    out = []
    for r in range(len(L)):
        P = parse_poly(L["poly_world_json"].iloc[r])
        src = 0
        if P is None:
            vr = V.loc[vidx[r]]
            P = np.array([[vr[f"corner{k}_x_um"], vr[f"corner{k}_y_um"]] for k in range(4)],
                         float)
            src = 1
        out.append(dict(poly=P, src=src, flagged=bool(L["flagged"].iloc[r]),
                        area=float(L["lumen_area_um2"].iloc[r])
                        if pd.notna(L["lumen_area_um2"].iloc[r]) else None,
                        vid=int(vidx[r])))
    return out


def densify(P, step=STEP):
    """Points along the closed ring, roughly `step` um apart."""
    Q = np.vstack([P, P[:1]])
    seg = np.diff(Q, axis=0)
    L = np.hypot(seg[:, 0], seg[:, 1])
    pts = []
    for i, li in enumerate(L):
        n = max(int(np.ceil(li / step)), 1)
        t = np.arange(n)[:, None] / n
        pts.append(Q[i] + t * seg[i])
    return np.vstack(pts)


def signed_distance(xy, polys):
    """Distance to the nearest lumen boundary; negative for cells inside a lumen."""
    bpts = np.vstack([densify(p["poly"]) for p in polys])
    owner = np.concatenate([np.full(len(densify(p["poly"])), i)
                            for i, p in enumerate(polys)])
    tree = cKDTree(bpts)
    d, j = tree.query(xy, k=1)
    nearest = owner[j]
    inside = np.zeros(len(xy), bool)
    cell_tree = cKDTree(xy)
    for i, p in enumerate(polys):
        P = p["poly"]
        lo, hi = P.min(0), P.max(0)
        c, r = (lo + hi) / 2, np.hypot(*(hi - lo)) / 2 + 1.0
        cand = cell_tree.query_ball_point(c, r)
        if not cand:
            continue
        cand = np.asarray(cand)
        m = MplPath(P).contains_points(xy[cand])
        inside[cand[m]] = True
    return np.where(inside, -d, d), nearest, inside


def auroc(score, label):
    """Rank-based AUROC with tie handling. label is boolean."""
    n1 = int(label.sum())
    n0 = int((~label).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = pd.Series(score).rank(method="average").to_numpy()
    return float((r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main():
    global NEAR, FAR
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--near", type=float, default=NEAR)
    ap.add_argument("--far", type=float, default=FAR)
    a = ap.parse_args()
    NEAR, FAR = a.near, a.far
    assert FAR > NEAR
    if a.out is None:
        a.out = ("/data1/tanseyw/projects/fanj2/alarmist_luad/"
                 f"motif_vessel_alignment_ring{int(NEAR)}.json")

    f = h5py.File(H5, "r")
    obs = f["obs"]
    sec_cats, sec_codes = read_cat(obs, "sample_id")
    ct_cats, ct_codes = read_cat(obs, "cell_type")
    tl_cats, tl_codes = read_cat(obs, "tissue_layer")
    xy = f["obsm"]["spatial"][:]
    N = len(sec_codes)
    log(f"{N:,} cells")

    vasc_codes = {ct_cats.index(t) for t in VASCULAR if t in ct_cats}
    is_vasc = np.isin(ct_codes, list(vasc_codes))
    log(f"vascular cell types {VASCULAR} -> {is_vasc.sum():,} cells "
        f"({100*is_vasc.mean():.1f}%)")

    load = np.empty((NM, N), np.float32)
    pos = np.empty((NM, N), bool)
    for k in range(NM):
        load[k] = obs[f"motif_{k}_loading"][:]
        c, cc = read_cat(obs, f"motif_{k}_state")
        pos[k] = cc == c.index("positive")
    f.close()
    log("loadings + states read")

    res = {"sections": {}, "params": dict(bands=BAND_NAMES, near=NEAR, far=FAR,
                                          step=STEP, vascular=VASCULAR)}
    for si, sec in enumerate(sec_cats):
        m = sec_codes == si
        P = xy[m]
        polys = load_polys(sec)
        t0 = time.time()
        d, nearest, inside = signed_distance(P, polys)
        log(f"{sec}: {m.sum():,} cells, {len(polys)} lumens, "
            f"distances in {time.time()-t0:.1f}s  "
            f"(inside {inside.sum():,}; median d {np.median(d):.0f} um)")

        band_idx = np.full(len(P), -1, np.int8)
        for bi, (lo, hi) in enumerate(BANDS):
            band_idx[(d > lo) & (d <= hi)] = bi
        near = (d > 0) & (d <= NEAR)
        far = d > FAR

        sub = {}
        sub["n"] = int(m.sum())
        sub["n_lumens"] = len(polys)
        sub["n_inside"] = int(inside.sum())
        sub["band_n"] = {BAND_NAMES[bi]: int((band_idx == bi).sum()) for bi in range(len(BANDS))}
        sub["band_vasc_frac"] = {BAND_NAMES[bi]: float(is_vasc[m][band_idx == bi].mean())
                                 if (band_idx == bi).any() else None
                                 for bi in range(len(BANDS))}
        # tissue-layer composition of each band -- the confound to look at first
        tl_m = tl_codes[m]
        sub["band_tissue"] = {
            BAND_NAMES[bi]: {tl_cats[c]: int((tl_m[band_idx == bi] == c).sum())
                             for c in range(len(tl_cats))}
            for bi in range(len(BANDS))}

        nv = ~is_vasc[m]                      # non-vascular cells of this section
        motifs = []
        for k in range(NM):
            pk = pos[k][m]
            lk = load[k][m]
            base = float(pk.mean())
            rec = dict(k=k, frac_pos=base)
            for tag, keep in (("all", np.ones(len(P), bool)), ("nonvasc", nv)):
                nearK = near & keep
                farK = far & keep
                bk = pk[keep]
                b = float(bk.mean()) if keep.any() else float("nan")
                e25 = float(pk[nearK].mean() / b) if nearK.any() and b > 0 else float("nan")
                near10 = (d > 0) & (d <= NEAR / 2.5) & keep
                e10 = float(pk[near10].mean() / b) if near10.any() and b > 0 else float("nan")
                sel = nearK | farK
                au = auroc(lk[sel], nearK[sel]) if sel.any() else float("nan")
                spec = float((pk & nearK).sum() / max(int(pk[keep].sum()), 1))
                rec[tag] = dict(
                    base=b, enr10=e10, enr25=e25, auroc=au, spec25=spec,
                    n_near=int(nearK.sum()), n_far=int(farK.sum()),
                    band_pos={BAND_NAMES[bi]: (float(pk[(band_idx == bi) & keep].mean())
                                               if ((band_idx == bi) & keep).any() else None)
                              for bi in range(len(BANDS))},
                    band_load={BAND_NAMES[bi]: (float(lk[(band_idx == bi) & keep].mean())
                                                if ((band_idx == bi) & keep).any() else None)
                               for bi in range(len(BANDS))})
            # per-vessel: is this lumen's 0-25 um ring enriched over the section base?
            hits = 0, 0
            ring_fracs = []
            for vi in range(len(polys)):
                r = near & (nearest == vi)
                if r.sum() < 20:
                    continue
                ring_fracs.append(float(pk[r].mean()))
            ring_fracs = np.asarray(ring_fracs)
            rec["n_vessels_scored"] = int(len(ring_fracs))
            rec["vessel_ring_median"] = float(np.median(ring_fracs)) if len(ring_fracs) else None
            rec["vessel_frac_enriched"] = (float((ring_fracs > base).mean())
                                           if len(ring_fracs) else None)
            motifs.append(rec)
        sub["motifs"] = motifs
        res["sections"][sec] = sub

    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    log(f"wrote {a.out}")

    # ---- console summary: rank motifs by non-vascular AUROC, per section --------------
    for sec in sec_cats:
        s = res["sections"][sec]
        rows = sorted(s["motifs"], key=lambda r: -(r["nonvasc"]["auroc"]
                                                   if r["nonvasc"]["auroc"] == r["nonvasc"]["auroc"]
                                                   else 0))
        print(f"\n=== {sec}  n={s['n']:,}  lumens={s['n_lumens']}  "
              f"near(<=25um)={sum(s['band_n'][b] for b in ('0-10','10-25')):,} ===")
        print(f"{'k':>3} {'%pos':>6} | {'AUROC':>6} {'enr10':>6} {'enr25':>6} {'spec25':>7} "
              f"{'vesEnr':>7} | {'AUROC':>6} {'enr25':>6}   (left = non-vascular cells only, "
              f"right = all cells)")
        for r in rows[:8]:
            nvr, alr = r["nonvasc"], r["all"]
            print(f"{r['k']:>3} {100*r['frac_pos']:>5.1f}% | {nvr['auroc']:>6.3f} "
                  f"{nvr['enr10']:>6.2f} {nvr['enr25']:>6.2f} {100*nvr['spec25']:>6.1f}% "
                  f"{100*(r['vessel_frac_enriched'] or 0):>6.0f}% | {alr['auroc']:>6.3f} "
                  f"{alr['enr25']:>6.2f}")


if __name__ == "__main__":
    main()
