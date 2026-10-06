#!/usr/bin/env python
"""Pick the H&E-overlay cell-type palette by distance from the real stain colours.

Rather than excluding a hue band, this measures the stain itself and keeps every candidate
colour far from it in OKLab:

1. sample tissue pixels from the H&E of every crop the figure draws, balanced across the four
   section images, dropping near-white glass and very dark pixels (LIGHT_RANGE);
2. k-means (K clusters) on those samples -- drawn as swatches so the pink / purple / red
   stain clusters are visible;
3. candidate grid in OKLCH (L, C, H below), out-of-gamut candidates dropped;
4. score each candidate by its OKLab distance to its NN-th nearest sampled pixel (KD-tree);
   the NN-th rather than the 1st so a few stray pixels cannot dominate. Candidates scoring
   above THRESH survive;
5. survivors plotted as hue vs lightness, which shows which hue ranges are usable.

Steps 1-4 only decide WHICH HUES ARE ALLOWED. The colours themselves follow a fixed rule:
chroma is CHROMA for everything (hue does the work), the 19 cell types are laid evenly along
the allowed hues in group order -- within-group step STEP, between-group GAP_FACTOR x STEP --
and lightness cycles LIGHTNESS inside each group so neighbours differ. Each colour is then
checked against the stain (>= THRESH); a failure is repaired by moving its LIGHTNESS ONLY,
never its hue.

Writes he_overlay_palette.py (one dict, hand-editable), plus he_stain_clusters,
he_palette_candidates, he_palette_swatches and he_palette.csv.

Usage:
    export PYTHONNOUSERSITE=1
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/research/luad_revision/make_he_palette.py [--thresh 0.12]
"""
from __future__ import annotations

import argparse
import itertools
import json
import logging
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from coloraide import Color
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_state_palette import srgb_to_oklab, oklab_to_srgb, _M1, _M2
from motif_vessel_alignment import load_polys, log

logging.getLogger("fontTools").setLevel(logging.WARNING)   # subsetting chatter on svg save

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = "/data1/tanseyw/projects/fanj2/alarmist_luad/figures_motif_specific_vessels"
MM = 1 / 25.4

GROUPS = [
    ("epithelial", ["Tumor_epi", "Alveolar_epi", "Airway_epi"]),
    ("vascular/stromal", ["Vas_Endo", "Lym_Endo", "Pericyte", "SMC", "Fibro"]),
    ("lymphoid", ["T", "NK", "B", "Plasma"]),
    ("myeloid", ["pDC", "cDC", "Langhans_cell", "Macro", "Myeloid_other", "Neutro", "Mast"]),
]
N_SAMPLE = 50_000
CHROMA_FRAC = 0.90          # chroma per colour: this fraction of the max in-gamut C at its L,h
LIGHTNESS = [0.55, 0.75, 0.65]
YELLOW_BAND = (70.0, 125.0)         # yellow/olive hues go brown when dark...
YELLOW_LIGHTNESS = [0.72, 0.90, 0.81]   # ...so the same alternation, shifted up, applies there
L_SCAN = np.arange(0.40, 0.9001, 0.02)  # lightnesses tried when deciding if a hue is usable
GAP_FACTOR = 2.0            # between-group hue gap, in units of the within-group step
LIGHT_RANGE = (0.25, 0.92)      # OKLab L kept: drops near-white glass and very dark nuclei
K = 8
L_GRID = np.arange(0.35, 0.8501, 0.05)
C_GRID = np.arange(0.10, 0.2201, 0.02)
H_GRID = np.arange(0, 359, 2.0)
NN = 50
THRESH = 0.12
SEED = 0

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{OUT}/{stem}.{ext}", dpi=400 if ext == "png" else None)
    plt.close(fig)
    print(f"  {stem}.{{png,pdf,svg}}")


def hexof(lab):
    return Color("oklab", list(lab)).convert("srgb").fit().to_string(hex=True)


def sample_pixels(n=N_SAMPLE, seed=SEED):
    """tissue pixels in OKLab, balanced across the section images"""
    from plot_vessel_he_morphology import he_crop, FOV_MIN, FOV_MAX, CONTEXT_UM
    E = pd.read_csv(f"{OUT}/example_vessels.csv")
    POLY = {s: {f"{s}__{p['vid']}": p["poly"] for p in load_polys(s)} for s in set(E.section)}
    rng = np.random.default_rng(seed)
    secs = sorted(E.section.unique())
    per_section = n // len(secs)
    out = []
    for s in secs:
        rows = E[E.section == s]
        per_crop = per_section // len(rows)
        got = []
        for _, r in rows.iterrows():
            poly = POLY[s][r.vessel_id]
            fov = float(np.clip(np.ptp(poly, 0).max() + 2 * CONTEXT_UM, FOV_MIN, FOV_MAX))
            c = poly.mean(0)
            img, _, _ = he_crop(s, c - fov / 2, c + fov / 2)
            P = img.reshape(-1, 3).astype(float) / 255
            lab = srgb_to_oklab(P)
            keep = (lab[:, 0] > LIGHT_RANGE[0]) & (lab[:, 0] < LIGHT_RANGE[1])
            lab = lab[keep]
            take = rng.choice(len(lab), min(per_crop, len(lab)), replace=False)
            got.append(lab[take])
        got = np.vstack(got)
        out.append(got)
        log(f"  {s}: {len(got):,} pixels from {len(rows)} crops")
    return np.vstack(out)


def stain_clusters(lab, k=K, seed=SEED):
    km = KMeans(n_clusters=k, random_state=seed, n_init=10).fit(lab)
    order = np.argsort(-np.bincount(km.labels_, minlength=k))
    cen, size = km.cluster_centers_[order], np.bincount(km.labels_, minlength=k)[order]
    fig, ax = plt.subplots(figsize=(120 * MM, 34 * MM))
    for i, (c, n) in enumerate(zip(cen, size)):
        h = np.degrees(np.arctan2(c[2], c[1])) % 360
        ax.add_patch(plt.Rectangle((i, 0), .86, 1, facecolor=hexof(c), edgecolor="#4d4d4d",
                                   lw=.3))
        ax.text(i + .43, -.12, f"{hexof(c)}\nL {c[0]:.2f}  h {h:.0f}°\n{n / len(lab):.0%}",
                ha="center", va="top", fontsize=5.5, linespacing=1.3)
    ax.set_xlim(-.2, k); ax.set_ylim(-1.5, 1.15)
    ax.axis("off")
    ax.set_title(f"H&E tissue colours, k-means k={k} on {len(lab):,} sampled pixels",
                 loc="left", pad=4)
    save(fig, "he_stain_clusters")
    return cen, size


def candidates():
    """OKLCH grid, in-gamut only; returns (lch, lab)"""
    lch = np.array([[L, C, h] for L in L_GRID for C in C_GRID for h in H_GRID])
    lab = np.column_stack([lch[:, 0], lch[:, 1] * np.cos(np.deg2rad(lch[:, 2])),
                           lch[:, 1] * np.sin(np.deg2rad(lch[:, 2]))])
    lin = ((lab @ np.linalg.inv(_M2).T) ** 3) @ np.linalg.inv(_M1).T
    ok = (lin > -1e-4).all(1) & (lin < 1 + 1e-4).all(1)
    return lch[ok], lab[ok]


def in_gamut(lab):
    lin = ((np.asarray(lab, float) @ np.linalg.inv(_M2).T) ** 3) @ np.linalg.inv(_M1).T
    return bool(np.all(lin > -1e-4) and np.all(lin < 1 + 1e-4))


def max_chroma(L, h, hi=0.45, iters=24):
    """largest in-gamut chroma at this lightness and hue"""
    t = np.deg2rad(h)
    lo = 0.0
    for _ in range(iters):
        mid = (lo + hi) / 2
        if in_gamut([L, mid * np.cos(t), mid * np.sin(t)]):
            lo = mid
        else:
            hi = mid
    return lo


def lch_at(L, h):
    """the colour this rule gives at (L, h): chroma = CHROMA_FRAC of the in-gamut maximum"""
    return [float(L), float(CHROMA_FRAC * max_chroma(L, h)), float(h)]


def lightness_cycle(h):
    return YELLOW_LIGHTNESS if YELLOW_BAND[0] <= h <= YELLOW_BAND[1] else LIGHTNESS


def hue_scan(tree, thresh, nn=NN):
    """per hue: the best stain distance any lightness can reach under the chroma rule"""
    rows = []
    for h in H_GRID:
        for L in L_SCAN:
            lch = lch_at(L, h)
            if lch[1] < 0.04:                      # no usable chroma here
                continue
            rows.append((h, L, lch[1], stain_distance(lch, tree, nn)))
    R = np.array(rows)
    allowed = np.zeros(len(H_GRID), bool)
    for i, h in enumerate(H_GRID):
        m = R[:, 0] == h
        allowed[i] = m.any() and R[m, 3].max() >= thresh
    return R, allowed


def assign(allowed, thresh):
    """19 colours evenly along the allowed hues, group order, lightness alternating"""
    idx = np.flatnonzero(~allowed)
    gaps = [(len(H_GRID) + idx[(i + 1) % len(idx)] - idx[i]) % len(H_GRID) for i in range(len(idx))]
    start = (idx[int(np.argmax(gaps))] + 1) % len(H_GRID)
    ring = [H_GRID[(start + i) % len(H_GRID)] for i in range(len(H_GRID))
            if allowed[(start + i) % len(H_GRID)]]
    n_in = sum(len(g) - 1 for _, g in GROUPS)
    step = (len(ring) - 1) / (n_in + GAP_FACTOR * (len(GROUPS) - 1))
    log(f"  {len(ring)} allowed hue steps ({ring[0]:.0f}-{ring[-1]:.0f} deg walking the ring)"
        f" | within-group {step * 2:.1f} deg, between-group {step * 2 * GAP_FACTOR:.1f} deg")
    out, pos, li = {}, 0.0, 0
    for gi, (gname, members) in enumerate(GROUPS):
        for j, ct in enumerate(members):
            h = ring[int(round(pos))]
            out[ct] = lch_at(lightness_cycle(h)[li % len(LIGHTNESS)], h)
            li += 1
            if j < len(members) - 1:
                pos += step
        pos += GAP_FACTOR * step
    return out, ring


def stain_distance(lch, tree, nn=NN):
    """distance from the gamut-mapped colour to its nn-th nearest sampled H&E pixel"""
    rgb = Color("oklch", list(lch)).convert("srgb").fit()
    lab = srgb_to_oklab(np.array(rgb[:3], float))
    return float(tree.query(lab, k=nn)[0][-1])


def repair_lightness(out, tree, thresh):
    """a colour too close to the stain moves in LIGHTNESS ONLY (chroma follows the rule);
    among passing lightnesses take the one furthest from the rest of the palette"""
    for ct, lch in out.items():
        d = stain_distance(lch, tree)
        if d >= thresh:
            continue
        h = lch[2]
        floor = 0.70 if YELLOW_BAND[0] <= h <= YELLOW_BAND[1] else 0.35
        grid = np.arange(floor, 0.9001, 0.01)
        olab = np.array([Color("oklch", list(v)).convert("oklab")[:3]
                         for k, v in out.items() if k != ct], float)
        best = None
        for L in grid:
            cand = lch_at(L, h)
            dd = stain_distance(cand, tree)
            if dd < thresh:
                continue
            lab = np.array(Color("oklch", cand).convert("oklab")[:3], float)
            key = (round(float(np.linalg.norm(olab - lab, axis=1).min()), 3), -abs(L - lch[0]))
            if best is None or key > best[0]:
                best = (key, cand, dd)
        if best is None:
            log(f"  repair: {ct} at {d:.3f} -- no lightness clears {thresh}, left as is")
            continue
        _, cand, dd = best
        log(f"  repair: {ct} at {d:.3f} from the stain -> lightness {lch[0]:.2f} to "
            f"{cand[0]:.2f}, chroma {lch[1]:.3f} to {cand[1]:.3f} ({dd:.3f})")
        out[ct] = cand
    return out


def candidate_plot(R, out, thresh):
    ok = R[:, 3] >= thresh
    fig, ax = plt.subplots(figsize=(150 * MM, 52 * MM))
    ax.scatter(R[ok, 0], R[ok, 1], s=2, lw=0,
               c=[hexof(Color("oklch", [L, C, h]).convert("oklab")[:3])
                  for h, L, C, _ in R[ok]])
    for ct, l in out.items():
        ax.plot(l[2], l[0], "o", ms=5, mfc="none", mec="#1a1a1a", mew=.7)
        ax.annotate(ct, (l[2], l[0]), xytext=(0, 5), textcoords="offset points", fontsize=4.6,
                    ha="center", rotation=90)
    ax.set(xlim=(0, 360), ylim=(0.35, 0.95), xlabel="OKLCH hue (deg)", ylabel="lightness",
           xticks=range(0, 361, 30))
    ax.set_title(f"reachable at {CHROMA_FRAC:.0%} of max chroma and at least {thresh} from the "
                 f"H&E pixels ({ok.sum():,} of {len(R):,}); circles = chosen", loc="left", pad=4)
    save(fig, "he_palette_candidates")


def swatches(hexes, bg, previous=None, stem="he_palette_swatches"):
    rows = [("new", hexes)] + ([("previous", previous)] if previous else [])
    fig, axes = plt.subplots(2, 1, figsize=(170 * MM, 44 + 28 * len(rows)))
    fig.set_size_inches(170 * MM, (34 + 26 * len(rows)) * MM * 2 / 2)
    for ax, back, lab in zip(axes, (bg, "#ffffff"), (f"mean H&E tissue {bg}", "white")):
        ax.set_facecolor(back)
        fg = "#1a1a1a" if back == "#ffffff" else "#ffffff"
        for ri, (rname, pal) in enumerate(rows):
            y = -ri * 1.15
            x = 0
            for gname, members in GROUPS:
                for ct in members:
                    ax.plot(x, y, "o", ms=10, color=pal[ct], mec="none")
                    if ri == len(rows) - 1:
                        ax.text(x, y - .55, ct, rotation=90, ha="center", va="top",
                                fontsize=5.5, color=fg)
                    x += 1.35
                if ri == 0:
                    ax.text(x - 1.35 * len(members) / 2 - .35, y + .5, gname, ha="center",
                            va="bottom", fontsize=6, color=fg)
                x += 1.5
            ax.text(-1.0, y, rname, ha="right", va="center", fontsize=6, color=fg)
        ax.set_xlim(-5.0, x); ax.set_ylim(-1.15 * (len(rows) - 1) - 3.4, 1.15)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(lab, loc="left", fontsize=7)
    save(fig, stem)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--thresh", type=float, default=THRESH)
    ap.add_argument("--k", type=int, default=K)
    ap.add_argument("--nn", type=int, default=NN)
    ap.add_argument("--dry-run", action="store_true", help="measure only; write nothing")
    a = ap.parse_args()

    he = sample_pixels()
    log(f"{len(he):,} tissue pixels sampled (L kept in {LIGHT_RANGE})")
    stain_clusters(he, a.k)
    bg = hexof(he.mean(0))
    tree = cKDTree(he)

    R, allowed = hue_scan(tree, a.thresh, a.nn)
    log(f"hue scan at {CHROMA_FRAC:.0%} of max chroma: {allowed.sum()} of {len(H_GRID)} hues "
        f"reach {a.thresh}")
    out, ring = assign(allowed, a.thresh)
    out = repair_lightness(out, tree, a.thresh)
    hexes = {ct: Color("oklch", list(v)).convert("srgb").fit().to_string(hex=True)
             for ct, v in out.items()}
    dist = {ct: stain_distance(v, tree) for ct, v in out.items()}

    pairs = sorted((Color(hexes[x]).delta_e(hexes[y], method="ok"), x, y)
                   for x, y in itertools.combinations(hexes, 2))
    log(f"closest pair {pairs[0][0]:.3f} ({pairs[0][1]}/{pairs[0][2]}); "
        f"closest to the stain {min(dist.values()):.3f}; "
        f"chroma {min(v[1] for v in out.values()):.3f}-{max(v[1] for v in out.values()):.3f}")
    log("all pairs below 0.08: " + (", ".join(f"{x}/{y} {d:.3f}" for d, x, y in pairs if d < .08)
                                    or "none"))
    if a.dry_run:
        return

    prev_path = f"{HERE}/he_palette_previous.json"
    previous = json.load(open(prev_path)) if os.path.exists(prev_path) else None
    candidate_plot(R, out, a.thresh)
    swatches(hexes, bg, previous)

    pd.DataFrame([dict(cell_type=ct, group=gn, hex=hexes[ct], L=round(out[ct][0], 3),
                       C=round(out[ct][1], 3), h=round(out[ct][2], 1),
                       dist_to_he=round(dist[ct], 3))
                  for gn, g in GROUPS for ct in g]).to_csv(f"{OUT}/he_palette.csv", index=False)
    with open(f"{HERE}/he_overlay_palette.py", "w") as f:
        f.write('"""Cell-type colours for the H&E overlay figure -- edit this dict to tweak.\n\n'
                "Generated by make_he_palette.py: hue set by group along the hues that can clear\n"
                f"{a.thresh} from the measured H&E stain, chroma {CHROMA_FRAC:.0%} of the sRGB "
                f"maximum at\neach colour's lightness and hue, lightness alternating "
                f"{LIGHTNESS} ({YELLOW_LIGHTNESS} over\nthe yellow/olive band "
                f"{YELLOW_BAND[0]:.0f}-{YELLOW_BAND[1]:.0f} deg, where dark reads brown).\n"
                '"""\n\n'
                f"BACKGROUND = {bg!r}   # mean sampled H&E tissue pixel\n\n"
                "HE_CELLTYPE_COLORS = {\n")
        for gn, g in GROUPS:
            f.write(f"    # {gn}\n")
            for ct in g:
                f.write(f"    {ct!r}: {hexes[ct]!r},\n")
        f.write("}\n")
    print("  wrote he_overlay_palette.py, he_palette.csv")
    print(pd.read_csv(f"{OUT}/he_palette.csv").to_string(index=False))


if __name__ == "__main__":
    main()
