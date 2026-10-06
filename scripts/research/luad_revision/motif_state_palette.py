#!/usr/bin/env python
"""A shared 16-state colour map for the four vessel motifs, mixed in OKLab.

The four motif states are binary but not exclusive, so every cell (or vessel) falls in one of
2**4 = 16 states: all-negative, 4 singles, 6 doubles, 4 triples, 1 quadruple. Each motif gets
an anchor colour; a multi-positive state is the average of its anchors **in OKLab**, which
mixes perceptually rather than in sRGB; the all-negative state is a very light grey so it
falls back into the page.

State codes are bit-packed in MOTIFS order: bit j is set when MOTIFS[j] is positive, so the
code of a cell positive for MOTIFS[0] and MOTIFS[2] is 0b0101 = 5. `state_code` builds them
from a boolean matrix, `state_palette` returns {code: hex}, `draw_state_legend` draws the
16-swatch membership legend. Import these anywhere the same colours are needed (the section
maps and the UpSet plot both do).

OKLab conversion follows Bjorn Ottosson's matrices (2020); round trip is exact to < 1e-6.
"""
from __future__ import annotations

import numpy as np
from matplotlib.colors import to_hex, to_rgb
from matplotlib.patches import Rectangle

MOTIFS = [1, 23, 24, 10]
NAME = {1: "SMC vascular stabilization motif", 23: "Vascular homeostasis motif",
        24: "Healthy alveolar motif", 10: "Tumor vasculature motif"}
SHORT = {1: "SMC stabilization", 23: "Vascular homeostasis", 24: "Healthy alveolar",
         10: "Tumor vasculature"}
# Anchors chosen for separation, not inherited. The old set put m1 (#e69f00 orange) and m24
# (#d55e00 burnt orange) in the same hue family, so their mixes collapsed into one brown band.
# These four sit far apart in OKLab hue and were found by random + local search maximising the
# closest of the 16 state colours, subject to: every state in gamut, anchor lightness inside
# validate_palette.py's light band, chroma >= 0.10, and every anchor pair clearing both the
# normal floor (>= 15) and the colour-vision floor (min(protan, deutan) >= 6).
# Result: closest pair of the 16 is 0.083 in OKLab, against 0.054 for the old anchors, and the
# four anchors themselves come back RESULT: PASS from validate_palette.py (worst normal pair
# 24.8; one WARN, red-orange vs green under deuteranopia at 6.2).
ANCHOR = {1: "#d94527", 23: "#2b86d8", 24: "#74a54b", 10: "#e98aea"}
NEG = "#ededed"             # all-negative: very light grey, meant to disappear

_M1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                [0.2119034982, 0.6806995451, 0.1073969566],
                [0.0883024619, 0.2817188376, 0.6299787005]])
_M2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                [1.9779984951, -2.4285922050, 0.4505937099],
                [0.0259040371, 0.7827717662, -0.8086757660]])


def _srgb_to_linear(c):
    c = np.asarray(c, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _linear_to_srgb(c):
    c = np.asarray(c, float)
    return np.where(c <= 0.0031308, c * 12.92, 1.055 * np.clip(c, 0, None) ** (1 / 2.4) - .055)


def srgb_to_oklab(rgb):
    lms = _srgb_to_linear(rgb) @ _M1.T
    return np.cbrt(lms) @ _M2.T


def oklab_to_srgb(lab):
    lms = (np.asarray(lab, float) @ np.linalg.inv(_M2).T) ** 3
    return np.clip(_linear_to_srgb(lms @ np.linalg.inv(_M1).T), 0, 1)


def state_code(M):
    """M: n x len(MOTIFS) boolean -> integer state code per row"""
    M = np.asarray(M, bool)
    return (M * (1 << np.arange(M.shape[1]))).sum(1).astype(np.int64)


def members(code):
    """the motifs positive in this state, in MOTIFS order"""
    return [k for j, k in enumerate(MOTIFS) if code >> j & 1]


def state_label(code, short=True):
    ms = members(code)
    if not ms:
        return "none of the four"
    d = SHORT if short else NAME
    return " + ".join(d[k] for k in ms)


LSTEP = 0.13    # OKLab lightness removed per motif beyond the first (see state_palette)


def state_palette(anchor=None, neg=NEG, chroma="average", lstep=LSTEP):
    """{state code: hex}: singles keep their anchor, multi-positives are the OKLab mean.

    Two knobs on top of the plain mean, both measured with palette_separation():
      lstep   lightness taken off per extra motif. Averaging alone cannot separate all 16
              states -- the mean of two anchors can equal the mean of four (exactly so for
              symmetric anchors) -- and with these anchors it leaves two states 0.032 apart
              in OKLab, under the ~0.04 people can tell apart. lstep 0.13 lifts the closest
              pair to 0.083 with the current anchors and makes "more motifs" read as
              "darker". lstep=0 is the plain recipe.
      chroma  'restore' rescales the mixed a/b back to the mean anchor chroma so mixes stay
              saturated; measured worse here (0.031) than plain 'average', so not the default.
    """
    anchor = anchor or ANCHOR
    lab = {k: srgb_to_oklab(to_rgb(anchor[k])) for k in MOTIFS}
    out = {0: neg}
    for code in range(1, 1 << len(MOTIFS)):
        ms = members(code)
        mix = np.mean([lab[k] for k in ms], axis=0)
        if chroma == "restore" and len(ms) > 1:
            c_mix = np.hypot(*mix[1:])
            c_target = float(np.mean([np.hypot(*lab[k][1:]) for k in ms]))
            if c_mix > 1e-9:
                mix = np.r_[mix[0], mix[1:] * (c_target / c_mix)]
        mix = np.r_[mix[0] - lstep * (len(ms) - 1), mix[1:]]
        out[code] = to_hex(oklab_to_srgb(mix))
    return out


def palette_separation(pal):
    """min OKLab distance between any two state colours (higher = easier to tell apart)"""
    codes = sorted(pal)
    lab = np.array([srgb_to_oklab(to_rgb(pal[c])) for c in codes])
    d = np.linalg.norm(lab[:, None] - lab[None], axis=-1)
    np.fill_diagonal(d, np.inf)
    i, j = np.unravel_index(np.argmin(d), d.shape)
    return float(d.min()), (codes[i], codes[j]), d


def draw_state_legend(ax, pal, counts=None, name_fs=6, swatch_fs=5.5):
    """16 swatches with a membership matrix underneath, in a dedicated axes.

    counts: optional {code: n} printed under each swatch.
    """
    codes = sorted(pal, key=lambda c: (bin(c).count("1"), c))
    k = len(MOTIFS)
    for x, code in enumerate(codes):
        ax.add_patch(Rectangle((x - .38, k - .30), .76, .62, facecolor=pal[code],
                               edgecolor="#4d4d4d", lw=.3, clip_on=False))
        if counts is not None:
            v = counts.get(code, 0)
            ax.text(x, k + .40, f"{v/1000:.0f}k" if v >= 10_000 else f"{v:,}", ha="left",
                    va="bottom", fontsize=swatch_fs, color="#444444", rotation=90)
        for j in range(k):
            y = k - 1 - j - .55
            on = bool(code >> j & 1)
            ax.plot(x, y, "o", ms=4.2, color=ANCHOR[MOTIFS[j]] if on else "#dcdcdc",
                    clip_on=False)
        ys = [k - 1 - j - .55 for j in range(k) if code >> j & 1]
        if len(ys) > 1:
            ax.plot([x, x], [min(ys), max(ys)], color="#4d4d4d", lw=.9, zorder=1)
    ax.set_xlim(-.7, len(codes) - .3)
    ax.set_ylim(-.9, k + .45)
    ax.set_yticks([k - 1 - j - .55 for j in range(k)],
                  [SHORT[MOTIFS[j]] for j in range(k)], fontsize=name_fs)
    ax.set_xticks([])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
