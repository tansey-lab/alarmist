#!/usr/bin/env python
"""Figures: where each comparator ranks ALARMIST's motif-10 / motif-24 LRIs.

Reads `luad_motif_ranks_long.csv` from `luad_motif_rank_gather.py` and draws, per
motif and per level, one row per method split into an AIS violin and a LUAD violin,
with one point per (target LRI, section).  ALARMIST sits in its own top panel on its
own axis, because its points are motif ON-fractions and not ranks -- putting them on
a shared axis would assert an equivalence that does not hold.

Two x-scales are produced for each figure:

    pct     top-percentile, 100 * (1 - (rank-1)/N).  100 % = the method's top hit.
            Comparable across methods; this is the one to read.
    rank    raw rank on a log10 axis, 1 = best, drawn on the right so "better" is the
            same direction in both.  Preserves "how many places down", but each method
            has its own N, so the shapes are not comparable across rows.

Entries the method never tested (stLearn drops every multi-subunit complex; the
CellChat DB tier has no TBXAS1_TBXA2R) are drawn in a grey band OUTSIDE the axis and
excluded from the violin -- "not tested" is a different claim from "ranked last".
Entries the method tested but did not report sit at the axis floor as open circles.

Cell-type resolution differs per method and the row labels say which:
CellChat is the only one that resolves a direction; stLearn's pair is undirected;
LIANA+ inflow names the sender only.  See the gather script's docstring.

Usage
-----
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/comparators/luad_motif_rank_plot.py [--in-dir DIR] [--out-dir DIR]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/comparators
from _common.plotting import apply_publication_style, save_all_formats  # noqa: E402

HERE = Path(__file__).resolve().parent

# dataviz categorical slots 1 and 2.  Validated with the skill's own
# scripts/validate_palette.py: all six checks PASS on the all-pairs list
# (CVD dE 24.7 protan, normal-vision dE 33.6, contrast >= 3:1 vs surface).
STAGE_COLOR = {"AIS": "#2a78d6", "LUAD": "#eb6834"}
GREY, INK, MUTED, RULE = "#9a9a95", "#1a1a19", "#5c5c58", "#e0e0dc"
BAND = "#f2f2ef"

# Methods that resolve at least one cell-type axis first, then the LR-only ones.
METHOD_ORDER = ["CellChat", "stLearn", "LIANA+ (inflow)", "SpatialDM", "CytoSignal"]
CT_RESOLUTION = {
    "CellChat": "sender -> receiver",
    "stLearn": "undirected pair",
    "LIANA+ (inflow)": "sender only",
}
LEVEL_TITLE = {"lri_ct": "LRI x cell types", "lri": "LRI only"}
LEVEL_SUB = {
    "lri_ct": ("only CellChat resolves a direction; stLearn's cell-type pair is "
               "undirected by construction and LIANA+ inflow names the sender only"),
    "lri": "ligand-receptor pair, cell types collapsed -- every method contributes",
}
MOTIF_TITLE = {
    10: "Motif 10 — tumour vasculature (Tumor_epi ↔ Vas_Endo / Pericyte)",
    24: "Motif 24 — healthy alveolar vasculature (Alveolar_epi / Vas_Endo / Fibro)",
}
MARKERS = {"P17": "o", "P21": "^", "P17+P21": "s"}
CENSOR_FRAC = 0.15          # width of the out-of-axis band, as a fraction of the axis


class Scale:
    """Everything the panel needs to know about the chosen x-scale."""

    def __init__(self, xmode, ranked_vals):
        self.xmode = xmode
        if xmode == "pct":
            self.col = "pct_rank"
            self.best, self.worst = 100.0, 0.0
            self.band_outer = -CENSOR_FRAC * 100
            self.band_inner = -0.8
            self.log = False
        else:
            self.col = "rank"
            top = float(np.nanmax(ranked_vals)) if len(ranked_vals) else 10.0
            self.best, self.worst = 0.82, max(10.0, top * 1.35)
            self.band_inner = self.worst
            self.band_outer = 10 ** (np.log10(self.worst) + CENSOR_FRAC
                                     * np.log10(self.worst / self.best))
            self.log = True
        self.band_mid = (10 ** np.mean(np.log10([self.band_inner, self.band_outer]))
                         if self.log else (self.band_inner + self.band_outer) / 2)

    def apply(self, ax):
        if self.log:
            ax.set_xscale("log")
        ax.set_xlim(self.band_outer, self.best)     # "better" is always to the RIGHT
        if self.xmode == "pct":
            ax.set_xticks([0, 25, 50, 75, 100])
            ax.set_xticklabels(["0 %", "25 %", "50 %", "75 %", "100 %"])
            ax.set_xlabel("top-percentile within the method's own ranking "
                          "(100 % = its top-ranked interaction)")
        else:
            ax.set_xlabel("rank within the method's own ranking (1 = its top hit; "
                          "log scale — N differs per method, so widths are not comparable)")


def _violin(ax, y, vals, color, width=0.36, clip=None):
    vals = np.asarray([v for v in vals if np.isfinite(v)], dtype=float)
    if len(vals) == 0:
        return
    if len(vals) >= 8 and np.ptp(vals) > 0:
        parts = ax.violinplot([vals], positions=[y], orientation="horizontal",
                              widths=width * 2, showextrema=False, showmedians=False)
        for b in parts["bodies"]:
            b.set_facecolor(color)
            b.set_alpha(0.15)
            b.set_edgecolor(color)
            b.set_linewidth(0.5)
            if clip is not None:                    # keep the KDE tails inside the axis
                v = b.get_paths()[0].vertices
                v[:, 0] = np.clip(v[:, 0], min(clip), max(clip))
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    ax.plot([q1, q3], [y, y], color=color, lw=2.8, solid_capstyle="butt",
            alpha=0.95, zorder=4)
    ax.plot([med], [y], marker="|", color="white", ms=6.5, mew=1.5, zorder=6)


def _panel(ax, sub, scale):
    methods = [m for m in METHOD_ORDER if m in set(sub.method)]
    scale.apply(ax)
    clip = (min(scale.best, scale.worst), max(scale.best, scale.worst))

    yticks, ylabels, method_mid = [], [], {}
    y = 0.0
    for m in methods:
        msub = sub[sub.method == m]
        first_y = y
        for stage in ("LUAD", "AIS"):               # AIS ends up on top
            rows = msub[msub.stage == stage]
            if rows.empty:
                y += 1.0
                continue
            colour = STAGE_COLOR[stage]
            rng = np.random.default_rng(abs(hash((m, stage, scale.xmode))) % (2 ** 32))

            ranked = rows[rows.status == "ranked"]
            _violin(ax, y, ranked[scale.col].astype(float).values, colour, clip=clip)
            for pat, grp in ranked.groupby("patient", observed=True):
                v = grp[scale.col].astype(float).values
                ax.scatter(v, np.full(len(v), y) + rng.uniform(-0.18, 0.18, len(v)),
                           s=11, marker=MARKERS.get(pat, "o"), facecolor=colour,
                           edgecolor="white", linewidth=0.35, alpha=0.85, zorder=5)

            zero = rows[rows.status == "zero"]
            if len(zero):
                ax.scatter(np.full(len(zero), scale.worst),
                           np.full(len(zero), y) + rng.uniform(-0.18, 0.18, len(zero)),
                           s=14, marker="o", facecolor="none", edgecolor=colour,
                           linewidth=0.65, alpha=0.9, zorder=5, clip_on=False)

            nt = rows[rows.status == "not_tested"]
            if len(nt):
                ax.scatter(np.full(len(nt), scale.band_mid),
                           np.full(len(nt), y) + rng.uniform(-0.15, 0.15, len(nt)),
                           s=14, marker="x", color=GREY, linewidth=0.75,
                           zorder=5, clip_on=False)

            if (rows.status == "no_data").all():
                ax.text(np.mean(clip) if not scale.log
                        else 10 ** np.mean(np.log10(clip)), y,
                        "run did not finish", ha="center", va="center",
                        fontsize=5.6, color=GREY, style="italic")

            yticks.append(y)
            ylabels.append(stage)
            y += 1.0
        method_mid[m] = (first_y + y - 1.0) / 2
        y += 0.62

    ax.set_ylim(y - 0.62 - 0.55, -0.55)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.tick_params(axis="y", length=0, pad=1.5)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color=RULE, lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    ax.axvspan(scale.band_outer, scale.band_inner, color=BAND, zorder=-1, lw=0)
    ax.axvline(scale.band_inner, color="#d2d2ce", lw=0.7, zorder=1)
    ax.annotate("not\ntested", xy=(scale.band_mid, 0), xycoords=("data", "axes fraction"),
                xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
                fontsize=5.0, color=GREY, linespacing=1.05)
    return methods, method_mid


def _method_labels(ax, methods, method_mid, sub, level):
    for m in methods:
        msub = sub[sub.method == m]
        n_rank, n_slot = int((msub.status == "ranked").sum()), len(msub)
        tag = CT_RESOLUTION.get(m, "") if level == "lri_ct" else ""
        ax.annotate(m, xy=(0, method_mid[m]), xycoords=("axes fraction", "data"),
                    xytext=(-31, 5.5 if tag else 3.5), textcoords="offset points",
                    ha="right", va="center", fontsize=7.2, color=INK)
        ax.annotate(f"{n_rank}/{n_slot} ranked", xy=(0, method_mid[m]),
                    xycoords=("axes fraction", "data"),
                    xytext=(-31, -3.0 if tag else -5.0), textcoords="offset points",
                    ha="right", va="center", fontsize=5.2, color=MUTED)
        if tag:
            ax.annotate(tag, xy=(0, method_mid[m]), xycoords=("axes fraction", "data"),
                        xytext=(-31, -10.5), textcoords="offset points",
                        ha="right", va="center", fontsize=5.0, color=GREY, style="italic")


def _alarmist_panel(ax, onfrac, motif):
    rows = onfrac[onfrac.motif == motif].dropna(subset=["on_fraction"])
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0 %", "25 %", "50 %", "75 %", "100 %"])
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("ALARMIST motif ON-fraction of cells — its own axis, NOT a rank",
                  fontsize=6.2, labelpad=3)
    for y, stage in enumerate(("LUAD", "AIS")):
        r = rows[rows.stage == stage].sort_values("on_fraction")
        if not len(r):
            ax.text(50, y, f"motif-{motif} ON-fractions not yet supplied",
                    ha="center", va="center", fontsize=5.8, color=GREY, style="italic")
            continue
        if len(r) > 1:
            ax.plot([r.on_fraction.min(), r.on_fraction.max()], [y, y],
                    color=STAGE_COLOR[stage], lw=1.4, alpha=0.45, zorder=3)
        for i, (_, s) in enumerate(r.iterrows()):
            ax.scatter([s.on_fraction], [y], s=26, marker=MARKERS.get(s.patient, "o"),
                       facecolor=STAGE_COLOR[stage], edgecolor="white",
                       linewidth=0.5, zorder=5)
            ax.annotate(f"{s.patient} {s.on_fraction:.1f}%", (s.on_fraction, y),
                        xytext=(0, 7 if i % 2 == 0 else -11), textcoords="offset points",
                        ha="center", fontsize=5.2, color=MUTED)
    ax.set_ylim(1.75, -0.75)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["LUAD", "AIS"])
    ax.tick_params(axis="y", length=0, pad=1.5)
    for s in ("bottom", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color=RULE, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.annotate("ALARMIST", xy=(0, 0.5), xycoords=("axes fraction", "data"),
                xytext=(-31, 0), textcoords="offset points", ha="right", va="center",
                fontsize=7.6, color=INK, weight="bold")


def make_figure(df, onfrac, motif, level, xmode, out_stem):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    sub = df[(df.motif == motif) & (df.level == level)]
    n_methods = len([m for m in METHOD_ORDER if m in set(sub.method)])

    # Inch-exact layout: gridspec's hspace is relative to mean axes height, which
    # opens a large hole between a 0.7 in panel and a 3.6 in one.
    HEAD, A_H, GAP, FOOT = 1.05, 0.70, 0.34, 0.92
    body_h = 0.72 * n_methods
    h = HEAD + A_H + GAP + body_h + FOOT
    w, L, R = 7.4, 0.235, 0.985
    fig = plt.figure(figsize=(w, h))
    ax_a = fig.add_axes([L, (FOOT + body_h + GAP) / h, R - L, A_H / h])
    ax_m = fig.add_axes([L, FOOT / h, R - L, body_h / h])

    _alarmist_panel(ax_a, onfrac, motif)
    scale = Scale(xmode, sub.loc[sub.status == "ranked", "rank"].astype(float).values)
    methods, method_mid = _panel(ax_m, sub, scale)
    _method_labels(ax_m, methods, method_mid, sub, level)

    fig.text(0.010, 1 - 0.22 / h, MOTIF_TITLE[motif], ha="left", va="baseline",
             fontsize=8.8, weight="bold", color=INK)
    fig.text(0.010, 1 - 0.40 / h,
             f"Top {sub.alarmist_rank.nunique()} LRIs by factor_lrnorm, located in each "
             f"comparator's own ranking   |   level: {LEVEL_TITLE[level]}",
             ha="left", va="baseline", fontsize=6.6, color=MUTED)
    fig.text(0.010, 1 - 0.55 / h, LEVEL_SUB[level], ha="left", va="baseline",
             fontsize=5.8, color=GREY)

    handles = [
        Patch(facecolor=STAGE_COLOR["AIS"], alpha=0.35,
              edgecolor=STAGE_COLOR["AIS"], label="AIS (precursor)"),
        Patch(facecolor=STAGE_COLOR["LUAD"], alpha=0.35,
              edgecolor=STAGE_COLOR["LUAD"], label="LUAD (invasive)"),
        Line2D([], [], marker="o", ls="none", mfc=MUTED, mec="white", ms=4,
               label="patient P17"),
        Line2D([], [], marker="^", ls="none", mfc=MUTED, mec="white", ms=4,
               label="patient P21"),
        Line2D([], [], marker="s", ls="none", mfc=MUTED, mec="white", ms=4,
               label="P17+P21 pooled (CellChat only)"),
        Line2D([], [], marker="o", ls="none", mfc="none", mec=MUTED, ms=4.5,
               label="tested, not reported"),
        Line2D([], [], marker="x", ls="none", color=GREY, ms=4.5,
               label="not in the method's LR universe"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=5.8,
               frameon=False, bbox_to_anchor=(0.56, 0.004), handletextpad=0.5,
               columnspacing=1.5)
    save_all_formats(fig, out_stem, dpi=450, close=True, verbose=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    default = "/data1/tanseyw/projects/fanj2/results/_figures/luad_motif_ranks"
    ap.add_argument("--in-dir", default=default)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--onfrac", default=str(HERE / "alarmist_motif_onfrac_luad.csv"))
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir or in_dir / "figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_dir / "luad_motif_ranks_long.csv")
    onfrac = pd.read_csv(args.onfrac)

    apply_publication_style(**{
        "font.size": 6.4, "axes.labelsize": 6.4, "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.2, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.2, "xtick.minor.size": 1.2,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "axes.edgecolor": "#cfcfca", "text.color": INK,
        "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": INK,
    })

    for motif in (10, 24):
        for level in ("lri_ct", "lri"):
            for xmode in ("pct", "rank"):
                make_figure(df, onfrac, motif, level, xmode,
                            out_dir / f"motif{motif}_{level}_{xmode}")
    print(f"\n[plot] figures in {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
