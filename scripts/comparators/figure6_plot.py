#!/usr/bin/env python
"""Figure 6 — the CCC comparator benchmark, and its density-collinearity supplement.

Reads the tidy CSVs written by figure6_gather.py; computes nothing.

Fig 6   a  every comparator detects both arms of motif 1
        b  ... but in tested sets whose sizes differ 30-fold
        c  no comparator puts the two arms in one data-derived object
        d  why: at the cell they barely co-occur, at the 50 um patch they do
        e  and only two comparators ship a between-condition test at the punch

Supp    grade and cellularity are collinear in this 13-core TMA

Env: /Users/jiayifan/anaconda3/envs/bptf/bin/python
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common.plotting import apply_publication_style
from _common.plotting import save_all_formats as _save_all_formats

ROOT = Path("/Users/jiayifan/tansey_lab/alarmist")
DATA = ROOT / "results/comparators/_benchmark/figure6"
FIG = ROOT / "results/comparators/_benchmark/figure6"
FIG.mkdir(parents=True, exist_ok=True)

U = 1.5  # the layout unit, in inches

# Okabe-Ito, safe under the common colour-vision deficiencies.
INK = "#000000"
GREY = "#4A4A4A"
ACCENT = "#D55E00"   # ALARMIST
BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
LIGHT = "#DDDDDD"

OBJ_COLOR = {"flat": "#FFFFFF", "curated": SKY, "learned": GREEN}
OBJ_LABEL = {"flat": "no grouping object exists",
             "curated": "grouping is a database column",
             "learned": "grouping estimated from data"}


def set_style() -> None:
    """Figure-6 family style: the CLAUDE.md baseline plus this figure's 6 pt scale.

    Re-exported for `figure6_supp_density.py` and `figure7_plot_spatial.py`,
    which import it from here.
    """
    apply_publication_style(**{
        "font.size": 6,
        "axes.labelsize": 6,
        "axes.titlesize": 6,
        "xtick.labelsize": 5.5,
        "ytick.labelsize": 5.5,
        "legend.fontsize": 5.5,
        "axes.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "legend.frameon": False,
    })


def save_all_formats(fig: plt.Figure, stem: Path) -> None:
    """One saver, three formats, at this figure family's 450 dpi.

    Re-exported for `figure6_supp_density.py` and `figure7_plot_spatial.py`.
    """
    _save_all_formats(fig, stem, dpi=450, verbose=True)


def panel_letter(ax: plt.Axes, letter: str, dx: float = -0.16, dy: float = 1.06) -> None:
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="left")


# Row order is shared by panels a, b, c and e so the eye tracks one method across.
ORDER = ["COMMOT", "NICHES", "stLearn", "CytoSignal", "CellChat",
         "LIANA+", "SpatialDM", "ALARMIST"]

SHORT_OBJ = {
    ("CytoSignal", 0): "ranked table", ("CytoSignal", 1): "ranked table",
    ("stLearn", 0): "per-LR file", ("stLearn", 1): "per-LR file",
    ("SpatialDM", 0): "global_res row", ("SpatialDM", 1): "global_res row",
    ("COMMOT", 0): "pathway GRN", ("COMMOT", 1): "pathway ANNEXIN",
    ("NICHES", 0): "Vascular block", ("NICHES", 1): "mGAM block",
    ("CellChat", 0): "outgoing P3", ("CellChat", 1): "incoming P2",
    ("LIANA+", 0): "NMF F1", ("LIANA+", 1): "NMF F3",
    ("ALARMIST", 0): "motif 1", ("ALARMIST", 1): "motif 1",
}


def main() -> None:
    set_style()

    pa = pd.read_csv(DATA / "panel_a_recovery.csv")
    pb = pd.read_csv(DATA / "panel_b_denominators.csv").set_index("method")
    pc = pd.read_csv(DATA / "panel_c_objects.csv").set_index("method")
    pdd = pd.read_csv(DATA / "panel_d_unit.csv")
    pe = pd.read_csv(DATA / "panel_e_grade.csv").set_index("method")

    ypos = {m: i for i, m in enumerate(ORDER)}

    fig = plt.figure(figsize=(5 * U, 4.6 * U))
    gs = fig.add_gridspec(
        3, 3, height_ratios=[1.15, 1.25, 1.0], width_ratios=[2.4, 0.95, 1.6],
        hspace=0.78, wspace=0.6,
    )

    # ------------------------------------------------------------------ a
    ax = fig.add_subplot(gs[0, :2])
    for m in ORDER:
        sub = pa[pa["method"] == m]
        a1 = sub[sub["arm"] == "GRN->SORT1"].iloc[0]
        a2 = sub[sub["arm"] == "ANXA1->FPR1"].iloc[0]
        y = ypos[m]
        col = ACCENT if m == "ALARMIST" else GREY
        ax.plot([a1["percentile"], a2["percentile"]], [y, y], color=col, lw=0.8, zorder=1)
        # Per-core methods: the point is the median of the per-core percentiles and the
        # whisker their IQR. Hiding that spread is the main way this panel could mislead --
        # SpatialDM ranks GRN->SORT1 in the top 1% of one core and the bottom 8% of another.
        for a, off in ((a1, 0.19), (a2, -0.19)):
            if bool(a["per_core"]):
                ax.plot([a["pct_q25"], a["pct_q75"]], [y + off] * 2, color=col, lw=0.8,
                        solid_capstyle="butt", alpha=0.65, zorder=2)
                for q in ("pct_q25", "pct_q75"):
                    ax.plot([a[q]] * 2, [y + off - 0.075, y + off + 0.075], color=col,
                            lw=0.8, alpha=0.65, zorder=2)
        ax.scatter(a1["percentile"], y, s=22, facecolor=col, edgecolor=col, lw=0.8, zorder=3)
        ax.scatter(a2["percentile"], y, s=22, facecolor="white", edgecolor=col, lw=0.8, zorder=3)
        if bool(a1["per_core"]):
            lab = f"median of {int(a1['n_cores'])} cores"
            if int(a1["n_cores"]) != int(a2["n_cores"]):
                lab = f"median of {int(a1['n_cores'])} / {int(a2['n_cores'])} cores"
        else:
            lab = f"{a1['rank']:.0f}/{a1['n_tested']:.0f} · {a2['rank']:.0f}/{a2['n_tested']:.0f}"
            if m == "ALARMIST":
                lab += "   (within motif 1)"
        ax.text(73.0, y, lab, va="center", ha="left", fontsize=5, color=INK)
    ax.axvline(70, color=LIGHT, lw=0.8, zorder=0)
    ax.set_xlim(-2, 104)
    ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70])
    ax.set_ylim(-0.7, len(ORDER) - 0.3)
    # Only stLearn ships a ranked list. Everyone else was sorted by us, and for the two
    # marked with a double dagger the ranked QUANTITY is ours too, because the method
    # emits no per-interaction scalar. That is why they sit at the top of this panel.
    PROV_MARK = {"native": "", "sorted": " †", "derived": " ‡"}
    prov = {m: pa.loc[pa["method"] == m, "rank_provenance"].iloc[0] for m in ORDER}
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels([m + PROV_MARK[prov[m]] for m in ORDER])
    ax.set_xlabel("Rank of the interaction within that method's own tested set (percentile)")
    ax.text(0.0, -0.27,
            "unmarked = the method's own ranked output (stLearn only)   "
            "† ranked by us from the method's own per-interaction statistic\n"
            "‡ the method emits no per-interaction statistic, so the ranked quantity is ours too — "
            "COMMOT: summed transport mass · NICHES: fraction of cells detected · ALARMIST: motif-1 weight",
            transform=ax.transAxes, ha="left", va="top", fontsize=4.5, color=GREY,
            linespacing=1.5)
    ax.spines["left"].set_bounds(-0.7, len(ORDER) - 0.3)
    ax.spines["bottom"].set_bounds(0, 70)
    handles = [
        Line2D([], [], marker="o", ls="none", mfc=GREY, mec=GREY, ms=4, label="GRN→SORT1"),
        Line2D([], [], marker="o", ls="none", mfc="white", mec=GREY, ms=4, label="ANXA1→FPR1"),
        Line2D([], [], color=GREY, lw=0.8, alpha=0.65, label="IQR across cores"),
    ]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2,
              handletextpad=0.3, columnspacing=1.2)
    panel_letter(ax, "a", dx=-0.13)

    # ------------------------------------------------------------------ b
    axb = fig.add_subplot(gs[0, 2])
    vals = [pb.loc[m, "n_tested"] for m in ORDER]
    cols = [ACCENT if m == "ALARMIST" else GREY for m in ORDER]
    axb.barh(range(len(ORDER)), vals, color=cols, height=0.55, lw=0)
    for i, v in enumerate(vals):
        axb.text(v * 1.15, i, f"{v:,.0f}", va="center", ha="left", fontsize=5)
    axb.set_xscale("log")
    axb.set_xlim(30, 9000)
    axb.set_xticks([100, 1000])
    axb.set_xticklabels(["100", "1,000"])
    axb.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    axb.set_ylim(-0.7, len(ORDER) - 0.3)
    axb.set_yticks(range(len(ORDER)))
    axb.set_yticklabels([])
    axb.set_xlabel("Interactions tested (n)")
    axb.spines["left"].set_bounds(-0.7, len(ORDER) - 0.3)
    panel_letter(axb, "b", dx=-0.10)

    # ------------------------------------------------------------------ c
    axc = fig.add_subplot(gs[1, :2])
    BW, BH = 0.80, 0.52          # token box, in data coordinates
    XC = (0.0, 1.05)
    for m in ORDER:
        y = ypos[m]
        r = pc.loc[m]
        fc = OBJ_COLOR[r["obj_type"]]
        ec = GREY if r["obj_type"] == "flat" else fc
        if r["same_object"]:
            axc.plot(XC, [y, y], color=ACCENT, lw=1.6, zorder=1)
        for j, xx in enumerate(XC):
            axc.add_patch(Rectangle((xx - BW / 2, y - BH / 2), BW, BH, facecolor=fc,
                                    edgecolor=ec, lw=0.9, zorder=3))
            axc.text(xx, y, SHORT_OBJ[(m, j)], ha="center", va="center", fontsize=4.8,
                     color=INK, zorder=4)
        axc.text(1.62, y, "same object" if r["same_object"] else "different objects",
                 va="center", ha="left", fontsize=5,
                 color=ACCENT if r["same_object"] else GREY)
    axc.set_xlim(-0.52, 2.85)
    axc.set_ylim(-0.7, len(ORDER) - 0.3)
    axc.set_xticks(list(XC))
    axc.set_xticklabels(["GRN→SORT1", "ANXA1→FPR1"])
    axc.set_yticks(range(len(ORDER)))
    axc.set_yticklabels(ORDER)
    axc.spines["left"].set_bounds(-0.7, len(ORDER) - 0.3)
    axc.spines["bottom"].set_bounds(XC[0], XC[1])
    axc.set_xlabel("Which output object each arm lives in")
    hc = [Line2D([], [], marker="s", ls="none", mfc=OBJ_COLOR[k], mec=GREY if k == "flat" else OBJ_COLOR[k],
                 ms=5, label=OBJ_LABEL[k]) for k in ("flat", "curated", "learned")]
    axc.legend(handles=hc, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=1,
               handletextpad=0.3, labelspacing=0.25)
    panel_letter(axc, "c", dx=-0.13)

    # ------------------------------------------------------------------ d
    gsd = gs[1, 2].subgridspec(1, 2, wspace=1.15)
    dcols = [GREY, ACCENT]
    for j, (col, lab, fmt) in enumerate([
        ("both_pct", "Units carrying both arms (%)", "{:.2f}"),
        ("pearson", "Correlation between arms ($r$)", "{:.3f}"),
    ]):
        axd = fig.add_subplot(gsd[0, j])
        v = pdd[col].to_numpy()
        axd.bar(range(2), v, color=dcols, width=0.62, lw=0)
        for i, vv in enumerate(v):
            axd.text(i, vv, fmt.format(vv), ha="center", va="bottom", fontsize=5)
        axd.set_xticks(range(2))
        axd.set_xticklabels(["cell", "50 µm\npatch"], fontsize=5)
        axd.set_xlim(-0.6, 1.6)
        axd.set_ylim(0, max(v) * 1.32)
        axd.set_ylabel(lab, fontsize=5.2, labelpad=2)
        axd.spines["bottom"].set_bounds(-0.5, 1.5)
        if j == 0:
            panel_letter(axd, "d", dx=-0.60)

    # ------------------------------------------------------------------ e
    axe = fig.add_subplot(gs[2, :])
    NOTE = {"CellChat": "no punch-level test",
            "NICHES": "no p-value emitted",
            "stLearn": "no between-condition test",
            "COMMOT": "no between-condition test"}
    XMAX = 3.2
    for m in ORDER:
        y = ypos[m]
        r = pe.loc[m]
        col = ACCENT if m == "ALARMIST" else GREY
        if m in NOTE:
            axe.add_patch(Rectangle((0, y - 0.3), XMAX, 0.6, facecolor=LIGHT, lw=0, zorder=0))
            axe.text(0.09, y, NOTE[m], va="center", ha="left", fontsize=5, color=GREY,
                     style="italic", zorder=2)
        else:
            axe.scatter(-np.log10(r["p_arm1"]), y, s=28, facecolor=col, edgecolor=col,
                        lw=0.8, zorder=3)
        unit = r["replicate_unit"]
        txt = unit if unit == "none" else f"{unit},  n = {r['n']:,.0f}"
        axe.text(XMAX + 0.15, y, txt, va="center", ha="left", fontsize=5, color=INK,
                 fontweight="bold" if unit == "TMA punch" else "normal")
    axe.axvline(-np.log10(0.05), color=INK, lw=0.8, ls=(0, (3, 2)), zorder=4)
    axe.text(-np.log10(0.05), len(ORDER) - 0.32, "$p$ = 0.05", fontsize=5, ha="center", va="bottom")
    axe.set_xlim(0, 7.0)
    axe.set_ylim(-0.7, len(ORDER) - 0.25)
    axe.set_xticks([0, 1, 2, 3])
    axe.set_yticks(range(len(ORDER)))
    axe.set_yticklabels(ORDER)
    axe.set_xlabel("Grade association of the motif-1 arms, $-$log$_{10}$ $p$")
    axe.spines["left"].set_bounds(-0.7, len(ORDER) - 0.25)
    axe.spines["bottom"].set_bounds(0, XMAX)
    panel_letter(axe, "e", dx=-0.088)

    save_all_formats(fig, FIG / "figure6_comparator_benchmark")
    plt.close(fig)


if __name__ == "__main__":
    print("writing figures to", FIG)
    main()
