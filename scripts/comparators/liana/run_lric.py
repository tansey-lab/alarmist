#!/usr/bin/env python
"""LIANA+ ``lric`` + ``cross_pcf`` run PER TISSUE PUNCH (TMA core).

Why per punch and not whole slide
---------------------------------
Both ``li.mt.lric`` and ``li.mt.cross_pcf`` are *density-normalised* point-pattern
statistics: they build a ``cKDTree``, count neighbours in annuli, and divide the
observed count by the count expected under complete spatial randomness at the
observed density.  The density they use is ``n_points / bounding-box area``
(``_LRIC.py::_corrected_areas`` / ``CrossPCF._compute_pair`` /
``LRIC._compute_pair``), and the bounding box is the axis-aligned box of the
points handed to the call.

On a tissue microarray the *global* bounding box is mostly empty space between
cores.  Running whole-slide therefore understates the density and inflates every
``g(r)`` by the ratio (global bbox) / (sum of per-core bboxes).  The script
re-verifies that ratio on the data it is given and records it in the manifest.

Contrast: ``li.mt.bivariate`` / ``li.mt.inflow`` *were* safe whole-slide, because
their Gaussian kernel has finite support (28.2 um at our bandwidth) which is far
below the minimum inter-core distance, so no cross-core pair ever forms.  Density
normalisation is a different exposure and it is not protected by that argument.

Call contract
-------------
Parameters follow ``LRIC_tutorial.ipynb`` and were verified against the installed
``liana/method/sp/_LRIC.py`` (liana 1.8.1) before this script was written:
``annulus_width=25``, ``max_radius=200``, ``radius_step=25``, ``min_cells=50``,
``min_expressing=20`` (cell-type mode only), ``use_raw=False``,
``extend_first_annulus=True`` (default).

Realised bins are ``[0,50) [50,75) [75,100) ... [200,225)`` -- 8 bins, the FIRST
IS DOUBLE-WIDTH because ``extend_first_annulus`` merges the contact band, and
``max_radius`` is the INNER edge of the last bin so the true reach is 225 um.

``li.ut.spatial_neighbors`` is deliberately NOT called: neither function consumes
a connectivity graph, they build their own KD-tree, so the 13.1454 um bandwidth
used by the bivariate/inflow runs is irrelevant here and must not be harmonised
in.

Nothing is hardcoded; every path, column and parameter is a CLI flag.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import resource as _res
import subprocess
import sys
import threading
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import plotting as _plot


# ── small utilities ───────────────────────────────────────────────────────────

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-file", required=True)
    p.add_argument("--resource", required=True,
                   help="LR database csv with 'ligand'/'receptor' columns")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--sample-column", default="tma_id",
                   help="obs column defining the independent point patterns (punches)")
    p.add_argument("--cell-type-column", default="cell_type")
    p.add_argument("--condition-column", default="grade",
                   help="obs column with the per-sample condition (one value per sample)")
    p.add_argument("--samples", default=None,
                   help="comma-separated subset of sample ids; default = all")
    p.add_argument("--spatial-key", default="spatial")
    p.add_argument("--count-layer", default="counts",
                   help="layer holding raw counts; '' to use .X as-is")
    p.add_argument("--drop-var-pattern", default="_",
                   help="drop var_names containing this substring (control probes; "
                        "'_' is also LIANA's heteromeric-subunit separator)")
    p.add_argument("--min-genes", type=int, default=10)
    p.add_argument("--min-cells-gene", type=int, default=3)
    p.add_argument("--target-sum", type=float, default=1e4)
    # LRIC / cross-PCF tutorial parameters
    p.add_argument("--annulus-width", type=float, default=25.0)
    p.add_argument("--max-radius", type=float, default=200.0)
    p.add_argument("--radius-step", type=float, default=25.0)
    p.add_argument("--min-cells", type=int, default=50)
    p.add_argument("--min-expressing", type=int, default=20)
    p.add_argument("--n-angle-samples", type=int, default=360)
    p.add_argument("--no-extend-first-annulus", action="store_true")
    # reporting focus
    p.add_argument("--targets", default="GRN^SORT1,ANXA1^FPR1",
                   help="LR pairs (liana lr_sep '^') to report in detail")
    p.add_argument("--focus-directions", default="mGAM->MES-like,MES-like->mGAM",
                   help="sender->receiver directed cell-type pairs to report in detail")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--seed-note", default="both methods are deterministic; --seed only "
                                          "fixes numpy global state for reproducibility "
                                          "hygiene", help=argparse.SUPPRESS)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--whole-slide-check", action="store_true",
                   help="additionally run the SAME calls on all samples pooled into a "
                        "single point pattern, to measure empirically how much the "
                        "TMA's empty inter-core space inflates g(r). This is a CONTROL, "
                        "not a result: it is the thing this script exists to avoid.")
    return p.parse_args(argv)


class RssSampler:
    """Sample RSS in a background thread so the peak is attributable to the call."""

    def __init__(self, interval=0.01):
        import psutil
        self._proc = psutil.Process()
        self._interval = interval
        self._stop = threading.Event()
        self.base = self._proc.memory_info().rss
        self.peak = self.base

    def _loop(self):
        while not self._stop.is_set():
            self.peak = max(self.peak, self._proc.memory_info().rss)
            time.sleep(self._interval)

    def __enter__(self):
        self.base = self._proc.memory_info().rss
        self.peak = self.base
        self._t0 = time.time()
        self._th = threading.Thread(target=self._loop, daemon=True)
        self._th.start()
        return self

    def __exit__(self, *a):
        self.wall = time.time() - self._t0
        self._stop.set()
        self._th.join()
        self.delta = self.peak - self.base


def ru_maxrss_gb():
    v = _res.getrusage(_res.RUSAGE_SELF).ru_maxrss
    return v / 1e9 if sys.platform == "darwin" else v * 1024 / 1e9


def nan_stats(m):
    """(mean, max, argmax-index) per column of a (n_bins, n_pairs) matrix, NaN-safe."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean(m, axis=0)
        mx = np.nanmax(m, axis=0)
    allnan = np.isnan(m).all(axis=0)
    am = np.where(allnan, 0, np.nanargmax(np.nan_to_num(m, nan=-np.inf), axis=0))
    return mean, mx, am, allnan


def bbox_area(xy):
    return float((xy[:, 0].max() - xy[:, 0].min()) * (xy[:, 1].max() - xy[:, 1].min()))


# ── per-punch worker ──────────────────────────────────────────────────────────

def prepare_punch(adata_full, sid, args):
    """Subset to one punch and apply the project preprocessing recipe."""
    import scanpy as sc

    if sid is None:                       # pooled: every cell, one point pattern
        sub = adata_full.copy()
    else:
        sub = adata_full[
            adata_full.obs[args.sample_column].astype(str) == str(sid)].copy()
    n_raw = sub.n_obs
    if args.count_layer:
        sub.X = sub.layers[args.count_layer].copy()
    sub.layers.clear()
    sub.raw = None
    sub.uns.pop("log1p", None)
    n_var_raw = sub.n_vars
    if args.drop_var_pattern:
        keep = ~pd.Index(sub.var_names).str.contains(args.drop_var_pattern, regex=False)
        sub = sub[:, np.asarray(keep)].copy()
    n_dropped_var = n_var_raw - sub.n_vars
    sc.pp.filter_cells(sub, min_genes=args.min_genes)
    sc.pp.filter_genes(sub, min_cells=args.min_cells_gene)
    sc.pp.normalize_total(sub, target_sum=args.target_sum)
    sc.pp.log1p(sub)
    return sub, n_raw, n_dropped_var


def target_expression_support(sub, args, targets, focus_cts):
    """Per-cell-type expressing-cell counts for the target ligands/receptors.

    ``min_expressing`` in ``li.mt.lric`` is applied to ``(weights > 0).sum(0)``
    within each cell type, and ``_linear_transform`` preserves zeros, so this
    equals the raw nonzero count per cell type -- which lets us attribute every
    NaN target to a specific gate (cell type below ``min_cells``, gene lost to
    ``filter_genes``, or expressing cells below ``min_expressing``).
    """
    from scipy import sparse

    genes = sorted({g for t in targets for g in t.split("^")})
    ctv = sub.obs[args.cell_type_column].astype(str).values
    rows = []
    for ct in focus_cts:
        m = ctv == ct
        for g in genes:
            present = g in set(map(str, sub.var_names))
            if present:
                col = sub[:, g].X
                col = col.toarray().ravel() if sparse.issparse(col) else np.ravel(col)
                n_expr = int((col[m] > 0).sum())
            else:
                n_expr = 0
            rows.append(dict(cell_type=ct, gene=g, gene_present_after_filter=present,
                             n_cells=int(m.sum()), n_expressing=n_expr))
    return rows


def run_punch(sub, args, resource, outdir):
    """Run cross_pcf + agnostic LRIC + cell-type LRIC on one prepared punch."""
    import liana as li

    extend = not args.no_extend_first_annulus
    common = dict(spatial_key=args.spatial_key,
                  max_radius=args.max_radius,
                  radius_step=args.radius_step,
                  annulus_width=args.annulus_width,
                  extend_first_annulus=extend,
                  n_angle_samples=args.n_angle_samples)
    bench = []

    with RssSampler() as s:
        li.mt.cross_pcf(sub, groupby=args.cell_type_column,
                        min_cells=args.min_cells, key_added="cross_pcf",
                        verbose=False, **common)
    bench.append(dict(call="cross_pcf", wall_s=round(s.wall, 3),
                      peak_rss_gb=round(s.peak / 1e9, 3),
                      delta_rss_gb=round(s.delta / 1e9, 3)))

    with RssSampler() as s:
        li.mt.lric(sub, resource=resource, use_raw=False,
                   key_added="lric_ag", verbose=False, **common)
    bench.append(dict(call="lric_agnostic", wall_s=round(s.wall, 3),
                      peak_rss_gb=round(s.peak / 1e9, 3),
                      delta_rss_gb=round(s.delta / 1e9, 3)))

    with RssSampler() as s:
        li.mt.lric(sub, resource=resource, groupby=args.cell_type_column,
                   min_cells=args.min_cells, min_expressing=args.min_expressing,
                   use_raw=False, key_added="lric_ct", verbose=False, **common)
    bench.append(dict(call="lric_celltype", wall_s=round(s.wall, 3),
                      peak_rss_gb=round(s.peak / 1e9, 3),
                      delta_rss_gb=round(s.delta / 1e9, 3)))

    return sub.uns["cross_pcf"], sub.uns["lric_ag"], sub.uns["lric_ct"], bench


def persist_punch(sid, cpcf, ag, ct, radii_inner, radii_outer, outdir):
    """Write the per-punch result tables; returns the cell-type summary frame."""
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for (s, r), g in cpcf["results"].items():
        for b, val in enumerate(np.asarray(g, dtype=float)):
            rows.append((s, r, radii_inner[b], radii_outer[b], val))
    cpcf_df = pd.DataFrame(rows, columns=["sender", "receiver", "r_inner", "r_outer", "g"])
    cpcf_df.to_csv(outdir / "cross_pcf_long.csv", index=False)

    L = np.asarray(ag["lric"], dtype=float)
    ag_df = pd.DataFrame(L, index=[f"r{int(x)}" for x in radii_inner],
                         columns=ag["pair_names"]).T
    ag_df.index.name = "pair"
    ag_df.to_csv(outdir / "lric_agnostic_matrix.csv.gz", compression="gzip")

    mean, mx, am, allnan = nan_stats(L)
    ag_scores = (pd.DataFrame({"pair": ag["pair_names"], "mean_lric": mean,
                               "max_lric": mx, "peak_radius": radii_inner[am]})
                 .assign(max_minus_mean=lambda d: d["max_lric"] - d["mean_lric"])
                 .sort_values("mean_lric", ascending=False).reset_index(drop=True))
    ag_scores.to_csv(outdir / "lric_agnostic_pair_scores.csv", index=False)

    ct_keys = list(ct["results"].keys())
    stack = np.stack([np.asarray(ct["results"][k], dtype=np.float32) for k in ct_keys])
    np.savez_compressed(outdir / "lric_celltype.npz", stack=stack,
                        ct_pairs=np.array([f"{s}->{r}" for s, r in ct_keys]),
                        pair_names=np.array(ct["pair_names"]),
                        radii_inner=radii_inner, radii_outer=radii_outer)

    pn = np.array(ct["pair_names"])
    recs = []
    for (s, r), m in ct["results"].items():
        m = np.asarray(m, dtype=float)
        mean, mx, am, allnan = nan_stats(m)
        d = pd.DataFrame({"sender": s, "receiver": r, "pair": pn,
                          "mean_lric": mean, "max_lric": mx,
                          "peak_radius": radii_inner[am], "all_nan": allnan})
        # rank by mean among the LR pairs that survived min_expressing (1 = highest)
        valid = ~allnan
        d["n_valid_lr_pairs"] = int(valid.sum())
        rk = np.full(len(d), np.nan)
        order = np.argsort(-np.where(valid, np.nan_to_num(mean, nan=-np.inf), -np.inf))
        rk[order[: valid.sum()]] = np.arange(1, valid.sum() + 1)
        d["rank_by_mean"] = rk
        for b in range(m.shape[0]):
            d[f"g_bin{b}"] = m[b]
        recs.append(d)
    ct_df = pd.concat(recs, ignore_index=True)
    ct_df.to_csv(outdir / "lric_celltype_summary.csv.gz", index=False, compression="gzip")
    return cpcf_df, ag_scores, ct_df, L


def extract_targets(sid, cond, cpcf, ag, ct, ct_df, ag_scores, L, ct_counts,
                    targets, focus_dirs, radii_inner, radii_outer):
    """Pull the focus LR x direction curves (+ agnostic curves) out of one run."""
    agnostic_rows, target_rows = [], []
    for t in targets:
        if t in ag["pair_names"]:
            i = ag["pair_names"].index(t)
            rank_ag = int(ag_scores.index[ag_scores["pair"] == t][0]) + 1
            for b in range(L.shape[0]):
                agnostic_rows.append(dict(sample=sid, condition=cond, lr=t,
                                          r_inner=radii_inner[b],
                                          r_outer=radii_outer[b], g=L[b, i],
                                          rank_by_mean=rank_ag,
                                          n_lr_pairs=len(ag["pair_names"])))

    cpcf_lookup = {(s, r): np.asarray(g, dtype=float)
                   for (s, r), g in cpcf["results"].items()}
    for (snd, rcv) in focus_dirs:
        key = (snd, rcv)
        present_ct = key in ct["results"]
        present_pcf = key in cpcf_lookup
        for t in targets:
            base = dict(sample=sid, condition=cond, lr=t, sender=snd, receiver=rcv,
                        direction=f"{snd}->{rcv}",
                        n_sender_cells=int(ct_counts.get(snd, 0)),
                        n_receiver_cells=int(ct_counts.get(rcv, 0)),
                        celltypes_present=bool(present_ct))
            if not present_ct or t not in ct["pair_names"]:
                target_rows.append(dict(base, r_inner=np.nan, r_outer=np.nan,
                                        lric_g=np.nan, cross_pcf_g=np.nan,
                                        ratio=np.nan, agnostic_g=np.nan,
                                        rank_by_mean=np.nan, n_valid_lr_pairs=np.nan,
                                        lr_in_resource=t in ct["pair_names"]))
                continue
            j = ct["pair_names"].index(t)
            m = np.asarray(ct["results"][key], dtype=float)[:, j]
            pcf = cpcf_lookup[key] if present_pcf else np.full_like(m, np.nan)
            sel = ct_df[(ct_df["sender"] == snd) & (ct_df["receiver"] == rcv)
                        & (ct_df["pair"] == t)]
            rk = float(sel["rank_by_mean"].iloc[0]) if len(sel) else np.nan
            nv = float(sel["n_valid_lr_pairs"].iloc[0]) if len(sel) else np.nan
            agc = (L[:, ag["pair_names"].index(t)]
                   if t in ag["pair_names"] else np.full_like(m, np.nan))
            for b in range(len(m)):
                target_rows.append(dict(base, r_inner=radii_inner[b],
                                        r_outer=radii_outer[b],
                                        lric_g=m[b], cross_pcf_g=pcf[b],
                                        ratio=m[b] / pcf[b] if pcf[b] else np.nan,
                                        agnostic_g=agc[b], rank_by_mean=rk,
                                        n_valid_lr_pairs=nv, lr_in_resource=True))
    return target_rows, agnostic_rows


# ── plotting ──────────────────────────────────────────────────────────────────

def setup_style():
    import logging

    import matplotlib
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _plot.apply_publication_style(**{
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 150,
        "savefig.bbox": "tight",
    })
    return plt


def save_all_formats(fig, path_no_ext):
    _plot.save_all_formats(fig, path_no_ext, dpi=300)


COND_COLORS = {"high": "#c0392b", "low": "#2c6fbb"}


def _cond_color(c):
    return COND_COLORS.get(str(c), "#7f7f7f")


def plot_curves(long, figdir, value, ylabel, fname, hline=None):
    if not len(long):
        print(f"[plot] skipping {fname}: no data", flush=True)
        return
    plt = setup_style()
    feats = list(dict.fromkeys(zip(long["lr"], long["direction"])))
    ncol = len(feats)
    fig, axes = plt.subplots(1, ncol, figsize=(3.2 * ncol, 3.1), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, (lr, dirn) in zip(axes, feats):
        d = long[(long["lr"] == lr) & (long["direction"] == dirn)]
        for sid, g in d.groupby("sample"):
            g = g.sort_values("r_inner")
            ax.plot(g["r_mid"], g[value], color=_cond_color(g["condition"].iloc[0]),
                    alpha=0.45, lw=1.0)
        med = d.groupby("r_mid")[value].median()
        ax.plot(med.index, med.values, color="k", lw=2.4, label="median across punches")
        if hline is not None:
            ax.axhline(hline, color="0.4", ls="--", lw=1.0)
        ax.set_title(f"{lr}\n{dirn}", fontsize=9)
        ax.set_xlabel("distance r (um, bin midpoint)")
    axes[0].set_ylabel(ylabel)
    handles = [plt.Line2D([], [], color=COND_COLORS["high"], lw=1.5, label="high grade"),
               plt.Line2D([], [], color=COND_COLORS["low"], lw=1.5, label="low grade"),
               plt.Line2D([], [], color="k", lw=2.4, label="median")]
    axes[-1].legend(handles=handles, fontsize=7, frameon=False, loc="best")
    fig.tight_layout()
    save_all_formats(fig, Path(figdir) / fname)
    plt.close(fig)


def plot_condition(summary, figdir, value, ylabel, fname, hline=None):
    if not len(summary):
        print(f"[plot] skipping {fname}: no data", flush=True)
        return
    plt = setup_style()
    feats = list(dict.fromkeys(zip(summary["lr"], summary["direction"])))
    fig, axes = plt.subplots(1, len(feats), figsize=(2.3 * len(feats), 3.2), sharey=True)
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(0)
    for ax, (lr, dirn) in zip(axes, feats):
        d = summary[(summary["lr"] == lr) & (summary["direction"] == dirn)]
        conds = sorted(d["condition"].dropna().unique())
        for i, c in enumerate(conds):
            v = d[d["condition"] == c][value].dropna().values
            if len(v) == 0:
                continue
            ax.boxplot([v], positions=[i], widths=0.55, showfliers=False,
                       medianprops=dict(color="k"))
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=18,
                       color=_cond_color(c), zorder=3, edgecolor="none")
        ax.set_xticks(range(len(conds)))
        ax.set_xticklabels(conds, fontsize=8)
        if hline is not None:
            ax.axhline(hline, color="0.4", ls="--", lw=1.0)
        ax.set_title(f"{lr}\n{dirn}", fontsize=8)
    axes[0].set_ylabel(ylabel)
    fig.tight_layout()
    save_all_formats(fig, Path(figdir) / fname)
    plt.close(fig)


def plot_asymmetry(summary, figdir, fname):
    """Paired forward vs reverse mean g, one line per punch."""
    if not len(summary):
        print(f"[plot] skipping {fname}: no data", flush=True)
        return
    plt = setup_style()
    lrs = list(dict.fromkeys(summary["lr"]))
    fig, axes = plt.subplots(1, len(lrs), figsize=(2.6 * len(lrs), 3.2))
    axes = np.atleast_1d(axes)
    for ax, lr in zip(axes, lrs):
        d = summary[summary["lr"] == lr]
        dirs = list(dict.fromkeys(d["direction"]))
        piv = d.pivot_table(index=["sample", "condition"], columns="direction",
                            values="mean_g", observed=True)
        for (sid, cond), row in piv.iterrows():
            if row.isna().any():
                continue
            ax.plot(range(len(dirs)), [row[x] for x in dirs],
                    marker="o", ms=4, lw=1.0, color=_cond_color(cond), alpha=0.8)
        ax.set_xticks(range(len(dirs)))
        ax.set_xticklabels(dirs, fontsize=7, rotation=15)
        ax.axhline(1.0, color="0.4", ls="--", lw=1.0)
        ax.set_title(lr, fontsize=9)
        ax.set_ylabel("mean LRIC g(r) over bins")
    fig.tight_layout()
    save_all_formats(fig, Path(figdir) / fname)
    plt.close(fig)


REASON_COLORS = {
    "ok": "#2e7d32",
    "below_min_expressing": "#f39c12",
    "celltype_below_min_cells": "#c0392b",
    "gene_absent_after_filter_genes": "#7f5539",
    "unknown": "#999999",
}


def plot_availability(reasons, figdir, fname):
    """punch x feature grid, coloured by the gate that removed each cell."""
    if not len(reasons):
        return
    plt = setup_style()
    feats = list(dict.fromkeys(zip(reasons["lr"], reasons["direction"])))
    samples = list(dict.fromkeys(reasons["sample"]))
    fig, ax = plt.subplots(figsize=(0.62 * len(samples) + 3.0, 0.5 * len(feats) + 2.2))
    for j, (lr, dirn) in enumerate(feats):
        for i, sid in enumerate(samples):
            row = reasons[(reasons["lr"] == lr) & (reasons["direction"] == dirn)
                          & (reasons["sample"] == sid)]
            if not len(row):
                continue
            ax.add_patch(plt.Rectangle((i - 0.5, j - 0.5), 1, 1,
                                       color=REASON_COLORS.get(row["reason"].iloc[0],
                                                               "#999999"),
                                       ec="white", lw=1.2))
    ax.set_xlim(-0.5, len(samples) - 0.5)
    ax.set_ylim(-0.5, len(feats) - 0.5)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels([str(s) for s in samples], fontsize=8)
    ax.set_yticks(range(len(feats)))
    ax.set_yticklabels([f"{a}\n{b}" for a, b in feats], fontsize=7)
    ax.set_xlabel("punch")
    ax.invert_yaxis()
    for sp in ax.spines.values():
        sp.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=v, label=k)
               for k, v in REASON_COLORS.items() if k != "unknown"]
    ax.legend(handles=handles, fontsize=7, frameon=False,
              loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    save_all_formats(fig, Path(figdir) / fname)
    plt.close(fig)


def plot_wholeslide(long, ws, figdir, fname):
    """Per-punch g(r) with the pooled whole-slide curve overlaid (the control)."""
    if not len(long) or not len(ws):
        return
    plt = setup_style()
    feats = list(dict.fromkeys(zip(long["lr"], long["direction"])))
    fig, axes = plt.subplots(1, len(feats), figsize=(3.2 * len(feats), 3.1), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, (lr, dirn) in zip(axes, feats):
        d = long[(long["lr"] == lr) & (long["direction"] == dirn)]
        for _, g in d.groupby("sample"):
            g = g.sort_values("r_inner")
            ax.plot(g["r_mid"], g["lric_g"], color="0.7", lw=0.9)
        med = d.groupby("r_mid")["lric_g"].median()
        ax.plot(med.index, med.values, color="k", lw=2.4)
        w = ws[(ws["lr"] == lr) & (ws["direction"] == dirn)].sort_values("r_inner")
        if len(w):
            ax.plot((w["r_inner"] + w["r_outer"]) / 2, w["lric_g"],
                    color="#8e44ad", lw=2.4, ls="--")
        ax.axhline(1.0, color="0.4", ls=":", lw=1.0)
        ax.set_title(f"{lr}\n{dirn}", fontsize=9)
        ax.set_xlabel("distance r (um, bin midpoint)")
    axes[0].set_ylabel("LRIC g(r)")
    handles = [plt.Line2D([], [], color="0.7", lw=1.0, label="individual punches"),
               plt.Line2D([], [], color="k", lw=2.4, label="per-punch median"),
               plt.Line2D([], [], color="#8e44ad", lw=2.4, ls="--",
                          label="whole slide pooled (WRONG density)")]
    axes[-1].legend(handles=handles, fontsize=7, frameon=False, loc="best")
    fig.tight_layout()
    save_all_formats(fig, Path(figdir) / fname)
    plt.close(fig)


def plot_support(punch_df, figdir, ct, fname):
    plt = setup_style()
    d = punch_df.sort_values(f"n_{ct}", ascending=False)
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    ax.bar(np.arange(len(d)), d[f"n_{ct}"],
           color=[_cond_color(c) for c in d["condition"]])
    ax.axhline(d["min_cells"].iloc[0], color="k", ls="--", lw=1.2)
    ax.set_xticks(np.arange(len(d)))
    ax.set_xticklabels(d["sample"].astype(str), fontsize=8)
    ax.set_xlabel("punch")
    ax.set_ylabel(f"{ct} cells after QC")
    ax.set_yscale("log")
    ax.text(len(d) - 0.5, d["min_cells"].iloc[0] * 1.1,
            f"min_cells = {int(d['min_cells'].iloc[0])}", ha="right", fontsize=7)
    fig.tight_layout()
    save_all_formats(fig, Path(figdir) / fname)
    plt.close(fig)


# ── main ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    args = parse_args(argv)
    np.random.seed(args.seed)
    import scanpy as sc
    import liana as li

    t_start = time.time()
    out = Path(args.output_dir)
    (out / "punches").mkdir(parents=True, exist_ok=True)
    (out / "combined").mkdir(parents=True, exist_ok=True)
    figdir = out / "figures"

    targets = [t for t in args.targets.split(",") if t]
    focus_dirs = [tuple(x.split("->")) for x in args.focus_directions.split(",") if x]

    print(f"[load] {args.data_file}", flush=True)
    adata_full = sc.read_h5ad(args.data_file)
    obs = adata_full.obs
    sample_ids = sorted(obs[args.sample_column].astype(str).unique(),
                        key=lambda s: (len(s), s))
    if args.samples:
        want = set(args.samples.split(","))
        sample_ids = [s for s in sample_ids if s in want]
    cond_map = (obs.groupby(obs[args.sample_column].astype(str), observed=True)
                [args.condition_column].agg(lambda x: str(x.astype(str).iloc[0])).to_dict())

    # ── the density argument, verified on this dataset ────────────────────────
    xy_all = np.asarray(adata_full.obsm[args.spatial_key], dtype=float)
    global_bbox = bbox_area(xy_all)
    per_bbox = {}
    for sid in sample_ids:
        m = (obs[args.sample_column].astype(str) == sid).values
        per_bbox[sid] = bbox_area(xy_all[m])
    sum_bbox = float(sum(per_bbox.values()))
    density_inflation = global_bbox / sum_bbox
    print(f"[bbox] global={global_bbox:,.0f} um^2  sum_per_punch={sum_bbox:,.0f} um^2 "
          f"({sum_bbox / global_bbox:.4%})  whole-slide density understated "
          f"{density_inflation:.4f}x", flush=True)

    db = pd.read_csv(args.resource)
    resource = (db[["ligand", "receptor"]].astype(str)
                .drop_duplicates().reset_index(drop=True))
    print(f"[resource] {args.resource}: {len(resource)} unique LR pairs", flush=True)

    inner_true = np.arange(args.radius_step, args.max_radius + args.radius_step,
                           args.radius_step, dtype=float)
    radii_outer = inner_true + args.annulus_width
    radii_inner = inner_true.copy()
    if not args.no_extend_first_annulus:
        radii_inner[0] = 0.0

    punch_rows, target_rows, skip_rows, bench_rows, agnostic_rows = [], [], [], [], []
    support_rows = []
    focus_cts = list(dict.fromkeys([c for pair in focus_dirs for c in pair]))

    for sid in sample_ids:
        pdir = out / "punches" / f"punch_{sid}"
        print(f"\n=== punch {sid} ({cond_map.get(sid)}) ===", flush=True)
        sub, n_raw, n_dropped_var = prepare_punch(adata_full, sid, args)
        ct_counts = sub.obs[args.cell_type_column].astype(str).value_counts().to_dict()
        kept_cts = sorted([c for c, n in ct_counts.items() if n >= args.min_cells])
        for c, n in sorted(ct_counts.items()):
            skip_rows.append(dict(sample=sid, condition=cond_map.get(sid),
                                  cell_type=c, n_cells=int(n),
                                  kept=bool(n >= args.min_cells),
                                  min_cells=args.min_cells))
        print(f"  cells {n_raw} -> {sub.n_obs}, genes {sub.n_vars} "
              f"(dropped {n_dropped_var} '{args.drop_var_pattern}' probes); "
              f"cell types kept {len(kept_cts)}/{len(ct_counts)}", flush=True)

        for r in target_expression_support(sub, args, targets, focus_cts):
            support_rows.append(dict(sample=sid, condition=cond_map.get(sid), **r))

        cpcf, ag, ct, bench = run_punch(sub, args, resource, pdir)
        for b in bench:
            bench_rows.append(dict(sample=sid, n_cells=int(sub.n_obs), **b))
        cpcf_df, ag_scores, ct_df, L = persist_punch(sid, cpcf, ag, ct,
                                                     radii_inner, radii_outer, pdir)

        xy = np.asarray(sub.obsm[args.spatial_key], dtype=float)
        area = bbox_area(xy)
        prow = dict(sample=sid, condition=cond_map.get(sid), n_cells_raw=int(n_raw),
                    n_cells_used=int(sub.n_obs), n_genes=int(sub.n_vars),
                    n_lr_pairs_agnostic=len(ag["pair_names"]),
                    n_lr_pairs_ct=len(ct["pair_names"]),
                    n_cell_types_total=len(ct_counts),
                    n_cell_types_kept=len(kept_cts),
                    cell_types_kept="|".join(kept_cts),
                    cell_types_skipped="|".join(sorted(set(ct_counts) - set(kept_cts))),
                    bbox_area_um2=area, density_cells_per_um2=sub.n_obs / area,
                    min_cells=args.min_cells,
                    wall_s=round(sum(b["wall_s"] for b in bench), 2),
                    peak_rss_gb=max(b["peak_rss_gb"] for b in bench))
        for c in set(list(ct_counts) + [d for pair in focus_dirs for d in pair]):
            prow[f"n_{c}"] = int(ct_counts.get(c, 0))
        punch_rows.append(prow)

        t_rows, a_rows = extract_targets(sid, cond_map.get(sid), cpcf, ag, ct, ct_df,
                                         ag_scores, L, ct_counts, targets, focus_dirs,
                                         radii_inner, radii_outer)
        target_rows.extend(t_rows)
        agnostic_rows.extend(a_rows)

        with open(pdir / "punch_report.json", "w") as f:
            json.dump({"sample": sid, "condition": cond_map.get(sid),
                       "n_cells_raw": int(n_raw), "n_cells_used": int(sub.n_obs),
                       "n_genes": int(sub.n_vars),
                       "cell_type_counts": {k: int(v) for k, v in ct_counts.items()},
                       "cell_types_kept": kept_cts,
                       "n_lr_pairs_agnostic": len(ag["pair_names"]),
                       "n_lr_pairs_ct": len(ct["pair_names"]),
                       "radii_inner": radii_inner.tolist(),
                       "radii_outer": radii_outer.tolist(),
                       "bbox_area_um2": area, "bench": bench}, f, indent=2)
        del sub

    # ── whole-slide control (the thing per-punch analysis exists to avoid) ────
    ws_rows = []
    if args.whole_slide_check:
        print("\n=== whole-slide control (ALL cells pooled, one point pattern) ===",
              flush=True)
        wdir = out / "whole_slide_check"
        sub, n_raw, _ = prepare_punch(adata_full, None, args)
        ct_counts = sub.obs[args.cell_type_column].astype(str).value_counts().to_dict()
        cpcf, ag, ct, bench = run_punch(sub, args, resource, wdir)
        for b in bench:
            bench_rows.append(dict(sample="ALL_POOLED", n_cells=int(sub.n_obs), **b))
        _, ag_scores, ct_df, L = persist_punch("ALL_POOLED", cpcf, ag, ct,
                                               radii_inner, radii_outer, wdir)
        ws_rows, ws_agn = extract_targets("ALL_POOLED", "pooled", cpcf, ag, ct, ct_df,
                                          ag_scores, L, ct_counts, targets, focus_dirs,
                                          radii_inner, radii_outer)
        pd.DataFrame(ws_rows).to_csv(wdir / "targets_long.csv", index=False)
        pd.DataFrame(ws_agn).to_csv(wdir / "targets_agnostic_long.csv", index=False)
        del sub

    # ── combine ───────────────────────────────────────────────────────────────
    punch_df = pd.DataFrame(punch_rows)
    punch_df.to_csv(out / "combined" / "punch_summary.csv", index=False)
    pd.DataFrame(skip_rows).to_csv(out / "combined" / "celltype_kept_skipped.csv",
                                   index=False)
    sup = pd.DataFrame(support_rows)
    sup.to_csv(out / "combined" / "target_expression_support.csv", index=False)
    pd.DataFrame(bench_rows).to_csv(out / "combined" / "timings.csv", index=False)
    tgt = pd.DataFrame(target_rows)
    tgt["r_mid"] = (tgt["r_inner"] + tgt["r_outer"]) / 2
    tgt.to_csv(out / "combined" / "targets_long.csv", index=False)
    agn = pd.DataFrame(agnostic_rows)
    agn_summ = pd.DataFrame()
    if len(agn):
        agn["r_mid"] = (agn["r_inner"] + agn["r_outer"]) / 2
        agn.to_csv(out / "combined" / "targets_agnostic_long.csv", index=False)
        agn_summ = (agn.groupby(["sample", "condition", "lr"], observed=True)
                    .agg(mean_g=("g", "mean"), max_g=("g", "max"),
                         g_contact_bin=("g", "first"),
                         rank_by_mean=("rank_by_mean", "first"),
                         n_lr_pairs=("n_lr_pairs", "first")).reset_index())
        agn_summ["rank_pct"] = agn_summ["rank_by_mean"] / agn_summ["n_lr_pairs"]
        agn_summ.to_csv(out / "combined" / "agnostic_per_punch_summary.csv", index=False)

    # ── attribute every missing target to the gate that removed it ────────────
    sup_idx = sup.set_index(["sample", "cell_type", "gene"])
    reason_rows = []
    for _, r in (tgt.drop_duplicates(subset=["sample", "lr", "direction"])).iterrows():
        lig, rec = r["lr"].split("^")
        snd, rcv = r["sender"], r["receiver"]
        why, detail = "ok", ""
        if not np.isnan(r["lric_g"]):
            pass
        elif int(r["n_sender_cells"]) < args.min_cells or \
                int(r["n_receiver_cells"]) < args.min_cells:
            why = "celltype_below_min_cells"
            detail = (f"{snd}={int(r['n_sender_cells'])}, "
                      f"{rcv}={int(r['n_receiver_cells'])} vs min_cells={args.min_cells}")
        else:
            try:
                sl = sup_idx.loc[(r["sample"], snd, lig)]
                sr = sup_idx.loc[(r["sample"], rcv, rec)]
            except KeyError:
                sl = sr = None
            if sl is not None and (not sl["gene_present_after_filter"]
                                   or not sr["gene_present_after_filter"]):
                why = "gene_absent_after_filter_genes"
                detail = (f"{lig} present={bool(sl['gene_present_after_filter'])}, "
                          f"{rec} present={bool(sr['gene_present_after_filter'])}")
            elif sl is not None:
                below = []
                if sl["n_expressing"] < args.min_expressing:
                    below.append(f"sender {snd} {lig}+={int(sl['n_expressing'])}")
                if sr["n_expressing"] < args.min_expressing:
                    below.append(f"receiver {rcv} {rec}+={int(sr['n_expressing'])}")
                why = "below_min_expressing" if below else "unknown"
                detail = "; ".join(below) + f" (min_expressing={args.min_expressing})"
        reason_rows.append(dict(sample=r["sample"], condition=r["condition"],
                                lr=r["lr"], direction=r["direction"],
                                available=bool(not np.isnan(r["lric_g"])),
                                reason=why, detail=detail,
                                n_sender_cells=int(r["n_sender_cells"]),
                                n_receiver_cells=int(r["n_receiver_cells"])))
    reasons = pd.DataFrame(reason_rows)
    reasons.to_csv(out / "combined" / "target_availability.csv", index=False)

    # per-punch summary of each focus feature
    have = tgt.dropna(subset=["lric_g"])
    summ = (have.groupby(["sample", "condition", "lr", "direction"], observed=True)
            .agg(mean_g=("lric_g", "mean"), max_g=("lric_g", "max"),
                 min_g=("lric_g", "min"), g_contact_bin=("lric_g", "first"),
                 mean_cross_pcf=("cross_pcf_g", "mean"),
                 mean_ratio=("ratio", "mean"), median_ratio=("ratio", "median"),
                 min_ratio=("ratio", "min"), max_ratio=("ratio", "max"),
                 mean_agnostic=("agnostic_g", "mean"),
                 rank_by_mean=("rank_by_mean", "first"),
                 n_valid_lr_pairs=("n_valid_lr_pairs", "first"),
                 n_sender_cells=("n_sender_cells", "first"),
                 n_receiver_cells=("n_receiver_cells", "first"))
            .reset_index())
    summ["rank_pct"] = summ["rank_by_mean"] / summ["n_valid_lr_pairs"]
    summ.to_csv(out / "combined" / "per_punch_feature_summary.csv", index=False)

    # across-punch aggregate: per bin and overall
    def q(x, p):
        return np.nanpercentile(x, p) if len(x) else np.nan

    agg_bin = (have.groupby(["lr", "direction", "r_inner", "r_outer"], observed=True)
               .agg(n_punches=("lric_g", "size"),
                    median_g=("lric_g", "median"),
                    q25_g=("lric_g", lambda x: q(x, 25)),
                    q75_g=("lric_g", lambda x: q(x, 75)),
                    min_g=("lric_g", "min"), max_g=("lric_g", "max"),
                    median_cross_pcf=("cross_pcf_g", "median"),
                    median_ratio=("ratio", "median"),
                    q25_ratio=("ratio", lambda x: q(x, 25)),
                    q75_ratio=("ratio", lambda x: q(x, 75)))
               .reset_index())
    agg_bin.to_csv(out / "combined" / "aggregate_per_bin.csv", index=False)

    if ws_rows:
        ws = pd.DataFrame(ws_rows).dropna(subset=["lric_g"])
        cmp_ = agg_bin.merge(
            ws[["lr", "direction", "r_inner", "lric_g", "cross_pcf_g", "ratio"]]
            .rename(columns={"lric_g": "whole_slide_lric_g",
                             "cross_pcf_g": "whole_slide_cross_pcf_g",
                             "ratio": "whole_slide_ratio"}),
            on=["lr", "direction", "r_inner"], how="left")
        cmp_["lric_inflation_wholeslide_over_perpunch_median"] = (
            cmp_["whole_slide_lric_g"] / cmp_["median_g"])
        cmp_["crosspcf_inflation_wholeslide_over_perpunch_median"] = (
            cmp_["whole_slide_cross_pcf_g"] / cmp_["median_cross_pcf"])
        cmp_.to_csv(out / "combined" / "whole_slide_vs_perpunch.csv", index=False)
        print("\n[whole-slide control] median inflation of g(r) over the per-punch "
              f"median: LRIC "
              f"{cmp_['lric_inflation_wholeslide_over_perpunch_median'].median():.3f}x, "
              "cross-PCF "
              f"{cmp_['crosspcf_inflation_wholeslide_over_perpunch_median'].median():.3f}x "
              f"(bbox-derived expectation {density_inflation:.3f}x)", flush=True)

    agg_feat = (summ.groupby(["lr", "direction"], observed=True)
                .agg(n_punches=("mean_g", "size"),
                     median_mean_g=("mean_g", "median"),
                     q25=("mean_g", lambda x: q(x, 25)),
                     q75=("mean_g", lambda x: q(x, 75)),
                     min=("mean_g", "min"), max=("mean_g", "max"),
                     median_rank=("rank_by_mean", "median"),
                     median_n_valid=("n_valid_lr_pairs", "median"),
                     median_rank_pct=("rank_pct", "median"),
                     median_ratio=("median_ratio", "median"),
                     q25_ratio=("median_ratio", lambda x: q(x, 25)),
                     q75_ratio=("median_ratio", lambda x: q(x, 75)))
                .reset_index())
    agg_feat.to_csv(out / "combined" / "aggregate_per_feature.csv", index=False)

    # ── statistics: asymmetry (paired) and condition (unpaired) ───────────────
    from scipy import stats
    stat_rows = []
    for lr in dict.fromkeys(summ["lr"]):
        d = summ[summ["lr"] == lr]
        dirs = list(dict.fromkeys(d["direction"]))
        if len(dirs) == 2:
            piv = d.pivot_table(index="sample", columns="direction",
                                values="mean_g", observed=True).dropna()
            a, b = piv[dirs[0]].values, piv[dirs[1]].values
            if len(a) >= 3:
                w = stats.wilcoxon(a, b)
                stat_rows.append(dict(test="wilcoxon_paired_direction", lr=lr,
                                      group_a=dirs[0], group_b=dirs[1],
                                      n_a=len(a), n_b=len(b),
                                      median_a=float(np.median(a)),
                                      median_b=float(np.median(b)),
                                      median_diff=float(np.median(a - b)),
                                      n_a_gt_b=int((a > b).sum()),
                                      statistic=float(w.statistic),
                                      p_value=float(w.pvalue)))
        for dirn in dirs:
            dd = d[d["direction"] == dirn]
            groups = sorted(dd["condition"].dropna().unique())
            if len(groups) == 2:
                a = dd[dd["condition"] == groups[0]]["mean_g"].dropna().values
                b = dd[dd["condition"] == groups[1]]["mean_g"].dropna().values
                if len(a) >= 2 and len(b) >= 2:
                    u = stats.mannwhitneyu(a, b, alternative="two-sided")
                    stat_rows.append(dict(test="mannwhitney_condition_mean_g", lr=lr,
                                          direction=dirn, group_a=groups[0],
                                          group_b=groups[1], n_a=len(a), n_b=len(b),
                                          median_a=float(np.median(a)),
                                          median_b=float(np.median(b)),
                                          median_diff=float(np.median(a) - np.median(b)),
                                          statistic=float(u.statistic),
                                          p_value=float(u.pvalue)))
            # ratio vs 1
            v = dd["median_ratio"].dropna().values
            if len(v) >= 3:
                w = stats.wilcoxon(v - 1.0)
                stat_rows.append(dict(test="wilcoxon_ratio_vs_1", lr=lr, direction=dirn,
                                      n_a=len(v), median_a=float(np.median(v)),
                                      statistic=float(w.statistic),
                                      p_value=float(w.pvalue)))
        # condition effect on the asymmetry itself
        if len(dirs) == 2:
            piv = d.pivot_table(index=["sample", "condition"], columns="direction",
                                values="mean_g", observed=True).dropna().reset_index()
            piv["delta"] = piv[dirs[0]] - piv[dirs[1]]
            groups = sorted(piv["condition"].dropna().unique())
            if len(groups) == 2:
                a = piv[piv["condition"] == groups[0]]["delta"].values
                b = piv[piv["condition"] == groups[1]]["delta"].values
                if len(a) >= 2 and len(b) >= 2:
                    u = stats.mannwhitneyu(a, b, alternative="two-sided")
                    stat_rows.append(dict(test="mannwhitney_condition_asymmetry", lr=lr,
                                          direction=f"{dirs[0]} minus {dirs[1]}",
                                          group_a=groups[0], group_b=groups[1],
                                          n_a=len(a), n_b=len(b),
                                          median_a=float(np.median(a)),
                                          median_b=float(np.median(b)),
                                          median_diff=float(np.median(a) - np.median(b)),
                                          statistic=float(u.statistic),
                                          p_value=float(u.pvalue)))
    # condition effect on the cell-type-AGNOSTIC target curves; this is the only
    # condition contrast with a usable n on both sides when the focus cell types
    # are too rare in one arm to survive min_cells / min_expressing.
    for lr in dict.fromkeys(agn_summ["lr"]) if len(agn_summ) else []:
        d = agn_summ[agn_summ["lr"] == lr]
        groups = sorted(d["condition"].dropna().unique())
        if len(groups) == 2:
            a = d[d["condition"] == groups[0]]["mean_g"].dropna().values
            b = d[d["condition"] == groups[1]]["mean_g"].dropna().values
            if len(a) >= 2 and len(b) >= 2:
                u = stats.mannwhitneyu(a, b, alternative="two-sided")
                stat_rows.append(dict(test="mannwhitney_condition_agnostic_mean_g",
                                      lr=lr, direction="agnostic",
                                      group_a=groups[0], group_b=groups[1],
                                      n_a=len(a), n_b=len(b),
                                      median_a=float(np.median(a)),
                                      median_b=float(np.median(b)),
                                      median_diff=float(np.median(a) - np.median(b)),
                                      statistic=float(u.statistic),
                                      p_value=float(u.pvalue)))

    stats_df = pd.DataFrame(stat_rows)
    stats_df.to_csv(out / "combined" / "statistics.csv", index=False)

    # ── plots ─────────────────────────────────────────────────────────────────
    if not args.no_plots:
        plot_curves(have, figdir, "lric_g", "LRIC g(r)", "gr_curves_per_punch", hline=1.0)
        plot_curves(have, figdir, "cross_pcf_g", "cross-PCF g(r)",
                    "cross_pcf_curves_per_punch", hline=1.0)
        plot_curves(have, figdir, "ratio", "LRIC / cross-PCF", "ratio_vs_distance",
                    hline=1.0)
        plot_condition(summ, figdir, "mean_g", "mean LRIC g(r) over bins",
                       "grade_comparison_mean_g", hline=1.0)
        plot_condition(summ, figdir, "median_ratio", "median LRIC / cross-PCF",
                       "grade_comparison_ratio", hline=1.0)
        plot_asymmetry(summ, figdir, "direction_asymmetry_paired")
        plot_availability(reasons, figdir, "target_availability")
        if ws_rows:
            plot_wholeslide(have, pd.DataFrame(ws_rows).dropna(subset=["lric_g"]),
                            figdir, "whole_slide_vs_perpunch")
        if len(agn):
            a2 = agn.copy()
            a2["direction"] = "agnostic (all cells)"
            a2["lric_g"] = a2["g"]
            plot_curves(a2, figdir, "lric_g", "agnostic LRIC g(r)",
                        "gr_curves_agnostic", hline=1.0)
        if len(agn_summ):
            s2 = agn_summ.copy()
            s2["direction"] = "agnostic (all cells)"
            plot_condition(s2, figdir, "mean_g", "mean agnostic LRIC g(r)",
                           "grade_comparison_agnostic", hline=1.0)
        senders = [s for s, _ in focus_dirs]
        for ctname in dict.fromkeys(senders):
            if f"n_{ctname}" in punch_df.columns:
                plot_support(punch_df, figdir, ctname, f"support_{ctname}_per_punch")

    # ── manifest ──────────────────────────────────────────────────────────────
    manifest = {
        "script": str(Path(__file__).resolve()),
        "command": " ".join([sys.executable] + sys.argv),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "wall_s_total": round(time.time() - t_start, 1),
        "peak_rss_gb_process": round(ru_maxrss_gb(), 3),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "versions": {"liana": li.__version__, "scanpy": sc.__version__,
                     "numpy": np.__version__, "pandas": pd.__version__},
        "seed": args.seed,
        "data_file": args.data_file,
        "resource": args.resource,
        "n_resource_pairs": int(len(resource)),
        "params": {k: v for k, v in vars(args).items()},
        "realised_bins": [{"r_inner": float(a), "r_outer": float(b)}
                          for a, b in zip(radii_inner, radii_outer)],
        "bins_note": ("8 bins; the first is DOUBLE-WIDTH because "
                      "extend_first_annulus merges the [0, radius_step) contact band "
                      "into it. max_radius is the INNER edge of the last bin, so the "
                      "true reach is max_radius + annulus_width."),
        "why_per_punch": {
            "reason": ("li.mt.lric and li.mt.cross_pcf normalise observed neighbour "
                       "counts by the density expected under complete spatial "
                       "randomness, computed as n_points / bounding-box area "
                       "(_LRIC.py::_corrected_areas, CrossPCF._compute_pair, "
                       "LRIC._compute_pair). On a TMA the global bounding box is "
                       "mostly empty space between cores, so a whole-slide run "
                       "understates density and inflates every g(r)."),
            "global_bbox_um2": global_bbox,
            "sum_per_punch_bbox_um2": sum_bbox,
            "occupied_fraction": sum_bbox / global_bbox,
            "whole_slide_density_understated_by": density_inflation,
            "per_punch_bbox_um2": per_bbox,
            "empirically_measured": (
                {"median_lric_inflation_wholeslide_over_perpunch_median": float(
                    pd.read_csv(out / "combined" / "whole_slide_vs_perpunch.csv")
                    ["lric_inflation_wholeslide_over_perpunch_median"].median()),
                 "median_crosspcf_inflation_wholeslide_over_perpunch_median": float(
                     pd.read_csv(out / "combined" / "whole_slide_vs_perpunch.csv")
                     ["crosspcf_inflation_wholeslide_over_perpunch_median"].median())}
                if ws_rows else "not run (pass --whole-slide-check)"),
            "contrast_bivariate_inflow": ("bivariate/inflow WERE safe whole-slide: "
                                          "their Gaussian kernel support (28.2 um at "
                                          "bandwidth 13.1454, cutoff 0.1) is far below "
                                          "the minimum inter-punch distance, so zero "
                                          "cross-punch pairs form. Density "
                                          "normalisation is a different exposure and is "
                                          "NOT protected by that argument."),
            "spatial_neighbors_not_called": ("li.ut.spatial_neighbors is deliberately "
                                             "not called: neither lric nor cross_pcf "
                                             "consumes a connectivity graph, they build "
                                             "their own cKDTree, so the 13.1454 um "
                                             "bandwidth is irrelevant here."),
        },
        "n_samples": len(sample_ids),
        "samples": sample_ids,
        "conditions": {k: cond_map.get(k) for k in sample_ids},
    }
    try:
        manifest["git_commit"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[3],
            text=True).strip()
    except Exception as e:  # noqa: BLE001
        manifest["git_commit"] = f"unavailable: {e}"
    with open(out / "run_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    print("\n[done]", json.dumps({"wall_s": manifest["wall_s_total"],
                                  "peak_rss_gb": manifest["peak_rss_gb_process"],
                                  "n_punches": len(sample_ids)}), flush=True)
    print(agg_feat.to_string(index=False), flush=True)
    if len(stats_df):
        print(stats_df.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
