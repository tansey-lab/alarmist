#!/usr/bin/env python
"""Where do ALARMIST's motif-10 / motif-24 LRIs sit in each comparator's own ranking?

For every comparator method that has LUAD/AIS output on disk, and for every section,
this finds each target ligand-receptor interaction (optionally with its sender and
receiver cell types) inside that method's *own* ranked result table and records its
rank, the size of the ranked set, and a top-percentile.

Two levels are produced, because the methods do not all resolve cell types:

    lri     LR pair only                             -- every method
    lri_ct  LR pair x sender -> receiver, directed  -- CellChat only
            LR pair x {cell type, cell type}, UNDIRECTED -- stLearn
            LR pair x sender only                    -- LIANA+ inflow

Only CellChat resolves a *direction*.  stLearn's cell-type matrix looks directed but
is not: get_interaction_matrix (stlearn/tl/cci/het.py:228-276) accepts a spot pair if
EITHER orientation works and de-duplicates on an unordered edge key, so ligand-side and
receptor-side are interchangeable by construction -- the only thing making the matrix
asymmetric is which spot was called the LR hotspot, not signalling direction.  Its
statistic is therefore symmetrised, max(A[s,r], A[r,s]), and ranked over unordered
cell-type pairs.  Do not read a direction out of the stLearn rows.

Ranking rule, applied uniformly
-------------------------------
Rank inside the method's own exported result table for that section, by that
method's own headline statistic, descending, `min` ties.  N is the size of that
table.  `pct_rank = 100 * (1 - (rank - 1) / N)`, so 100 % is the top entry.

Cross-method rank *counts* are meaningless here -- METHODS.md's "The methods do NOT
share a spatial support" applies -- which is exactly why `pct_rank` exists and why
the figure uses it.

Three non-ranked outcomes are kept distinct, because "the method never tested this"
is not the same claim as "the method tested it and found nothing":

    ranked      found, has a rank
    zero        the LR is inside the method's tested universe, but this entry
                carries no signal (absent from a positives-only export, or below
                the method's own significance filter)
    not_tested  the LR is not in this method's tested universe at all -- e.g.
                stLearn dropped every multi-subunit complex, CellChat's DB tier
                has no TBXAS1_TBXA2R
    no_data     the section has not finished running

Usage
-----
    /data1/tanseyw/fanj2/envs/comp-liana/bin/python \
        scripts/comparators/luad_motif_rank_gather.py [--out-dir DIR] [--methods a,b]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = Path(os.environ.get("LUAD_RESULTS_DIR", "/data1/tanseyw/projects/fanj2/results"))
TIER = os.environ.get("LUAD_TIER", "cellchatdb2")

# Section order is ALARMIST's own row ordering; do not reorder (luad_config.sh:70).
SECTIONS = ["P17_AIS", "P17_LUAD", "P21_AIS", "P21_LUAD"]
SECTION_META = {s: dict(patient=s.split("_")[0], stage=s.split("_")[1]) for s in SECTIONS}

# ----------------------------------------------------------------------------- helpers


def _rank_table(df: pd.DataFrame, stat: str, ascending: bool = False) -> pd.DataFrame:
    """Add `rank` (1 = best, `min` ties) and `pct_rank` to an already-subset table."""
    out = df.copy()
    out["rank"] = out[stat].rank(ascending=ascending, method="min").astype(int)
    n = len(out)
    out["n_ranked"] = n
    out["pct_rank"] = 100.0 * (1.0 - (out["rank"] - 1) / n) if n else np.nan
    return out


def _row(target, method, level, section, *, status, statistic_name,
         statistic=np.nan, rank=np.nan, n_ranked=np.nan, pct_rank=np.nan,
         observed_key="", note=""):
    meta = SECTION_META.get(section, dict(patient="P17+P21", stage=section.split("(")[0]))
    return dict(
        motif=target.motif, alarmist_rank=target.rank_in_motif,
        sender=target.sender, receiver=target.receiver,
        ligand=target.ligand, receptor=target.receptor,
        lr=f"{target.ligand}|{target.receptor}",
        lr_ct=f"{target.sender}|{target.receiver}|{target.ligand}|{target.receptor}",
        pathway=target.pathway, factor_lrnorm=target.factor_lrnorm,
        method=method, level=level, section=section,
        patient=meta["patient"], stage=meta["stage"],
        statistic_name=statistic_name, statistic=statistic,
        rank=rank, n_ranked=n_ranked, pct_rank=pct_rank,
        status=status, observed_key=observed_key, note=note,
    )


# ----------------------------------------------------------------------------- stLearn
# lr_summary.csv arrives already sorted descending by n_spots_sig -- verified
# `s.n_spots_sig.is_monotonic_decreasing == True` -- which is stLearn's own ranking
# column (METHODS.md, stLearn / Data outputs).
#
# per_lr_cci_cell_type is ALREADY significance-masked -- analysis.py:751 does
# `sig_int_matrix[int_pvals > p_cutoff] = 0` -- so `count > 0` in that matrix already
# means p <= 0.05 and no second gate is needed (18,562 nonzero vs 18,450 at p<0.05 on
# P17_AIS).  Ranking by the count is stLearn's own cell-type statistic: get_int_df
# (pl/cci_plot_helpers.py:394) returns exactly this matrix, and it is what ccinet_plot,
# cci_map and lr_cci_map display.

def gather_stlearn(targets: pd.DataFrame) -> list[dict]:
    rows = []
    root = RESULTS / "stlearn" / "LUAD" / TIER
    requested = root.parent / f"{TIER}_lrs.txt"
    universe = set(requested.read_text().split()) if requested.exists() else set()

    for section in SECTIONS:
        d = root / section / "data"
        summ = d / "lr_summary.csv"
        if not summ.exists():
            for t in targets.itertuples():
                rows.append(_row(t, "stLearn", "lri", section, status="no_data",
                                 statistic_name="n_spots_sig"))
            continue

        s = pd.read_csv(summ, index_col=0)
        s = _rank_table(s.reset_index().rename(columns={"index": "lr"}), "n_spots_sig")
        s = s.set_index("lr")

        for t in targets.itertuples():
            key = f"{t.ligand}_{t.receptor}"
            if key in s.index:
                r = s.loc[key]
                rows.append(_row(t, "stLearn", "lri", section, status="ranked",
                                 statistic_name="n_spots_sig", statistic=float(r.n_spots_sig),
                                 rank=int(r["rank"]), n_ranked=int(r.n_ranked),
                                 pct_rank=float(r.pct_rank), observed_key=key))
            else:
                # stLearn's requested-LR list holds zero keys with >=3 tokens: every
                # multi-subunit CellChatDB complex was dropped before the run.
                complexed = "_" in t.receptor or "_" in t.ligand
                note = ("multi-subunit complex; stLearn's LR universe has no 3-token keys"
                        if complexed else "not in stLearn's expression-filtered LR set")
                status = "not_tested" if (complexed or key not in universe) else "zero"
                rows.append(_row(t, "stLearn", "lri", section, status=status,
                                 statistic_name="n_spots_sig", note=note))

        # ---- cell-type level, SYMMETRISED (see the module docstring)
        cci_dir, p_dir = d / "per_lr_cci_cell_type", d / "per_lr_cci_pvals_cell_type"
        if not cci_dir.is_dir():
            continue
        frames = []
        for f in sorted(glob.glob(str(cci_dir / "*.csv"))):
            lr = Path(f).stem
            c = pd.read_csv(f, index_col=0)
            # max(A[i,j], A[j,i]): stLearn cannot tell the two apart, so collapse them
            # into one unordered pair rather than reporting two pseudo-directed ranks.
            sym = np.maximum(c.values, c.values.T)
            sym = pd.DataFrame(sym, index=c.index, columns=c.columns)
            keep = np.triu(np.ones(sym.shape, dtype=bool))       # upper triangle + diagonal
            m = sym.where(keep).stack().rename("count").reset_index()
            m.columns = ["ct_a", "ct_b", "count"]
            m = m[m["count"] > 0]
            m["lr"] = lr
            frames.append(m)
        allsig = pd.concat(frames, ignore_index=True)
        allsig = _rank_table(allsig, "count")
        # a sorted "A||B" string, not a frozenset: pandas would unpack a frozenset
        # inside .loc as a list of column labels.
        allsig["pair"] = ["||".join(sorted((a, b)))
                          for a, b in zip(allsig.ct_a, allsig.ct_b)]
        idx = allsig.set_index(["lr", "pair"])

        lrs_with_signal = set(allsig.lr)
        for t in targets.itertuples():
            key = (f"{t.ligand}_{t.receptor}", "||".join(sorted((t.sender, t.receiver))))
            if key in idx.index:
                r = idx.loc[key]
                if isinstance(r, pd.DataFrame):
                    r = r.iloc[0]
                rows.append(_row(t, "stLearn", "lri_ct", section, status="ranked",
                                 statistic_name="n_sig spot pairs, symmetrised max",
                                 statistic=float(r["count"]), rank=int(r["rank"]),
                                 n_ranked=int(r.n_ranked), pct_rank=float(r.pct_rank),
                                 observed_key=f"{{{t.sender},{t.receiver}}}|{t.ligand}_{t.receptor}",
                                 note="UNDIRECTED -- stLearn's kernel cannot separate "
                                      "sender from receiver"))
            elif key[0] in lrs_with_signal:
                rows.append(_row(t, "stLearn", "lri_ct", section, status="zero",
                                 statistic_name="n_sig spot pairs, symmetrised max",
                                 note="LR tested; this cell-type pair not significant "
                                      "in either orientation"))
            else:
                complexed = "_" in t.receptor or "_" in t.ligand
                rows.append(_row(t, "stLearn", "lri_ct", section,
                                 status="not_tested" if complexed else "zero",
                                 statistic_name="n_sig spot pairs, symmetrised max",
                                 note="multi-subunit complex dropped by stLearn" if complexed
                                      else "LR significant in no cell-type pair"))
    return rows


# --------------------------------------------------------------------------- SpatialDM
# global_res.csv is the whole tested universe (1,692 rows, identical in all four
# sections).  `z` is the standardised global bivariate Moran's R; SpatialDM's own
# `selected` column is the FDR call.  All 1,692 carry a z, so the full table is
# ranked and `selected` is recorded alongside rather than used to truncate.
#
# There is no sender x receiver decomposition: the statistic is a single spatial
# cross-correlation over spots, not a directed cell-type quantity.

def gather_spatialdm(targets: pd.DataFrame) -> list[dict]:
    rows = []
    root = RESULTS / "spatialdm" / "LUAD" / TIER
    for section in SECTIONS:
        f = root / section / "data" / "global_res.csv"
        if not f.exists():
            for t in targets.itertuples():
                rows.append(_row(t, "SpatialDM", "lri", section, status="no_data",
                                 statistic_name="z"))
            continue
        g = pd.read_csv(f)
        g = _rank_table(g, "z").set_index("interaction_name")
        for t in targets.itertuples():
            key = f"{t.ligand}_{t.receptor}"
            if key in g.index:
                r = g.loc[key]
                rows.append(_row(t, "SpatialDM", "lri", section, status="ranked",
                                 statistic_name="z", statistic=float(r.z),
                                 rank=int(r["rank"]), n_ranked=int(r.n_ranked),
                                 pct_rank=float(r.pct_rank), observed_key=key,
                                 note=f"selected={r.selected}, fdr={r.fdr:.3g}"))
            else:
                rows.append(_row(t, "SpatialDM", "lri", section, status="not_tested",
                                 statistic_name="z",
                                 note="absent from global_res.csv"))
    return rows


# -------------------------------------------------------------------------- CytoSignal
# Three signif_summary files per section.  Measured on P17_AIS:
#   diffusion n contact = 0        -> each LRI belongs to exactly one signalling mode
#   contact n Raw       = 139      -> Raw is a strict SUBSET of contact, i.e. the same
#                                     contact-type LRIs scored on the un-imputed slot,
#                                     not a third mode.
# DO NOT POOL THE TWO MODES.  An LRI is assigned to exactly one mode by CellChatDB
# `signaling_type` (Cell-Cell Contact -> contact; Secreted / ECM-Receptor / Non-protein
# -> diffusion), and their denominators differ 5.3-fold (tested universe 912 diffusion
# vs 171 contact, constant across sections because the 5,101-gene panel is) with
# different cell-count scales -- a pooled rank would mostly encode mode membership.
# Rank WITHIN the mode.  Of our 35 target LRIs, 20 are diffusion and 15 are contact.
#
# The ranked set is `result.hq`, the tier CytoSignal itself reports: rows with
# n_hq > 0.  Rows with n_hq == 0 are present in the CSV (they had significant cells at
# p<0.05 but failed the reads.thresh/sig.thresh QC wholesale) and must NOT be ranked --
# e.g. EFNB2-EPHB1 is n_result=32, n_hq=0 in P17_AIS contact.  n_hq per section:
# diffusion 642/865/837/858, contact 146/166/169/166 (P17_AIS/P17_LUAD/P21_AIS/P21_LUAD).
# `Raw` is dropped: it is a strict SUBSET of the contact ids, a same-spot sensitivity
# view, not a third mode.
#
# CytoSignal scores are per receiver cell; the runner exports no cell-type breakdown,
# so there is no lri_ct level.

def gather_cytosignal(targets: pd.DataFrame) -> list[dict]:
    rows = []
    root = RESULTS / "cytosignal" / "LUAD" / TIER
    for section in SECTIONS:
        q = root / section / "quant"
        parts = []
        for mode in ("diffusion", "contact"):
            f = q / f"signif_summary_{mode}_Raw_smooth.csv"
            if f.exists():
                p = pd.read_csv(f)
                p["mode"] = mode
                parts.append(p)
        if not parts:
            for t in targets.itertuples():
                rows.append(_row(t, "CytoSignal", "lri", section, status="no_data",
                                 statistic_name="n_hq"))
            continue
        cs = pd.concat(parts, ignore_index=True)
        cs["ligand"] = cs.name.str.split(" - ").str[0].str.strip()
        cs["receptor"] = cs.name.str.split(" - ").str[1].str.strip()
        cs["lr"] = cs.ligand + "|" + cs.receptor
        tested = set(cs.lr)                       # everything CytoSignal scored at all
        # rank INSIDE each mode, over result.hq (n_hq > 0) only
        ranked = pd.concat([_rank_table(g[g.n_hq > 0], "n_hq")
                            for _, g in cs.groupby("mode")], ignore_index=True)
        idx = ranked.set_index("lr")
        for t in targets.itertuples():
            key = f"{t.ligand}|{t.receptor}"
            if key in idx.index:
                r = idx.loc[key]
                if isinstance(r, pd.DataFrame):
                    r = r.iloc[0]
                rows.append(_row(t, "CytoSignal", "lri", section, status="ranked",
                                 statistic_name="n_hq (significant cells), within mode",
                                 statistic=float(r.n_hq),
                                 rank=int(r["rank"]), n_ranked=int(r.n_ranked),
                                 pct_rank=float(r.pct_rank), observed_key=r["name"],
                                 note=f"mode={r['mode']}"))
            else:
                seen = key in tested
                rows.append(_row(t, "CytoSignal", "lri", section, status="zero",
                                 statistic_name="n_hq (significant cells), within mode",
                                 note="scored but n_hq == 0 (failed reads/sig QC)" if seen
                                      else "no significant cells in its signalling mode"))
    return rows


# ----------------------------------------------------------------------------- LIANA+
# `feature_var.csv` indexes the inflow feature space as  <sender>^<LIGAND>^<RECEPTOR>.
# METHODS.md, "Branch: inflow": inflow carries SENDER identity inside the feature
# (C_{j,s} is a hard sender-cell-type indicator), so the leading token is the sender
# and there is no receiver axis.  `run_inflow_downstream.py` -- which would have given
# the full source x target x LR table -- is deliberately NOT run for LUAD
# (05_liana.sh:89-101: it requires --region-col, and inflow here is already per section).
#
# `mean` is the mean inflow score over cells; ranking by it descending is the natural
# read of the feature table.  For the LR-only level the sender axis is collapsed with
# max(), i.e. "the strongest sender for that LR".

def gather_liana(targets: pd.DataFrame) -> list[dict]:
    rows = []
    root = RESULTS / "liana" / "LUAD"
    for section in SECTIONS:
        f = root / section / f"{TIER}_inflow" / "data" / "feature_var.csv"
        if not f.exists():
            for t in targets.itertuples():
                rows.append(_row(t, "LIANA+ (inflow)", "lri", section, status="no_data",
                                 statistic_name="mean inflow"))
            continue
        v = pd.read_csv(f, index_col=0)
        parts = v.index.to_series().str.split("^", expand=True)
        v["sender"], v["ligand"], v["receptor"] = parts[0], parts[1], parts[2]
        v["lr"] = v.ligand + "|" + v.receptor

        ct = _rank_table(v.reset_index(drop=True), "mean")
        ct_idx = ct.set_index(["lr", "sender"])
        lri = _rank_table(v.groupby("lr", as_index=False)["mean"].max(), "mean").set_index("lr")

        for t in targets.itertuples():
            key = f"{t.ligand}|{t.receptor}"
            if key in lri.index:
                r = lri.loc[key]
                rows.append(_row(t, "LIANA+ (inflow)", "lri", section, status="ranked",
                                 statistic_name="max mean inflow over senders",
                                 statistic=float(r["mean"]), rank=int(r["rank"]),
                                 n_ranked=int(r.n_ranked), pct_rank=float(r.pct_rank),
                                 observed_key=key))
            else:
                rows.append(_row(t, "LIANA+ (inflow)", "lri", section, status="zero",
                                 statistic_name="max mean inflow over senders",
                                 note="LR did not clear inflow's nz_prop=0.001 for any sender"))

            k2 = (key, t.sender)
            if k2 in ct_idx.index:
                r = ct_idx.loc[k2]
                if isinstance(r, pd.DataFrame):
                    r = r.iloc[0]
                rows.append(_row(t, "LIANA+ (inflow)", "lri_ct", section, status="ranked",
                                 statistic_name="mean inflow (sender x LR)",
                                 statistic=float(r["mean"]), rank=int(r["rank"]),
                                 n_ranked=int(r.n_ranked), pct_rank=float(r.pct_rank),
                                 observed_key=f"{t.sender}^{t.ligand}^{t.receptor}",
                                 note="SENDER ONLY -- inflow has no receiver axis"))
            else:
                rows.append(_row(t, "LIANA+ (inflow)", "lri_ct", section, status="zero",
                                 statistic_name="mean inflow (sender x LR)",
                                 note="SENDER ONLY -- this sender x LR feature not in the space"))
    return rows


# --------------------------------------------------------------------------- CellChat
# The LUAD run is POOLED BY STAGE, not per section: run.log reports
# "AIS: 475240 cells ... 2 samples" and "LUAD: 1200922 cells ... 2 samples", i.e. two
# objects with the two patients as meta$samples.  So CellChat contributes ONE value per
# condition, with no patient replicate.
#
# `<COND>_net_full.csv.gz` is the non-zero communication table (source x target x LR
# with prob > 0); `<COND>_net_significant.csv` is its pval < 0.05 subset and is what
# `subsetCommunication()` returns.  Ranking uses net_significant, because that is the
# table the authors' own code ranks: plot_cellchat_luad.R:278-279 does
#   agg <- aggregate(prob ~ interaction_name, data = subsetCommunication(obj), FUN = sum)
#   agg <- agg[order(-agg$prob), ]
# and writes it as <cond>_lr_ranked.csv.  The LR level reproduces exactly that sum; the
# triple level ranks the 3,737 significant links by prob.
#
# Denominators: 151 LR pairs have >= 1 significant link (and the net_full and
# net_significant interaction_name SETS are identical, so nothing with prob > 0 is
# non-significant at LR level).  The true TESTED universe is 709 over-expressed
# interactions (run.log; cellchat_io.R:121) -- the other 558 sit at prob = 0 and are
# dropped by array3_to_long(drop.zero=TRUE), which is why they read as `zero` here.
#
# Join on (ligand, receptor): CellChat's `interaction_name` is CellChat-native
# (VEGFA_VEGFR2, not VEGFA_KDR), but the `ligand` / `receptor` columns are exactly the
# CellChatDB v2 export naming ALARMIST uses -- VWF/ITGAV_ITGB3, VEGFD/FLT4_KDR,
# LAMC3/ITGA6_ITGB4 all match verbatim.

def gather_cellchat(targets: pd.DataFrame) -> list[dict]:
    rows = []
    q = RESULTS / "cellchat" / "LUAD" / TIER / "quant"
    for cond in ("AIS", "LUAD"):
        section = f"{cond}(pooled)"
        nf_p, db_p = q / f"{cond}_net_full.csv.gz", q / f"{cond}_db_used.csv"
        if not nf_p.exists():
            for t in targets.itertuples():
                for lvl in ("lri", "lri_ct"):
                    rows.append(_row(t, "CellChat", lvl, section, status="no_data",
                                     statistic_name="communication prob",
                                     note="pooled run still in progress"))
            continue
        nf = pd.read_csv(q / f"{cond}_net_significant.csv")
        db = pd.read_csv(db_p)
        universe = set(zip(db.ligand.astype(str), db.receptor.astype(str)))

        ct = _rank_table(nf, "prob")
        ct_idx = ct.set_index(["ligand", "receptor", "source", "target"])
        lri = _rank_table(nf.groupby(["ligand", "receptor"], as_index=False)["prob"].sum(),
                          "prob").set_index(["ligand", "receptor"])

        for t in targets.itertuples():
            lr = (t.ligand, t.receptor)
            in_db = lr in universe
            if lr in lri.index:
                r = lri.loc[lr]
                rows.append(_row(t, "CellChat", "lri", section, status="ranked",
                                 statistic_name="summed communication prob over significant links",
                                 statistic=float(r.prob), rank=int(r["rank"]),
                                 n_ranked=int(r.n_ranked), pct_rank=float(r.pct_rank),
                                 observed_key=f"{t.ligand}|{t.receptor}"))
            else:
                rows.append(_row(t, "CellChat", "lri", section,
                                 status="zero" if in_db else "not_tested",
                                 statistic_name="summed communication prob over significant links",
                                 note="no significant link (of 709 tested LRIs, 151 are rankable)" if in_db
                                      else "not in the cellchatdb2 DB tier"))

            k = (t.ligand, t.receptor, t.sender, t.receiver)
            if k in ct_idx.index:
                r = ct_idx.loc[k]
                if isinstance(r, pd.DataFrame):
                    r = r.iloc[0]
                rows.append(_row(t, "CellChat", "lri_ct", section, status="ranked",
                                 statistic_name="communication prob",
                                 statistic=float(r.prob), rank=int(r["rank"]),
                                 n_ranked=int(r.n_ranked), pct_rank=float(r.pct_rank),
                                 observed_key=f"{t.sender}|{t.receiver}|{t.ligand}|{t.receptor}",
                                 note=f"pval={r.pval:.3g}"))
            else:
                rows.append(_row(t, "CellChat", "lri_ct", section,
                                 status="zero" if in_db else "not_tested",
                                 statistic_name="communication prob",
                                 note="no significant link for this sender->receiver pair" if in_db
                                      else "not in the cellchatdb2 DB tier"))
    return rows


GATHERERS = {
    "cellchat": gather_cellchat,
    "stlearn": gather_stlearn,
    "spatialdm": gather_spatialdm,
    "cytosignal": gather_cytosignal,
    "liana": gather_liana,
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=str(RESULTS / "_figures" / "luad_motif_ranks"))
    ap.add_argument("--targets", default=str(HERE / "motif_targets_luad.csv"))
    ap.add_argument("--methods", default=",".join(GATHERERS))
    args = ap.parse_args()

    targets = pd.read_csv(args.targets)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for m in args.methods.split(","):
        m = m.strip()
        if not m:
            continue
        print(f"[gather] {m} ...", flush=True)
        rows.extend(GATHERERS[m](targets))

    df = pd.DataFrame(rows)
    df["section"] = pd.Categorical(df.section,
                                   categories=SECTIONS + ["AIS(pooled)", "LUAD(pooled)"],
                                   ordered=True)
    df = df.sort_values(["motif", "level", "method", "section", "alarmist_rank"])

    full = out_dir / "luad_motif_ranks_long.csv"
    df.to_csv(full, index=False)
    for lvl in ("lri", "lri_ct"):
        df[df.level == lvl].to_csv(out_dir / f"luad_motif_ranks_{lvl}.csv", index=False)

    # Wide, human-readable: one row per target LRI, one column per method x section.
    # Cells carry the rank, or the status token when there is no rank.
    def cell(r):
        if r.status == "ranked":
            return f"{int(r['rank'])}/{int(r.n_ranked)}"
        return {"zero": "ns", "not_tested": "-", "no_data": "?"}[r.status]

    for motif in sorted(df.motif.unique()):
        for lvl in ("lri", "lri_ct"):
            d = df[(df.motif == motif) & (df.level == lvl)].copy()
            if d.empty:
                continue
            d["cell"] = d.apply(cell, axis=1)
            idx = ["alarmist_rank", "lr"] if lvl == "lri" else ["alarmist_rank", "lr_ct"]
            for name, val in (("rank", "cell"), ("pct", "pct_rank")):
                w = d.pivot_table(index=idx + ["pathway", "factor_lrnorm"],
                                  columns=["method", "section"], values=val,
                                  aggfunc="first", observed=True)
                if name == "pct":
                    w = w.round(1)
                w = w.sort_index(level="alarmist_rank")
                w.to_csv(out_dir / f"wide_motif{motif}_{lvl}_{name}.csv")
    print("[gather] wide tables: wide_motif{10,24}_{lri,lri_ct}_{rank,pct}.csv"
          "   (cell = rank/N; 'ns' = tested but not reported; '-' = not in the "
          "method's LR universe; '?' = run did not finish)")

    print(f"\n[gather] wrote {full}  ({len(df)} rows)")
    print("\n=== coverage: status counts by method x level ===")
    print(pd.crosstab([df.method, df.level], df.status).to_string())
    print("\n=== ranked-set size N per method x level x section ===")
    n = (df[df.status == "ranked"]
         .groupby(["method", "level", "section"], observed=True)["n_ranked"].max())
    print(n.to_string())
    print("\n=== median top-percentile of the 20 target LRIs, by motif x method x stage ===")
    piv = (df[df.status == "ranked"]
           .groupby(["level", "motif", "method", "stage"], observed=True)["pct_rank"]
           .median().round(1).unstack("stage"))
    print(piv.to_string())

    (out_dir / "gather_manifest.json").write_text(json.dumps({
        "results_root": str(RESULTS), "tier": TIER, "sections": SECTIONS,
        "methods": args.methods, "targets": args.targets,
        "n_rows": int(len(df)), "python": sys.version.split()[0],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
