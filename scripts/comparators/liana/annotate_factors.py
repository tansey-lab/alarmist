#!/usr/bin/env python
"""Pathway-annotate EXISTING LIANA+ NMF factors. No factorisation is re-run.

Two independent annotation routes are computed for each branch:

  (A) PROGENy footprint route -- the authors' own recipe, from
      docs/notebooks/mofatalk.ipynb cells 62/63/65:

          net  = dc.op.progeny(organism='human', top=5000, thr_padj=0.25)
          lrgs = li.rs.generate_lr_geneset(lr_pairs, net, lr_sep='^')
                   .rename(columns={'interaction': 'target'})
          est, pval = dc.mt.mlm(loadings.T, lrgs, tmin=5)

      CAVEAT, stated up front: PROGENy is a *transcriptional footprint* resource -- its
      targets are genes whose EXPRESSION responds to pathway perturbation. An NMF factor
      over LR-pair features is not a gene expression profile, and generate_lr_geneset only
      transfers a gene-level weight onto the ligand and the receptor of a pair. The score
      is therefore "do the LR pairs this factor loads on happen to involve genes that are
      PROGENy targets of pathway P", which is a much weaker statement than pathway
      activity. Coverage is also poor (see run_manifest.json: pct_features_in_lr_progeny).
      Route (B) is the annotation we actually stand behind.

  (B) CellChatDB v2 'pathway_name' composition -- assumption-free. Every feature carries an
      exact curator-assigned pathway label, so a factor's pathway identity is read directly
      off its loading vector. Reported three ways:
        * loading-weighted composition  (fraction of a factor's total loading per pathway)
        * top-N feature pathway counts
        * dc.mt.ulm on the unweighted pathway membership net -> t-value, p, BH-FDR

Both NMF branches are annotated:
  nmf_bivariate  131  'LIG^REC'             features x 6 factors
  nmf_inflow    2704  'SENDER^LIG^REC'      features x 7 factors  (sender split off)

TRAP HONOURED: li.ut.get_variable_loadings re-sorts rows by |Factor1|, so NMF_H_loadings.csv
row order does NOT match nmf_WH.npz['features']. Everything here is joined on the feature
NAME held in the CSV index; nmf_WH.npz is never opened.

Usage:
  python annotate_factors.py --results-root results/comparators/liana/GBM \
      --cellchatdb data/LRdatabase/CellChatDBv2.0.human.csv \
      --out-dir results/comparators/liana/GBM/factor_annotation
"""
import argparse
import json
import os
import resource as _resource
import subprocess
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.multitest import multipletests

import decoupler as dc
import liana as li

# ----------------------------------------------------------------------------- args
p = argparse.ArgumentParser()
p.add_argument("--results-root", required=True,
               help="dir holding nmf_bivariate/ and nmf_inflow/")
p.add_argument("--cellchatdb", required=True,
               help="the exact LRI csv the LIANA runs used (ligand,receptor,pathway,...)")
p.add_argument("--out-dir", required=True)
p.add_argument("--branches", nargs="+", default=["bivariate", "inflow"])
p.add_argument("--progeny-top", type=int, default=5000, help="mofatalk.ipynb cell 62")
p.add_argument("--progeny-thr-padj", type=float, default=0.25, help="mofatalk.ipynb cell 62")
p.add_argument("--tmin", type=int, default=5, help="mofatalk.ipynb cell 65")
p.add_argument("--top-n", type=int, default=25, help="top features per factor for route (B)")
p.add_argument("--n-show", type=int, default=12, help="pathways per factor in the barplots")
p.add_argument("--seed", type=int, default=0)
p.add_argument("--focus-lr", nargs="+", default=["GRN^SORT1", "ANXA1^FPR1"],
               help="LR pairs to trace explicitly across factors")
a = p.parse_args()

np.random.seed(a.seed)
DATA = os.path.join(a.out_dir, "data")
PLOTS = os.path.join(a.out_dir, "plots")
for d in (DATA, PLOTS):
    os.makedirs(d, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common.plotting import apply_publication_style

apply_publication_style(**{"font.size": 8})


def savefig(fig, name):
    # NOTE: png only, not the png+pdf+svg house rule — these are working
    # annotation figures, not manuscript panels. Recorded in DEVIATIONS.md.
    fp = os.path.join(PLOTS, name + ".png")
    fig.savefig(fp, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log("  fig ->", fp)


def bh(pvals_df):
    """BH-FDR across pathways WITHIN each factor (row)."""
    out = pvals_df.copy().astype(float)
    for i in out.index:
        v = out.loc[i].values
        ok = np.isfinite(v)
        q = np.full_like(v, np.nan, dtype=float)
        if ok.sum():
            q[ok] = multipletests(v[ok], method="fdr_bh")[1]
        out.loc[i] = q
    return out


# ------------------------------------------------------------------ prior knowledge
log("loading CellChatDB:", a.cellchatdb)
db = pd.read_csv(a.cellchatdb)
for c in ("ligand", "receptor", "pathway"):
    assert c in db.columns, f"{c} missing from {a.cellchatdb}"
db["lr"] = db["ligand"] + "^" + db["receptor"]
log(f"  {len(db)} rows, {db['lr'].nunique()} unique LR pairs, {db['pathway'].nunique()} pathways")

# (A) PROGENy -> LR geneset. dc.op.progeny fetches from OmniPath: requires network.
log("fetching PROGENy from OmniPath (dc.op.progeny)")
progeny = dc.op.progeny(organism="human", top=a.progeny_top, thr_padj=a.progeny_thr_padj)
log(f"  PROGENy net {progeny.shape}, {progeny['source'].nunique()} pathways")
progeny.to_csv(os.path.join(DATA, "progeny_net.csv"), index=False)

# DEVIATION: mofatalk passes li.rs.select_resource('consensus'). Our LIANA runs used
# CellChatDB v2, so the LR universe must be CellChatDB or the geneset would not name the
# features we have to annotate.
lr_progeny = li.rs.generate_lr_geneset(
    db[["ligand", "receptor"]].drop_duplicates(), progeny, lr_sep="^"
).rename(columns={"interaction": "target"})
log(f"  lr_progeny {lr_progeny.shape}, {lr_progeny['source'].nunique()} pathways, "
    f"{lr_progeny['target'].nunique()} LR pairs")
lr_progeny.to_csv(os.path.join(DATA, "lr_progeny_geneset.csv"), index=False)

# (B) CellChatDB pathway membership net, unweighted (a pair belongs to its pathway).
cc_net = db[["pathway", "lr"]].drop_duplicates().rename(
    columns={"pathway": "source", "lr": "target"})
cc_net["weight"] = 1.0

manifest_branches = {}
report_rows = []


def annotate(branch, loadings, feat2lr, tag):
    """loadings: features x factors (index = raw feature name). feat2lr: name -> LR pair."""
    log(f"--- {tag}: {loadings.shape[0]} features x {loadings.shape[1]} factors")
    factors = list(loadings.columns)
    info = {"n_features": int(loadings.shape[0]), "n_factors": len(factors)}

    # --- collapse to LR-pair level (inflow: sum over senders; bivariate: identity).
    lr = pd.Series({f: feat2lr(f) for f in loadings.index}, name="lr")
    lr_load = loadings.copy()
    lr_load["__lr__"] = lr.reindex(lr_load.index).values     # aligned BY NAME, never zipped
    lr_load = lr_load.groupby("__lr__")[factors].sum()
    info["n_lr_pairs"] = int(lr_load.shape[0])
    lr_load.to_csv(os.path.join(DATA, f"{tag}_lrpair_loadings.csv"))

    covered = lr_load.index.isin(set(lr_progeny["target"]))
    info["pct_lr_pairs_in_lr_progeny"] = round(100 * float(covered.mean()), 2)
    info["n_lr_pairs_in_lr_progeny"] = int(covered.sum())
    log(f"  PROGENy coverage: {covered.sum()}/{len(covered)} LR pairs "
        f"({info['pct_lr_pairs_in_lr_progeny']}%)")

    # ================================ (A) PROGENy via MLM (authors' method) + ULM check
    prog_tables = {}
    for meth_name, meth in (("mlm", dc.mt.mlm), ("ulm", dc.mt.ulm)):
        est, pv = meth(lr_load.T, lr_progeny, tmin=a.tmin, verbose=False)
        fdr = bh(pv)
        long = (est.melt(ignore_index=False, value_name="score", var_name="pathway")
                   .reset_index(names="factor")
                   .merge(pv.melt(ignore_index=False, value_name="pval", var_name="pathway")
                            .reset_index(names="factor"), on=["factor", "pathway"])
                   .merge(fdr.melt(ignore_index=False, value_name="fdr", var_name="pathway")
                            .reset_index(names="factor"), on=["factor", "pathway"]))
        long["method"] = meth_name
        long.to_csv(os.path.join(DATA, f"{tag}_progeny_{meth_name}.csv"), index=False)
        prog_tables[meth_name] = (est, pv, fdr, long)
        log(f"  progeny/{meth_name}: {est.shape[1]} pathways passed tmin={a.tmin}; "
            f"{int((fdr < 0.05).values.sum())} factor-pathway hits at FDR<0.05")
    info["n_progeny_pathways_tested"] = int(prog_tables["mlm"][0].shape[1])
    info["n_progeny_hits_fdr05_mlm"] = int((prog_tables["mlm"][2] < 0.05).values.sum())

    est, pv, fdr, _ = prog_tables["mlm"]
    if est.shape[1]:
        fig, ax = plt.subplots(figsize=(0.55 * est.shape[1] + 2.4, 0.42 * len(factors) + 1.8))
        vmax = float(np.nanmax(np.abs(est.values))) or 1.0
        ann = np.where(fdr.values < 0.01, "**",
                       np.where(fdr.values < 0.05, "*",
                                np.where(pv.values < 0.05, "·", "")))
        sns.heatmap(est, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax, ax=ax,
                    annot=ann, fmt="", annot_kws={"size": 9},
                    cbar_kws={"label": "MLM t-value"}, linewidths=.5, linecolor="white")
        ax.set_title(f"{tag}: PROGENy enrichment of NMF factor loadings (MLM)\n"
                     f"** FDR<0.01, * FDR<0.05, · p<0.05 uncorrected  |  "
                     f"{info['n_lr_pairs_in_lr_progeny']}/{info['n_lr_pairs']} LR pairs annotated")
        ax.set_xlabel("PROGENy pathway"); ax.set_ylabel("")
        savefig(fig, f"{tag}_progeny_mlm_heatmap")

    # ================================ (B) CellChatDB pathway composition
    pw = db.drop_duplicates("lr").set_index("lr")["pathway"]
    lr_pw = pw.reindex(lr_load.index)                        # aligned BY NAME
    info["pct_lr_pairs_in_cellchatdb"] = round(100 * float(lr_pw.notna().mean()), 2)

    comp = lr_load.assign(__pw__=lr_pw.values).groupby("__pw__")[factors].sum()
    comp_frac = comp / comp.sum(axis=0)
    comp.to_csv(os.path.join(DATA, f"{tag}_cellchatdb_pathway_loading_sum.csv"))
    comp_frac.to_csv(os.path.join(DATA, f"{tag}_cellchatdb_pathway_fraction.csv"))

    # top-N feature composition
    tops = []
    for f in factors:
        s = loadings[f].sort_values(ascending=False).head(a.top_n)
        d = pd.DataFrame({"factor": f, "feature": s.index, "loading": s.values})
        d["lr"] = [feat2lr(x) for x in d["feature"]]
        if branch == "inflow":
            d["sender"] = [x.split("^")[0] for x in d["feature"]]
        d["pathway"] = pw.reindex(d["lr"]).values
        d["rank"] = np.arange(1, len(d) + 1)
        tops.append(d)
    tops = pd.concat(tops, ignore_index=True)
    tops.to_csv(os.path.join(DATA, f"{tag}_top{a.top_n}_features_annotated.csv"), index=False)
    topcount = (tops.groupby(["factor", "pathway"]).size().unstack(fill_value=0)
                    .reindex(factors).T)
    topcount.to_csv(os.path.join(DATA, f"{tag}_top{a.top_n}_pathway_counts.csv"))

    # ULM on unweighted CellChatDB pathway membership -> p / FDR for route (B)
    est_cc, pv_cc = dc.mt.ulm(lr_load.T, cc_net, tmin=a.tmin, verbose=False)
    fdr_cc = bh(pv_cc)
    long_cc = (est_cc.melt(ignore_index=False, value_name="score", var_name="pathway")
                     .reset_index(names="factor")
                     .merge(pv_cc.melt(ignore_index=False, value_name="pval", var_name="pathway")
                              .reset_index(names="factor"), on=["factor", "pathway"])
                     .merge(fdr_cc.melt(ignore_index=False, value_name="fdr", var_name="pathway")
                              .reset_index(names="factor"), on=["factor", "pathway"]))
    long_cc.sort_values(["factor", "pval"]).to_csv(
        os.path.join(DATA, f"{tag}_cellchatdb_pathway_ulm.csv"), index=False)
    info["n_cellchatdb_pathways_tested"] = int(est_cc.shape[1])
    info["n_cellchatdb_hits_fdr05_ulm"] = int((fdr_cc < 0.05).values.sum())
    log(f"  cellchatdb/ulm: {est_cc.shape[1]} pathways passed tmin={a.tmin}; "
        f"{info['n_cellchatdb_hits_fdr05_ulm']} hits at FDR<0.05")

    # heatmap: union of each factor's top pathways by loading fraction
    keep = sorted({pw_ for f in factors
                   for pw_ in comp_frac[f].sort_values(ascending=False).head(a.n_show).index})
    sub = comp_frac.loc[keep]
    order = sub.values.argmax(axis=1).argsort()
    sub = sub.iloc[order]
    fig, ax = plt.subplots(figsize=(0.42 * len(factors) + 3.2, 0.24 * len(sub) + 2.0))
    sns.heatmap(sub, cmap="magma_r", ax=ax, linewidths=.4, linecolor="white",
                cbar_kws={"label": "fraction of factor's total loading"})
    ax.set_title(f"{tag}: CellChatDB v2 pathway composition of NMF factors\n"
                 f"(loading-weighted; union of per-factor top {a.n_show})")
    ax.set_xlabel(""); ax.set_ylabel("CellChatDB pathway")
    savefig(fig, f"{tag}_cellchatdb_pathway_heatmap")

    # per-factor top-pathway barplot
    ncol = min(4, len(factors))
    nrow = int(np.ceil(len(factors) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.5 * ncol, 0.24 * a.n_show * nrow + 1.2))
    axes = np.atleast_1d(axes).ravel()
    for ax_, f in zip(axes, factors):
        s = comp_frac[f].sort_values(ascending=False).head(a.n_show)[::-1]
        ax_.barh(s.index, s.values, color="#3b6ea5")
        ax_.set_title(f, fontsize=9)
        ax_.set_xlabel("frac. of factor loading", fontsize=7)
        ax_.tick_params(labelsize=7)
    for ax_ in axes[len(factors):]:
        ax_.axis("off")
    fig.suptitle(f"{tag}: top CellChatDB pathways per factor (loading-weighted)", y=1.005)
    fig.tight_layout()
    savefig(fig, f"{tag}_cellchatdb_top_pathways_barplot")

    # --------------------------------------------- inflow only: sender composition
    if branch == "inflow":
        snd_lab = pd.Series({f: f.split("^")[0] for f in loadings.index})   # BY NAME
        snd_comp = loadings.assign(__s__=snd_lab.reindex(loadings.index).values) \
                           .groupby("__s__")[factors].sum()
        snd_frac = snd_comp / snd_comp.sum(axis=0)
        snd_frac.to_csv(os.path.join(DATA, f"{tag}_sender_fraction.csv"))
        fig, ax = plt.subplots(figsize=(0.5 * len(factors) + 3.4, 0.3 * len(snd_frac) + 1.8))
        sns.heatmap(snd_frac, cmap="magma_r", ax=ax, annot=True, fmt=".2f",
                    annot_kws={"size": 6}, linewidths=.4, linecolor="white",
                    cbar_kws={"label": "fraction of factor's total loading"})
        ax.set_title(f"{tag}: sender cell-type composition of NMF factors")
        ax.set_ylabel("sender cell type"); ax.set_xlabel("")
        savefig(fig, f"{tag}_sender_composition_heatmap")

    # --------------------------------------------- inflow only: sender-resolved PROGENy
    if branch == "inflow":
        senders = sorted({f.split("^")[0] for f in loadings.index})
        rows = []
        for snd in senders:
            m = [f for f in loadings.index if f.split("^")[0] == snd]
            sl = loadings.loc[m].copy()
            sl.index = [feat2lr(f) for f in m]
            sl = sl.groupby(level=0).sum()
            e, pvv = dc.mt.mlm(sl.T, lr_progeny, tmin=a.tmin, verbose=False)
            d = (e.melt(ignore_index=False, value_name="score", var_name="pathway")
                  .reset_index(names="factor")
                  .merge(pvv.melt(ignore_index=False, value_name="pval", var_name="pathway")
                            .reset_index(names="factor"), on=["factor", "pathway"]))
            d["sender"] = snd
            rows.append(d)
        snd_df = pd.concat(rows, ignore_index=True)
        snd_df["fdr"] = multipletests(snd_df["pval"], method="fdr_bh")[1]
        snd_df.to_csv(os.path.join(DATA, f"{tag}_progeny_mlm_by_sender.csv"), index=False)
        piv = snd_df.pivot_table(index=["sender", "factor"], columns="pathway", values="score")
        fig, ax = plt.subplots(figsize=(0.55 * piv.shape[1] + 3.6, 0.2 * piv.shape[0] + 2.0))
        v = float(np.nanmax(np.abs(piv.values))) or 1.0
        sns.heatmap(piv, cmap="RdBu_r", center=0, vmin=-v, vmax=v, ax=ax,
                    cbar_kws={"label": "MLM t-value"})
        ax.set_title(f"{tag}: PROGENy enrichment per sender x factor (MLM)")
        ax.set_ylabel("sender x factor"); ax.set_xlabel("PROGENy pathway")
        savefig(fig, f"{tag}_progeny_mlm_by_sender_heatmap")

    # --------------------------------------------- focus LR pairs
    foc = []
    for target in a.focus_lr:
        feats = [f for f in loadings.index if feat2lr(f) == target]
        for f in feats:
            for k in factors:
                foc.append({"lr": target, "feature": f, "factor": k,
                            "loading": float(loadings.loc[f, k]),
                            "frac_of_feature": float(loadings.loc[f, k] / loadings.loc[f].sum())
                                               if loadings.loc[f].sum() > 0 else np.nan,
                            "rank_in_factor": int(
                                (loadings[k] > loadings.loc[f, k]).sum() + 1),
                            "n_features": int(loadings.shape[0]),
                            "cellchatdb_pathway": pw.get(target, None)})
    foc = pd.DataFrame(foc)
    if len(foc):
        foc.to_csv(os.path.join(DATA, f"{tag}_focus_lr_loadings.csv"), index=False)
        wide = foc.pivot_table(index="feature", columns="factor", values="loading")
        fig, ax = plt.subplots(figsize=(0.5 * len(factors) + 3.2, 0.3 * len(wide) + 1.8))
        sns.heatmap(wide, cmap="viridis", ax=ax, annot=True, fmt=".2g",
                    annot_kws={"size": 6}, cbar_kws={"label": "NMF loading"})
        ax.set_title(f"{tag}: loadings of {', '.join(a.focus_lr)} across factors")
        ax.set_ylabel(""); ax.set_xlabel("")
        savefig(fig, f"{tag}_focus_lr_heatmap")

    # --------------------------------------------- per-factor call
    for f in factors:
        top_cc = comp_frac[f].sort_values(ascending=False).head(5)
        pg = prog_tables["mlm"][3]
        pgf = pg[pg["factor"] == f].sort_values("pval")
        report_rows.append({
            "branch": tag, "factor": f,
            "top_senders": ("; ".join(f"{i} ({v:.1%})" for i, v in
                            snd_frac[f].sort_values(ascending=False).head(3).items())
                            if branch == "inflow" else "n/a (no sender axis)"),
            "top_cellchatdb_pathways": "; ".join(
                f"{i} ({v:.1%})" for i, v in top_cc.items()),
            "top_feature": loadings[f].idxmax(),
            "top_feature_loading": round(float(loadings[f].max()), 3),
            "best_progeny": (f"{pgf.iloc[0]['pathway']} t={pgf.iloc[0]['score']:.2f} "
                             f"p={pgf.iloc[0]['pval']:.3g} FDR={pgf.iloc[0]['fdr']:.3g}")
                            if len(pgf) else "n/a",
        })
    return info


for branch in a.branches:
    src = os.path.join(a.results_root, f"nmf_{branch}", "data", "NMF_H_loadings.csv")
    log("reading", src)
    L = pd.read_csv(src, index_col=0)
    L.index = L.index.astype(str)
    assert L.index.is_unique, "duplicate feature names"
    nsep = {f.count("^") for f in L.index}
    if branch == "inflow":
        assert nsep == {2}, f"expected SENDER^LIG^REC, got sep counts {nsep}"
        f2lr = lambda x: "^".join(x.split("^")[1:])
    else:
        assert nsep == {1}, f"expected LIG^REC, got sep counts {nsep}"
        f2lr = lambda x: x
    manifest_branches[branch] = annotate(branch, L, f2lr, branch)
    manifest_branches[branch]["loadings_csv"] = src

pd.DataFrame(report_rows).to_csv(os.path.join(DATA, "factor_calls_summary.csv"), index=False)
log("\n" + pd.DataFrame(report_rows).to_string(max_colwidth=80))

sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                     capture_output=True, text=True).stdout.strip() or None
peak_rss_gb = _resource.getrusage(_resource.RUSAGE_SELF).ru_maxrss / 1024**3  # bytes on macOS
json.dump({
    "script": "annotate_factors.py",
    "purpose": "PROGENy + CellChatDB pathway annotation of pre-existing NMF factors",
    "refactorised": False,
    "inputs": {b: manifest_branches[b]["loadings_csv"] for b in manifest_branches},
    "cellchatdb": a.cellchatdb,
    "non_default_parameters": {
        "progeny_top": a.progeny_top, "progeny_thr_padj": a.progeny_thr_padj,
        "tmin": a.tmin, "top_n": a.top_n, "n_show": a.n_show,
        "lr_sep": "^", "progeny_organism": "human",
        "lr_universe": "CellChatDBv2.0.human (DEVIATION: mofatalk uses select_resource('consensus'))",
        "methods": ["dc.mt.mlm (authors)", "dc.mt.ulm (robustness / CellChatDB net)"],
        "fdr": "Benjamini-Hochberg across pathways within each factor",
    },
    "seed": a.seed,
    "branches": manifest_branches,
    "versions": {"decoupler": dc.__version__, "liana": li.__version__,
                 "pandas": pd.__version__, "numpy": np.__version__},
    "git_sha": sha,
    "wall_min": round((time.time() - t0) / 60, 2),
    "peak_rss_gb": round(peak_rss_gb, 3),
}, open(os.path.join(a.out_dir, "run_manifest.json"), "w"), indent=2)
log(f"DONE in {(time.time()-t0)/60:.2f} min -> {a.out_dir}")
