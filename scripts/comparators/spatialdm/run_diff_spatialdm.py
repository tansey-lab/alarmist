#!/usr/bin/env python
"""SpatialDM native cross-condition differential test (`spatialdm.diff_utils`).

Tutorial: /Users/jiayifan/tansey_lab/SpatialDM/tutorial/differential_test_intestine.ipynb
Call contract: scripts/comparators/spatialdm/NOTES.md   Deviations: DEVIATIONS.md

Consumes the per-split objects already written by `run_spatialdm.py`
(`<run-dir>/<split>/data/spatialdm.h5ad`) and runs the authors' pipeline unmodified:

    concat_obj -> differential_test -> group_differential_pairs
    -> differential_dendrogram / differential_volcano / clustermap / compute_pathway+dot_path

NOTHING is recomputed: every z-score, p-value and `selected` flag comes from the persisted
per-core objects.

Two things are patched, both recorded in DEVIATIONS.md:

  1. `diff_utils.concat_db()` downloads **CellChatDB v1** from figshare and then does
     `geneInter.loc[ligand.index]`. Our per-core `geneInter` was injected from CellChatDB
     **v2** (`attach_cellchat_v2` in run_spatialdm.py), so its interaction names are not in
     the v1 index and the call dies with a KeyError. We monkeypatch `concat_db` to build the
     union frame from the samples' own v2 `geneInter` instead -- same three return values,
     same de-duplication, same `.loc[ligand.index]` reindex. This also removes the network
     dependency. `concat_obj` itself runs untouched.
  2. `dot_path` is called as `pl.dot_path(concat, '<cond>_specific', cut_off=2, ...)` in the
     tutorial, but installed v0.3.1 has signature `dot_path(pathway_res, cut_off=1, ...)` --
     the tutorial call would bind the string to `cut_off`. Version gap: we call
     `compute_pathway(dic={'<c1>_specific': ..., '<c2>_specific': ...})` first and pass the
     resulting frame, which is what 0.3.1 expects and yields one panel per condition.

Usage:
  python run_diff_spatialdm.py --run-dir results/comparators/spatialdm/GBM/cellchatdb2 \
      --out-dir results/comparators/spatialdm/GBM/cellchatdb2/differential_grade \
      --condition-col grade --c1 high --c2 low
"""
import argparse, json, os, resource, subprocess, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import spatialdm as sdm
import spatialdm.plottings as pl
from spatialdm import diff_utils
from spatialdm.diff_utils import concat_obj, differential_test, group_differential_pairs

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common.plotting import apply_publication_style

apply_publication_style(**{"figure.dpi": 110})

p = argparse.ArgumentParser()
p.add_argument("--run-dir", required=True, help="dir holding the per-split run (one subdir per split)")
p.add_argument("--out-dir", required=True)
p.add_argument("--condition-col", default="grade", help="obs column carrying the condition label")
p.add_argument("--c1", required=True, help="condition coded 1 in the OLS design")
p.add_argument("--c2", required=True, help="condition coded 0")
p.add_argument("--splits", default="", help="comma list; default = every split with a saved h5ad")
p.add_argument("--requested-lrs", default="GRN_SORT1,ANXA1_FPR1")
p.add_argument("--fdr", action="store_true", default=False,
               help="concat_obj(fdr=): use per-core fdr instead of z_pval in p_df. Tutorial: False")
p.add_argument("--diff-quantile1", type=float, default=0.7)
p.add_argument("--diff-quantile2", type=float, default=0.3)
p.add_argument("--fdr-co", type=float, default=0.1)
p.add_argument("--path-cut-off", type=int, default=2, help="dot_path min pairs per pathway (tutorial: 2)")
p.add_argument("--density-col", default="median_within_radius",
               help="per_split_summary.csv column used for the density sensitivity design")
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()
np.random.seed(a.seed)

OUT = a.out_dir
DATA, PLOTS, REQ = (os.path.join(OUT, "data"), os.path.join(OUT, "plots"),
                    os.path.join(OUT, "plots", "requested"))
for d in (DATA, PLOTS, REQ):
    os.makedirs(d, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
FAILED, SKIPPED = [], []


# ----------------------------------------------------------------- figure saver
def _has_content(fig):
    for ax in fig.get_axes():
        if (ax.collections or ax.images or ax.lines or ax.patches or ax.texts
                or len(ax.get_children()) > 10):
            return True
    return False


def save_all(stem, fig=None, d=None):
    """png + pdf + svg through one saver (CLAUDE.md); refuse blank figures."""
    fig = fig if fig is not None else plt.gcf()
    d = d or PLOTS
    if not fig.get_axes() or not _has_content(fig):
        plt.close(fig)
        raise RuntimeError("figure is BLANK -- wrong save mechanism for this function")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(d, f"{stem}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def guard(fn, nm):
    try:
        fn(); log(f"  ok: {nm}")
    except Exception as e:
        plt.close("all"); FAILED.append(f"{nm}: {type(e).__name__}: {e}")
        log(f"  FAILED {nm}: {type(e).__name__}: {e}")


# ------------------------------------------------------ patch 1: v2-aware concat_db
def concat_db_v2(adatas, species):
    """Authors' concat_db with the figshare CellChatDB v1 download replaced by the samples'
    own injected v2 geneInter. Identical structure, identical de-duplication."""
    ligand = pd.concat([s.uns["ligand"] for s in adatas], axis=1)
    receptor = pd.concat([s.uns["receptor"] for s in adatas], axis=1)
    ligand = ligand[~ligand.index.duplicated()]
    receptor = receptor[~receptor.index.duplicated()]
    gi = pd.concat([s.uns["geneInter"] for s in adatas], axis=0)
    gi = gi[~gi.index.duplicated()]
    geneInter = gi.loc[ligand.index]
    # v1-column shim. run_spatialdm.py put `interaction_name` in geneInter as a column, but
    # anndata wrote it out AS THE INDEX (`uns/geneInter` attrs carry `_index: interaction_name`)
    # and dropped the column, so after read_spatialdm_h5ad it is gone. compute_pathway does
    # `geneInter.groupby('pathway_name').interaction_name` and chord_LR does
    # `geneInter.interaction_name_2.str.split('-')` -- both need the columns back.
    if "interaction_name" not in geneInter.columns:
        geneInter["interaction_name"] = geneInter.index.astype(str)
    if "interaction_name_2" not in geneInter.columns:
        L = ligand.apply(lambda r: "_".join(pd.unique(r.dropna().astype(str))), axis=1)
        R = receptor.apply(lambda r: "_".join(pd.unique(r.dropna().astype(str))), axis=1)
        geneInter["interaction_name_2"] = (L.reindex(geneInter.index).astype(str) + "-" +
                                           R.reindex(geneInter.index).astype(str))
    return ligand, receptor, geneInter


diff_utils.concat_db = concat_db_v2

# ------------------------------------------------------------------- load samples
splits = ([s.strip() for s in a.splits.split(",") if s.strip()] or
          sorted([d for d in os.listdir(a.run_dir)
                  if os.path.exists(os.path.join(a.run_dir, d, "data", "spatialdm.h5ad"))],
                 key=lambda x: int(x) if x.isdigit() else 0))
log(f"splits with a persisted object: {splits}")

samples, names, cond_label, n_cells = [], [], [], []
for s in splits:
    sub = sdm.read_spatialdm_h5ad(os.path.join(a.run_dir, s, "data", "spatialdm.h5ad"))
    lab = sub.obs[a.condition_col].astype(str).unique()
    if len(lab) != 1:
        raise SystemExit(f"ERROR: split {s} is not homogeneous in '{a.condition_col}': {lab}")
    # Light copy. `anndata.concat` drops obsp/.raw anyway and nothing in the differential path
    # reads X/obsm; we keep X + spatial so `cdata` stays a usable object, and drop the source
    # h5ad's scanpy baggage (pca/umap/neighbors/rank_genes_groups) and the per-core
    # obsm['celltypes'] dummy frame, whose columns differ between cores and would break concat.
    light = ad.AnnData(X=sub.X, obs=sub.obs[[a.condition_col, "cell_type"]].astype(str).copy(),
                       var=sub.var.copy())
    light.obs["split"] = s
    light.obsm["spatial"] = np.asarray(sub.obsm["spatial"], dtype=float)
    for k in ("ligand", "receptor", "geneInter", "global_res", "global_I", "global_stat",
              "num_pairs", "mean", "single_cell"):
        if k in sub.uns:
            light.uns[k] = sub.uns[k]
    samples.append(light); names.append(s); cond_label.append(lab[0]); n_cells.append(sub.shape[0])
    log(f"  {s}: {sub.shape[0]:,} cells | {a.condition_col}={lab[0]} | "
        f"{len(sub.uns['global_res'])} pairs, {int(sub.uns['global_res'].selected.sum())} selected")
    del sub

cond_label = np.array(cond_label)
bad = set(cond_label) - {a.c1, a.c2}
if bad:
    raise SystemExit(f"ERROR: unexpected {a.condition_col} level(s) {bad}; expected {a.c1}/{a.c2}")
conditions = (cond_label == a.c1).astype(int)
log(f"design: {a.c1}=1 -> {[n for n,c in zip(names,conditions) if c==1]} | "
    f"{a.c2}=0 -> {[n for n,c in zip(names,conditions) if c==0]}")

# --------------------------------------------------------------- concat_obj (authors')
concat = concat_obj(samples, names, "human", "z-score", fdr=a.fdr)
p_df, z_df, tf_df = concat.uns["p_df"], concat.uns["zscore_df"], concat.uns["tf_df"]
log(f"concat: {concat.shape[0]:,} cells x {p_df.shape[0]} union LR pairs across {len(names)} splits")
p_df.to_csv(os.path.join(DATA, "p_df.csv"))
z_df.to_csv(os.path.join(DATA, "zscore_df.csv"))
tf_df.to_csv(os.path.join(DATA, "tf_df.csv"))

# how much of the union is zero-filled? concat_obj substitutes z=0 / p=1 / selected=False for a
# pair a core never tested, so an untested pair is modelled as a *null* pair, not as missing.
tested = pd.DataFrame({n: p_df.index.isin(s.uns["global_res"].index)
                       for n, s in zip(names, samples)}, index=p_df.index)
tested.to_csv(os.path.join(DATA, "tested_mask.csv"))
log(f"union pairs {p_df.shape[0]} | tested in all {len(names)} splits: {int(tested.all(1).sum())} "
    f"| zero-filled cells in zscore_df: {int((~tested).sum().sum())} "
    f"({100*(~tested).values.mean():.1f}%)")

# ---------------------------------------------------------- differential_test (authors')
subset = list(names)
differential_test(concat, subset, conditions)
group_differential_pairs(concat, a.c1, a.c2, diff_quantile1=a.diff_quantile1,
                         diff_quantile2=a.diff_quantile2, fdr_co=a.fdr_co)
c1_spec = list(concat.uns[f"{a.c1}_specific"]); c2_spec = list(concat.uns[f"{a.c2}_specific"])
c1_only = list(concat.uns[f"{a.c1}_only"]);     c2_only = list(concat.uns[f"{a.c2}_only"])
# snapshot NOW: the density sensitivity at the end re-runs differential_test on the same object
# and overwrites uns['p_val'/'diff'/'diff_fdr'].
N_SIG = int((concat.uns["diff_fdr"] < a.fdr_co).sum())
MIN_FDR = float(np.min(concat.uns["diff_fdr"]))
log(f"differential: {N_SIG} pairs at diff FDR < {a.fdr_co} (min diff FDR = {MIN_FDR:.4g}) "
    f"| {a.c1}_specific={len(c1_spec)} ({a.c1}_only={len(c1_only)}) "
    f"| {a.c2}_specific={len(c2_spec)} ({a.c2}_only={len(c2_only)})")

# --------------------------------------------------------------------- results table
gi = concat.uns["geneInter"]
np_origin = set()
_db = os.environ.get("SPATIALDM_DB_CSV", "data/LRdatabase/CellChatDBv2.0.human.csv")
if os.path.exists(_db):
    _d = pd.read_csv(_db)
    _d["interaction_name"] = _d["ligand"].astype(str) + "_" + _d["receptor"].astype(str)
    np_origin = set(_d.loc[_d.signaling_type == "Non-protein Signaling", "interaction_name"])

m1, m0 = np.array(subset)[conditions == 1], np.array(subset)[conditions == 0]
res = pd.DataFrame({
    "diff": concat.uns["diff"], "p_val": concat.uns["p_val"], "diff_fdr": concat.uns["diff_fdr"],
}, index=p_df.index)
res[f"mean_z_{a.c1}"] = z_df[m1].mean(1)
res[f"mean_z_{a.c2}"] = z_df[m0].mean(1)
res[f"n_tested_{a.c1}"] = tested[m1].sum(1);      res[f"n_tested_{a.c2}"] = tested[m0].sum(1)
res[f"n_selected_{a.c1}"] = tf_df[m1].sum(1);     res[f"n_selected_{a.c2}"] = tf_df[m0].sum(1)
res["n_selected_total"] = tf_df.sum(1)
# concat_obj fills z=0 for a pair a split never tested, so a pair that is merely UNTESTABLE in
# one arm (min_cell filter) separates the design perfectly and gets a tiny likelihood-ratio p.
# `group_differential_pairs` guards `_specific` with tf_df.sum(1) in 1..n_sub-1, but the raw
# diff_fdr ranking and the volcano's point cloud do not. Flag it so nobody reads a testability
# difference as a signalling difference.
res["tested_in_all_splits"] = tested.all(1)
res["confounded_by_testability"] = (~tested.all(1)) & (tf_df.sum(1) == 0)
res[f"{a.c1}_specific"] = res.index.isin(c1_spec); res[f"{a.c2}_specific"] = res.index.isin(c2_spec)
res[f"{a.c1}_only"] = res.index.isin(c1_only);     res[f"{a.c2}_only"] = res.index.isin(c2_only)
res["annotation"] = gi["annotation"].reindex(res.index)
res["pathway_name"] = gi["pathway_name"].reindex(res.index)
res["non_protein_origin"] = res.index.isin(np_origin)
res = res.sort_values("diff_fdr")
res.to_csv(os.path.join(DATA, "differential_results.csv"))
for nm, ls in ((f"{a.c1}_specific", c1_spec), (f"{a.c2}_specific", c2_spec),
               (f"{a.c1}_only", c1_only), (f"{a.c2}_only", c2_only)):
    pd.Series(ls, name="interaction_name").to_csv(os.path.join(DATA, f"{nm}.csv"), index=False)

# requested LRs -- reported whatever their rank, per the benchmark's standing request
requested = [s.strip() for s in a.requested_lrs.split(",") if s.strip()]
rq = res.reindex([r for r in requested if r in res.index])
missing = [r for r in requested if r not in res.index]
rq.to_csv(os.path.join(DATA, "requested_lr_differential.csv"))
if missing:
    log(f"  requested LR(s) absent from the union DB entirely: {missing}")

# ------------------------------------------------------------------------- plots
# tutorial cell 20: overview clustermap of 1 - p over every pair x sample
def _cm_all():
    g = sns.clustermap(1 - p_df, figsize=(7, 9), cmap="viridis",
                       xticklabels=True, yticklabels=False)
    g.ax_heatmap.set_xlabel("split"); g.ax_heatmap.set_ylabel("LR pair")
    save_all("clustermap_1_minus_p_all_pairs", g.figure)
guard(_cm_all, "clustermap_1_minus_p_all_pairs")

# tutorial cell 30
def _dend():
    g = pl.differential_dendrogram(concat)
    save_all("differential_dendrogram", g.figure)
guard(_dend, "differential_dendrogram")

# tutorial cell 31
def _volc():
    plt.figure(figsize=(6, 5))
    pl.differential_volcano(concat, legend=[f"{a.c1} specific", f"{a.c2} specific"])
    save_all("differential_volcano")
guard(_volc, "differential_volcano")

# tutorial cell 33 -- same volcano with the standing-request pairs highlighted
def _volc_req():
    plt.figure(figsize=(6, 5))
    pl.differential_volcano(concat, pairs=list(rq.index),
                            legend=[f"{a.c1} specific", f"{a.c2} specific"])
    save_all("differential_volcano_requested", d=REQ)
guard(_volc_req, "differential_volcano_requested")

# tutorial cells 35-36, via the 0.3.1 signature (see module docstring, patch 2).
# `dot_path` visualises the pathway enrichment OF THE `_specific` SETS, so if
# group_differential_pairs returns nothing there is no input -- that is a result, not a failure.
_dic = {k: pd.Index(v) for k, v in ((f"{a.c1}_specific", c1_spec), (f"{a.c2}_specific", c2_spec)) if len(v)}


def _paths():
    dic = _dic
    pr = sdm.compute_pathway(concat, dic=dic)
    pr.to_csv(os.path.join(DATA, "pathway_enrichment_specific.csv"))
    before = set(plt.get_fignums())
    pl.dot_path(pr, cut_off=a.path_cut_off, figsize=(6, 10))
    made = [n for n in plt.get_fignums() if n not in before]
    pats = list(pr[pr.selected >= a.path_cut_off].pattern.unique())
    if not made:
        raise RuntimeError(f"no pathway passed cut_off={a.path_cut_off}")
    for n, pat in zip(made, pats):
        save_all(f"dot_path_{pat}", plt.figure(n))
    log(f"    dot_path panels: {pats}")


if _dic:
    guard(_paths, "dot_path")
else:
    SKIPPED.append(f"dot_path / compute_pathway: both `_specific` sets are empty "
                   f"(no pair reached diff FDR < {a.fdr_co}, and none is selected in every "
                   f"{a.c1} split while absent from every {a.c2} split) -- nothing to enrich")
    log(f"  skipped dot_path: {SKIPPED[-1]}")

# ------------------------------------------- sensitivity: density instead of grade
# The z-score null derives from W, so power scales with neighbourhood size, and core density
# correlates with grade (DEVIATIONS.md: r = 0.659, p = 0.014). Re-running the AUTHORS' OWN test
# with a density-median split instead of grade shows how separable the two designs are here.
sens = {}
try:
    ps = pd.read_csv(os.path.join(a.run_dir, "per_split_summary.csv")).set_index("split")
    dens = ps.loc[[int(n) if n.isdigit() else n for n in names], a.density_col].values
    k = int(conditions.sum())
    cond_d = np.zeros(len(names), int); cond_d[np.argsort(-dens)[:k]] = 1
    swap = [n for n, g, d in zip(names, conditions, cond_d) if g != d]
    log(f"density design: dense=1 -> {[n for n,c in zip(names,cond_d) if c==1]} | "
        f"differs from the {a.condition_col} design at split(s) {swap}")
    differential_test(concat, subset, cond_d)          # overwrites uns; all plots already saved
    group_differential_pairs(concat, "dense", "sparse", diff_quantile1=a.diff_quantile1,
                             diff_quantile2=a.diff_quantile2, fdr_co=a.fdr_co)
    d1, d2 = set(concat.uns["dense_specific"]), set(concat.uns["sparse_specific"])
    jac = lambda x, y: len(x & y) / max(len(x | y), 1)
    sens = dict(density_col=a.density_col, dense_splits=[n for n, c in zip(names, cond_d) if c == 1],
                designs_differ_at=swap, n_dense_specific=len(d1), n_sparse_specific=len(d2),
                jaccard_c1_vs_dense=round(jac(set(c1_spec), d1), 3),
                jaccard_c2_vs_sparse=round(jac(set(c2_spec), d2), 3))
    pd.DataFrame({"diff": concat.uns["diff"], "p_val": concat.uns["p_val"],
                  "diff_fdr": concat.uns["diff_fdr"],
                  "dense_specific": p_df.index.isin(d1),
                  "sparse_specific": p_df.index.isin(d2)}, index=p_df.index
                 ).sort_values("diff_fdr").to_csv(
        os.path.join(DATA, "sensitivity_density_design.csv"))
    log(f"sensitivity: dense_specific={len(d1)} sparse_specific={len(d2)} | "
        f"Jaccard vs {a.c1}/{a.c2} = {sens['jaccard_c1_vs_dense']}/{sens['jaccard_c2_vs_sparse']}")
except Exception as e:
    log(f"  FAILED density sensitivity: {type(e).__name__}: {e}")
    FAILED.append(f"density_sensitivity: {type(e).__name__}: {e}")

# ---------------------------------------------------------------------- manifest
try:
    sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                  stderr=subprocess.DEVNULL).decode().strip()
except Exception:
    sha = None
json.dump({
    "script": "run_diff_spatialdm.py", "dataset": "GBM", "tier": os.path.basename(
        os.path.dirname(os.path.normpath(OUT)) if os.path.basename(os.path.normpath(OUT)) else OUT),
    "run_dir": a.run_dir, "out_dir": OUT,
    "condition_col": a.condition_col, "c1": a.c1, "c2": a.c2,
    "splits": names, "conditions": conditions.tolist(), "n_cells_per_split": n_cells,
    "n_cells_total": int(sum(n_cells)), "n_union_pairs": int(p_df.shape[0]),
    "n_pairs_tested_in_all_splits": int(tested.all(1).sum()),
    "frac_zero_filled": round(float((~tested).values.mean()), 4),
    "concat_obj_fdr": a.fdr, "diff_quantile1": a.diff_quantile1,
    "diff_quantile2": a.diff_quantile2, "fdr_co": a.fdr_co, "path_cut_off": a.path_cut_off,
    "n_diff_fdr_sig": N_SIG, "min_diff_fdr": round(MIN_FDR, 6),
    "n_raw_p_lt_0.05": int((res.p_val < 0.05).sum()),
    "n_raw_p_lt_0.05_expected_by_chance": round(0.05 * len(res), 1),
    f"n_{a.c1}_specific": len(c1_spec), f"n_{a.c2}_specific": len(c2_spec),
    f"n_{a.c1}_only": len(c1_only), f"n_{a.c2}_only": len(c2_only),
    "requested_lrs": requested, "requested_absent_from_union": missing,
    "density_sensitivity": sens,
    "patched": ["diff_utils.concat_db -> v2 union from samples' own geneInter",
                "dot_path called with the installed 0.3.1 signature, not the tutorial's"],
    "spatialdm": sdm.__version__, "anndata": ad.__version__, "seed": a.seed, "git_sha": sha,
    "failed": FAILED, "skipped_no_input": SKIPPED,
    "wall_min": round((time.time() - t0) / 60, 2),
    "peak_rss_gb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2),
}, open(os.path.join(OUT, "run_manifest.json"), "w"), indent=2)

log(f"DONE in {(time.time()-t0)/60:.1f} min -> {OUT} | {len(FAILED)} failure(s)")
if not rq.empty:
    print("\nRequested LRs:")
    print(rq[["diff", "p_val", "diff_fdr", f"mean_z_{a.c1}", f"mean_z_{a.c2}",
              f"n_selected_{a.c1}", f"n_selected_{a.c2}",
              f"{a.c1}_specific", f"{a.c2}_specific"]].to_string())
COLS = ["diff", "diff_fdr", f"n_tested_{a.c1}", f"n_tested_{a.c2}",
        f"n_selected_{a.c1}", f"n_selected_{a.c2}", "annotation", "pathway_name",
        "non_protein_origin"]
print("\nTop 15 by differential FDR (RAW -- includes testability artefacts):")
print(res.head(15)[COLS].to_string())
interp = res[res.tested_in_all_splits & (res.n_selected_total >= 1)]
interp.to_csv(os.path.join(DATA, "differential_results_interpretable.csv"))
print(f"\nTop 15 among the {len(interp)} pairs tested in ALL splits and selected somewhere:")
print(interp.head(15)[COLS].to_string())
