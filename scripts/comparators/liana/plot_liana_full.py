#!/usr/bin/env python
"""Full visualisation suite for the LIANA+ inflow + NMF runs.

Covers the visualisations the authors demonstrate that were missed on the first pass:
  inflow_score.ipynb : cell-type map, per-interaction spatial, ligand & receptor gene maps,
                       sender-cell-type mask, sc.pl.violin by receiver cell type,
                       li.pl.feature_by_group, li.mt.compute_global_specificity + heatmap
  bivariate.ipynb    : li.pl.connectivity
plus FACTOR IDENTITY panels, which answer "what IS this factor" -- top loadings
(<cell_type>^<lig>^<rec>), sender composition, receiver composition, spatial map.

Usage: python plot_liana_full.py --inflow-dir DIR --nmf-dir DIR --h5ad PATH --out-dir DIR
"""
import argparse, json, os, time, warnings
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import liana as li
warnings.filterwarnings("ignore")

p = argparse.ArgumentParser()
p.add_argument("--input", default="inflow", choices=["inflow", "bivariate"])
p.add_argument("--inflow-dir", required=True)
p.add_argument("--nmf-dir", required=True)
p.add_argument("--h5ad", required=True)
p.add_argument("--out-dir", required=True)
p.add_argument("--count-layer", default="counts")
p.add_argument("--cell-type-col", default="cell_type")
p.add_argument("--n-top-interactions", type=int, default=5)
p.add_argument("--global-specificity", action="store_true", default=False,
               help="run li.mt.compute_global_specificity ONLY if the canonical copy written by "
                    "run_inflow_downstream.py is absent; see --reuse-global-specificity")
p.add_argument("--reuse-global-specificity", action=argparse.BooleanOptionalAction, default=True,
               help="prefer <inflow-dir>/data/global_interactions.csv (n_perms=1000, seed=1337) "
                    "over recomputing here. Recomputing produced a SECOND, unseeded n_perms=100 "
                    "copy that disagreed with the canonical one on 393 rows at p<0.05.")
p.add_argument("--n-perms", type=int, default=1000)
p.add_argument("--seed", type=int, default=1337,
               help="seed for compute_global_specificity; the earlier pass left it unset")
p.add_argument("--required-lrs", default="GRN^SORT1,ANXA1^FPR1",
               help="project standing rule (comparator-benchmark SKILL.md): always plot these "
                    "alongside the method's own top hits, whatever their rank, into a SEPARATE "
                    "'requested/' directory so they are never mistaken for the method's ranking")
p.add_argument("--bandwidth", type=float, default=13.1454,
               help="must match the run being plotted; used by li.pl.connectivity")
p.add_argument("--cutoff", type=float, default=0.1)
a = p.parse_args()
HAS_SENDER = a.input == "inflow"   # inflow features are '<cell_type>^<lig>^<rec>'

OUT = a.out_dir
D_INT = os.path.join(OUT, "interactions"); D_FAC = os.path.join(OUT, "factors")
D_GLB = os.path.join(OUT, "global"); D_REQ = os.path.join(D_INT, "requested")
for d in (D_INT, D_FAC, D_GLB, D_REQ): os.makedirs(d, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
BLANK = []


def _is_blank(fig):
    """A figure with no drawn content. li.pl.connectivity used to slip through here: it draws
    with PLOTNINE, so plt.gcf() was an empty matplotlib canvas and 6.4 KB of white PNG was
    written with no exception raised. Same detector as run_inflow_downstream.py::_save."""
    if fig is None or not fig.axes:
        return True
    return not any(ax.has_data() or ax.get_images() or ax.collections or ax.patches or ax.texts
                   for ax in fig.axes)


def save(fig, path):
    fig = fig if fig is not None else (plt.gcf() if plt.get_fignums() else None)
    if _is_blank(fig):
        plt.close("all"); BLANK.append(os.path.basename(path))
        log(f"  BLANK, not written: {os.path.basename(path)}")
        return
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close("all")


def save_gg(gg, path, width=6, height=6):
    """plotnine ggplots have their own .save(); they never touch the matplotlib canvas."""
    gg.save(path, dpi=180, width=width, height=height, verbose=False)


def guard(fn, name):
    try: fn(); log(f"  ok: {name}")
    except Exception as e: plt.close("all"); log(f"  FAILED {name}: {type(e).__name__}: {e}")

# ------------------------------------------------------------------ load
log("loading inflow lrdata + NMF")
meta = pd.read_csv(os.path.join(a.inflow_dir, "data", "cell_meta.csv"))
if HAS_SENDER:
    lrdata = sc.read_h5ad(os.path.join(a.inflow_dir, "data", "inflow_lrdata.h5ad"))
else:
    _z = np.load(os.path.join(a.inflow_dir, "data", "local_scores.npz"), allow_pickle=True)
    lrdata = ad.AnnData(X=np.asarray(_z["values"], dtype=np.float32),
                        obs=pd.DataFrame(index=_z["cells"].astype(str)),
                        var=pd.DataFrame(index=_z["pairs"].astype(str)))
lrdata.obs["cell_type"] = pd.Categorical(meta["cell_type"].astype(str).values)
if "spatial" not in lrdata.obsm:
    lrdata.obsm["spatial"] = meta[["x", "y"]].values.astype(float)
xy = np.asarray(lrdata.obsm["spatial"], float)

z = np.load(os.path.join(a.nmf_dir, "data", "nmf_WH.npz"), allow_pickle=True)
W, H, feats = z["W"], z["H"], z["features"].astype(str)
loadings = pd.read_csv(os.path.join(a.nmf_dir, "data", "NMF_H_loadings.csv"), index_col=0)

# --------------------------------------------------------------- ROW-ORDER CONTRACT
# li.ut.get_variable_loadings SORTS its rows by |Factor1| descending, so NMF_H_loadings.csv is
# NOT in nmf_WH.npz['features'] order. Every annotation vector below (sender_of, lr_of) is derived
# from `feats`, and every value vector from `loadings` -- so masking one by the other silently
# pairs a feature's NAME with a different feature's LOADING. That is what corrupted the SENDER and
# LR-pair axes of the earlier plots_full trees (0/2704 positions agreed on inflow, 1/131 on
# bivariate). Reindex once, here, and assert; do not rely on the two happening to align.
_missing = set(feats) - set(loadings.index.astype(str))
if _missing:
    raise SystemExit(f"STOP: {len(_missing)} npz features absent from NMF_H_loadings.csv "
                     f"(e.g. {sorted(_missing)[:3]}) -- the two files are not from the same run.")
loadings.index = loadings.index.astype(str)
loadings = loadings.loc[list(feats)]
assert (np.asarray(loadings.index) == feats).all(), "loadings/feats still misaligned"
factors = list(loadings.columns)
log(f"  lrdata {lrdata.shape} | W {W.shape} | H {H.shape} | {len(factors)} factors")
log(f"  loadings reindexed onto nmf_WH.npz feature order ({len(feats)} features)")

# lrdata carries the UNFILTERED feature space; the NMF ran on whatever survived run_nmf.py's
# cross-punch filter. Keep both, but never let a plot mix them.
NMF_FEATS = set(feats)
if lrdata.shape[1] != len(feats):
    log(f"  note: lrdata has {lrdata.shape[1]} features, NMF used {len(feats)} "
        f"(run_nmf.py cross-punch filter) -- top interactions are chosen from the NMF set")

# raw expression for ligand/receptor gene maps (tutorial plots these next to the interaction)
adata = sc.read_h5ad(a.h5ad)
if a.count_layer in adata.layers: adata.X = adata.layers[a.count_layer].copy()
adata.layers.clear(); adata.raw = None; adata.uns.pop("log1p", None)
adata = adata[:, [g for g in adata.var_names.astype(str) if "_" not in g]].copy()
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
adata = adata[lrdata.obs_names].copy()
adata.obsm["spatial"] = xy
adata.obs["cell_type"] = lrdata.obs["cell_type"].values
log(f"  expression aligned: {adata.shape}")

ASPECT = (xy[:, 1].max() - xy[:, 1].min()) / (xy[:, 0].max() - xy[:, 0].min())   # ~2.2, tall
DOT = 0.25          # 100k cells: keep dots sub-pixel-ish so structure is visible, not blobs
def panel_size(ncol, width_per=3.0):
    """Figure size that RESPECTS the tissue aspect ratio -- never squash to a pancake."""
    return (width_per * ncol, width_per * ASPECT)
def scat(ax, v, title, cmap="magma", q=0.99):
    v = np.asarray(v).ravel()
    ax.scatter(xy[:, 0], xy[:, 1], c=v, s=DOT, cmap=cmap, linewidths=0, rasterized=True,
               vmax=np.quantile(v, q) if np.any(v > 0) else None)
    ax.set_title(title, fontsize=8); ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(xy[:, 0].min(), xy[:, 0].max()); ax.set_ylim(xy[:, 1].min(), xy[:, 1].max())

# ----------------------------------------------------- 0. cell type map
def _ct():
    cats = list(lrdata.obs["cell_type"].cat.categories)
    fig, ax = plt.subplots(figsize=(4.2, 4.2 * ASPECT))
    for i, c in enumerate(cats):
        mk = (lrdata.obs["cell_type"] == c).values
        ax.scatter(xy[mk, 0], xy[mk, 1], s=DOT, linewidths=0, rasterized=True,
                   label=f"{c} ({mk.sum():,})",
                   color=plt.cm.tab20(i % 20))
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title("cell types")
    ax.legend(markerscale=12, fontsize=7, loc="center left", bbox_to_anchor=(1, .5))
    save(fig, os.path.join(D_GLB, "cell_types.png"))
guard(_ct, "cell_types")

# --------- 1. per top interaction: inflow + ligand + receptor + sender mask
X = lrdata.X.toarray() if hasattr(lrdata.X, "toarray") else np.asarray(lrdata.X)
_col = "total_inflow" if HAS_SENDER else "total_local_cosine"   # bivariate sums a COSINE, not inflow
strength = pd.Series(np.asarray(X.sum(0)).ravel(), index=lrdata.var_names).sort_values(ascending=False)
pd.DataFrame({_col: strength,
              "in_nmf": [f in NMF_FEATS for f in strength.index]}
             ).to_csv(os.path.join(D_GLB, "interaction_total_inflow.csv"))
# Rank within the NMF's own feature space -- this tree describes the NMF fit, and picking a top
# interaction that the cross-punch filter had already removed would be incoherent.
top_int = [f for f in strength.index if f in NMF_FEATS][:a.n_top_interactions]
log(f"top interactions by {_col} (within the NMF feature set): {top_int}")

# ---- the project's standing required-LR rule (comparator-benchmark SKILL.md) ----
# Always render the two ALARMIST motif-1 arms, whatever their rank, into a SEPARATE directory.
required_lrs = [s.strip() for s in a.required_lrs.split(",") if s.strip()]
req_feats, req_report = [], {}
_lr_of_feat = {f: ("^".join(f.split("^")[1:]) if HAS_SENDER else f) for f in feats}
for lr in required_lrs:
    hits = [f for f in feats if _lr_of_feat[f] == lr]
    if not hits:
        req_report[lr] = None
        log(f"  REQUIRED {lr}: ABSENT from the NMF feature set "
            f"(filtered upstream, or not in the LR resource) -- that is itself a result")
        continue
    hits = sorted(hits, key=lambda f: -float(strength.get(f, 0.0)))
    ranks = {f: int(list(strength.index).index(f)) + 1 for f in hits}
    req_report[lr] = {"n_features": len(hits), "best_feature": hits[0],
                      "best_rank_of": int(len(strength)), "ranks": ranks}
    req_feats += hits
    log(f"  REQUIRED {lr}: {len(hits)} feature(s), best {hits[0]} "
        f"at rank {ranks[hits[0]]}/{len(strength)} by {_col}")
# plot at most the strongest few per LR, but never zero
req_plot = []
for lr in required_lrs:
    r = req_report.get(lr)
    if r:
        req_plot += sorted(r["ranks"], key=lambda f: r["ranks"][f])[:3]

def _tag(inter):
    return inter.replace("^", "-").replace("/", "-")


def _interaction_panel(inter, dest=None):
    sender, lig, rec = inter.split("^")
    fig, ax = plt.subplots(1, 4, figsize=panel_size(4, 2.8))
    scat(ax[0], X[:, list(lrdata.var_names).index(inter)], f"inflow\n{inter}")
    for k, g in enumerate((lig, rec), start=1):
        sub = [x for x in g.split("_") if x in adata.var_names]
        v = np.asarray(adata[:, sub].X.mean(1)).ravel() if sub else np.zeros(len(xy))
        scat(ax[k], v, f"{'ligand' if k==1 else 'receptor'} {g}" + ("" if sub else " (absent)"),
             cmap="viridis")
    m = (lrdata.obs["cell_type"] == sender).values.astype(float)
    ax[3].scatter(xy[~m.astype(bool), 0], xy[~m.astype(bool), 1], s=DOT, c="#dddddd",
                  linewidths=0, rasterized=True)
    ax[3].scatter(xy[m.astype(bool), 0], xy[m.astype(bool), 1], s=DOT * 2, c="crimson",
                  linewidths=0, rasterized=True)
    ax[3].set_xlim(xy[:, 0].min(), xy[:, 0].max()); ax[3].set_ylim(xy[:, 1].min(), xy[:, 1].max())
    ax[3].set_title(f"sender: {sender} (n={int(m.sum()):,})", fontsize=8)
    ax[3].set_aspect("equal"); ax[3].axis("off")
    fig.suptitle(inter, fontsize=11)
    save(fig, os.path.join(dest or D_INT, f"panel_{_tag(inter)}.png"))


# -------- 2. violin: which RECEIVER cell types get this interaction (tutorial)
def _violin(inter, dest=None):
    fig, ax = plt.subplots(figsize=(11, 4))
    sc.pl.violin(lrdata, groupby="cell_type", keys=inter, size=0.4, rotation=90, ax=ax, show=False)
    ax.set_title(f"{inter} — {'inflow' if HAS_SENDER else 'local cosine'} by RECEIVER cell type",
                 fontsize=9)
    save(fig, os.path.join(dest or D_INT, f"violin_{_tag(inter)}.png"))


# --------------------------- 3. li.pl.feature_by_group (tutorial function)
def _fbg(inter, dest=None):
    cats = list(lrdata.obs["cell_type"].cat.categories)[:4]
    # returns a matplotlib figure (unlike connectivity/dotplot, which return plotnine ggplots)
    _r = li.pl.feature_by_group(adata=lrdata, spatial_key="spatial", feature=inter,
                                groupby="cell_type", percentile_scaling=(1, 97),
                                labels=cats, show_counts=False, normalize=True,
                                figure_size=(10, 8))
    _f = _r[0] if isinstance(_r, tuple) else _r
    save(_f, os.path.join(dest or D_INT, f"feature_by_group_{_tag(inter)}.png"))


def _interaction_suite(inters, dest, label):
    for it in inters:
        if HAS_SENDER:
            guard(lambda it=it: _interaction_panel(it, dest), f"{label} panel {it}")
        guard(lambda it=it: _violin(it, dest), f"{label} violin {it}")
        # was previously drawn for top_int[0] only -- now one per interaction, as the tutorial does
        guard(lambda it=it: _fbg(it, dest), f"{label} feature_by_group {it}")


_interaction_suite(top_int, D_INT, "top")
if req_plot:
    _interaction_suite(req_plot, D_REQ, "REQUIRED")
    pd.DataFrame([{"required_lr": lr, **({"status": "absent"} if req_report[lr] is None else
                  {"status": "present", "n_features": req_report[lr]["n_features"],
                   "best_feature": req_report[lr]["best_feature"],
                   "best_rank": req_report[lr]["ranks"][req_report[lr]["best_feature"]],
                   "of_n_features": req_report[lr]["best_rank_of"]})}
                  for lr in required_lrs]).to_csv(
        os.path.join(D_REQ, "requested_lr_ranks.csv"), index=False)
else:
    log("  no required LR present in the NMF feature set -- 'requested/' left empty")

# --------------------------------- 4. li.pl.connectivity (bivariate tutorial)
def _conn():
    tmp = adata.copy()
    li.ut.spatial_neighbors(tmp, bandwidth=a.bandwidth, cutoff=a.cutoff, kernel="gaussian",
                            set_diag=True)
    # li.pl.connectivity draws with PLOTNINE. Without return_fig=True the matplotlib canvas stays
    # empty and plt.gcf() writes a blank white PNG with no exception -- the bug that produced the
    # two 6.4 KB blanks in the earlier plots_full trees. run_inflow_downstream.py:167 does it right.
    gg = li.pl.connectivity(tmp, idx=int(np.argmin(np.abs(xy - xy.mean(0)).sum(1))),
                            size=0.6, figure_size=(6, 6), spatial_key="spatial", return_fig=True)
    save_gg(gg, os.path.join(D_GLB, "connectivity.png"))
guard(_conn, "connectivity")

# ============================ FACTOR IDENTITY PANELS ======================
# what IS each factor: top loadings, sender composition, receiver composition, spatial map
sender_of = (np.array([f.split("^")[0] for f in feats]) if HAS_SENDER
             else np.array(["(none)"] * len(feats)))
lr_of = (np.array(["^".join(f.split("^")[1:]) for f in feats]) if HAS_SENDER
         else feats.copy())
ct = lrdata.obs["cell_type"].values

send_mat = pd.DataFrame(index=sorted(set(sender_of)), columns=factors, dtype=float)
recv_mat = pd.DataFrame(index=list(lrdata.obs["cell_type"].cat.categories), columns=factors, dtype=float)
for j, f in enumerate(factors):
    h = loadings[f].values
    for s in send_mat.index: send_mat.loc[s, f] = h[sender_of == s].sum()
    for r in recv_mat.index: recv_mat.loc[r, f] = W[ct == r, j].mean()
# bivariate features carry no sender, so send_mat degenerates to a single '(none)' row whose
# values are just the per-factor column sums. Writing it invited the reader to interpret a
# one-row "SENDER composition" that does not exist -- emit it only where it means something.
if HAS_SENDER:
    send_mat.to_csv(os.path.join(D_FAC, "factor_by_SENDER_celltype.csv"))
recv_mat.to_csv(os.path.join(D_FAC, "factor_by_RECEIVER_celltype.csv"))

def _factor_panel(j, f):
    top = loadings[f].sort_values(ascending=False).head(10)
    ncell = 4 if HAS_SENDER else 3          # bivariate has no sender axis -> omit that panel
    fig = plt.figure(figsize=(15 if HAS_SENDER else 11.5, 4.0 * max(1.0, ASPECT * 0.62)))
    gs = fig.add_gridspec(1, ncell,
                          width_ratios=[2.0, 1, 1, 1.25] if HAS_SENDER else [2.0, 1, 1.25])
    ax0 = fig.add_subplot(gs[0]); ax0.barh(range(len(top))[::-1], top.values, color="#4C72B0")
    ax0.set_yticks(range(len(top))[::-1]); ax0.set_yticklabels(top.index, fontsize=7)
    ax0.set_xlabel("NMF_H loading"); ax0.set_title(f"{f}: top-10 features", fontsize=9)
    s = send_mat[f].sort_values(ascending=False)
    if HAS_SENDER:
        ax1 = fig.add_subplot(gs[1])
        ax1.barh(range(len(s))[::-1], s.values, color="#DD8452")
        ax1.set_yticks(range(len(s))[::-1]); ax1.set_yticklabels(s.index, fontsize=7)
        ax1.set_title("SENDER composition\n(sum of loadings)", fontsize=8)
    ax2 = fig.add_subplot(gs[2 if HAS_SENDER else 1]); r = recv_mat[f].sort_values(ascending=False)
    ax2.barh(range(len(r))[::-1], r.values, color="#55A868")
    ax2.set_yticks(range(len(r))[::-1]); ax2.set_yticklabels(r.index, fontsize=7)
    ax2.set_title("RECEIVER composition\n(mean NMF_W)", fontsize=8)
    ax3 = fig.add_subplot(gs[3 if HAS_SENDER else 2]); scat(ax3, W[:, j], "spatial (NMF_W)")
    lbl = (f"top sender: {s.index[0]}   " if HAS_SENDER else "")
    fig.suptitle(f"{f}  —  {lbl}top receiver: {r.index[0]}", fontsize=11)
    fig.tight_layout()
    save(fig, os.path.join(D_FAC, f"identity_{f}.png"))
for j, f in enumerate(factors): guard(lambda j=j, f=f: _factor_panel(j, f), f"identity {f}")

def _heat(df, name, title):
    fig, ax = plt.subplots(figsize=(1.0 * len(df.columns) + 3, 0.42 * len(df) + 2.2))
    sns.heatmap(df.astype(float), cmap="magma", ax=ax, cbar_kws=dict(shrink=.6))
    ax.set_title(title, fontsize=10); fig.tight_layout()
    save(fig, os.path.join(D_FAC, name))
if HAS_SENDER:
    guard(lambda: _heat(send_mat, "heatmap_factor_by_SENDER.png",
                        "NMF factor x SENDER cell type (summed H loadings)"), "sender heatmap")
else:
    log("  skip sender heatmap: bivariate features carry no sender cell type")
guard(lambda: _heat(recv_mat, "heatmap_factor_by_RECEIVER.png",
                    "NMF factor x RECEIVER cell type (mean W)"), "receiver heatmap")

# factor x LR-pair (collapsed over sender) -- which LRIs define each factor
lr_mat = pd.DataFrame(0.0, index=sorted(set(lr_of)), columns=factors)
for f in factors:
    h = loadings[f].values
    for lr in lr_mat.index: lr_mat.loc[lr, f] = h[lr_of == lr].sum()
lr_mat.to_csv(os.path.join(D_FAC, "factor_by_LRpair.csv"))
def _lrheat():
    keep = lr_mat.loc[lr_mat.max(1).sort_values(ascending=False).index[:30]]
    _heat(keep, "heatmap_factor_by_LRpair_top30.png",
          "NMF factor x LR pair (top 30, summed over senders)")
guard(_lrheat, "LR-pair heatmap")

# ============ ALARMIST-style summaries (our additions; see DEVIATIONS.md) ============
# Mimics of alarmist.plotting.motif_plots.plot_celltype_communication_by_motif and
# plot_top_lri_interactions_dot, applied to NMF factors instead of BPTF motifs. These are
# OUR plots -- LIANA demonstrates neither. No new analysis: both read NMF_W / NMF_H only.

def _celltype_communication_by_factor():
    """One sender x receiver heatmap per factor (mimics plot_celltype_communication_by_motif).

    An NMF factor is rank-1, so its (sender, receiver) contribution is exactly the outer product
    of the sender loading profile (from H, via the feature's sender tag) and the receiver profile
    (from W, aggregated over each cell's own type). Stated here because it means the heatmap
    carries no interaction beyond those two margins -- unlike ALARMIST, whose factor values are
    indexed by (celltype1, celltype2) directly.
    """
    if not HAS_SENDER:
        raise RuntimeError("bivariate features carry no sender cell type -- not applicable")
    ncol = min(4, len(factors)); nrow = int(np.ceil(len(factors) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.4 * nrow), squeeze=False)
    axes = axes.ravel()
    for j, f in enumerate(factors):
        M = np.outer(recv_mat[f].values.astype(float), send_mat[f].values.astype(float))
        df = pd.DataFrame(M, index=recv_mat.index, columns=send_mat.index)
        sns.heatmap(df, cmap="Blues", ax=axes[j], cbar_kws=dict(shrink=.6))
        axes[j].set_title(f, fontsize=9)
        axes[j].set_xlabel("sender", fontsize=8); axes[j].set_ylabel("receiver", fontsize=8)
        axes[j].tick_params(labelsize=6)
    for k in range(len(factors), len(axes)): axes[k].axis("off")
    fig.suptitle("cell-type communication per NMF factor  (receiver x sender, rank-1 outer product)",
                 fontsize=11)
    fig.tight_layout()
    save(fig, os.path.join(D_FAC, "celltype_communication_by_factor.png"))
guard(_celltype_communication_by_factor, "celltype_communication_by_factor")

def _top_lri_dot():
    """Lollipop of the top LR pairs per factor (mimics plot_top_lri_interactions_dot).

    bivariate: bare LR pairs, single colour. inflow: dot coloured by SENDER cell type.
    """
    top_n = 20
    ncol = min(3, len(factors)); nrow = int(np.ceil(len(factors) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 0.30 * top_n * nrow + 1.2),
                             squeeze=False)
    axes = axes.ravel()
    senders = sorted(set(sender_of)) if HAS_SENDER else []
    cmap = {c: plt.cm.tab20(i % 20) for i, c in enumerate(senders)}
    for j, f in enumerate(factors):
        top = loadings[f].sort_values(ascending=False).head(top_n)[::-1]
        ax = axes[j]; ypos = np.arange(len(top))
        cols = ([cmap[t.split("^")[0]] for t in top.index] if HAS_SENDER else "#4C72B0")
        ax.hlines(ypos, 0, top.values, color="#BBBBBB", lw=1.1, zorder=1)
        ax.scatter(top.values, ypos, s=42, c=cols, zorder=2, edgecolors="white", linewidths=.5)
        lbl = ([t.split("^", 1)[1] for t in top.index] if HAS_SENDER else list(top.index))
        ax.set_yticks(ypos); ax.set_yticklabels(lbl, fontsize=6.5)
        ax.set_xlabel("NMF_H loading", fontsize=8); ax.set_title(f, fontsize=9)
        ax.grid(axis="x", alpha=.25)
    for k in range(len(factors), len(axes)): axes[k].axis("off")
    if HAS_SENDER:
        from matplotlib.lines import Line2D
        fig.legend(handles=[Line2D([0], [0], marker="o", ls="", mfc=cmap[c], mec="white",
                                   label=c, ms=7) for c in senders],
                   loc="lower center", ncol=min(9, len(senders)), fontsize=7,
                   frameon=False, bbox_to_anchor=(0.5, -0.012))
    fig.suptitle(f"top {top_n} LR interactions per NMF factor"
                 + ("  (dot colour = sender cell type)" if HAS_SENDER else "  (bare LR pairs)"),
                 fontsize=11)
    fig.tight_layout(rect=(0, 0.03 if HAS_SENDER else 0, 1, 0.98))
    save(fig, os.path.join(D_FAC, "top_lri_dot_by_factor.png"))
guard(_top_lri_dot, "top_lri_dot_by_factor")

# ------------------------- global specificity (sender x receiver per LRI)
# PREFER the canonical table written by run_inflow_downstream.py (n_perms=1000, seed=1337).
# Recomputing here previously wrote a SECOND global_interactions.csv at an unseeded n_perms=100
# that disagreed with the canonical one on 393 rows at p<0.05 -- two files, same name, different
# statistics, neither recording its n_perms.
GS_SRC, gi = None, None
_canon = os.path.join(a.inflow_dir, "data", "global_interactions.csv")
if HAS_SENDER and a.reuse_global_specificity and os.path.exists(_canon):
    gi = pd.read_csv(_canon)
    GS_SRC = f"reused {os.path.relpath(_canon)} (n_perms=1000, seed=1337)"
    log(f"  global specificity: {GS_SRC}, {gi.shape[0]:,} rows")
elif a.global_specificity and HAS_SENDER:
    def _compute_gs():
        global gi, GS_SRC
        li.mt.compute_global_specificity(lrdata, groupby="cell_type", n_perms=a.n_perms,
                                         seed=a.seed, use_raw=False, verbose=True)
        gi = lrdata.uns["global_interactions"].copy()
        gi.to_csv(os.path.join(D_GLB, "global_interactions.csv"), index=False)
        GS_SRC = f"computed here (n_perms={a.n_perms}, seed={a.seed})"
        log(f"  global_interactions columns: {list(gi.columns)}")
    guard(_compute_gs, "global_specificity")

if gi is not None:
    def _st_heat(inter, dest):
        parts = inter.split("^")
        lig, rec = (parts[1], parts[2]) if HAS_SENDER else (parts[0], parts[1])
        sub = gi[(gi["ligand_complex"] == lig) & (gi["receptor_complex"] == rec)]
        if not len(sub):
            log(f"  no global_interactions rows for {lig}->{rec}"); return
        val = next(c for c in ("lr_mean", "specificity") if c in sub)
        piv = sub.pivot_table(index="source", columns="target", values=val)
        ann = (sub.pivot_table(index="source", columns="target", values="pval", aggfunc="min")
                  .reindex(index=piv.index, columns=piv.columns)
                  .map(lambda x: "" if pd.isna(x) else ("**" if x < 0.01 else ("*" if x < 0.05 else "")))
               if "pval" in sub else None)
        fig, ax = plt.subplots(figsize=(1 + 0.62 * piv.shape[1], 1 + 0.62 * piv.shape[0]))
        sns.heatmap(piv, annot=ann, fmt="", cmap="magma", square=True, ax=ax,
                    cbar_kws={"label": val}, annot_kws={"size": 8})
        ax.set_title(f"{lig} -> {rec}: sender x receiver", fontsize=9)
        ax.set_xlabel("target (receiver)"); ax.set_ylabel("source (sender)")
        fig.tight_layout()
        save(fig, os.path.join(dest, f"sender_receiver_{_tag(inter)}.png"))
    guard(lambda: _st_heat(top_int[0], D_GLB), "sender_receiver heatmap (top)")
    for it in req_plot:
        guard(lambda it=it: _st_heat(it, D_REQ), f"sender_receiver heatmap {it}")

# count real artefacts only -- a bare walk counts macOS .DS_Store and inflates n_files by one
# per directory, which then gets quoted as a figure count
n_png = sum(1 for _, _, fs in os.walk(OUT) for f in fs if f.endswith(".png"))
n = sum(1 for _, _, fs in os.walk(OUT) for f in fs if not f.startswith("."))
try:
    import subprocess
    _sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                          text=True).stdout.strip() or None
except Exception:
    _sha = None
json.dump({"script": "plot_liana_full.py", "input": a.input, "dataset": "GBM", "tier": "cellchatdb2",
           "inflow_dir": a.inflow_dir, "nmf_dir": a.nmf_dir, "h5ad": a.h5ad, "out_dir": OUT,
           "n_factors": len(factors), "n_features_nmf": int(len(feats)),
           "n_features_lrdata": int(lrdata.shape[1]), "n_cells": int(lrdata.shape[0]),
           "top_interactions": top_int, "required_lrs": required_lrs,
           "required_lr_report": req_report, "required_features_plotted": req_plot,
           "global_specificity": GS_SRC, "bandwidth": a.bandwidth, "cutoff": a.cutoff,
           "seed": a.seed, "blank_figures_suppressed": BLANK,
           "loadings_reindexed_to_npz_order": True,
           "n_files": n, "n_png": n_png, "liana": li.__version__, "git_sha": _sha,
           "wall_min": round((time.time() - t0) / 60, 1)},
          open(os.path.join(OUT, "plot_manifest.json"), "w"), indent=2, default=str)
log(f"DONE in {(time.time()-t0)/60:.1f} min -> {OUT} ({n} files); "
    f"{len(BLANK)} blank figure(s) suppressed")
