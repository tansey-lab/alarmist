#!/usr/bin/env python
"""LIANA+ inflow -- the DOWNSTREAM half of inflow_score.ipynb that run_inflow.py omitted.

run_inflow.py stopped after `li.mt.inflow` and wrote no figures at all (it created plots/ and
left it empty). The tutorial does not stop there: it goes on to compute source->target
specificity, a region-level variant, cell-type proximity, and a spatially-constrained
rank_aggregate, and it plots each. This script performs those steps and every figure the
tutorial draws.

Tutorial cell -> what is reproduced here
  11/43  sc.pl.embedding(cell_type [, region])            -> global/spatial_celltype*.png
  19/20  li.ut.query_bandwidth + chosen-bandwidth vline   -> global/bandwidth_query.png
  24     li.pl.connectivity                               -> global/connectivity_idx*.png
  45     sc.pl.embedding(interaction) + (ligand,receptor) -> interactions/<LR>_*.png
  46     sender cell-type indicator embedding             -> interactions/<LR>_sender_*.png
  48     sc.pl.violin(groupby=cell_type)                  -> interactions/<LR>_violin.png
  51     li.pl.feature_by_group                           -> interactions/<LR>_feature_by_group.png
  54-57  li.mt.compute_global_specificity + heatmap       -> interactions/<LR>_source_target.png
  59     li.pl.dotplot(uns_key='global_interactions')     -> global/dotplot_global.png
  62/64  compute_global_specificity(groupby=region)       -> global/region_*.png
  66-71  cell_type::region 'niche' dotplot               -> global/niche_dotplot.png
  75-78  li.ut.spatial_pair_proximity + heatmap           -> global/pair_proximity.png
  83-88  li.mt.rank_aggregate(spatial) + dotplot/heatmap  -> global/rank_aggregate*.png

Deviations (also recorded in DEVIATIONS.md):
  * The tutorial's `major_brain_region` has no counterpart here. `--region-col grade`
    (high/low) is used as the tissue-class analogue. obs also carries `motif` and `patch_id`,
    which are ALARMIST outputs -- these are NEVER read by this script.
  * Interactions plotted = the method's own top-N by lr_mean, PLUS the two required
    LRs (GRN^SORT1, ANXA1^FPR1) regardless of rank.

Usage: python run_inflow_downstream.py --in-dir DIR --h5ad PATH --db PATH [options]
"""
import argparse, json, os, time, warnings
import numpy as np, pandas as pd
import scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import liana as li

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

p = argparse.ArgumentParser()
p.add_argument("--in-dir", required=True, help="cellchatdb2_inflow/ from run_inflow.py")
p.add_argument("--h5ad", required=True, help="source h5ad, to rebuild gene-level adata")
p.add_argument("--db", required=True,
               help="'consensus' = LIANA's own resource (the `default` tier); else a path to a "
                    "CSV with ['ligand','receptor'] columns. Must match the DB the upstream "
                    "run_inflow.py used -- rank_aggregate rescores from genes, so a mismatch "
                    "would silently score a different LR space than the inflow output.")
p.add_argument("--cell-type-col", default="cell_type")
p.add_argument("--region-col", default="grade",
               help="analogue of the tutorial's major_brain_region. NEVER 'motif'/'patch_id' "
                    "-- those are ALARMIST outputs.")
p.add_argument("--count-layer", default="counts")
p.add_argument("--bandwidth", type=float, default=13.1454)
p.add_argument("--cutoff", type=float, default=0.1)
p.add_argument("--n-perms", type=int, default=1000)
p.add_argument("--top-n", type=int, default=6)
p.add_argument("--required-lrs", default="GRN^SORT1,ANXA1^FPR1")
p.add_argument("--proximity-bandwidth", type=float, default=None,
               help="tutorial uses 100 um for proximity (a coarser scale than the inflow "
                    "kernel, deliberately). Default = 100.")
p.add_argument("--seed", type=int, default=0)
p.add_argument("--skip-rank-aggregate", action="store_true")
a = p.parse_args()
np.random.seed(a.seed)
PROX_BW = 100.0 if a.proximity_bandwidth is None else a.proximity_bandwidth

OUT = a.in_dir
PLOTS = os.path.join(OUT, "plots")
G, I = os.path.join(PLOTS, "global"), os.path.join(PLOTS, "interactions")
DATA = os.path.join(OUT, "data")
for d in (G, I, DATA): os.makedirs(d, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
SAVED, SKIPPED = [], []


def _save(fig, path, dpi=180):
    """Save, but refuse blank figures -- the SpatialDM lesson (143/156 empty PNGs)."""
    if fig is None:
        SKIPPED.append((os.path.basename(path), "no figure object")); return
    has = any(ax.has_data() or ax.get_images() or ax.collections or ax.patches or ax.texts
              for ax in fig.axes)
    if not fig.axes or not has:
        plt.close(fig); SKIPPED.append((os.path.basename(path), "blank")); return
    fig.savefig(path, dpi=dpi, bbox_inches="tight"); plt.close(fig)
    SAVED.append(os.path.basename(path))


def _cur(path, dpi=180):
    _save(plt.gcf() if plt.get_fignums() else None, path, dpi)


def _fs(ncol=1, per=4.0):
    """Figure size honouring the TMA's tall aspect so maps are not squashed."""
    return (per * ncol, per * ASPECT)


# ------------------------------------------------------------------ load lrdata
log("loading inflow lrdata")
lrdata = sc.read_h5ad(os.path.join(DATA, "inflow_lrdata.h5ad"))
lrdata.raw = None
for k in ("neighbors", "pca", "umap", "rank_genes_groups"):
    lrdata.uns.pop(k, None)
for k in ("X_pca", "X_umap"):
    lrdata.obsm.pop(k, None)
xy = np.asarray(lrdata.obsm["spatial"], float)
ASPECT = float(np.ptp(xy[:, 1]) / np.ptp(xy[:, 0]))
log(f"  {lrdata.shape[0]} cells x {lrdata.shape[1]} features | aspect {ASPECT:.2f}")

CT, REG = a.cell_type_col, a.region_col
for c in (CT, REG):
    if c not in lrdata.obs:
        raise SystemExit(f"STOP: obs has no '{c}'")
for banned in ("motif", "patch_id"):
    if banned in (CT, REG):
        raise SystemExit(f"STOP: '{banned}' is an ALARMIST output and must not be used.")

# ------------------------------------------- rebuild gene-level adata (same recipe)
log("rebuilding gene-level adata (same QC as run_inflow.py)")
adata = sc.read_h5ad(a.h5ad)
if a.count_layer != "X":
    adata.X = adata.layers[a.count_layer].copy()
adata.layers.clear(); adata.raw = None; adata.uns.pop("log1p", None)
adata = adata[:, [g for g in adata.var_names.astype(str) if "_" not in g]].copy()
adata.obs["cell_type"] = adata.obs[a.cell_type_col].astype(str).astype("category")
adata.obsm["spatial"] = np.asarray(adata.obsm["spatial"], float)
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
gm = pd.read_csv(os.path.join(DATA, "gene_moranI.csv"), index_col=0)
svgs = gm.index[(gm["pval_norm_fdr_bh"] < 0.05) & (gm["I"] > 0.01)]
adata = adata[:, [g for g in svgs if g in adata.var_names]].copy()
adata = adata[[c for c in lrdata.obs_names if c in adata.obs_names]].copy()
assert adata.shape[0] == lrdata.shape[0], f"cell mismatch {adata.shape[0]} vs {lrdata.shape[0]}"
li.ut.spatial_neighbors(adata, bandwidth=a.bandwidth, cutoff=a.cutoff, kernel="gaussian",
                        set_diag=True)
log(f"  adata {adata.shape} | graph nnz {adata.obsp['spatial_connectivities'].nnz:,}")

# ============================================================ GLOBAL (cells 11/43)
log("[cell 11/43] spatial embeddings")
for col, cmap in [(CT, None), (REG, None)]:
    sc.pl.embedding(lrdata, basis="spatial", color=col, s=3, frameon=False, show=False)
    f = plt.gcf(); f.set_size_inches(*_fs(1, 4.2))
    for ax in f.axes: ax.set_aspect("equal")
    _cur(os.path.join(G, f"spatial_{col}.png"))

# ------------------------------------------------------- (cells 19/20) bandwidth
log("[cell 19/20] query_bandwidth")
try:
    import plotnine as p9
    plot, dfb = li.ut.query_bandwidth(coordinates=xy, start=5, end=120, interval_n=40)
    dfb.to_csv(os.path.join(DATA, "bandwidth_query.csv"), index=False)
    # query_bandwidth's x-axis is a HARD query RADIUS (see liana/utils/query_bandwidth.py; its
    # column is literally spelled 'bandwith'), not the gaussian sigma. Marking sigma=13.1454 on it
    # compares two different quantities -- it reads as ~2.8 neighbours against the realised median
    # of 14. The comparable point is the kernel's SUPPORT radius, sigma*sqrt(-2 ln cutoff).
    R = a.bandwidth * np.sqrt(-2 * np.log(a.cutoff))
    pl = (plot + p9.geom_vline(xintercept=float(R), linetype="dashed", color="red")
          + p9.annotate("text", x=float(R), y=dfb.neighbours.max(),
                        label=f"support radius = {R:.1f} um   (gaussian sigma = {a.bandwidth})",
                        angle=90, va="bottom", ha="right", color="red"))
    pl.save(os.path.join(G, "bandwidth_query.png"), dpi=180, width=7, height=4.5, verbose=False)
    SAVED.append("bandwidth_query.png")
except Exception as e:
    log(f"  bandwidth_query FAILED: {type(e).__name__}: {e}"); SKIPPED.append(("bandwidth_query.png", str(e)[:60]))

# ---------------------------------------------------------- (cell 24) connectivity
log("[cell 24] connectivity")
for idx in (int(lrdata.shape[0] // 2), int(lrdata.shape[0] // 4)):
    try:
        fig = li.pl.connectivity(adata, idx=idx, size=0.01, figure_size=(6, 5),
                                 spatial_key="spatial", return_fig=True)
        # returns a plotnine ggplot, NOT a matplotlib Figure -> it has no .axes
        fig.save(os.path.join(G, f"connectivity_idx{idx}.png"), dpi=180, width=6, height=5,
                 verbose=False)
        SAVED.append(f"connectivity_idx{idx}.png")
    except Exception as e:
        plt.close("all"); log(f"  connectivity idx={idx} FAILED: {e}")

# ============================================ (cells 54-55) global specificity
log(f"[cell 54] compute_global_specificity(groupby='{CT}', n_perms={a.n_perms})")
li.mt.compute_global_specificity(lrdata, groupby=CT, n_perms=a.n_perms, use_raw=False,
                                 seed=1337, verbose=False)
gi = lrdata.uns["global_interactions"].copy()
gi.to_csv(os.path.join(DATA, "global_interactions.csv"), index=False)
log(f"  global_interactions: {gi.shape[0]:,} source-target-LR rows, "
    f"{(gi.pval < 0.05).sum():,} with p<0.05")

gi["lr"] = gi["ligand_complex"] + "^" + gi["receptor_complex"]
rank = (gi[gi.pval < 0.05].groupby("lr")["lr_mean"].max()
        .sort_values(ascending=False))
required = [s for s in a.required_lrs.split(",") if s]
top_lrs = [l for l in rank.index[:a.top_n]]
for r in required:
    if r not in top_lrs:
        top_lrs.append(r)
log(f"  plotting {len(top_lrs)} LRs: top-{a.top_n} by lr_mean + required {required}")
# NAMING: this is NOT a plain ranking by lr_mean. It is the MAXIMUM lr_mean over the
# source->target pairs that reached p<0.05, so an LR with no significant pair is absent entirely
# (616 of 633 unique LRs survive). Spell that out in the columns rather than in the filename,
# which downstream readers already reference.
(rank.reset_index()
     .rename(columns={"lr_mean": "max_lr_mean_over_signif_pairs"})
     .merge(gi[gi.pval < 0.05].groupby("lr").size().rename("n_signif_pairs").reset_index(),
            on="lr", how="left")
 ).to_csv(os.path.join(DATA, "lr_ranking_by_lr_mean.csv"), index=False)
# companion: the unconditional ranking, so "absent" and "present but not significant" are
# distinguishable without recomputing from the 41k-row table
(gi.groupby("lr")["lr_mean"].max().sort_values(ascending=False)
   .rename("max_lr_mean_all_pairs").reset_index()
 ).to_csv(os.path.join(DATA, "lr_ranking_all_pairs.csv"), index=False)
for r in required:
    pos = list(rank.index).index(r) + 1 if r in rank.index else None
    log(f"    {r}: rank {pos} of {len(rank)}" if pos else f"    {r}: not significant anywhere")

# --------------------------------------------------------- (cell 59) global dotplot
log("[cell 59] dotplot(uns_key='global_interactions')")
try:
    thr = float(gi["lr_mean"].quantile(0.999))
    fig = li.pl.dotplot(adata=lrdata, colour="lr_mean", size="pval", inverse_size=True,
                        uns_key="global_interactions", top_n=25, orderby="lr_mean",
                        orderby_ascending=False, filter_fun=lambda x: x["lr_mean"] > thr,
                        figure_size=(9, 8), return_fig=True)
    fig.save(os.path.join(G, "dotplot_global.png"), dpi=180, width=9, height=8, verbose=False)
    SAVED.append("dotplot_global.png")
except Exception as e:
    plt.close("all"); log(f"  dotplot FAILED: {type(e).__name__}: {e}")

# ================================================= (cells 62/64) region specificity
log(f"[cell 62] compute_global_specificity(groupby='{REG}', uns_key='region_global_interactions')")
sub = lrdata[~lrdata.obs[REG].isna()].copy()
li.mt.compute_global_specificity(sub, groupby=REG, n_perms=a.n_perms, use_raw=False,
                                 seed=1337, uns_key="region_global_interactions", verbose=False)
rgi = sub.uns["region_global_interactions"].copy()
rgi.to_csv(os.path.join(DATA, "region_global_interactions.csv"), index=False)
log(f"  region_global_interactions: {rgi.shape[0]:,} rows")

# =========================================== (cells 66-71) cell_type::region niche
log(f"[cell 66-68] group_labels = {CT}::{REG}")
sub.obs["group_labels"] = sub.obs[CT].astype(str) + "::" + sub.obs[REG].astype(str)
li.mt.compute_global_specificity(sub, groupby="group_labels", n_perms=a.n_perms,
                                 use_raw=False, seed=1337, uns_key="niche_interactions",
                                 verbose=False)
ni = sub.uns["niche_interactions"].copy()
ni[["target", "niche"]] = ni["target"].str.split("::", expand=True)
ni.to_csv(os.path.join(DATA, "niche_interactions.csv"), index=False)


def _heat(df, title, path, value="lr_mean"):
    if df.empty:
        SKIPPED.append((os.path.basename(path), "no rows")); return
    piv = df.pivot_table(index="source", columns="target", values=value)
    ann = (df.pivot_table(index="source", columns="target", values="pval", aggfunc="min")
             .reindex(index=piv.index, columns=piv.columns)
             .map(lambda x: "" if pd.isna(x) else ("**" if x < 0.01 else ("*" if x < 0.05 else ""))))
    n = max(piv.shape)
    fig, ax = plt.subplots(figsize=(1 + 0.62 * piv.shape[1], 1 + 0.62 * piv.shape[0]))
    sns.heatmap(piv, annot=ann, fmt="", square=True, cmap="viridis", ax=ax,
                cbar_kws={"label": value}, annot_kws={"size": 8})
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("target (receiver)"); ax.set_ylabel("source (sender)")
    fig.tight_layout(); _save(fig, path)


# ================================================ (cells 45-57) per-interaction
log("[cells 45-57] per-interaction figures")
for lr in top_lrs:
    lig, rec = lr.split("^", 1)
    tag = lr.replace("^", "-").replace("_", "+")
    feats = [f for f in lrdata.var_names.astype(str) if f.endswith("^" + lr)]
    if not feats:
        log(f"  {lr}: no feature -- skipped"); SKIPPED.append((tag, "no feature")); continue
    star = " [REQUIRED]" if lr in required else ""
    log(f"  {lr}{star}: {len(feats)} sender-typed features")

    # (cell 45a) inflow score map, one panel per sender cell type
    ncol = 3; nrow = int(np.ceil(len(feats) / ncol))
    sc.pl.embedding(lrdata, basis="spatial", color=feats, s=2, cmap="magma", ncols=ncol,
                    frameon=False, show=False, vmin=0, vmax="p99.5", sort_order=True)
    f = plt.gcf(); f.set_size_inches(3.4 * ncol, 3.4 * ASPECT * nrow / 1.6)
    for ax in f.axes: ax.set_aspect("equal")
    _cur(os.path.join(I, f"{tag}_inflow_by_sender.png"))

    # (cell 45b) the ligand and receptor genes themselves
    genes = [g for g in set(lig.split("_") + rec.split("_")) if g in adata.var_names]
    if genes:
        genes = sorted(genes)
        nc = min(3, len(genes)); nr = int(np.ceil(len(genes) / nc))
        sc.pl.embedding(adata, basis="spatial", color=genes, s=2, cmap="magma", ncols=nc,
                        use_raw=False, frameon=False, show=False, vmin=0, vmax="p99.5",
                        sort_order=True)
        f = plt.gcf(); f.set_size_inches(3.4 * nc, 3.4 * ASPECT * nr / 1.6)
        for ax in f.axes: ax.set_aspect("equal")
        _cur(os.path.join(I, f"{tag}_genes.png"))

    # (cell 46) sender cell-type indicator
    # var['mean'] is written by li.mt.inflow; avoids materialising a sparse column view
    top_feat = (lrdata.var.loc[feats, "mean"].astype(float).idxmax()
                if "mean" in lrdata.var.columns else feats[0])
    sender = top_feat.split("^")[0]
    key = f"__sender__{sender}"
    adata.obs[key] = (adata.obs["cell_type"].astype(str) == sender).astype(int)
    sc.pl.embedding(adata, basis="spatial", color=key, s=2, cmap="magma_r", frameon=False,
                    show=False, title=f"sender: {sender}")
    f = plt.gcf(); f.set_size_inches(*_fs(1, 4.2))
    for ax in f.axes: ax.set_aspect("equal")
    _cur(os.path.join(I, f"{tag}_sender_{sender.replace('/','-')}.png"))
    del adata.obs[key]

    # (cell 48) violin by cell type
    try:
        fig, ax = plt.subplots(figsize=(max(6, 0.9 * lrdata.obs[CT].nunique()), 4))
        sc.pl.violin(lrdata, groupby=CT, keys=top_feat, size=0.4, rotation=90, ax=ax,
                     show=False)
        ax.set_title(f"{top_feat}", fontsize=9); fig.tight_layout()
        _save(fig, os.path.join(I, f"{tag}_violin.png"))
    except Exception as e:
        plt.close("all"); log(f"    violin FAILED: {e}")

    # (cell 51) feature_by_group
    try:
        labs = list(pd.Series(lrdata.obs[CT].astype(str)).value_counts().index[:4])
        # returns (fig, ax) -- a matplotlib tuple, unlike connectivity/dotplot which
        # return plotnine ggplots. The three plotting entry points disagree on return type.
        fig, _ax = li.pl.feature_by_group(adata=lrdata, spatial_key="spatial", feature=top_feat,
                                          groupby=CT, percentile_scaling=(1, 97), labels=labs,
                                          show_counts=False, normalize=True, figure_size=(10, 8))
        _save(fig, os.path.join(I, f"{tag}_feature_by_group.png"), dpi=170)
    except Exception as e:
        plt.close("all"); log(f"    feature_by_group FAILED: {type(e).__name__}: {e}")

    # (cell 57) source x target heatmap
    _heat(gi[(gi.ligand_complex == lig) & (gi.receptor_complex == rec)],
          f"{lig} -> {rec}   (cell type)", os.path.join(I, f"{tag}_source_target.png"))
    # (cell 64) same, region level
    _heat(rgi[(rgi.ligand_complex == lig) & (rgi.receptor_complex == rec)],
          f"{lig} -> {rec}   ({REG})", os.path.join(I, f"{tag}_region.png"))

    # (cell 70/71) niche dotplot
    d = ni[(ni.ligand_complex == lig) & (ni.receptor_complex == rec)].copy()
    d = d[d.pval < 0.05]
    if not d.empty:
        d["interaction"] = d["source"] + "->" + d["target"]
        fig, ax = plt.subplots(figsize=(3.0 + 0.9 * d["niche"].nunique(),
                                        2.2 + 0.30 * d["interaction"].nunique()))
        sizes = 18 + 90 * (-np.log10(np.clip(d["pval"], 1e-4, 1))) / 4
        sca = ax.scatter(d["niche"], d["interaction"], c=d["lr_mean"], s=sizes,
                         cmap="viridis", edgecolor="black", linewidth=0.4, alpha=0.9)
        for _, row in d.iterrows():
            s = "**" if row["pval"] < 0.01 else "*"
            ax.text(row["niche"], row["interaction"], s, fontsize=7, ha="center", va="center")
        fig.colorbar(sca, ax=ax, label="lr_mean")
        ax.set_title(f"{lig} -> {rec} by {CT}::{REG}  (p<0.05; * p<0.05, ** p<0.01)", fontsize=9)
        ax.tick_params(axis="x", rotation=45); ax.grid(alpha=.25)
        fig.tight_layout(); _save(fig, os.path.join(I, f"{tag}_niche_dotplot.png"))
    else:
        SKIPPED.append((f"{tag}_niche_dotplot.png", "no p<0.05 rows"))

# ------------------------------------------------- (cells 75-78) pair proximity
log(f"[cell 75-78] spatial_pair_proximity(bandwidth={PROX_BW})")
try:
    prox = li.ut.spatial_pair_proximity(adata, groupby=CT, kernel="gaussian",
                                        bandwidth=PROX_BW, spatial_key="spatial", verbose=False)
    prox.to_csv(os.path.join(DATA, "pair_proximity.csv"), index=False)
    piv = prox.pivot(index="source", columns="target", values="proximity")
    ann = (prox.pivot(index="source", columns="target", values="interacting")
               .fillna(0).astype(int).replace({1: "o", 0: ""}))
    fig, ax = plt.subplots(figsize=(1 + 0.66 * piv.shape[1], 1 + 0.66 * piv.shape[0]))
    sns.heatmap(piv, annot=ann, fmt="", cmap="cividis", square=True, ax=ax,
                cbar_kws={"label": "proximity"}, annot_kws={"size": 8})
    ax.set_title(f"cell-type pair proximity (bandwidth={PROX_BW} um; o = interacting)",
                 fontsize=9)
    fig.tight_layout(); _save(fig, os.path.join(G, "pair_proximity.png"))
except Exception as e:
    plt.close("all"); log(f"  pair_proximity FAILED: {type(e).__name__}: {e}")

# ------------------------------------------- (cells 83-88) spatial rank_aggregate
if not a.skip_rank_aggregate:
    log(f"[cell 83] rank_aggregate(spatial, bandwidth={PROX_BW}) -- this is the slow step")
    try:
        if a.db == "consensus":
            resource = li.rs.select_resource("consensus")[["ligand", "receptor"]]
        else:
            resource = pd.read_csv(a.db)[["ligand", "receptor"]]
        resource = resource.astype(str).drop_duplicates().reset_index(drop=True)
        li.mt.rank_aggregate(adata, groupby=CT, resource=resource, expr_prop=0.01,
                             use_raw=False, verbose=False, n_perms=a.n_perms, seed=1337,
                             spatial_key="spatial",
                             spatial_kwargs={"kernel": "gaussian", "bandwidth": PROX_BW})
        lres = adata.uns["liana_res"].copy()
        lres.to_csv(os.path.join(DATA, "rank_aggregate_liana_res.csv"), index=False)
        log(f"  liana_res: {lres.shape[0]:,} rows")
        fig = li.pl.dotplot(adata=adata, colour="lr_means", size="magnitude_rank",
                            inverse_size=True, top_n=25, orderby="magnitude_rank",
                            orderby_ascending=True, figure_size=(9, 8), return_fig=True)
        fig.save(os.path.join(G, "rank_aggregate_dotplot.png"), dpi=180, width=9, height=8,
                 verbose=False)
        SAVED.append("rank_aggregate_dotplot.png")
        for lr in top_lrs:
            lig, rec = lr.split("^", 1)
            d = lres[(lres.ligand_complex == lig) & (lres.receptor_complex == rec)]
            if d.empty: continue
            tag = lr.replace("^", "-").replace("_", "+")
            piv = d.pivot_table(index="source", columns="target", values="lr_means")
            ann = (d.pivot_table(index="source", columns="target", values="specificity_rank",
                                 aggfunc="min").reindex(index=piv.index, columns=piv.columns)
                    .map(lambda x: "" if pd.isna(x) else ("**" if x < 0.01 else ("*" if x < 0.05 else ""))))
            fig, ax = plt.subplots(figsize=(1 + 0.62 * piv.shape[1], 1 + 0.62 * piv.shape[0]))
            sns.heatmap(piv, annot=ann, fmt="", square=True, cmap="viridis", ax=ax,
                        cbar_kws={"label": "lr_means"}, annot_kws={"size": 8})
            ax.set_title(f"{lig} -> {rec}  (rank_aggregate)", fontsize=10)
            fig.tight_layout(); _save(fig, os.path.join(G, f"rank_aggregate_{tag}.png"))
    except Exception as e:
        plt.close("all"); log(f"  rank_aggregate FAILED: {type(e).__name__}: {e}")

# ------------------------------------------------------------------- persist
lrdata.uns["region_global_interactions"] = rgi
lrdata.uns["niche_interactions"] = ni
lrdata.write(os.path.join(DATA, "inflow_lrdata.h5ad"))

try:
    import subprocess
    _sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True).stdout.strip() or None
except Exception:
    _sha = None
json.dump({"script": "run_inflow_downstream.py", "n_perms": a.n_perms, "seed": a.seed,
           "specificity_seed": 1337,
           "dataset": "GBM", "tier": "cellchatdb2", "h5ad": a.h5ad, "db": os.path.abspath(a.db),
           "cell_type_col": a.cell_type_col,
           "region_col": REG, "proximity_bandwidth": PROX_BW,
           "bandwidth": a.bandwidth, "cutoff": a.cutoff,
           "n_cells": int(lrdata.shape[0]), "n_features": int(lrdata.shape[1]),
           "n_global_interactions": int(gi.shape[0]),
           "n_signif_p05": int((gi.pval < 0.05).sum()),
           # NOT a differential test: compute_global_specificity permutes labels across CELLS, so
           # these p-values are cell-level (n = n_cells) and are pseudoreplicated with respect to
           # the 13 TMA cores. See DEVIATIONS.md.
           "replicate_unit": "cell (NOT tma_id core) -- see DEVIATIONS.md",
           "lrs_plotted": top_lrs, "required_lrs": required,
           "n_figures_saved": len(SAVED), "n_skipped": len(SKIPPED),
           "skipped": SKIPPED, "liana": li.__version__, "git_sha": _sha,
           "wall_min": round((time.time() - t0) / 60, 1)},
          open(os.path.join(OUT, "downstream_manifest.json"), "w"), indent=2, default=str)

log("")
log(f"SAVED {len(SAVED)} figures; SKIPPED {len(SKIPPED)}")
for n, why in SKIPPED: log(f"  skip {n}: {why}")
log(f"DONE in {(time.time()-t0)/60:.1f} min -> {PLOTS}")
