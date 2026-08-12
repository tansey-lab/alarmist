#!/usr/bin/env python
"""COMMOT downstream impact: ct.tl.communication_impact on the GBM cores.

The tutorial's downstream chain (visium-mouse_brain.ipynb) is

    df_deg, df_yhat = ct.tl.communication_deg_detection(adata, database_name, pathway_name, summary='receiver')
    df_deg_clus, df_yhat_clus = ct.tl.communication_deg_clustering(df_deg, df_yhat, deg_clustering_res=0.4)
    top_de_genes = ct.pl.plot_communication_dependent_genes(df_deg_clus, df_yhat_clus,
                       top_ngene_per_cluster=5, return_genes=True)
    df_impact = ct.tl.communication_impact(adata, database_name, pathway_name=..., tree_combined=True,
                    method='treebased_score', tree_ntrees=100, tree_repeat=100, tree_method='rf',
                    ds_genes=top_de_genes, bg_genes=500, normalize=True)
    ct.pl.plot_communication_impact(df_impact, summary='receiver', ...)

**Step 1 is unavailable on this machine.** `communication_deg_detection` imports rpy2 + anndata2ri
and calls R's tradeSeq; the docstring pins tradeSeq 1.0.1 / R 3.6.3 (2020). Neither rpy2 nor
anndata2ri is installed, and none of the four R installs on this box has tradeSeq (R 4.4.2 at
/usr/local/bin, R 4.3.3 in comp-{cellchat,niches,cytosignal}). So `ds_genes` cannot come from the
tutorial's source. `communication_impact`'s own docstring sanctions the substitute: "A list of
genes ... for example, the highly variable genes." We therefore run two ds_genes variants:

  --ds-source hvg       top-N highly variable genes  (the docstring's substitute; method-faithful)
  --ds-source alarmist  top-N genes from an ALARMIST motif impact table, so the two methods are
                        asked about the SAME genes and the answers are directly comparable

Everything else follows the tutorial's argument values verbatim.

NOTE: `communication_impact` reads `adata.raw` (`adata_bg = adata.raw.to_adata()`), but
run_commot.py dropped `.raw` before writing adata_commot.h5ad to halve its size. This script
rebuilds `.raw` from the source h5ad's layers['counts'] -- which is exactly what `.raw` held
during the run -- so no OT re-run is needed.

Usage: python run_commot_impact.py --run-dir DIR --pathway GRN --ds-source hvg [--cores 13]
"""
import argparse, json, os, time
import numpy as np, pandas as pd, h5py
import anndata as ad, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import commot as ct

ROOT = "/Users/jiayifan/tansey_lab/alarmist"
p = argparse.ArgumentParser()
p.add_argument("--run-dir", default=f"{ROOT}/results/comparators/commot/GBM/cellchatdb2")
p.add_argument("--source-h5ad", default=f"{ROOT}/data/xenium_mm_final_cell_id.h5ad")
p.add_argument("--out-dir", default=f"{ROOT}/results/comparators/commot/GBM/impact")
p.add_argument("--cores", default=None, help="comma-separated; default = all cores with an adata")
p.add_argument("--pathway", default="GRN")
p.add_argument("--db-name", default="cellchat")
p.add_argument("--ds-source", default="hvg", choices=["hvg", "alarmist"])
p.add_argument("--alarmist-impact",
               default=f"{ROOT}/results/GBM/impact/motif_1_celltype_mGAM_de_results.csv")
p.add_argument("--n-ds-genes", type=int, default=30)
p.add_argument("--bg-genes", type=int, default=500)      # tutorial value
p.add_argument("--method", default="treebased_score")     # tutorial value
p.add_argument("--tree-ntrees", type=int, default=100)    # tutorial value
p.add_argument("--tree-repeat", type=int, default=100)    # tutorial value
p.add_argument("--tree-method", default="rf")             # tutorial value
p.add_argument("--tree-combined", action="store_true", default=True)  # tutorial value
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()
np.random.seed(a.seed)

OUT = os.path.join(a.out_dir, f"{a.pathway}_{a.ds_source}")
os.makedirs(OUT, exist_ok=True)
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
t0 = time.time()

cores = (a.cores.split(",") if a.cores else
         sorted([d for d in os.listdir(a.run_dir)
                 if os.path.exists(os.path.join(a.run_dir, d, "adata_commot.h5ad"))],
                key=lambda x: int(x)))
log(f"cores: {cores} | pathway={a.pathway} | ds_source={a.ds_source} | method={a.method}")

# ---- source raw counts, read once (dense 100197 x 5119 in this file) -------------------------
src = h5py.File(a.source_h5ad, "r")
_ix = src["obs"].attrs.get("_index", "_index")
_ix = _ix.decode() if isinstance(_ix, bytes) else str(_ix)
src_names = [x.decode() if isinstance(x, bytes) else str(x) for x in src["obs"][_ix][:]]
row_of = {c: i for i, c in enumerate(src_names)}
var_names = [x.decode() if isinstance(x, bytes) else str(x) for x in src["var"]["_index"][:]]
log(f"source h5ad: {len(src_names)} cells x {len(var_names)} genes")

# ---- the downstream gene set --------------------------------------------------------------
ds_from_alarmist = None
if a.ds_source == "alarmist":
    d = pd.read_csv(a.alarmist_impact)
    if "qval" in d.columns:
        d = d[d.qval < 0.05]
    d = d.reindex(d.logFC.abs().sort_values(ascending=False).index)
    ds_from_alarmist = [g for g in d.gene.tolist() if g in set(var_names)][:a.n_ds_genes]
    log(f"ALARMIST ds_genes ({len(ds_from_alarmist)}): {', '.join(ds_from_alarmist[:10])}...")

summary_rows, failures = [], []
for c in cores:
    ts = time.time()
    ap = os.path.join(a.run_dir, c, "adata_commot.h5ad")
    A = ad.read_h5ad(ap)
    # is this pathway present in this core's LR set at all?
    info = A.uns[f"commot-{a.db_name}-info"]["df_ligrec"]
    pw_col = info.columns[2]
    n_pairs = int((info[pw_col].astype(str) == a.pathway).sum())
    if n_pairs == 0:
        log(f"--- core {c}: pathway '{a.pathway}' absent from its 671 pairs, SKIPPED")
        summary_rows.append(dict(core=c, status="pathway_absent")); continue

    # rebuild .raw (raw counts) -- run_commot.py dropped it to save space
    idx = np.array([row_of[x] for x in A.obs_names])
    order = np.argsort(idx)
    counts = np.empty((len(idx), len(var_names)), dtype=np.float32)
    counts[order] = src["layers/counts"][np.sort(idx), :]
    raw = ad.AnnData(X=counts, obs=pd.DataFrame(index=A.obs_names),
                     var=pd.DataFrame(index=var_names))
    import scipy.sparse as sp
    raw.X = sp.csr_matrix(raw.X)          # .X.toarray() is called downstream; keep it sparse-typed
    A.raw = raw
    # communication_impact runs normalize_total+log1p on adata.raw.to_adata(). That is correct --
    # verified the rebuilt .raw holds genuine integer counts (max 128) -- but the stale
    # uns['log1p'] key left by the main run makes scanpy warn "already log-transformed". Same trap
    # as run_commot.py; popping it is behaviour-neutral and removes a misleading warning.
    A.uns.pop("log1p", None)
    log(f"--- core {c}: {A.shape[0]} cells | pathway '{a.pathway}' has {n_pairs} LR pairs | .raw rebuilt")

    if a.ds_source == "hvg":
        bg = raw.copy()
        sc.pp.normalize_total(bg, inplace=True); sc.pp.log1p(bg)
        sc.pp.highly_variable_genes(bg, n_top_genes=a.n_ds_genes)
        ds_genes = list(bg.var_names[bg.var.highly_variable])
    else:
        ds_genes = list(ds_from_alarmist)
    ds_genes = [g for g in ds_genes if g in set(A.raw.var_names)]
    log(f"    ds_genes = {len(ds_genes)}  |  bg_genes = {a.bg_genes}")

    try:
        df_impact = ct.tl.communication_impact(
            A, database_name=a.db_name, pathway_name=a.pathway,
            tree_combined=a.tree_combined, method=a.method,
            tree_ntrees=a.tree_ntrees, tree_repeat=a.tree_repeat, tree_method=a.tree_method,
            ds_genes=ds_genes, bg_genes=a.bg_genes, normalize=True)
    except Exception as e:
        log(f"    communication_impact FAILED: {type(e).__name__}: {e}")
        failures.append(dict(core=c, step="communication_impact", err=f"{type(e).__name__}: {e}"))
        summary_rows.append(dict(core=c, status="impact_failed")); del A; continue

    df_impact.to_csv(os.path.join(OUT, f"impact_core{int(c):02d}.csv"))
    # plot_communication_impact clusters the communication rows, so cluster_knn must be smaller
    # than the number of rows in a summary (= n_pairs + 1, the pair(s) plus the pathway aggregate).
    # The tutorial hits this too and passes cluster_knn=2 "here we only have two LR pairs in PSAP"
    # (3 rows). A single-pair pathway like GRN has 2 rows and needs cluster_knn=1. Fall back rather
    # than skip -- the figure is producible, it just needs the right knob.
    for smry in ("sender", "receiver"):
        try:
            done = False
            for knn in (5, 2, 1):          # 5 = function default, 2 = the tutorial's PSAP value
                try:
                    ct.pl.plot_communication_impact(
                        df_impact, summary=smry, top_ngene=30, top_ncomm=5, colormap="coolwarm",
                        font_scale=1.2, linewidth=0, cluster_knn=knn,
                        show_gene_names=True, show_comm_names=True,
                        filename=os.path.join(OUT, f"impact_core{int(c):02d}_{smry}.png"))
                    plt.close("all"); done = True
                    log(f"    plot_communication_impact({smry}): OK at cluster_knn={knn}")
                    break
                except IndexError:
                    plt.close("all"); continue
            if not done:
                raise RuntimeError("plot failed at cluster_knn 5, 2 and 1")
        except Exception as e:
            plt.close("all")
            log(f"    plot_communication_impact({smry}) FAILED: {type(e).__name__}: {e}")
            failures.append(dict(core=c, step=f"plot_{smry}", err=f"{type(e).__name__}: {e}"))

    rt = (time.time() - ts) / 60
    r_rows = [i for i in df_impact.index if i.startswith("r-")]
    s_rows = [i for i in df_impact.index if i.startswith("s-")]
    gcols = [x for x in df_impact.columns if x != "average"]
    summary_rows.append(dict(
        core=c, status="ok", n_cells=A.shape[0], n_pairs_in_pathway=n_pairs,
        n_ds_genes=len(ds_genes), impact_shape=f"{df_impact.shape[0]}x{df_impact.shape[1]}",
        max_receiver_impact=float(df_impact.loc[r_rows, gcols].values.max()),
        max_sender_impact=float(df_impact.loc[s_rows, gcols].values.max()),
        mean_receiver_impact=float(df_impact.loc[r_rows, gcols].values.mean()),
        top_receiver_gene=str(df_impact.loc[r_rows, gcols].max(axis=0).idxmax()),
        runtime_min=round(rt, 2)))
    log(f"    df_impact {df_impact.shape} | done in {rt:.1f} min")
    del A

src.close()
pd.DataFrame(summary_rows).to_csv(os.path.join(OUT, "impact_summary.csv"), index=False)
json.dump({"script": "run_commot_impact.py", "run_dir": a.run_dir, "pathway": a.pathway,
           "ds_source": a.ds_source, "n_ds_genes": a.n_ds_genes, "bg_genes": a.bg_genes,
           "method": a.method, "tree_ntrees": a.tree_ntrees, "tree_repeat": a.tree_repeat,
           "tree_method": a.tree_method, "tree_combined": a.tree_combined,
           "alarmist_impact": (a.alarmist_impact if a.ds_source == "alarmist" else None),
           "ds_genes_alarmist": ds_from_alarmist, "cores": cores, "seed": a.seed,
           "deg_detection_run": False,
           "deg_detection_blocked_by": "rpy2 + anndata2ri absent; tradeSeq absent from all four "
                                       "R installs (4.4.2, 3x 4.3.3); docstring pins tradeSeq "
                                       "1.0.1 / R 3.6.3",
           "failures": failures, "wall_min": round((time.time() - t0) / 60, 1)},
          open(os.path.join(OUT, "run_manifest.json"), "w"), indent=2)
log(f"ALL DONE in {(time.time()-t0)/60:.1f} min -> {OUT}")
print(pd.DataFrame(summary_rows).to_string(index=False))
