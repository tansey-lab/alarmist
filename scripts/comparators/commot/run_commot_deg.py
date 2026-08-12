#!/usr/bin/env python
"""COMMOT's full downstream chain on the GBM cores — the authors' workflow, verbatim.

    (1) ct.tl.communication_deg_detection      -> df_deg, df_yhat      [needs R: tradeSeq]
    (2) ct.tl.communication_deg_clustering     -> df_deg_clus, df_yhat_clus
    (3) ct.pl.plot_communication_dependent_genes(return_genes=True) -> top_de_genes
    (4) ct.tl.communication_impact(ds_genes=top_de_genes)           -> df_impact
    (5) ct.pl.plot_communication_impact
    (6) the tutorial's 3-panel example figure (received signal | a negative DE gene | a positive one)

Source of truth: /Users/jiayifan/tansey_lab/COMMOT/docs/notebooks/visium-mouse_brain.ipynb,
"Downstream analysis". Setup, traps and runtime expectations: SETUP_tradeSeq.md (read it first).

REQUIRES the `comp-commot-r` env (rpy2 + anndata2ri + R tradeSeq + clusterExperiment). It will not
run in `comp-commot`. ***This script has never been executed*** — it was written against the
installed source (`commot/tools/_downstream_analysis.py:32-183`), not against a successful run.
Treat the first invocation as a test; `--cores 2 --n-var-genes 500` is the cheap one.

Reads each core's persisted `adata_commot.h5ad`, so **no optimal-transport re-run** is needed;
the transport plans come off disk.

Usage: python run_commot_deg.py [--cores 2] --pathway GRN --n-var-genes 2000
"""
import argparse, json, os, pickle, sys, time
import numpy as np, pandas as pd, h5py
import scipy.sparse as sp
import anndata as ad, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Trap 1 of SETUP_tradeSeq.md, pre-empted ---------------------------------------------------
# commot calls `ro.numpy2ri.activate()` / `ro.pandas2ri.activate()` (_downstream_analysis.py:76-78).
# Since rpy2 3.5 those are submodules that must be imported before they exist as attributes of
# rpy2.robjects. commot relies on something else having imported them. Import them here so the
# attribute access resolves -- this touches OUR process, never the installed package.
try:
    import rpy2.robjects as _ro                      # noqa: F401
    from rpy2.robjects import numpy2ri as _n2r, pandas2ri as _p2r   # noqa: F401
except Exception as _e:                              # pragma: no cover
    sys.exit(f"ERROR: rpy2 not importable ({type(_e).__name__}: {_e}).\n"
             f"This script needs the comp-commot-r env — see SETUP_tradeSeq.md section 1.")
import commot as ct

ROOT = "/Users/jiayifan/tansey_lab/alarmist"
p = argparse.ArgumentParser()
p.add_argument("--run-dir", default=f"{ROOT}/results/comparators/commot/GBM/cellchatdb2")
p.add_argument("--source-h5ad", default=f"{ROOT}/data/xenium_mm_final_cell_id.h5ad")
p.add_argument("--out-dir", default=f"{ROOT}/results/comparators/commot/GBM/deg")
p.add_argument("--cores", default=None, help="comma-separated; default = every core, smallest first")
p.add_argument("--pathway", default="GRN")
p.add_argument("--db-name", default="cellchat")
p.add_argument("--summary", default="receiver", choices=["sender", "receiver"])  # tutorial value
p.add_argument("--n-var-genes", type=int, default=2000,
               help="Genes handed to tradeSeq. LEAVE THIS SET. With None, commot falls back to "
                    "highly_variable_genes(min_mean=0.0125, max_mean=3, min_disp=0.5) -- Seurat "
                    "defaults tuned for whole-transcriptome data, unpredictable on a 5,119-gene "
                    "targeted panel, and this number sets how long fitGAM runs.")
p.add_argument("--nknots", type=int, default=6)            # tutorial/function default
p.add_argument("--n-points", type=int, default=50)         # tutorial/function default
p.add_argument("--deg-pvalue-cutoff", type=float, default=0.05)
p.add_argument("--deg-clustering-res", type=float, default=0.4)   # tutorial value
p.add_argument("--top-ngene-per-cluster", type=int, default=5)    # tutorial value
p.add_argument("--bg-genes", type=int, default=500)        # tutorial value
p.add_argument("--tree-ntrees", type=int, default=100)     # tutorial value
p.add_argument("--tree-repeat", type=int, default=100)     # tutorial value
p.add_argument("--tree-method", default="rf")              # tutorial value
p.add_argument("--skip-deg", action="store_true",
               help="reuse deg_core<NN>.pkl if present -- step 1 is the expensive one")
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()
np.random.seed(a.seed)

OUT = os.path.join(a.out_dir, a.pathway); os.makedirs(OUT, exist_ok=True)
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
T0 = time.time()
SUM_ABRV = "s" if a.summary == "sender" else "r"

cores = (a.cores.split(",") if a.cores else
         [d for d in os.listdir(a.run_dir)
          if os.path.exists(os.path.join(a.run_dir, d, "adata_commot.h5ad"))])
# smallest first: fitGAM cost tracks cell count, so failures surface cheaply
sizes = {}
for c in cores:
    with h5py.File(os.path.join(a.run_dir, c, "adata_commot.h5ad"), "r") as f:
        sizes[c] = int(f["obsm/spatial"].shape[0])
cores = sorted(cores, key=lambda c: sizes[c])
log(f"cores (smallest first): {[(c, sizes[c]) for c in cores]}")
log(f"pathway={a.pathway} summary={a.summary} n_var_genes={a.n_var_genes} nknots={a.nknots}")

# ---- source raw counts, opened once ---------------------------------------------------------
src = h5py.File(a.source_h5ad, "r")
_ix = src["obs"].attrs.get("_index", "_index")
_ix = _ix.decode() if isinstance(_ix, bytes) else str(_ix)
src_names = [x.decode() if isinstance(x, bytes) else str(x) for x in src["obs"][_ix][:]]
row_of = {c: i for i, c in enumerate(src_names)}
var_names = [x.decode() if isinstance(x, bytes) else str(x) for x in src["var"]["_index"][:]]

summary_rows, failures = [], []
for c in cores:
    cid = f"{int(c):02d}"
    tc = time.time()
    A = ad.read_h5ad(os.path.join(a.run_dir, c, "adata_commot.h5ad"))

    # is the requested pathway present in this core's LR set?
    info = A.uns[f"commot-{a.db_name}-info"]["df_ligrec"]
    n_pairs = int((info[info.columns[2]].astype(str) == a.pathway).sum())
    if n_pairs == 0:
        log(f"--- core {c}: pathway '{a.pathway}' absent, SKIPPED")
        summary_rows.append(dict(core=c, status="pathway_absent")); del A; continue
    sum_key = f"commot-{a.db_name}-sum-{a.summary}"
    col = f"{SUM_ABRV}-{a.pathway}"
    if col not in A.obsm[sum_key].columns:
        log(f"--- core {c}: '{col}' missing from obsm['{sum_key}'], SKIPPED")
        summary_rows.append(dict(core=c, status="summary_column_missing")); del A; continue

    # Trap 2: run_commot.py cleared .layers and dropped .raw to halve the file. Both are needed --
    # deg_detection reads layers['counts'] (NOTE: the docstring says 'count', the code at :139 says
    # 'counts'; the code wins), communication_impact reads .raw. Rebuild from the source h5ad.
    idx = np.array([row_of[x] for x in A.obs_names])
    order = np.argsort(idx)
    counts = np.empty((len(idx), len(var_names)), dtype=np.float32)
    counts[order] = src["layers/counts"][np.sort(idx), :]
    counts = sp.csr_matrix(counts)
    A.layers["counts"] = counts
    A.raw = ad.AnnData(X=counts, obs=pd.DataFrame(index=A.obs_names),
                       var=pd.DataFrame(index=var_names))
    A.uns.pop("log1p", None)      # stale key -> spurious "already log-transformed" warning
    log(f"--- core {c}: {A.shape[0]} cells | '{a.pathway}' has {n_pairs} LR pair(s) | "
        f"layers['counts'] + .raw rebuilt")

    rec = dict(core=c, status="ok", n_cells=int(A.shape[0]), n_pairs_in_pathway=n_pairs)
    pkl = os.path.join(OUT, f"deg_core{cid}.pkl")

    # ---- (1) tradeSeq association test ------------------------------------------------------
    if a.skip_deg and os.path.exists(pkl):
        with open(pkl, "rb") as fh:
            d = pickle.load(fh)
        df_deg, df_yhat = d["df_deg"], d["df_yhat"]
        log(f"    (1) reused {os.path.basename(pkl)}: df_deg {df_deg.shape}, df_yhat {df_yhat.shape}")
    else:
        try:
            t = time.time()
            df_deg, df_yhat = ct.tl.communication_deg_detection(
                A, database_name=a.db_name, pathway_name=a.pathway, summary=a.summary,
                n_var_genes=a.n_var_genes, nknots=a.nknots, n_points=a.n_points,
                deg_pvalue_cutoff=a.deg_pvalue_cutoff)
            rec["deg_min"] = round((time.time() - t) / 60, 2)
            log(f"    (1) deg_detection OK in {rec['deg_min']:.1f} min | "
                f"df_deg {df_deg.shape}, df_yhat {df_yhat.shape}")
        except Exception as e:
            log(f"    (1) deg_detection FAILED: {type(e).__name__}: {e}")
            log("        -> see SETUP_tradeSeq.md section 3 to tell a version gap from a data problem")
            failures.append(dict(core=c, step="deg_detection", err=f"{type(e).__name__}: {e}"))
            summary_rows.append(dict(core=c, status="deg_failed")); del A; continue
        with open(pkl, "wb") as fh:                       # the tutorial's own checkpoint pattern
            pickle.dump({"df_deg": df_deg, "df_yhat": df_yhat}, fh,
                        protocol=pickle.HIGHEST_PROTOCOL)

    df_deg.to_csv(os.path.join(OUT, f"df_deg_core{cid}.csv"))
    df_yhat.to_csv(os.path.join(OUT, f"df_yhat_core{cid}.csv"))
    # COMMOT does not FDR-correct; add BH so nobody reads raw p-values as significance
    if "pvalue" in df_deg.columns:
        pv = df_deg["pvalue"].to_numpy()
        o = np.argsort(pv); n = len(pv)
        q = np.empty(n); q[o] = np.minimum.accumulate((pv[o] * n / (np.arange(n) + 1))[::-1])[::-1]
        # write the annotated copy to disk but keep `df_deg` itself untouched -- it is handed to
        # communication_deg_clustering next and should reach it exactly as commot produced it
        df_deg.assign(qvalue_BH=np.clip(q, 0, 1)).to_csv(
            os.path.join(OUT, f"df_deg_core{cid}.csv"))
        rec["n_deg_q05"] = int((np.clip(q, 0, 1) < 0.05).sum())
        rec["n_deg_p05"] = int((pv < 0.05).sum())
        log(f"    (1) {rec['n_deg_p05']} genes p<0.05 | {rec['n_deg_q05']} survive BH q<0.05")

    # ---- (2) cluster the fitted expression patterns -----------------------------------------
    try:
        df_deg_clus, df_yhat_clus = ct.tl.communication_deg_clustering(
            df_deg, df_yhat, deg_clustering_res=a.deg_clustering_res)
        df_deg_clus.to_csv(os.path.join(OUT, f"df_deg_clustered_core{cid}.csv"))
        df_yhat_clus.to_csv(os.path.join(OUT, f"df_yhat_clustered_core{cid}.csv"))
        log(f"    (2) deg_clustering OK: {df_deg_clus.shape}")
    except Exception as e:
        log(f"    (2) deg_clustering FAILED: {type(e).__name__}: {e}")
        failures.append(dict(core=c, step="deg_clustering", err=f"{type(e).__name__}: {e}"))
        summary_rows.append({**rec, "status": "clustering_failed"}); del A; continue

    # ---- (3) heatmap + the gene list the whole chain exists to produce -----------------------
    top_de_genes = []
    try:
        top_de_genes = ct.pl.plot_communication_dependent_genes(
            df_deg_clus, df_yhat_clus, top_ngene_per_cluster=a.top_ngene_per_cluster,
            font_scale=1.2, filename=os.path.join(OUT, f"heatmap_deg_core{cid}.png"),
            return_genes=True)
        plt.close("all")
        top_de_genes = [str(g) for g in list(top_de_genes)]
        with open(os.path.join(OUT, f"top_de_genes_core{cid}.txt"), "w") as fh:
            fh.write("\n".join(top_de_genes) + "\n")
        rec["n_top_de_genes"] = len(top_de_genes)
        log(f"    (3) top_de_genes ({len(top_de_genes)}): {', '.join(top_de_genes[:10])}...")
    except Exception as e:
        plt.close("all")
        log(f"    (3) plot_communication_dependent_genes FAILED: {type(e).__name__}: {e}")
        failures.append(dict(core=c, step="dependent_genes", err=f"{type(e).__name__}: {e}"))

    # ---- (4) impact of the signal on those genes --------------------------------------------
    if top_de_genes:
        try:
            t = time.time()
            df_impact = ct.tl.communication_impact(
                A, database_name=a.db_name, pathway_name=a.pathway, tree_combined=True,
                method="treebased_score", tree_ntrees=a.tree_ntrees, tree_repeat=a.tree_repeat,
                tree_method=a.tree_method, ds_genes=top_de_genes, bg_genes=a.bg_genes,
                normalize=True)
            df_impact.to_csv(os.path.join(OUT, f"impact_core{cid}.csv"))
            rec["impact_min"] = round((time.time() - t) / 60, 2)
            rec["impact_shape"] = f"{df_impact.shape[0]}x{df_impact.shape[1]}"
            g = [x for x in df_impact.columns if x != "average"]
            for pre, nm in (("r-", "receiver"), ("s-", "sender")):
                rows = [i for i in df_impact.index if i.startswith(pre)]
                if rows:
                    # NB: the null of this score is 0.5, not 0 (importance percentile vs bg genes)
                    rec[f"median_{nm}_impact"] = float(np.median(df_impact.loc[rows, g].values))
            log(f"    (4) impact {df_impact.shape} in {rec['impact_min']:.1f} min | "
                f"median receiver {rec.get('median_receiver_impact', float('nan')):.3f} vs "
                f"sender {rec.get('median_sender_impact', float('nan')):.3f}  (null = 0.5)")

            # ---- (5) impact heatmap; cluster_knn must be < rows per summary (= n_pairs + 1).
            # The tutorial passes cluster_knn=2 for PSAP (2 pairs -> 3 rows); a single-pair
            # pathway like GRN has 2 rows and needs 1. Fall back rather than skip.
            for smry in ("sender", "receiver"):
                done = False
                for knn in (5, 2, 1):
                    try:
                        ct.pl.plot_communication_impact(
                            df_impact, summary=smry, top_ngene=30, top_ncomm=5,
                            colormap="coolwarm", font_scale=1.2, linewidth=0, cluster_knn=knn,
                            show_gene_names=True, show_comm_names=True,
                            filename=os.path.join(OUT, f"impact_core{cid}_{smry}.png"))
                        plt.close("all"); done = True
                        log(f"    (5) plot_communication_impact({smry}) OK at cluster_knn={knn}")
                        break
                    except IndexError:
                        plt.close("all"); continue
                if not done:
                    failures.append(dict(core=c, step=f"plot_impact_{smry}",
                                         err="IndexError at cluster_knn 5, 2 and 1"))
        except Exception as e:
            log(f"    (4) communication_impact FAILED: {type(e).__name__}: {e}")
            failures.append(dict(core=c, step="impact", err=f"{type(e).__name__}: {e}"))

    # ---- (6) the tutorial's 3-panel example figure -------------------------------------------
    # Tutorial hand-picks one negative (Ctxn1) and one positive (Gpr37) DE gene. Choose them
    # programmatically instead: the largest fall / rise of the fitted curve along the signal axis.
    try:
        # Use df_yhat (step 1's yhatScaled: genes x n_points, all numeric) rather than the
        # clustered frame, which may carry a numeric cluster column that would corrupt iloc[:, -1].
        num = df_yhat.select_dtypes(include=[np.number])
        trend = (num.iloc[:, -1] - num.iloc[:, 0]).dropna()
        keep = set(A.var_names) & (set(top_de_genes) if top_de_genes else set(A.var_names))
        cand = [g for g in trend.index if str(g) in keep]
        if not cand:
            raise ValueError("no fitted gene is present in adata.var_names")
        trend = trend.loc[cand]
        gneg, gpos = str(trend.idxmin()), str(trend.idxmax())
        X = A.obsm["spatial"]
        sig = A.obsm[sum_key][col].values
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        for k, (vals, title) in enumerate((
                (sig, f"Amount of {a.summary} signal ({a.pathway})"),
                (np.asarray(A[:, gneg].X.todense()).ravel() if sp.issparse(A.X)
                 else np.asarray(A[:, gneg].X).ravel(), f"Negative DE gene ({gneg})"),
                (np.asarray(A[:, gpos].X.todense()).ravel() if sp.issparse(A.X)
                 else np.asarray(A[:, gpos].X).ravel(), f"Positive DE gene ({gpos})"))):
            o = np.argsort(vals)
            ax[k].scatter(X[o, 0], X[o, 1], c=np.asarray(vals)[o], cmap="coolwarm", s=4)
            ax[k].set_title(title, fontsize=10); ax[k].set_aspect("equal"); ax[k].axis("off")
        fig.savefig(os.path.join(OUT, f"examples_core{cid}.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)
        rec["example_neg_gene"], rec["example_pos_gene"] = gneg, gpos
        log(f"    (6) examples: negative {gneg}, positive {gpos}")
    except Exception as e:
        plt.close("all")
        log(f"    (6) example figure FAILED: {type(e).__name__}: {e}")
        failures.append(dict(core=c, step="examples", err=f"{type(e).__name__}: {e}"))

    rec["runtime_min"] = round((time.time() - tc) / 60, 2)
    summary_rows.append(rec)
    log(f"    core {c} done in {rec['runtime_min']:.1f} min")
    del A
    pd.DataFrame(summary_rows).to_csv(os.path.join(OUT, "deg_summary.csv"), index=False)

src.close()
pd.DataFrame(summary_rows).to_csv(os.path.join(OUT, "deg_summary.csv"), index=False)

versions = {"python": sys.version.split()[0], "commot": "0.0.3", "numpy": np.__version__,
            "pandas": pd.__version__, "scanpy": sc.__version__, "anndata": ad.__version__}
try:
    import rpy2, anndata2ri, rpy2.robjects as ro
    versions.update({"rpy2": rpy2.__version__, "anndata2ri": anndata2ri.__version__,
                     "R": str(ro.r('R.version.string')[0]),
                     "tradeSeq": str(ro.r('as.character(packageVersion("tradeSeq"))')[0]),
                     "clusterExperiment": str(ro.r('as.character(packageVersion("clusterExperiment"))')[0])})
except Exception as e:
    versions["r_stack"] = f"unavailable: {type(e).__name__}: {e}"

json.dump({"script": "run_commot_deg.py", "run_dir": a.run_dir, "pathway": a.pathway,
           "summary": a.summary, "n_var_genes": a.n_var_genes, "nknots": a.nknots,
           "n_points": a.n_points, "deg_pvalue_cutoff": a.deg_pvalue_cutoff,
           "deg_clustering_res": a.deg_clustering_res,
           "top_ngene_per_cluster": a.top_ngene_per_cluster, "bg_genes": a.bg_genes,
           "tree_ntrees": a.tree_ntrees, "tree_repeat": a.tree_repeat,
           "tree_method": a.tree_method, "tree_combined": True, "seed": a.seed,
           "cores": cores, "versions": versions, "failures": failures,
           "note_impact_null": "communication_impact scores are importance percentiles against "
                               "~bg_genes background genes; the null is 0.5, not 0",
           "wall_min": round((time.time() - T0) / 60, 1)},
          open(os.path.join(OUT, "run_manifest.json"), "w"), indent=2)
log(f"ALL DONE in {(time.time()-T0)/60:.1f} min -> {OUT}")
print(pd.DataFrame(summary_rows).to_string(index=False))
