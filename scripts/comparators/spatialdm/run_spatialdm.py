#!/usr/bin/env python
"""SpatialDM on a spatial h5ad, following the authors' tutorials.

Tutorials: /Users/jiayifan/tansey_lab/SpatialDM/tutorial/{melanoma,differential_test_intestine}.ipynb
Call contract: scripts/comparators/spatialdm/NOTES.md   Deviations: DEVIATIONS.md

Usage: python run_spatialdm.py --h5ad PATH --out-dir DIR [options]
"""
import argparse, json, os, sys, time
# bokeh's chord export shells out to geckodriver/firefox, which live in the env's bin/ -- and that
# bin/ is NOT on PATH when the interpreter is invoked by absolute path, which is how this repo has
# to run (conda activate is broken here; see CLAUDE.md). Without this, every chord plot fails.
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", "")
from itertools import zip_longest
import numpy as np, pandas as pd
import anndata as ad, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import spatialdm as sdm
import spatialdm.plottings as pl
# bokeh/holoviews (used by chord_celltype / chord_LR) writes a default output_file named
# after the CALLING SCRIPT into the script's own directory -- i.e. a 600 KB .html dropped into
# scripts/, violating the code/output separation rule. Pin it to a temp path.
try:
    import tempfile as _tf
    from bokeh.io import output_file as _bokeh_output_file
    _bokeh_output_file(os.path.join(_tf.gettempdir(), "bokeh_scratch.html"))
except Exception:
    pass

p = argparse.ArgumentParser()
p.add_argument("--h5ad", required=True)
p.add_argument("--out-dir", required=True)
p.add_argument("--cell-type-col", default="cell_type")
p.add_argument("--count-layer", default="counts", help="layer with RAW counts -> goes to .raw")
p.add_argument("--l", type=float, default=75.0,
               help="RBF length scale in COORDINATE UNITS. 75 reproduces the intestine "
                    "tutorial's absolute 135um effective radius (see NOTES.md)")
p.add_argument("--cutoff", type=float, default=0.2)
p.add_argument("--single-cell", action="store_true", default=True)
p.add_argument("--min-cell", type=int, default=3)
p.add_argument("--n-perm", type=int, default=1000)
p.add_argument("--nproc", type=int, default=1)
p.add_argument("--db", default="builtin",
               help="'builtin' = CellChatDB v1 via extract_lr (see --datahost); "
                    "else path to CellChatDBv2.0.human.csv")
p.add_argument("--datahost", default="package", choices=["package", "builtin"],
               help="extract_lr's datahost, and its name is INVERTED from what it does: "
                    "'package' reads the files shipped inside the installed wheel "
                    "(spatialdm/datasets/LR_data/), while 'builtin' -- the function's own "
                    "DEFAULT, and what the tutorial gets -- DOWNLOADS from figshare on every "
                    "call. We default to 'package': same CellChatDB v1, but offline, pinned to "
                    "the installed version, and not re-fetched once per split.")
p.add_argument("--requested-lrs", default="GRN_SORT1,ANXA1_FPR1")
p.add_argument("--n-top-plots", type=int, default=3)
p.add_argument("--seed", type=int, default=0)
p.add_argument("--split-col", default=None,
               help="obs column to run SEPARATELY per level (e.g. tma_id). NOT primarily about "
                    "cross-piece weights -- on the GBM TMA the cores are >=222.9 um apart and the "
                    "cutoff zeroes every cross-core weight anyway. The real reasons: "
                    "spatialdm_global returns ONE Moran's R per pair per object (pooling gives a "
                    "tissue-wide number, not one per piece); diff_utils.concat_obj needs a LIST of "
                    "separately fitted objects; and extract_lr's min_cell filter should adapt per "
                    "piece. See DEVIATIONS.md.")
p.add_argument("--n-neighbors", default="auto",
               help="size of the large RBF graph. 'auto' = (max cells within the effective "
                    "radius in THAT split) + 1, so `cutoff` truncates and the kNN cap never "
                    "binds -- the regime both tutorials operate in. INDEPENDENT of "
                    "--n-nearest-neighbors (rbfweight only derives it when n_neighbors is None).")
p.add_argument("--n-nearest-neighbors", type=int, default=6,
               help="separate SMALL graph for adjacent signalling; NOT coupled to --n-neighbors")
a = p.parse_args()
np.random.seed(a.seed)

OUT = a.out_dir; os.makedirs(OUT, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
R = a.l * np.sqrt(-2 * np.log(a.cutoff))   # effective radius; see NOTES.md -- never use eff_dist

# ------------------------------------------------------------------ load once
log("reading", a.h5ad)
adata = sc.read_h5ad(a.h5ad)
log(f"loaded {adata.shape[0]} x {adata.shape[1]}")

# SpatialDM INPUT CONTRACT: log-transformed .X AND raw counts in .raw. Both are used --
# extract_lr/global Moran read .X, but LOCAL spot selection reads .raw
# (utils.spot_selection_matrix: raw_norm = adata.raw.to_adata(), max-normalise, sc.pp.scale).
# Leaving .raw as the log matrix does not error; it silently changes the local calls.
if a.count_layer != "X":
    if a.count_layer not in adata.layers:
        raise SystemExit(f"ERROR: layer '{a.count_layer}' absent; have {list(adata.layers)}")
    adata.raw = ad.AnnData(X=adata.layers[a.count_layer].copy(),
                           obs=adata.obs.copy(), var=adata.var.copy())
    log(f".raw <- layers['{a.count_layer}'] (raw counts; used by local spot selection)")
adata.layers.clear()
log(f".X max = {adata.X.max():.2f} (expect log-scale; tutorial melanoma is 9.98)")

if a.cell_type_col not in adata.obs:
    raise SystemExit(f"ERROR: '{a.cell_type_col}' not in obs")
adata.obs["cell_type"] = adata.obs[a.cell_type_col].astype(str).astype("category")
adata.obsm["spatial"] = np.asarray(adata.obsm["spatial"], dtype=float)

splits = ([(str(v), np.asarray(adata.obs[a.split_col] == v)) for v in
           sorted(adata.obs[a.split_col].unique())] if a.split_col
          else [("all", np.ones(adata.shape[0], bool))])
log(f"splits: {len(splits)} ({a.split_col or 'none'}) | l={a.l} cutoff={a.cutoff} "
    f"-> effective radius {R:.1f} coord units")

# ------------------------------------------------------- LR db loader (v2 tier)
def attach_cellchat_v2(sub, csv):
    """Replicate extract_lr() with CellChatDB v2.

    extract_lr has NO custom-DB hook (`datahost` only selects builtin vs figshare, both v1), so
    we build uns['ligand'/'receptor'/'geneInter'/'num_pairs'] ourselves with identical logic:
    subunits split and filtered to present genes; a pair kept only if BOTH sides retain >=1 gene
    AND >= min_cell cells express each side's mean; NaN-padded frames indexed by interaction.
    LOSSLESS -- SpatialDM handles multi-subunit complexes natively, so no v2 row is dropped for
    being a complex (contrast stLearn, which must drop 57.5%).
    """
    db = pd.read_csv(csv)
    db["interaction_name"] = db["ligand"].astype(str) + "_" + db["receptor"].astype(str)
    db = db.drop_duplicates("interaction_name").set_index("interaction_name")
    genes = set(sub.var_names.astype(str))
    sp = lambda s: [g for g in str(s).split("_") if g]
    lig = [np.array([g for g in sp(v) if g in genes], dtype=object) for v in db["ligand"]]
    rec = [np.array([g for g in sp(v) if g in genes], dtype=object) for v in db["receptor"]]
    sub.uns["mean"] = "algebra"
    keep = []
    for i in range(len(lig)):
        if len(lig[i]) and len(rec[i]):
            mL = np.asarray(sub[:, list(lig[i])].X.mean(axis=1)).ravel()
            mR = np.asarray(sub[:, list(rec[i])].X.mean(axis=1)).ravel()
            keep.append(bool((mL > 0).sum() >= a.min_cell and (mR > 0).sum() >= a.min_cell))
        else:
            keep.append(False)
    keep = np.asarray(keep)
    cols = [c for c in ("pathway", "signaling_type", "version") if c in db.columns]
    gi = db.loc[db.index[keep], cols].rename(
        columns={"signaling_type": "annotation", "pathway": "pathway_name"})
    if "annotation" not in gi: gi["annotation"] = "Secreted Signaling"
    # CellChatDB v2 introduced a FOURTH signaling_type, 'Non-protein Signaling' (994/3233 rows),
    # that SpatialDM's v1-era code does not know about. globle_st_compute builds `st` as
    #   hstack(repeat(nm0, #{ECM-Receptor, Cell-Cell Contact}), repeat(nm, #{Secreted Signaling}))
    # so those rows get no variance term and `st` comes out short -> IndexError against the
    # length-N idx_use mask. It is also internally inconsistent beforehand, since
    # n_short_lri = (annotation != 'Secreted Signaling') counts them as SHORT-range while st
    # allocates them nothing. We remap them to 'Secreted Signaling' (long-range RBF weight),
    # matching this benchmark's CytoSignal mapping (Secreted/ECM/Non-protein -> diffusion) and
    # the fact that non-protein ligands are diffusible neurotransmitters/metabolites.
    n_np = int((gi["annotation"] == "Non-protein Signaling").sum())
    if n_np:
        gi["annotation"] = gi["annotation"].replace("Non-protein Signaling", "Secreted Signaling")
        print(f"  remapped {n_np} 'Non-protein Signaling' pairs -> 'Secreted Signaling' "
              f"(SpatialDB v1 vocabulary; see DEVIATIONS.md)", flush=True)
    # sort_values puts Cell-Cell Contact, ECM-Receptor first and Secreted last -- which is exactly
    # the block order globle_st_compute assumes (nm0 block then nm block) and that n_short_lri indexes.
    gi = gi.sort_values("annotation")
    assert set(gi["annotation"]) <= {"Cell-Cell Contact", "ECM-Receptor", "Secreted Signaling"}, \
        f"unhandled annotation(s) for SpatialDM: {set(gi['annotation'])}"
    pos = {k: i for i, k in enumerate(db.index)}
    sel = [pos[k] for k in gi.index]
    L = pd.DataFrame.from_records(zip_longest(*pd.Series([lig[i] for i in sel]).values)).transpose()
    L.columns = [f"Ligand{i}" for i in range(L.shape[1])]; L.index = gi.index
    Rf = pd.DataFrame.from_records(zip_longest(*pd.Series([rec[i] for i in sel]).values)).transpose()
    Rf.columns = [f"Receptor{i}" for i in range(Rf.shape[1])]; Rf.index = gi.index
    sub.uns["ligand"], sub.uns["receptor"] = L, Rf
    # SpatialDM's own plots expect CellChatDB v1 column names: compute_pathway reads
    # `.interaction_name` and chord_LR reads `.interaction_name_2` ("LIG-REC"). Our v2 frame
    # keeps the name in the index, so add both columns here.
    gi["interaction_name"] = gi.index.astype(str)
    gi["interaction_name_2"] = (L.apply(lambda r: "_".join(r.dropna().astype(str)), axis=1)
                                .reindex(gi.index).astype(str) + "-" +
                                Rf.apply(lambda r: "_".join(r.dropna().astype(str)), axis=1)
                                .reindex(gi.index).astype(str))
    sub.uns["geneInter"], sub.uns["num_pairs"] = gi, len(gi.index)

# total DB rows before any per-split min_cell filtering -- recorded in the manifest so the two
# tiers can be compared on equal footing.
if a.db == "builtin":
    if a.datahost == "package":
        # same file extract_lr's datahost='package' branch reads, resolved by module path rather
        # than pkg_resources (setuptools >= 81 no longer ships it).
        import spatialdm.datasets as _ds
        db_rows_total = int(len(pd.read_csv(os.path.join(
            os.path.dirname(_ds.__file__), "LR_data",
            "human-interaction_input_CellChatDB.csv.gz"), index_col=0, compression="gzip")))
    else:
        db_rows_total = None          # figshare copy, counted at extract_lr time only
else:
    db_rows_total = int(len(pd.read_csv(a.db)))
log(f"LR database: {a.db}"
    + (f" (datahost={a.datahost})" if a.db == "builtin" else "")
    + (f", {db_rows_total} rows before filtering" if db_rows_total else ""))

# --------------------------------------------------------------- per-split run
requested = [s.strip() for s in a.requested_lrs.split(",") if s.strip()]
summary = []
for name, mask in splits:
    ts = time.time()
    sub = adata[mask].copy()
    sdir = os.path.join(OUT, name if a.split_col else "")
    PLOTS, REQ, DATA = (os.path.join(sdir, "plots"), os.path.join(sdir, "plots", "requested"),
                        os.path.join(sdir, "data"))
    for d in (PLOTS, REQ, DATA): os.makedirs(d, exist_ok=True)
    xy = sub.obsm["spatial"]
    log(f"--- split {name}: {sub.shape[0]} cells ---")

    # n_neighbors so that `cutoff` truncates, never the kNN cap. Both tutorials operate this way
    # (melanoma ~14.6 units inside its radius, intestine ~5.7, vs a 186 default cap), so letting
    # the cap bind would be a structural departure, not a preservation of defaults.
    within = np.array([len(x) for x in cKDTree(xy).query_ball_point(xy, R)]) - 1
    if a.n_neighbors == "auto":
        n_nb = int(min(int(within.max()) + 1, sub.shape[0] - 1))
    else:
        n_nb = int(min(int(a.n_neighbors), sub.shape[0] - 1))
    bind = float((within > n_nb).mean())
    log(f"  n_neighbors={n_nb} (max within {R:.0f} units = {int(within.max())}, "
        f"median {int(np.median(within))}) | cap binds for {100*bind:.3f}% of cells")
    if bind > 0.01:
        raise SystemExit(f"STOP: kNN cap binds for {100*bind:.2f}% of cells in split {name} "
                         f"(>1%). Raise --n-neighbors.")

    sdm.weight_matrix(sub, l=a.l, cutoff=a.cutoff, n_neighbors=n_nb,
                      n_nearest_neighbors=a.n_nearest_neighbors, single_cell=a.single_cell)
    # rbfweight zeroes sub-cutoff entries with rbf_d[rows, cols] = 0 and never calls
    # eliminate_zeros(), so stored nnz stays n_cells * n_neighbors regardless of cutoff. Nothing
    # downstream uses .nonzero()/.nnz on W (all uses are .multiply/.sum/@), so compacting is safe.
    stored_before = sub.obsp["weight"].nnz
    sub.obsp["weight"].eliminate_zeros()
    sub.obsp["nearest_neighbors"].eliminate_zeros()
    stored = sub.obsp["weight"].nnz
    log(f"  W stored nnz {stored_before:,} -> {stored:,} after eliminate_zeros")

    if a.db == "builtin":
        sdm.extract_lr(sub, "human", min_cell=a.min_cell, datahost=a.datahost)
        db_name = f"CellChatDB v1 (extract_lr datahost={a.datahost})"
    else:
        attach_cellchat_v2(sub, a.db); db_name = os.path.basename(a.db)
    log(f"  LR db '{db_name}': {sub.uns['num_pairs']} valid pairs")

    sdm.spatialdm_global(sub, n_perm=a.n_perm, method="z-score", nproc=a.nproc)
    sdm.sig_pairs(sub, method="z-score", fdr=True, threshold=0.1)
    n_sel = int(sub.uns["global_res"].selected.sum())
    sdm.spatialdm_local(sub, n_perm=a.n_perm, method="z-score", nproc=a.nproc)
    sdm.sig_spots(sub, method="z-score", fdr=False, threshold=0.1)   # tutorial passes fdr=False
    log(f"  globally significant (FDR<0.1): {n_sel} / {sub.uns['num_pairs']}")

    gr = sub.uns["global_res"].copy(); gr.to_csv(os.path.join(DATA, "global_res.csv"))
    sub.uns["selected_spots"].to_csv(os.path.join(DATA, "selected_spots.csv.gz"), compression="gzip")
    pd.Series(sub.uns["local_stat"]["n_spots"], name="n_spots").to_csv(
        os.path.join(DATA, "local_n_spots.csv"))
    # local_z / local_z_p / local_perm_p live at uns TOP LEVEL (main.py:232-233, utils.py:246,266);
    # uns['local_stat'] holds only local_I, local_I_R, n_spots and the permutation arrays.
    # These are (pairs x cells) floats -- hundreds of MB per core as text -- so store float32 npz.
    for k in ("local_z_p", "local_perm_p", "local_z"):
        v = sub.uns.get(k, sub.uns["local_stat"].get(k))
        if v is None: continue
        arr = v.values if hasattr(v, "values") else np.asarray(v)
        np.savez_compressed(os.path.join(DATA, f"{k}.npz"), values=arr.astype(np.float32),
                            pairs=np.asarray(getattr(v, "index", np.arange(arr.shape[0]))).astype(str),
                            cells=np.asarray(getattr(v, "columns", np.arange(arr.shape[1]))).astype(str))
        log(f"  saved {k} {arr.shape} (float32 npz)")
    pd.DataFrame({"cell": sub.obs_names, "x": xy[:, 0], "y": xy[:, 1],
                  "cell_type": sub.obs["cell_type"].astype(str)}).to_csv(
        os.path.join(DATA, "cell_meta.csv"), index=False)
    sub.uns["ligand"].to_csv(os.path.join(DATA, "lr_ligand_subunits.csv"))
    sub.uns["receptor"].to_csv(os.path.join(DATA, "lr_receptor_subunits.csv"))

    ranked = gr.sort_values("fdr").index.tolist()
    top = [x for x in ranked if gr.loc[x, "selected"]][:a.n_top_plots]

    def savefig(nm, d): plt.gcf().savefig(os.path.join(d, nm + ".png"), dpi=200,
                                          bbox_inches="tight"); plt.close("all")
    def try_plot(fn, nm, d=None):
        d = d or PLOTS
        try: fn(); savefig(nm, d)
        except Exception as e: plt.close("all"); log(f"  plot FAILED {nm}: {type(e).__name__}: {e}")
    if top:
        try_plot(lambda: pl.global_plot(sub, pairs=top[:4], figsize=(6, 5), cmap="RdGy_r",
                                        vmin=-1.5, vmax=2), "global_plot")
        try_plot(lambda: pl.plot_pairs(sub, top, marker="s"), "plot_pairs_top")
    req_rank = {}
    for lr in requested:
        if lr in gr.index:
            rk = ranked.index(lr) + 1
            req_rank[lr] = (rk, bool(gr.loc[lr, "selected"]), float(gr.loc[lr, "fdr"]))
            try_plot(lambda lr=lr: pl.plot_pairs(sub, [lr], marker="s"), f"rank{rk}_{lr}", REQ)
        else:
            req_rank[lr] = (None, False, None)

    # ---- completeness pass: the plots that attribute a pair to CELL TYPES ----
    focus = top[:2] + [x for x in requested if x in gr.index]
    try:
        sub.obsm["celltypes"] = pd.get_dummies(sub.obs["cell_type"]).astype(float)
    except Exception as e:
        log(f"  celltypes obsm failed: {e}")
    for lr in focus:
        try_plot(lambda lr=lr: pl.ligand_ct(sub, lr), f"ligand_ct_{lr}")
        try_plot(lambda lr=lr: pl.receptor_ct(sub, lr), f"receptor_ct_{lr}")
        try_plot(lambda lr=lr: pl.chord_celltype(sub, pairs=[lr], ncol=1, min_quantile=0.01),
                 f"chord_celltype_{lr}")
    cts = list(sub.obs["cell_type"].cat.categories)
    try_plot(lambda: pl.chord_LR(sub, senders=cts, receivers=cts, min_quantile=0.5), "chord_LR_all")
    # pathway enrichment of the selected pairs (tutorial: compute_pathway -> dot_path)
    def _pathway():
        sel = [x for x in ranked if gr.loc[x, "selected"]]
        res = sdm.compute_pathway(sub, dic={f"core{name}": pd.Index(sel)})
        res.to_csv(os.path.join(DATA, "pathway_enrichment.csv"))
        pl.dot_path(res, figsize=(6, 10))
    try_plot(_pathway, "dot_path")

    # ---- PERSIST THE OBJECT via SpatialDM's own writer, so replots never need a re-run ----
    try:
        out_ad = sub.copy()
        ls = out_ad.uns.get("local_stat", {})
        for k in ("local_permI", "local_permI_R"):
            if k in ls: ls[k] = "None"
        sdm.write_spatialdm_h5ad(out_ad, filename=os.path.join(DATA, "spatialdm.h5ad"))
        log("  saved spatialdm.h5ad (re-readable via sdm.read_spatialdm_h5ad)")
    except Exception as e:
        log(f"  WARNING: write_spatialdm_h5ad failed: {type(e).__name__}: {e}")

    rt = (time.time() - ts) / 60
    row = dict(split=name, n_cells=int(sub.shape[0]), n_neighbors=n_nb,
               max_within_radius=int(within.max()), median_within_radius=int(np.median(within)),
               cap_bind_frac=bind, stored_nnz_before=int(stored_before), stored_nnz=int(stored),
               n_pairs_valid=int(sub.uns["num_pairs"]), n_significant=n_sel, runtime_min=round(rt, 2))
    for lr, (rk, sel_, fdr) in req_rank.items():
        row[f"{lr}_rank"] = rk; row[f"{lr}_selected"] = sel_; row[f"{lr}_fdr"] = fdr
    summary.append(row)
    log(f"  split {name} done in {rt:.1f} min")

pd.DataFrame(summary).to_csv(os.path.join(OUT, "per_split_summary.csv"), index=False)
try:
    import subprocess, resource
    sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                  stderr=subprocess.DEVNULL).decode().strip()
    rss = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2)
except Exception:
    sha, rss = None, None
json.dump({"script": "run_spatialdm.py", "h5ad": a.h5ad, "split_col": a.split_col,
           "n_splits": len(splits), "l": a.l, "cutoff": a.cutoff,
           "effective_radius": round(float(R), 1), "n_neighbors": a.n_neighbors,
           "n_nearest_neighbors": a.n_nearest_neighbors, "single_cell": a.single_cell,
           "min_cell": a.min_cell, "n_perm": a.n_perm,
           "db": a.db, "datahost": a.datahost if a.db == "builtin" else None,
           "db_rows_total": db_rows_total, "cell_type_column": a.cell_type_col,
           "count_layer": a.count_layer, "n_cells_total": int(adata.shape[0]),
           "n_cells_per_split": [int(r["n_cells"]) for r in summary],
           "n_genes": int(adata.shape[1]),
           "spatialdm": sdm.__version__, "seed": a.seed, "git_sha": sha, "peak_rss_gb": rss,
           "wall_min": round((time.time() - t0) / 60, 1)},
          open(os.path.join(OUT, "run_manifest.json"), "w"), indent=2)
log(f"ALL DONE in {(time.time()-t0)/60:.1f} min -> {OUT}")
print(pd.DataFrame(summary).to_string(index=False))
