#!/usr/bin/env python
"""COMMOT on a spatial h5ad, following the authors' Visium tutorial.

Tutorials: /Users/jiayifan/tansey_lab/COMMOT/docs/notebooks/{Basic_usage,visium-mouse_brain}.ipynb
Call contract: scripts/comparators/commot/NOTES.md
Deviations:    scripts/comparators/commot/DEVIATIONS.md

NOTE ON SOURCE PROVENANCE: the local clone /Users/jiayifan/tansey_lab/COMMOT is AHEAD of the
installed release (commot 0.0.3 from PyPI) -- its `_optimal_transport/_cot.py` carries two
post-0.0.3 upstream commits (np.Inf fix; sparse cost-matrix support) that the INSTALLED package
does not have. `tools/_spatial_communication.py` is identical between the two. Every call
signature below was verified against the INSTALLED package.

Usage: python run_commot.py --h5ad PATH --out-dir DIR [options]
"""
import argparse, gc, json, os, subprocess, sys, time
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import commot as ct

# --- RECORDED VERSION GAP: ct.pl.plot_cluster_communication_dotplot is UNRUNNABLE here --------
# commot 0.0.3's dotplot was written against matplotlib <3.9 and seaborn <0.13; installed here
# are matplotlib 3.10.9 and seaborn 0.13.2. Two independent breakages, verified by running it:
#   1. `plotting/_plotting.py:788` iterates `g.legend.legendHandles` -- an alias matplotlib
#      REMOVED in 3.9 (now `legend_handles`).  -> AttributeError
#   2. past that, the handles seaborn 0.13 returns are `Line2D`, not the `PathCollection` the
#      code assumes, so `set_edgecolor` does not exist. -> AttributeError
# Restoring (1) is a one-line alias, but (2) needs the function's internals rewritten, which
# would stop being "use the method as its authors wrote it". Recorded, NOT patched. The
# information the dotplot would show is fully persisted as cluster_comm_*.csv / cluster_pval_*.csv,
# and `ours_dotplot` below draws an equivalent from those -- labelled `ours_` so it is never
# mistaken for the method's own figure.

p = argparse.ArgumentParser()
p.add_argument("--h5ad", required=True)
p.add_argument("--out-dir", required=True)
p.add_argument("--cell-type-col", default="cell_type")
p.add_argument("--count-layer", default="counts")
p.add_argument("--db", default="builtin",
               help="'builtin' = COMMOT's bundled CellChatDB (v1); else path to CellChatDBv2.0.human.csv")
p.add_argument("--signaling-types", default="Secreted Signaling,Non-protein Signaling",
               help="CellChatDB signaling_type values to keep. COMMOT applies ONE dis_thr to every "
                    "pair, so only DIFFUSIBLE categories belong here (the tutorial uses Secreted "
                    "only; v2's Non-protein is also diffusible and did not exist in v1).")
p.add_argument("--dis-thr", type=float, default=365.0,
               help="distance threshold in COORDINATE UNITS. 365 reproduces the tutorial's "
                    "physical constraint: its dis_thr=500 on Visium coords = 365um (see NOTES.md)")
p.add_argument("--filter-criteria", default="min_cell", choices=["min_cell", "min_cell_pct"],
               help="COMMOT offers both; the tutorial calls min_cell_pct=0.05, which is 5%% of "
                    "Visium SPOTS (10-30 cells each). At single-cell resolution that keeps only "
                    "0.9-1.8%% of pairs vs the tutorial's 20.9%%, so we use the function's other "
                    "default, min_cell=100 -- an absolute floor that is stable across cores.")
p.add_argument("--filter-scope", default="per_split", choices=["per_split", "global"],
               help="WHERE the expression filter is evaluated. 'per_split' (the 2026-08-01 run) "
                    "applies min_cell inside each split, so the pair set differs per split and "
                    "sparse splits can end up with ZERO pairs. 'global' evaluates it once on the "
                    "whole object and hands the SAME pair set to every split -- pair sets become "
                    "comparable across splits and no split is dropped. The tutorial is a single "
                    "section, so it does not distinguish the two.")
p.add_argument("--min-cell", type=int, default=100)
p.add_argument("--min-cell-pct", type=float, default=0.05)
p.add_argument("--n-permutations", type=int, default=100)
p.add_argument("--split-col", default=None)
p.add_argument("--max-cells-per-split", type=int, default=0,
               help="0 = no limit. COMMOT builds a DENSE NxN distance matrix; guard for OOM.")
p.add_argument("--requested-lrs", default="GRN_SORT1,ANXA1_FPR1")
p.add_argument("--n-top-lr", type=int, default=3,
               help="how many of the method's own top LR PAIRS to draw (pathway aggregates are "
                    "ranked and drawn separately -- see --n-top-pathways)")
p.add_argument("--n-top-pathways", type=int, default=3,
               help="how many top PATHWAYS (by total received signal) get cluster_communication "
                    "permutation p-values and the native cluster-level plots")
p.add_argument("--save-adata", action="store_true",
               help="persist the per-split AnnData INCLUDING obsp (the cell x cell transport "
                    "plans). Without this, every COMMOT downstream tool needs a full OT re-run.")
p.add_argument("--no-native-plots", action="store_true",
               help="skip the ct.pl.* figures (vector fields, cluster network, dotplots)")
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()
np.random.seed(a.seed)

OUT = a.out_dir; os.makedirs(OUT, exist_ok=True)
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
DBN = "cellchat"
CLUST = "cell_type"

def rss_gb():
    try:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 ** 3)  # macOS: bytes
    except Exception:
        return float("nan")

def git_sha():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return None

# ------------------------------------------------------------------ load once
log("reading", a.h5ad)
adata_all = sc.read_h5ad(a.h5ad)
n_cells_in = adata_all.shape[0]
log(f"loaded {adata_all.shape[0]} x {adata_all.shape[1]}")

# COMMOT's stated input requirement: "non-negative values that reasonably reflect the abundancy
# of signaling molecules". The tutorial's recipe is exactly: adata.raw = adata;
# sc.pp.normalize_total; sc.pp.log1p -- applied to RAW COUNTS. Our .X is already log-normalised
# by an unknown recipe, so we start from layers['counts'] and reproduce the tutorial's two steps.
if a.count_layer != "X":
    if a.count_layer not in adata_all.layers:
        raise SystemExit(f"ERROR: layer '{a.count_layer}' absent; have {list(adata_all.layers)}")
    adata_all.X = adata_all.layers[a.count_layer].copy()
    log(f"X <- layers['{a.count_layer}'] (raw counts)")
adata_all.layers.clear()
if adata_all.raw is not None: adata_all.raw = None
adata_all.raw = adata_all
adata_all.uns.pop('log1p', None)   # stale key -> spurious warning
sc.pp.normalize_total(adata_all, inplace=True)
sc.pp.log1p(adata_all)
log(f"normalize_total + log1p done (tutorial recipe); X max = {adata_all.X.max():.2f}")

if a.cell_type_col not in adata_all.obs:
    raise SystemExit(f"ERROR: '{a.cell_type_col}' not in obs")
adata_all.obs[CLUST] = adata_all.obs[a.cell_type_col].astype(str).astype("category")
adata_all.obsm["spatial"] = np.asarray(adata_all.obsm["spatial"], dtype=float)

# ------------------------------------------------------------- LR database
if a.db == "builtin":
    df_lig = ct.pp.ligand_receptor_database(species="human", database="CellChat",
                                            signaling_type="Secreted Signaling")
    db_label = "CellChatDB v1 (bundled), Secreted Signaling"
else:
    raw = pd.read_csv(a.db)
    keep = [s.strip() for s in a.signaling_types.split(",") if s.strip()]
    sub = raw[raw["signaling_type"].isin(keep)].copy()
    # COMMOT's df_ligrec is exactly 3 columns: ligand, receptor, pathway -- with heteromeric
    # subunits joined by '_'. Our CellChatDB v2 CSV already uses that encoding, so this is a
    # direct, LOSSLESS handover (heteromeric=True + heteromeric_rule='min' aggregates subunits).
    df_lig = pd.DataFrame({0: sub["ligand"].astype(str), 1: sub["receptor"].astype(str),
                           2: sub["pathway"].astype(str)}).reset_index(drop=True)
    df_lig = df_lig.drop_duplicates(subset=[0, 1]).reset_index(drop=True)
    db_label = f"{os.path.basename(a.db)} [{'+'.join(keep)}]"
log(f"LR database: {db_label} -> {df_lig.shape[0]} pairs before expression filtering")

crit = f"min_cell={a.min_cell}" if a.filter_criteria == "min_cell" else f"min_cell_pct={a.min_cell_pct}"

def run_filter(obj):
    return ct.pp.filter_lr_database(df_lig, obj, heteromeric=True, heteromeric_delimiter="_",
                                    heteromeric_rule="min", filter_criteria=a.filter_criteria,
                                    min_cell=a.min_cell, min_cell_pct=a.min_cell_pct)

# --filter-scope global: evaluate the expression filter ONCE on the whole object and hand the
# same pair set to every split. Makes pair sets comparable across splits and stops sparse splits
# from being dropped entirely. See DEVIATIONS.md.
df_global = None
if a.filter_scope == "global":
    df_global = run_filter(adata_all)
    log(f"GLOBAL filter_lr_database({crit}) on all {n_cells_in} cells: "
        f"{df_lig.shape[0]} -> {df_global.shape[0]} pairs, reused for every split")
    df_global.to_csv(os.path.join(OUT, "lr_pairs_global.csv"), index=False)

splits = ([(str(v), np.asarray(adata_all.obs[a.split_col] == v)) for v in
           sorted(adata_all.obs[a.split_col].unique())] if a.split_col
          else [("all", np.ones(adata_all.shape[0], bool))])
# smallest first: COMMOT's dense N x N distance matrix means the big cores are the OOM risk
splits.sort(key=lambda t: t[1].sum())
log(f"splits: {len(splits)} ({a.split_col or 'none'}), smallest first | dis_thr={a.dis_thr} "
    f"| filter_scope={a.filter_scope}")

requested = [s.strip() for s in a.requested_lrs.split(",") if s.strip()]
summary = []
for name, mask in splits:
    n = int(mask.sum())
    if a.max_cells_per_split and n > a.max_cells_per_split:
        log(f"--- split {name}: {n} cells > --max-cells-per-split, SKIPPED ---")
        summary.append(dict(split=name, n_cells=n, status="skipped_too_large")); continue
    ts = time.time()
    sub = adata_all[mask].copy()
    sdir = os.path.join(OUT, name if a.split_col else "")
    PLOTS, DATA = os.path.join(sdir, "plots"), os.path.join(sdir, "data")
    for d in (PLOTS, DATA): os.makedirs(d, exist_ok=True)
    dense_gb = n * n * 8 / 1e9
    log(f"--- split {name}: {n} cells (dense NxN distance matrix = {dense_gb:.2f} GB) ---")

    df_f = df_global if df_global is not None else run_filter(sub)
    log(f"  filter_lr_database({crit}, scope={a.filter_scope}): {df_lig.shape[0]} -> "
        f"{df_f.shape[0]} pairs ({100*df_f.shape[0]/max(1,df_lig.shape[0]):.1f}% retained; "
        f"tutorial retained 20.9%)")
    if df_f.shape[0] == 0:
        summary.append(dict(split=name, n_cells=n, status="no_pairs_after_filter")); continue

    ct.tl.spatial_communication(sub, database_name=DBN, df_ligrec=df_f, dis_thr=a.dis_thr,
                                heteromeric=True, heteromeric_rule="min",
                                heteromeric_delimiter="_", pathway_sum=True)
    log(f"  spatial_communication done | peak RSS {rss_gb():.1f} GB")

    snd = sub.obsm[f"commot-{DBN}-sum-sender"]
    rcv = sub.obsm[f"commot-{DBN}-sum-receiver"]
    snd.to_csv(os.path.join(DATA, "sum_sender.csv.gz"), compression="gzip")
    rcv.to_csv(os.path.join(DATA, "sum_receiver.csv.gz"), compression="gzip")
    df_f.to_csv(os.path.join(DATA, "lr_pairs_used.csv"), index=False)
    pd.DataFrame({"cell": sub.obs_names, "x": sub.obsm["spatial"][:, 0],
                  "y": sub.obsm["spatial"][:, 1],
                  "cell_type": sub.obs[CLUST].astype(str)}).to_csv(
        os.path.join(DATA, "cell_meta.csv"), index=False)
    obsp_keys = [k for k in sub.obsp if k.startswith(f"commot-{DBN}-")]
    obsp_nnz = int(sum(sub.obsp[k].nnz for k in obsp_keys)) if obsp_keys else 0
    log(f"  obsp: {len(obsp_keys)} cell x cell transport plans, {obsp_nnz:,} total nonzeros")

    # Rank by total received signal (COMMOT has NO p-value at the LR-pair level).
    # pathway_sum=True adds PATHWAY columns alongside the per-pair ones, so they must be
    # separated -- ranking a pair against a list that also contains pathway aggregates would
    # inflate the denominator and shift every rank.
    pair_tags = {f"{r[0]}-{r[1]}" for r in df_f.itertuples(index=False)}
    pw_tags = set(df_f[2].astype(str)) - pair_tags
    lr_cols = [c for c in rcv.columns if c.startswith("r-") and c[2:] in pair_tags]
    pw_cols = [c for c in rcv.columns if c.startswith("r-") and c[2:] in pw_tags]
    strength = rcv[lr_cols].sum(axis=0).sort_values(ascending=False)
    strength.rename("total_received").to_csv(os.path.join(DATA, "lr_total_received.csv"))
    pw_strength = (rcv[pw_cols].sum(axis=0).sort_values(ascending=False) if pw_cols
                   else pd.Series(dtype=float))
    if pw_cols:
        pw_strength.rename("total_received").to_csv(
            os.path.join(DATA, "pathway_total_received.csv"))
    ranked = [c[2:] for c in strength.index]
    log(f"  ranked {len(ranked)} LR pairs (+{len(pw_cols)} pathway aggregates kept separate)")

    # ---- cell-type level communication with permutation p-values -------------------------
    # For the top-N PATHWAYS *by total received signal* (the pre-2026-08-10 run selected these
    # in column order, i.e. alphabetically, so its p-values were for near-irrelevant pathways),
    # plus each REQUESTED LR pair -- cluster_communication accepts lr_pair= directly, which is
    # the only significance test COMMOT can give the requested pairs.
    top_pw = [c[2:] for c in pw_strength.index][:a.n_top_pathways]
    tag_of = {f"{r[0]}_{r[1]}": f"{r[0]}-{r[1]}" for r in df_f.itertuples(index=False)}
    lr_of = {f"{r[0]}_{r[1]}": (r[0], r[1]) for r in df_f.itertuples(index=False)}
    log(f"  top {len(top_pw)} pathways by total received signal: {', '.join(top_pw)}")

    cluster_uns, cluster_labels = [], {}
    def do_cluster(label, **kw):
        try:
            ct.tl.cluster_communication(sub, database_name=DBN, clustering=CLUST,
                                        n_permutations=a.n_permutations, random_seed=a.seed, **kw)
            key = (f"commot_cluster-{CLUST}-{DBN}-{kw['pathway_name']}" if "pathway_name" in kw
                   else f"commot_cluster-{CLUST}-{DBN}-{kw['lr_pair'][0]}-{kw['lr_pair'][1]}")
            if key in sub.uns:
                cm = sub.uns[key]["communication_matrix"]
                pv = sub.uns[key]["communication_pvalue"]
                cm.to_csv(os.path.join(DATA, f"cluster_comm_{label}.csv"))
                pv.to_csv(os.path.join(DATA, f"cluster_pval_{label}.csv"))
                cluster_uns.append(key); cluster_labels[label] = (cm, pv)
            else:
                log(f"  cluster_communication({label}): uns key '{key}' absent")
        except Exception as e:
            log(f"  cluster_communication({label}) FAILED: {type(e).__name__}: {e}")

    for pw in top_pw:
        do_cluster(pw, pathway_name=pw)
    for lr in requested:
        if lr in lr_of:
            do_cluster(lr, lr_pair=lr_of[lr])

    # ---- figures --------------------------------------------------------------------------
    def try_plot(fn, nm):
        try:
            fn(); plt.gcf().savefig(os.path.join(PLOTS, nm + ".png"), dpi=200, bbox_inches="tight")
            plt.close("all")
        except Exception as e:
            plt.close("all"); log(f"  plot FAILED {nm}: {type(e).__name__}: {e}")

    def try_native(fn, nm):
        try:
            fn(); plt.close("all")
        except Exception as e:
            plt.close("all"); log(f"  native plot FAILED {nm}: {type(e).__name__}: {e}")

    def ours_dotplot(res, p_cut=0.05):
        """Stand-in for the unrunnable ct.pl.plot_cluster_communication_dotplot.

        Rows = the keys we permutation-tested (top pathways + requested LR pairs); columns =
        sender->receiver cell-type pairs significant in at least one key. Dot area = -log10(p),
        colour = communication strength. Drawn from the same matrices written to
        cluster_comm_*.csv / cluster_pval_*.csv, so it is fully reproducible from disk.
        """
        if not res: raise ValueError("no cluster_communication results")
        cols, cells = [], {}
        for lab, (cm, pv) in res.items():
            for s in cm.index:
                for r in cm.columns:
                    if pv.loc[s, r] < p_cut and cm.loc[s, r] > 0:
                        cols.append(f"{s}->{r}")
                        cells[(lab, f"{s}->{r}")] = (cm.loc[s, r], pv.loc[s, r])
        cols = sorted(set(cols))
        if not cols: raise ValueError(f"no cell-type pair reaches p<{p_cut}")
        labs = list(res)
        fig, ax = plt.subplots(figsize=(max(6, 0.34 * len(cols) + 3), 0.55 * len(labs) + 2.2))
        xs, ys, ss, cs = [], [], [], []
        for i, lab in enumerate(labs):
            for j, c in enumerate(cols):
                if (lab, c) in cells:
                    v, pval = cells[(lab, c)]
                    xs.append(j); ys.append(i); cs.append(v)
                    ss.append(12 + 60 * min(-np.log10(max(pval, 1e-3)) / 3.0, 1.0))
        sc_ = ax.scatter(xs, ys, s=ss, c=cs, cmap="cool", edgecolors="none")
        ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=6)
        ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=8)
        ax.set_title(f"cell-type communication (p<{p_cut}); size=-log10(p), colour=strength",
                     fontsize=9)
        ax.margins(0.02)
        fig.colorbar(sc_, ax=ax, shrink=0.6, label="communication strength")
        fig.tight_layout()

    def sender_receiver_map(tag, nm):
        pts = sub.obsm["spatial"]
        s, r = snd.get(f"s-{tag}"), rcv.get(f"r-{tag}")
        if s is None or r is None: raise KeyError(tag)
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
        ax[0].scatter(pts[:, 0], pts[:, 1], c=s, s=2, cmap="Blues"); ax[0].set_title(f"{tag} sender")
        ax[1].scatter(pts[:, 0], pts[:, 1], c=r, s=2, cmap="Reds");  ax[1].set_title(f"{tag} receiver")
        for x in ax: x.set_aspect("equal"); x.axis("off")

    # (a) the method's own top LR PAIRS -- these are genuine pairs now, not pathway aggregates
    for tag in ranked[:a.n_top_lr]:
        try_plot(lambda tag=tag: sender_receiver_map(tag, tag), f"signal_{tag}")
    # (a2) the top PATHWAY aggregates, labelled as such so they cannot be mistaken for pairs
    for pw in top_pw:
        try_plot(lambda pw=pw: sender_receiver_map(pw, pw), f"pathway_{pw}")
    # (b) the standing requested LRIs. COMMOT names a pair '<ligand>-<receptor>'; complexes keep
    # their internal '_'. Resolve "LIG_REC" against the actual pairs used, not by string-munging.
    req_rank = {}
    for lr in requested:
        t = tag_of.get(lr)
        if t is not None and t in ranked:
            rk = ranked.index(t) + 1
            req_rank[lr] = rk
            log(f"  requested {lr}: rank {rk}/{len(ranked)} by total received signal")
            try_plot(lambda t=t: sender_receiver_map(t, t), f"requested_rank{rk}_{lr}")
        else:
            req_rank[lr] = None
            why = "filtered out by expression" if lr in {f"{x}_{y}" for x, y in zip(df_lig[0], df_lig[1])} \
                  else "absent from this signaling-type subset"
            log(f"  REQUESTED {lr}: not tested -- {why}")

    # ---- native ct.pl.* figures (none of these were produced by the 2026-08-01 run) ---------
    if not a.no_native_plots:
        # vector fields: communication_direction -> plot_cell_communication. background='summary'
        # is the INSTALLED default and needs no H&E image (the earlier "needs background='image'"
        # justification was wrong against the source).
        for pw in top_pw:
            for smry in ("sender", "receiver"):
                try_native(lambda pw=pw, smry=smry: (
                    ct.tl.communication_direction(sub, database_name=DBN, pathway_name=pw, k=5),
                    ct.pl.plot_cell_communication(
                        sub, database_name=DBN, pathway_name=pw, plot_method="cell",
                        background="summary", summary=smry, ndsize=1, scale=1.0,
                        filename=os.path.join(PLOTS, f"native_vf_{smry}_pathway_{pw}.png"))),
                    f"native_vf_{smry}_pathway_{pw}")
        for lr in requested:
            if lr not in lr_of: continue
            for smry in ("sender", "receiver"):
                try_native(lambda lr=lr, smry=smry: (
                    ct.tl.communication_direction(sub, database_name=DBN, lr_pair=lr_of[lr], k=5),
                    ct.pl.plot_cell_communication(
                        sub, database_name=DBN, lr_pair=lr_of[lr], plot_method="cell",
                        background="summary", summary=smry, ndsize=1, scale=1.0,
                        filename=os.path.join(PLOTS, f"native_vf_{smry}_requested_{lr}.png"))),
                    f"native_vf_{smry}_requested_{lr}")
        # cluster-level network + dotplot (nx_node_pos='cluster' needs cluster_position first;
        # the network route goes through networkx.nx_agraph.to_agraph -> needs pygraphviz)
        if cluster_uns:
            try_native(lambda: ct.tl.cluster_position(sub, clustering=CLUST), "cluster_position")
            try_native(lambda: ct.pl.plot_cluster_communication_network(
                sub, uns_names=cluster_uns, clustering=CLUST,
                filename=os.path.join(PLOTS, "native_cluster_network.pdf")),
                "native_cluster_network")
            # ct.pl.plot_cluster_communication_dotplot is unrunnable on this matplotlib/seaborn
            # (see the version-gap note at the top). Equivalent drawn from the persisted CSVs.
            try_plot(lambda: ours_dotplot(cluster_labels), "ours_dotplot_top")

    # ---- persist the AnnData WITH obsp so downstream tools never need another OT run --------
    adata_path = None
    if a.save_adata:
        try:
            sub.raw = None          # reconstructible from --h5ad; halves the file
            adata_path = os.path.join(sdir, "adata_commot.h5ad")
            sub.write_h5ad(adata_path, compression="gzip")
            log(f"  saved AnnData (+obsp) -> {adata_path} "
                f"({os.path.getsize(adata_path)/1e9:.2f} GB)")
        except Exception as e:
            adata_path = None
            log(f"  save_adata FAILED: {type(e).__name__}: {e}")

    rt = (time.time() - ts) / 60
    summary.append(dict(split=name, n_cells=n, status="ok", dense_dist_gb=round(dense_gb, 2),
                        n_pairs_input=int(df_lig.shape[0]), n_pairs_used=int(df_f.shape[0]),
                        n_obsp=len(obsp_keys), obsp_nnz=obsp_nnz,
                        top_pathways=";".join(top_pw),
                        n_cluster_tests=len(cluster_uns),
                        adata_gb=(round(os.path.getsize(adata_path)/1e9, 2) if adata_path else None),
                        # NOTE: getrusage is process-wide and MONOTONIC, so this is a running max
                        # over all splits so far, not this split's own peak.
                        peak_rss_gb_running_max=round(rss_gb(), 2), runtime_min=round(rt, 2),
                        **{f"{k}_rank": v for k, v in req_rank.items()}))
    log(f"  split {name} done in {rt:.1f} min")
    del sub; gc.collect()
    pd.DataFrame(summary).to_csv(os.path.join(OUT, "per_split_summary.csv"), index=False)

pd.DataFrame(summary).to_csv(os.path.join(OUT, "per_split_summary.csv"), index=False)

import scipy, sklearn
json.dump({
    "script": "run_commot.py", "git_sha": git_sha(), "h5ad": a.h5ad,
    "cell_type_col": a.cell_type_col, "count_layer": a.count_layer,
    "n_cells_in": int(n_cells_in), "n_cells_analysed": int(sum(
        d["n_cells"] for d in summary if d.get("status") == "ok")),
    "n_splits_ok": sum(1 for d in summary if d.get("status") == "ok"),
    "n_splits_total": len(summary),
    "db": db_label, "signaling_types": a.signaling_types,
    "n_pairs_input": int(df_lig.shape[0]),
    "dis_thr": a.dis_thr,
    # BOTH filter parameters recorded, plus which criterion was actually in force -- the
    # pre-2026-08-10 manifest recorded only min_cell_pct, which was NOT the one used.
    "filter_criteria": a.filter_criteria, "filter_scope": a.filter_scope,
    "min_cell": a.min_cell, "min_cell_pct": a.min_cell_pct,
    "n_pairs_global": (int(df_global.shape[0]) if df_global is not None else None),
    "n_permutations": a.n_permutations, "cluster_random_seed": a.seed,
    "n_top_lr": a.n_top_lr, "n_top_pathways": a.n_top_pathways,
    "requested_lrs": requested, "split_col": a.split_col,
    "save_adata": bool(a.save_adata), "native_plots": not a.no_native_plots,
    "seed": a.seed,
    "versions": {"python": sys.version.split()[0], "commot": "0.0.3", "numpy": np.__version__,
                 "scipy": scipy.__version__, "pandas": pd.__version__,
                 "scanpy": sc.__version__, "anndata": ad.__version__,
                 "sklearn": sklearn.__version__},
    "peak_rss_gb": round(rss_gb(), 2),
    "wall_min": round((time.time() - t0) / 60, 1)},
    open(os.path.join(OUT, "run_manifest.json"), "w"), indent=2)
log(f"ALL DONE in {(time.time()-t0)/60:.1f} min -> {OUT}")
print(pd.DataFrame(summary).to_string(index=False))
