#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cpdb_module_scores_kw_dunn.py

Compute CPDB module scores (ligand / receptor / union) across cell types using:
1) Single pass of scanpy.tl.score_genes on ALL CELLS (to keep scores comparable).
2) Distribution plots (violin/box) grouped by cell_type.
3) Statistics:
   - Global multi-group difference: Kruskal–Wallis
   - Post-hoc pairwise: Dunn test with FDR BH correction
4) A small clustered heatmap of per–cell_type summary (mean or median).

NO single-gene or per-entry comparisons are performed.

Dependencies:
  - scanpy, anndata, numpy, pandas, seaborn, matplotlib, scipy
  - scikit-posthocs  (pip install scikit-posthocs)

Usage example:
python cpdb_module_scores_kw_dunn.py \
  --h5ad /path/to/data.h5ad \
  --celltype-col cell_type \
  --cellphone /path/to/cellphonedb/interactions.csv \
  --outdir /path/to/outdir \
  --plot-type violin \
  --summary-stat mean \
  --random-state 0
"""

import os
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import anndata as ad
import scanpy as sc
from scipy.stats import kruskal
import scikit_posthocs as sp  # Dunn posthoc

# ----------------------------
# Utilities
# ----------------------------
def _split_multi(s: str):
    """Split a CPDB entry (multi-subunit joined by '_') into a list of gene symbols."""
    if pd.isna(s) or s == "":
        return []
    return [g.strip() for g in str(s).split("_") if g.strip()]

def _ensure_outdir(p: str):
    """Create directory if it does not exist."""
    Path(p).mkdir(parents=True, exist_ok=True)

def _resolve_cpdb_path(p: str) -> Path:
    """Resolve CPDB interactions CSV file; allow passing a directory with common filenames."""
    p = Path(p)
    if p.is_file():
        return p
    for name in ["interactions.csv", "interactions_filtered.csv", "cellphonedb_interactions.csv"]:
        cand = p / name
        if cand.exists():
            return cand
    raise FileNotFoundError("No CellPhoneDB interactions CSV found (tried common filenames).")

def _read_cpdb_gene_sets(cp_csv: Path):
    """Read CPDB interactions and construct ligand, receptor, and union gene lists (unique symbols)."""
    df = pd.read_csv(cp_csv)
    colmap = {c.lower(): c for c in df.columns}
    if "ligand" not in colmap or "receptor" not in colmap:
        raise ValueError("CellPhoneDB CSV must contain 'ligand' and 'receptor' columns.")
    lig_col, rec_col = colmap["ligand"], colmap["receptor"]

    lig_lists = df[lig_col].apply(_split_multi)
    rec_lists = df[rec_col].apply(_split_multi)

    lig_genes = sorted(set(g for lst in lig_lists for g in lst))
    rec_genes = sorted(set(g for lst in rec_lists for g in lst))
    union_genes = sorted(set(lig_genes).union(rec_genes))
    return lig_genes, rec_genes, union_genes

def _score_one(
    adata: ad.AnnData,
    genes: list,
    score_name: str,
    random_state: int = 0,
    n_bins: int = 25,
    ctrl_size: int = 50,
):
    """
    Run scanpy.tl.score_genes on ALL CELLS ONCE for a given gene set.
    ctrl_size is fixed to an int (e.g., 50) to be compatible with legacy scanpy.
    """
    genes_present = [g for g in genes if g in adata.var_names]
    if len(genes_present) == 0:
        adata.obs[score_name] = np.nan
        return genes_present

    sc.tl.score_genes(
        adata,
        gene_list=genes_present,
        score_name=score_name,
        use_raw=False,
        random_state=random_state,
        n_bins=n_bins,
        ctrl_size=ctrl_size,  # must be int for legacy versions
    )
    return genes_present

def _kruskal_and_dunn(df: pd.DataFrame, group_col: str, value_col: str, p_adjust: str = "fdr_bh"):
    """
    Perform Kruskal–Wallis test across groups and Dunn post-hoc comparisons with FDR BH.
    Returns:
      - H (float): Kruskal–Wallis statistic
      - p (float): Kruskal–Wallis p-value
      - dunn_long (DataFrame): long-form pairwise adjusted p-values [group1, group2, p_adj]
    Notes:
      * Uses 'stack' to build long table from the Dunn matrix to avoid index-name pitfalls.
      * Casts labels to str to avoid 'Unordered Categoricals' comparison issues.
    """
    # Kruskal–Wallis: use observed=False to silence future warnings and keep current behavior
    groups = [sub[value_col].dropna().values for _, sub in df.groupby(group_col, observed=False)]
    if len(groups) >= 2 and all(len(g) > 0 for g in groups):
        H, p = kruskal(*groups)
    else:
        H, p = np.nan, np.nan

    # Dunn posthoc with FDR BH
    if df[group_col].nunique() >= 2 and df[value_col].notna().sum() > 0:
        # Returns a square DataFrame with row/col labels = group labels
        dunn_mat = sp.posthoc_dunn(df, val_col=value_col, group_col=group_col, p_adjust=p_adjust)

        # Ensure both axes are strings (avoid categorical comparison issues)
        dunn_mat.index = dunn_mat.index.astype(str)
        dunn_mat.columns = dunn_mat.columns.astype(str)

        # Build long form via 'stack'; drop diagonal (NaNs or zeros not guaranteed)
        dunn_long = (
            dunn_mat.where(~np.eye(len(dunn_mat), dtype=bool))
                    .stack()
                    .reset_index()
        )
        dunn_long.columns = ["group1", "group2", "p_adj"]

        # Keep only upper triangle using a stable alphabetical order (or customize if needed)
        order = {k: i for i, k in enumerate(sorted(dunn_mat.index))}
        dunn_long = dunn_long[
            dunn_long["group1"].map(order) < dunn_long["group2"].map(order)
        ].reset_index(drop=True)
    else:
        dunn_long = pd.DataFrame(columns=["group1", "group2", "p_adj"])

    return H, p, dunn_long

# ----------------------------
# Main
# ----------------------------
def main(args):
    _ensure_outdir(args.outdir)
    figs_dir = os.path.join(args.outdir, "figs")
    _ensure_outdir(figs_dir)

    # Load AnnData
    adata = ad.read_h5ad(args.h5ad)
    if args.celltype_col not in adata.obs.columns:
        raise ValueError(f"{args.celltype_col} not found in adata.obs")

    # Resolve CPDB gene sets
    cp_csv = _resolve_cpdb_path(args.cellphone)
    lig_genes, rec_genes, union_genes = _read_cpdb_gene_sets(cp_csv)

    # Score on ALL CELLS ONCE (keeps scores on same scale across cell types)
    print("Scoring: ligand...")
    _ = _score_one(adata, lig_genes, "cpdb_ligand_score", random_state=args.random_state)
    print("Scoring: receptor...")
    _ = _score_one(adata, rec_genes, "cpdb_receptor_score", random_state=args.random_state)
    print("Scoring: union...")
    _ = _score_one(adata, union_genes, "cpdb_union_score", random_state=args.random_state)

    # Save obs with scores for reproducibility
    adata.obs.to_csv(os.path.join(args.outdir, "cell_scores_obs.csv"))

    # Tidy long table: one row per cell per score type
    score_cols = ["cpdb_ligand_score", "cpdb_receptor_score", "cpdb_union_score"]
    long = (
        adata.obs[[args.celltype_col] + score_cols]
        .melt(id_vars=[args.celltype_col], var_name="score_type", value_name="score")
        .dropna(subset=["score"])
    )

    # Statistics: Kruskal–Wallis + Dunn(FDR BH), and distribution plots
    kw_rows = []
    dunn_all = []
    for s in score_cols:
        df = long[long["score_type"] == s].copy()

        # Stats
        H, p, dunn_long = _kruskal_and_dunn(df, args.celltype_col, "score", p_adjust="fdr_bh")
        kw_rows.append({"score_type": s, "kw_H": H, "kw_pvalue": p})
        if not dunn_long.empty:
            dunn_long.insert(0, "score_type", s)
            dunn_all.append(dunn_long)

        # Plot distributions with global p-value annotated
        plt.figure(figsize=(max(7, 0.7 * df[args.celltype_col].nunique()), 5))
        if args.plot_type == "box":
            sns.boxplot(data=df, x=args.celltype_col, y="score", showfliers=False)
        else:
            sns.violinplot(data=df, x=args.celltype_col, y="score", cut=0, inner="box")
        sns.stripplot(data=df, x=args.celltype_col, y="score", alpha=0.25, size=2)

        title = s.replace("cpdb_", "").replace("_score", "").capitalize()
        if np.isfinite(p):
            plt.title(f"CPDB {title} module score  |  Kruskal–Wallis p = {p:.2e}")
        else:
            plt.title(f"CPDB {title} module score")

        plt.xlabel("cell_type")
        plt.ylabel("Module score (scanpy.tl.score_genes)")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        base = f"{s}_{args.plot_type}"
        plt.savefig(os.path.join(figs_dir, f"{base}.png"), dpi=300, bbox_inches="tight")
        plt.savefig(os.path.join(figs_dir, f"{base}.svg"), bbox_inches="tight")
        plt.close()

    # Save stats
    kw_df = pd.DataFrame(kw_rows)
    kw_df.to_csv(os.path.join(args.outdir, "kruskal_wallis_overall.csv"), index=False)
    if dunn_all:
        dunn_df = pd.concat(dunn_all, ignore_index=True)
        dunn_df.to_csv(os.path.join(args.outdir, "dunn_pairwise_fdrbh.csv"), index=False)

    # Per–cell_type summary heatmap (mean or median of scores)
    agg_func = np.mean if args.summary_stat == "mean" else np.median
    mat = (
        long.pivot_table(
            index=args.celltype_col, columns="score_type", values="score",
            aggfunc=agg_func, observed=False  # keep current pandas behavior
        ).fillna(0.0)
    )

    # Column-wise z-score for visualization
    mat_z = (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-9)
    cg = sns.clustermap(
        mat_z, cmap="vlag", center=0,
        figsize=(6, max(4, 0.6 * mat_z.shape[0])),
        cbar_kws={"label": f"Z-scored {args.summary_stat} score"},
    )
    out_png = os.path.join(figs_dir, f"celltype_x_cpdb_scores_{args.summary_stat}_clustermap.png")
    out_svg = os.path.join(figs_dir, f"celltype_x_cpdb_scores_{args.summary_stat}_clustermap.svg")
    cg.savefig(out_png, dpi=300, bbox_inches="tight")
    cg.savefig(out_svg, bbox_inches="tight")
    plt.close("all")

    print("✅ Done.")
    print(f"Figures: {figs_dir}")
    print(f"KW table: {os.path.join(args.outdir, 'kruskal_wallis_overall.csv')}")
    if os.path.exists(os.path.join(args.outdir, 'dunn_pairwise_fdrbh.csv')):
        print(f"Dunn pairwise (FDR BH): {os.path.join(args.outdir, 'dunn_pairwise_fdrbh.csv')}")

# ----------------------------
# CLI
# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CPDB module scores across cell types with Kruskal–Wallis + Dunn(FDR BH) (no single-gene comparisons)."
    )
    parser.add_argument("--h5ad", required=True, help="Path to AnnData .h5ad file.")
    parser.add_argument("--celltype-col", required=True, help="Column in adata.obs with cell type labels.")
    parser.add_argument("--cellphone", required=True, help="Path to CellPhoneDB interactions CSV or a directory containing it.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--plot-type", choices=["violin", "box"], default="violin", help="Distribution plot style.")
    parser.add_argument("--summary-stat", choices=["mean", "median"], default="mean", help="Heatmap aggregation per cell type.")
    parser.add_argument("--random-state", type=int, default=0, help="Random state for score_genes control genes.")
    args = parser.parse_args()
    main(args)
