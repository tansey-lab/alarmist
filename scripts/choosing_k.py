#!/usr/bin/env python3
"""
Ward clustering with spatial connectivity and model selection via:
- Inertia (SSE)
- Gaussian BIC (PI-style)
- Gaussian AICc (PI-style simplified)
- Mixture IC (BIC/AICc blend)
- Calinski–Harabasz (CH) Index

Also computes & plots the Ward **merging cost** curve by fitting the full
hierarchical tree once (n_clusters=None, compute_distances=True).

Outputs (in --outdir):
- ward_ic_scores.csv  (k, inertia_sse, bic, aicc, mixture_ic, ch_index)
- seg_labels_k{k}.npy and seg_labels_k{k}.csv  (cached labels for every k)
- best_k_{k}_labels.csv  (with cell_id if available)
- ward_ic_curves.png    (Inertia/BIC/AICc/Mixture IC)
- ward_ch_curve.png     (CH index vs k; higher is better)
- ward_merging_cost.csv (step, n_clusters, merge_cost)
- ward_merging_cost_curve.png

Example:
python ward_ic_runner.py \
  --u-cell-npy /path/to/cell_loadings.npy \
  --cell-lri-dir /path/to/results/single_cell_neighborhood_lri \
  --graph radius --radius 30.0 \
  --k-min 2 --k-max 20 \
  --bic-proportion 0.5 \
  --outdir ward_ic_outputs_radius30
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances_argmin_min
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from scipy.sparse.csgraph import connected_components


# ------------------------------ IO helpers ------------------------------

def load_inputs(u_cell_npy: str,
                cell_lri_dir: str,
                metadata_csv: str):
    """
    Load U_cell and metadata_df (must contain x_coord, y_coord; optionally cell_id).
    Priority for metadata:
      1) --metadata-csv (CSV file with columns)
      2) --cell-lri-dir via neighborhood_lri_analysis.load_cell_lri_results
    """
    U_cell = np.load(u_cell_npy)
    if metadata_csv is not None:
        metadata_df = pd.read_csv(metadata_csv)
    else:
        if cell_lri_dir is None:
            raise ValueError("Either --metadata-csv or --cell-lri-dir must be provided.")
        # Lazy import only if used
        from neighborhood_lri_analysis import load_cell_lri_results
        results = load_cell_lri_results(cell_lri_dir)
        metadata_df = results["cell_metadata_df"]
    # Basic checks
    if len(metadata_df) != U_cell.shape[0]:
        raise ValueError(f"metadata_df rows ({len(metadata_df)}) != U_cell rows ({U_cell.shape[0]}).")
    for col in ["x_coord", "y_coord"]:
        if col not in metadata_df.columns:
            raise ValueError(f"metadata_df is missing required column '{col}'.")
    return U_cell, metadata_df


# ------------------------------ scoring (PI-style + CH) ------------------------------

def _sse_from_labels(X: np.ndarray, labels: np.ndarray) -> float:
    """Sum of squared distances to cluster centers (a.k.a. inertia/SSE)."""
    centers = np.vstack([X[labels == k].mean(axis=0) for k in np.unique(labels)])
    _, dists = pairwise_distances_argmin_min(X, centers)
    return float((dists ** 2).sum())

def gaussian_bic_from_sse(sse: float, n_samples: int, k: int) -> float:
    """PI’s simplified BIC: k*log(n) + 2*sse"""
    return k * np.log(n_samples) + 2.0 * sse

def gaussian_aicc_from_sse(sse: float, n_samples: int, k: int) -> float:
    """PI’s simplified AICc: 2*k + 2*sse (same as PI's code, small-sample correction omitted)"""
    return 2.0 * k + 2.0 * sse

def gaussian_mixture_ic_from_sse(sse: float, n_samples: int, k: int, bic_proportion: float = 0.5) -> tuple[float, float, float]:
    """Weighted mixture of BIC and AICc (default 50/50), returns (mixture, bic, aicc)."""
    bic = gaussian_bic_from_sse(sse, n_samples, k)
    aicc = gaussian_aicc_from_sse(sse, n_samples, k)
    mix = bic_proportion * bic + (1.0 - bic_proportion) * aicc
    return mix, bic, aicc

def ch_index_from_labels(X: np.ndarray, labels: np.ndarray, precomputed_total_ss: float = None) -> float:
    """
    Calinski–Harabasz index:
      CH(k) = (B_k / (k-1)) / (W_k / (n-k))
    where W_k = within-cluster SSE, B_k = total_ss - W_k,
    total_ss = sum ||x - global_mean||^2.
    """
    n = X.shape[0]
    k = len(np.unique(labels))
    if k < 2 or k >= n:
        return np.nan
    # W_k
    w_ss = _sse_from_labels(X, labels)
    # total_ss
    if precomputed_total_ss is None:
        xbar = X.mean(axis=0, keepdims=True)
        total_ss = float(((X - xbar) ** 2).sum())
    else:
        total_ss = float(precomputed_total_ss)
    b_ss = total_ss - w_ss
    # CH
    return (b_ss / (k - 1)) / (w_ss / (n - k))


# ------------------------------ merging cost (Ward distances) ------------------------------

def compute_merging_cost_curve(X: np.ndarray, connectivity, outdir: str):
    """
    Fit full Ward tree once (n_clusters=None, distance_threshold=0) with compute_distances=True
    to extract the sequence of merge 'distances' == Ward's linkage costs (increase in within-SS).
    Saves CSV and a plot.
    """
    csv_path = os.path.join(outdir, "ward_merging_cost.csv")
    png_path = os.path.join(outdir, "ward_merging_cost_curve.png")

    # Fit full tree once
    try:
        model = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=0.0,
            linkage="ward",
            connectivity=connectivity,
            compute_distances=True,  # sklearn >= 0.22 supports this flag
        )
    except TypeError:
        # Fallback (older sklearn): compute_distances not available
        model = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=0.0,
            linkage="ward",
            connectivity=connectivity,
        )

    model.fit(X)

    if not hasattr(model, "distances_"):
        # If we can't access distances, save a note and return quietly
        pd.DataFrame(
            {"note": ["distances_ not available in this sklearn version; upgrade to use merging curve."]}
        ).to_csv(csv_path, index=False)
        return

    # distances_: length n-1; step t merges two clusters -> new cost
    # When there are n samples, the resulting number of clusters after step t is: n - t
    n = X.shape[0]
    steps = np.arange(1, len(model.distances_) + 1)
    n_clusters = n - steps
    merge_cost = model.distances_.astype(float)

    df = pd.DataFrame({"step": steps, "n_clusters": n_clusters, "merge_cost": merge_cost})
    df.to_csv(csv_path, index=False)

    # Plot as merge_cost vs n_clusters (descending k)
    plt.figure(figsize=(7.5, 5))
    plt.plot(df["n_clusters"], df["merge_cost"], marker=".", linestyle="-")
    plt.gca().invert_xaxis()  # show from large k -> small k
    plt.xlabel("Number of clusters (k)")
    plt.ylabel("Ward merge cost")
    plt.title("Merging cost curve (Ward linkage distances)")
    plt.tight_layout()
    plt.savefig(png_path, dpi=300)
    plt.close()


# ------------------------------ main runner ------------------------------

def ward_grid_with_metrics(
    X: np.ndarray,
    connectivity,
    k_min: int,
    k_max: int,
    bic_proportion: float,
    outdir: str,
    cell_ids: np.ndarray,
    overwrite: bool = False,
):
    """
    Run Ward for k in [k_min, k_max], compute inertia/BIC/AICc/mixture/CH,
    cache per-k labels, save scores CSV and curves plots, and return best-k & labels (by mixture IC).
    """
    os.makedirs(outdir, exist_ok=True)
    n = X.shape[0]
    rows = []

    # Precompute total sum of squares for CH
    xbar = X.mean(axis=0, keepdims=True)
    total_ss = float(((X - xbar) ** 2).sum())

    best = {"k": None, "ic": np.inf, "labels": None}

    for k in range(k_min, k_max + 1):
        labels_npy = os.path.join(outdir, f"seg_labels_k{k}.npy")
        labels_csv = os.path.join(outdir, f"seg_labels_k{k}.csv")

        if os.path.exists(labels_npy) and not overwrite:
            labels = np.load(labels_npy)
        else:
            ward = AgglomerativeClustering(n_clusters=k, linkage="ward", connectivity=connectivity)
            labels = ward.fit_predict(X)
            np.save(labels_npy, labels)

        # Save per-k CSV (with cell_id if present)
        if (not os.path.exists(labels_csv)) or overwrite:
            if cell_ids is not None:
                pd.DataFrame({"cell_id": cell_ids, "segment_id": labels}).to_csv(labels_csv, index=False)
            else:
                pd.DataFrame({"segment_id": labels}).to_csv(labels_csv, index=False)

        # Scores
        sse = _sse_from_labels(X, labels)
        mix, bic, aicc = gaussian_mixture_ic_from_sse(sse, n_samples=n, k=k, bic_proportion=bic_proportion)
        ch = ch_index_from_labels(X, labels, precomputed_total_ss=total_ss)

        rows.append({
            "k": k,
            "inertia_sse": sse,
            "bic": bic,
            "aicc": aicc,
            "mixture_ic": mix,
            "ch_index": ch,
        })

        if mix < best["ic"]:
            best.update({"k": k, "ic": mix, "labels": labels})

    # Save scores CSV
    scores_df = pd.DataFrame(rows).sort_values("k")
    scores_csv = os.path.join(outdir, "ward_ic_scores.csv")
    scores_df.to_csv(scores_csv, index=False)

    # Save best-k labels (CSV)
    best_labels_csv = os.path.join(outdir, f"best_k_{best['k']}_labels.csv")
    if cell_ids is not None:
        pd.DataFrame({"cell_id": cell_ids, "segment_id": best["labels"]}).to_csv(best_labels_csv, index=False)
    else:
        pd.DataFrame({"segment_id": best["labels"]}).to_csv(best_labels_csv, index=False)

    # Plot curves: Inertia/BIC/AICc/Mixture IC
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    ax.plot(scores_df["k"], scores_df["inertia_sse"], marker="o", label="Inertia (SSE)")
    ax.plot(scores_df["k"], scores_df["bic"], marker="o", label="Gaussian BIC")
    ax.plot(scores_df["k"], scores_df["aicc"], marker="o", label="Gaussian AICc")
    ax.plot(scores_df["k"], scores_df["mixture_ic"], marker="o", linewidth=2.6,
            label=f"Mixture IC (BIC weight = {bic_proportion:.2f})")
    ybest = scores_df.loc[scores_df["k"] == best["k"], "mixture_ic"].values[0]
    ax.scatter([best["k"]], [ybest], color="red", zorder=5, label=f"Best k (mixture) = {best['k']}")
    ax.set_xlabel("Number of clusters (k)")
    ax.set_ylabel("Score (lower is better)")
    ax.set_title("Ward model selection: Inertia vs IC")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "ward_ic_curves.png"), dpi=300)
    plt.close(fig)

    # Plot CH index separately (higher is better)
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(scores_df["k"], scores_df["ch_index"], marker="o")
    # mark max CH (ignoring NaNs)
    valid = scores_df["ch_index"].replace([np.inf, -np.inf], np.nan).dropna()
    if not valid.empty:
        k_ch_best = int(scores_df.loc[scores_df["ch_index"] == valid.max(), "k"].iloc[0])
        y_ch_best = float(valid.max())
        ax.scatter([k_ch_best], [y_ch_best], color="orange", zorder=5, label=f"Best k (CH) = {k_ch_best}")
    ax.set_xlabel("Number of clusters (k)")
    ax.set_ylabel("CH index (higher is better)")
    ax.set_title("Calinski–Harabasz (CH) index vs k")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "ward_ch_curve.png"), dpi=300)
    plt.close(fig)

    return best["k"], best["labels"], scores_df


# ------------------------------ graph builder ------------------------------

def build_connectivity(coords: np.ndarray,
                       graph: str,
                       n_neighbors: int,
                       radius: float,
                       include_self: bool = False):
    """
    Build a sparse connectivity graph from 2D coords.
    graph: "knn" or "radius"
    """
    if graph == "knn":
        print("[Graph] using kNN graph with k =", n_neighbors)
        G = kneighbors_graph(coords, n_neighbors=n_neighbors, mode="connectivity", include_self=include_self)
    elif graph == "radius":
        print("[Graph] using radius graph with radius =", radius)
        G = radius_neighbors_graph(coords, radius=radius, mode="connectivity", include_self=include_self)
    else:
        raise ValueError("--graph must be 'knn' or 'radius'")
    n_comp, _ = connected_components(G)
    print(f"[Graph] connected components: {n_comp}")
    return G


# ------------------------------ CLI ------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Ward clustering with BIC/AICc/CH model selection + merging-cost curve.")
    p.add_argument("--u-cell-npy", required=False, default='/Users/jiayifan/Desktop/Lab/alarmist/results/single_cell_bptf/cell_loadings.npy',
                   help="Path to U_cell .npy (n_cells x K).")
    p.add_argument("--cell-lri-dir", required=False, default='/Users/jiayifan/Desktop/Lab/alarmist/results/single_cell_neighborhood_lri',
                   help="Directory for load_cell_lri_results (to get metadata_df).")
    p.add_argument("--metadata-csv", default=None,
                   help="Alternative to --cell-lri-dir: CSV with at least ['x_coord','y_coord'] and optionally 'cell_id'.")
    p.add_argument("--graph", choices=["knn", "radius"], default="radius", help="Graph type for spatial connectivity.")
    p.add_argument("--n-neighbors", type=int, default=15, help="k for kNN graph (if --graph knn).")
    p.add_argument("--radius", type=float, default=50.0, help="Radius for radius graph (if --graph radius).")
    p.add_argument("--k-min", type=int, default=2, help="Minimum number of clusters.")
    p.add_argument("--k-max", type=int, default=20, help="Maximum number of clusters.")
    p.add_argument("--bic-proportion", type=float, default=0.5,
                   help="Weight of BIC in the mixture IC (the rest is AICc).")
    p.add_argument("--outdir", default='/Users/jiayifan/Desktop/Lab/alarmist/results/single_cell_ward', help="Output directory to save results.")
    p.add_argument("--overwrite", action="store_true", help="Recompute and overwrite cached labels/CSVs.")
    return p.parse_args()


def main():
    args = parse_args()

    print("[IO] Loading inputs…")
    U_cell, metadata_df = load_inputs(args.u_cell_npy, args.cell_lri_dir, args.metadata_csv)
    X = U_cell.astype(float)  # (n_cells x K)
    coords = metadata_df[["x_coord", "y_coord"]].to_numpy()

    # Optional: standardize features (can help Ward with Euclidean distance)
    # X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-8)

    print("[Graph] Building spatial connectivity…")
    G = build_connectivity(coords, graph=args.graph, n_neighbors=args.n_neighbors,
                           radius=args.radius, include_self=False)

    os.makedirs(args.outdir, exist_ok=True)
    cell_ids = metadata_df["cell_id"].to_numpy() if "cell_id" in metadata_df.columns else None

    print("[Ward] Fitting full tree once for merging-cost curve…")
    compute_merging_cost_curve(X, G, args.outdir)

    print("[Ward] Running k-grid & computing all criteria (Inertia/BIC/AICc/Mixture/CH)…")
    best_k, best_labels, scores_df = ward_grid_with_metrics(
        X=X,
        connectivity=G,
        k_min=args.k_min,
        k_max=args.k_max,
        bic_proportion=args.bic_proportion,
        outdir=args.outdir,
        cell_ids=cell_ids,
        overwrite=args.overwrite,
    )

    print(f"[Done] Best k by Mixture IC = {best_k}")
    print(f"[Saved] Scores CSV: {os.path.join(args.outdir, 'ward_ic_scores.csv')}")
    print(f"[Saved] Curves: {os.path.join(args.outdir, 'ward_ic_curves.png')} and {os.path.join(args.outdir, 'ward_ch_curve.png')}")
    print(f"[Saved] Merging-cost CSV/plot: {os.path.join(args.outdir, 'ward_merging_cost.csv')} / ward_merging_cost_curve.png")
    print(f"[Saved] Best labels: {os.path.join(args.outdir, f'best_k_{best_k}_labels.csv')}")
    print(f"[Cache] Per-k labels saved as seg_labels_k*.npy and seg_labels_k*.csv in {args.outdir}")


if __name__ == "__main__":
    main()
