#!/usr/bin/env python3
"""
Build neighborhood membership (cell -> neighborhood_id) using an LRI-guided random walk.
FIXED: Always use adata.obs.index instead of cell_id to avoid duplicates.
"""

# ----------------------------------------------------------------------
# Environment: avoid multi-OpenMP conflicts on macOS/conda stacks
# ----------------------------------------------------------------------
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

# ----------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------
import argparse
from typing import List, Tuple
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse as sp
from numba import njit
import warnings

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def _to_list_genes(s: str) -> List[str]:
    """Split gene string into list. Treat underscores as multi-subunit separators."""
    s = str(s).strip()
    if "_" in s:
        return [g.strip() for g in s.split("_") if g.strip()]
    return [s]

def _prepare_lr_table(lr_df: pd.DataFrame, adata: ad.AnnData) -> pd.DataFrame:
    """
    Keep only LR rows whose all ligand & receptor genes are present in adata.
    Adds: lig_list, rec_list
    """
    df = lr_df.copy()
    assert {"ligand", "receptor"}.issubset(df.columns), "LR CSV must have 'ligand' and 'receptor' columns"
    df["lig_list"] = df["ligand"].apply(_to_list_genes)
    df["rec_list"] = df["receptor"].apply(_to_list_genes)

    gene_set = set(map(str, adata.var_names))
    keep = []
    for _, row in df.iterrows():
        lig_ok = all(g in gene_set for g in row["lig_list"])
        rec_ok = all(g in gene_set for g in row["rec_list"])
        keep.append(lig_ok and rec_ok)
    return df.loc[keep].reset_index(drop=True)

# ----------------------------------------------------------------------
# Expression masks (ligand AND, receptor AND)
# ----------------------------------------------------------------------
def compute_lr_masks_optimized(
    adata: ad.AnnData, 
    lr_df: pd.DataFrame,
    expr_min: int = 0
) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Optimized computation of ligand/receptor expression masks.
    - Ligand: ALL genes must be expressed
    - Receptor: ALL genes must be expressed
    
    Returns:
        lig_can_send: (n_cells, n_LR) bool
        rec_can_recv: (n_cells, n_LR) bool
        lr_df_filtered: filtered LR dataframe
    """
    lr_df_filtered = _prepare_lr_table(lr_df, adata)
    if len(lr_df_filtered) == 0:
        return np.zeros((adata.n_obs, 0), dtype=bool), np.zeros((adata.n_obs, 0), dtype=bool), lr_df_filtered
    
    n_cells = adata.n_obs
    n_lr = len(lr_df_filtered)
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}
    
    # index lists per LR
    all_genes_needed = set()
    lig_idx_per_lr, rec_idx_per_lr = [], []
    for _, row in lr_df_filtered.iterrows():
        lig_idx = [gene_to_idx[g] for g in row["lig_list"] if g in gene_to_idx]
        rec_idx = [gene_to_idx[g] for g in row["rec_list"] if g in gene_to_idx]
        # if anything missing, mark empty (no one can satisfy)
        if len(lig_idx) != len(row["lig_list"]) or len(rec_idx) != len(row["rec_list"]):
            lig_idx, rec_idx = [], []
        lig_idx_per_lr.append(lig_idx)
        rec_idx_per_lr.append(rec_idx)
        all_genes_needed.update(lig_idx)
        all_genes_needed.update(rec_idx)

    if len(all_genes_needed) == 0:
        warnings.warn("No valid genes found in LR database after filtering")
        return np.zeros((n_cells, 0), dtype=bool), np.zeros((n_cells, 0), dtype=bool), lr_df_filtered
    
    # one-shot extraction to dense bool
    gene_list = sorted(all_genes_needed)
    print(f"Extracting expression for {len(gene_list)} unique genes...")
    if sp.issparse(adata.X):
        X_bool = adata.X[:, gene_list].toarray() > expr_min
    else:
        X_bool = np.asarray(adata.X[:, gene_list]) > expr_min
    g2loc = {g: i for i, g in enumerate(gene_list)}

    lig_can_send = np.zeros((n_cells, n_lr), dtype=bool)
    rec_can_recv = np.zeros((n_cells, n_lr), dtype=bool)

    print(f"Computing masks for {n_lr} LR pairs...")
    for j, (L_idx, R_idx) in enumerate(zip(lig_idx_per_lr, rec_idx_per_lr)):
        if L_idx:
            cols = [g2loc[g] for g in L_idx]
            lig_can_send[:, j] = np.all(X_bool[:, cols], axis=1)
        if R_idx:
            cols = [g2loc[g] for g in R_idx]
            rec_can_recv[:, j] = np.all(X_bool[:, cols], axis=1)
        if (j + 1) % 200 == 0:
            print(f"  Processed {j + 1}/{n_lr} LR pairs...")
    return lig_can_send, rec_can_recv, lr_df_filtered

# ----------------------------------------------------------------------
# CSR/CSC adjacency builders
# ----------------------------------------------------------------------
def build_sender_csr(lig_can_send: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build CSR-like adjacency for senders.
    Returns:
        indices: concatenated LR indices
        indptr: length n_cells+1, row slice is indices[indptr[i]:indptr[i+1]]
    """
    rows, cols = np.nonzero(lig_can_send)
    n_cells = lig_can_send.shape[0]
    indptr = np.zeros(n_cells + 1, dtype=np.int64)
    np.add.at(indptr, rows + 1, 1)
    np.cumsum(indptr, out=indptr)
    indices = np.empty(len(cols), dtype=np.int32)
    counter = indptr.copy()
    for r, c in zip(rows, cols):
        pos = counter[r]
        indices[pos] = c
        counter[r] += 1
    return indices, indptr

def build_receiver_csc(rec_can_recv: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build CSC-like adjacency for receivers.
    Returns:
        indices: concatenated cell indices
        indptr: length n_lr+1, column slice is indices[indptr[j]:indptr[j+1]]
    """
    rows, cols = np.nonzero(rec_can_recv)
    n_lr = rec_can_recv.shape[1]
    indptr = np.zeros(n_lr + 1, dtype=np.int64)
    np.add.at(indptr, cols + 1, 1)
    np.cumsum(indptr, out=indptr)
    indices = np.empty(len(rows), dtype=np.int32)
    counter = indptr.copy()
    for r, c in zip(rows, cols):
        pos = counter[c]
        indices[pos] = r
        counter[c] += 1
    return indices, indptr

# ----------------------------------------------------------------------
# Numba kernel using CSR/CSC adjacency
# ----------------------------------------------------------------------
@njit(cache=True, parallel=False)
def build_neighborhoods_fast_csr(
    sender_idx: np.ndarray, sender_ptr: np.ndarray,   # CSR on cells -> LRs
    recv_idx: np.ndarray, recv_ptr: np.ndarray,       # CSC on LRs   -> cells
    neighborhood_size: int = 20,
    seed: int = 0
) -> Tuple[np.ndarray, int]:
    """
    Random walk over CSR/CSC adjacency.
    - Global no-replacement for cells.
    - LR sampling is with replacement.
    """
    np.random.seed(seed)
    n_cells = sender_ptr.size - 1
    assignments = np.full(n_cells, -1, dtype=np.int32)
    nb_id = 0

    # cells that can send at least one LR
    can_send = np.zeros(n_cells, dtype=np.bool_)
    for i in range(n_cells):
        can_send[i] = (sender_ptr[i+1] > sender_ptr[i])
    sender_cells = np.where(can_send)[0]
    # shuffle seeds
    for i in range(sender_cells.size - 1, 0, -1):
        j = np.random.randint(i + 1)
        tmp = sender_cells[i]
        sender_cells[i] = sender_cells[j]
        sender_cells[j] = tmp

    for s in range(sender_cells.size):
        seed_cell = sender_cells[s]
        if assignments[seed_cell] != -1:
            continue

        # start new neighborhood
        nb_cells = np.zeros(neighborhood_size, dtype=np.int32)
        nb_size = 0
        nb_cells[nb_size] = seed_cell
        assignments[seed_cell] = nb_id
        nb_size += 1

        current = seed_cell
        attempts = 0
        max_attempts = neighborhood_size * 20

        while nb_size < neighborhood_size and attempts < max_attempts:
            attempts += 1

            # sender LR slice
            s_beg = sender_ptr[current]
            s_end = sender_ptr[current + 1]
            if s_end == s_beg:
                # pick another sender inside the neighborhood
                found_sender = False
                if nb_size > 0:
                    start = np.random.randint(nb_size)
                    for k in range(nb_size):
                        cand = nb_cells[(start + k) % nb_size]
                        if sender_ptr[cand + 1] > sender_ptr[cand]:
                            current = cand
                            found_sender = True
                            break
                if not found_sender:
                    break
                continue

            # try different LRs for this sender
            tried = np.zeros(s_end - s_beg, dtype=np.bool_)
            advanced = False
            
            while (not advanced) and (not np.all(tried)):
                untried = np.where(~tried)[0]
                pick = untried[np.random.randint(untried.size)]
                lr_idx = sender_idx[s_beg + pick]
                tried[pick] = True

                # receiver slice
                r_beg = recv_ptr[lr_idx]
                r_end = recv_ptr[lr_idx + 1]
                if r_end == r_beg:
                    continue

                # sample a few receivers to find an unassigned one quickly
                local_trials = min(16, r_end - r_beg)
                for _ in range(local_trials):
                    ridx = r_beg + np.random.randint(r_end - r_beg)
                    nxt = recv_idx[ridx]
                    # Check: not current cell AND not already assigned
                    if nxt != current and assignments[nxt] == -1:
                        nb_cells[nb_size] = nxt
                        assignments[nxt] = nb_id
                        nb_size += 1
                        current = nxt
                        advanced = True
                        break
                
                # Break outer LR loop if we've advanced
                if advanced:
                    break

            if not advanced:
                # fallback: pick another sender in this neighborhood
                found_sender = False
                if nb_size > 0:
                    start = np.random.randint(nb_size)
                    for k in range(nb_size):
                        cand = nb_cells[(start + k) % nb_size]
                        if sender_ptr[cand + 1] > sender_ptr[cand]:
                            current = cand
                            found_sender = True
                            break
                if not found_sender:
                    break

        nb_id += 1

    # assign remaining cells as singleton neighborhoods
    for i in range(n_cells):
        if assignments[i] == -1:
            assignments[i] = nb_id
            nb_id += 1

    return assignments, nb_id

def merge_small_neighborhoods(
    membership_df: pd.DataFrame,
    target_size: int = 20,
    max_size: int = 30  # Allow some flexibility in merging
) -> pd.DataFrame:
    """
    Post-process neighborhoods to merge those smaller than target_size.
    
    Args:
        membership_df: DataFrame with 'cell' and 'neighborhood_id' columns
        target_size: Target neighborhood size (default: 20)
        max_size: Maximum allowed size after merging (default: 30)
    
    Returns:
        Updated membership_df with merged neighborhoods
    """
    df = membership_df.copy()
    
    # Calculate neighborhood sizes
    nb_sizes = df.groupby('neighborhood_id').size().to_dict()
    
    # Identify small and large neighborhoods
    small_nbs = [nb_id for nb_id, size in nb_sizes.items() if size < target_size]
    acceptable_nbs = [nb_id for nb_id, size in nb_sizes.items() if size >= target_size]
    
    if not small_nbs:
        print("No small neighborhoods to merge")
        return df
    
    print(f"\nMerging {len(small_nbs)} small neighborhoods (size < {target_size})...")
    
    # Sort small neighborhoods by size (largest first for better packing)
    small_nbs.sort(key=lambda x: nb_sizes[x], reverse=True)
    
    # Track merged neighborhoods
    merge_mapping = {}
    new_nb_id = max(df['neighborhood_id']) + 1
    
    # Strategy 1: Merge small neighborhoods together
    current_merger = []
    current_size = 0
    
    for nb_id in small_nbs:
        size = nb_sizes[nb_id]
        
        if current_size + size <= max_size:
            # Add to current merger group
            current_merger.append(nb_id)
            current_size += size
            
            # Check if we've reached a good size
            if current_size >= target_size:
                # Assign all neighborhoods in current merger to new ID
                for old_id in current_merger:
                    merge_mapping[old_id] = new_nb_id
                new_nb_id += 1
                current_merger = []
                current_size = 0
        else:
            # Current merger is complete, start new one
            if current_merger:
                for old_id in current_merger:
                    merge_mapping[old_id] = new_nb_id
                new_nb_id += 1
            current_merger = [nb_id]
            current_size = size
    
    # Handle remaining neighborhoods in current_merger
    if current_merger:
        if current_size >= target_size:
            # Create new neighborhood
            for old_id in current_merger:
                merge_mapping[old_id] = new_nb_id
            new_nb_id += 1
        else:
            # Try to merge with existing acceptable neighborhoods
            for old_id in current_merger:
                # Find an acceptable neighborhood that won't exceed max_size
                merged = False
                for accept_id in acceptable_nbs:
                    if nb_sizes[accept_id] + nb_sizes[old_id] <= max_size:
                        merge_mapping[old_id] = accept_id
                        nb_sizes[accept_id] += nb_sizes[old_id]  # Update size
                        merged = True
                        break
                
                if not merged:
                    # If can't merge, create singleton or small group
                    merge_mapping[old_id] = new_nb_id
                    new_nb_id += 1
    
    # Apply the merge mapping
    df['neighborhood_id'] = df['neighborhood_id'].map(lambda x: merge_mapping.get(x, x))
    
    # Renumber neighborhoods to be continuous from 0
    unique_nbs = sorted(df['neighborhood_id'].unique())
    renumber_map = {old_id: new_id for new_id, old_id in enumerate(unique_nbs)}
    df['neighborhood_id'] = df['neighborhood_id'].map(renumber_map)
    
    # Print statistics
    new_nb_sizes = df.groupby('neighborhood_id').size()
    print(f"After merging: {len(new_nb_sizes)} neighborhoods")
    print(f"  Size distribution: mean={new_nb_sizes.mean():.1f}, std={new_nb_sizes.std():.1f}")
    print(f"  Min={new_nb_sizes.min()}, Max={new_nb_sizes.max()}")
    print(f"  Target size ({target_size}): {int(np.sum(new_nb_sizes == target_size))} neighborhoods")
    print(f"  Small neighborhoods (< {target_size}): {int(np.sum(new_nb_sizes < target_size))}")
    
    # Show updated size distribution
    size_counts = new_nb_sizes.value_counts().sort_index()
    if len(size_counts) <= 15:
        print("\nUpdated size distribution:")
        for size, count in size_counts.items():
            print(f"  Size {size:2d}: {count:4d} neighborhoods")
    
    return df

# ----------------------------------------------------------------------
# Main wrapper - FIXED to use obs.index
# ----------------------------------------------------------------------
def build_neighborhoods_optimized(
    adata: ad.AnnData,
    lr_df: pd.DataFrame,
    neighborhood_size: int = 20,
    expr_min: int = 0,
    random_seed: int = 0
) -> Tuple[pd.DataFrame, int]:
    """
    Optimized neighborhood construction for large datasets.
    ALWAYS uses adata.obs.index for cell identifiers.
    """
    print("Computing LR expression masks (optimized)...")
    lig_can_send, rec_can_recv, lr_df_filtered = compute_lr_masks_optimized(
        adata, lr_df, expr_min
    )
    if lr_df_filtered.empty or lig_can_send.shape[1] == 0:
        print("WARNING: No valid LR pairs after filtering")
        # Use obs.index directly
        cell_ids = adata.obs.index.astype(str).values
        df = pd.DataFrame({"cell": cell_ids, "neighborhood_id": np.arange(adata.n_obs, dtype=int)})
        return df, adata.n_obs

    # stats
    n_senders = int(np.sum(np.any(lig_can_send, axis=1)))
    n_receivers = int(np.sum(np.any(rec_can_recv, axis=1)))
    print(f"  {n_senders}/{adata.n_obs} cells can send at least one ligand")
    print(f"  {n_receivers}/{adata.n_obs} cells can receive at least one signal")

    # adjacency
    print("Building sender/receiver adjacency (CSR/CSC)...")
    sender_idx, sender_ptr = build_sender_csr(lig_can_send)
    recv_idx, recv_ptr = build_receiver_csc(rec_can_recv)

    # random walk
    print("Building neighborhoods (Numba-accelerated over CSR/CSC)...")
    assignments, n_neighborhoods = build_neighborhoods_fast_csr(
        sender_idx, sender_ptr, recv_idx, recv_ptr,
        neighborhood_size=neighborhood_size,
        seed=random_seed
    )

    # VALIDATION
    print("\nValidating assignments...")
    
    # Check 1: All cells should be assigned
    assert len(assignments) == adata.n_obs, f"Assignment length mismatch"
    assert np.all(assignments >= 0), "Some cells have no assignment"
    assert np.all(assignments < n_neighborhoods), "Invalid neighborhood IDs"
    print(f"  ✓ All {adata.n_obs} cells assigned to neighborhoods")
    
    # Use obs.index for output - this guarantees unique identifiers
    cell_ids = adata.obs.index.astype(str).values
    membership_df = pd.DataFrame({"cell": cell_ids, "neighborhood_id": assignments})
    
    # Check 2: Verify no duplicate cells in output
    duplicated = membership_df[membership_df.duplicated(subset=['cell'], keep=False)]
    if not duplicated.empty:
        print(f"  ✗ WARNING: Found duplicate cell IDs in obs.index!")
        print(f"    Number of duplicated entries: {len(duplicated)}")
        print("    First 10 duplicates:")
        print(duplicated.head(10))
    else:
        print(f"  ✓ No duplicate cell IDs in output")

    # Neighborhood statistics
    nb_sizes = membership_df.groupby("neighborhood_id").size()
    print(f"\nCreated {n_neighborhoods} neighborhoods:")
    print(f"  Size distribution: mean={nb_sizes.mean():.1f}, std={nb_sizes.std():.1f}")
    print(f"  Min={nb_sizes.min()}, Max={nb_sizes.max()}")
    print(f"  Target size ({neighborhood_size}): {int(np.sum(nb_sizes == neighborhood_size))} neighborhoods")
    print(f"  Singleton neighborhoods: {int(np.sum(nb_sizes == 1))}")
    
    # Show size distribution if not too many unique sizes
    size_counts = nb_sizes.value_counts().sort_index()
    if len(size_counts) <= 15:
        print("\nDetailed size distribution:")
        for size, count in size_counts.items():
            print(f"  Size {size:2d}: {count:4d} neighborhoods")

    return membership_df, n_neighborhoods




# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Build cell neighborhoods using LRI-guided random walk",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--data-file", required=True, help="Path to .h5ad file")
    ap.add_argument("--lr-csv", default="data/LRdatabase/CellPhoneDBv5.0.human.csv", help="Path to LR database CSV")
    ap.add_argument("--output-dir", default="results/random_walk", help="Output directory for membership CSV")
    ap.add_argument("--neighborhood-size", type=int, default=20, help="Target cells per neighborhood (default: 20)")
    ap.add_argument("--expr-min", type=int, default=0, help="Expression threshold: count > expr_min => expressed (default: 0)")
    ap.add_argument("--random-seed", type=int, default=0, help="Random seed for reproducibility (default: 0)")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("Loading AnnData...")
    adata = ad.read_h5ad(args.data_file)
    print(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    # Check for duplicate obs.index (which shouldn't happen but let's be safe)
    if adata.obs.index.duplicated().any():
        n_dups = adata.obs.index.duplicated().sum()
        print(f"  ⚠️  WARNING: Found {n_dups} duplicate entries in obs.index!")
        print("     This may cause issues - consider reindexing your data.")

    print("\nLoading LR database...")
    lr_df = pd.read_csv(args.lr_csv)
    print(f"  Loaded {len(lr_df):,} LR pairs")

    print("\n" + "=" * 60)
    membership_df, n_nb = build_neighborhoods_optimized(
        adata, lr_df,
        neighborhood_size=args.neighborhood_size,
        expr_min=args.expr_min,
        random_seed=args.random_seed
    )

    # Add merging step
    membership_df = merge_small_neighborhoods(
        membership_df, 
        target_size=args.neighborhood_size,
        max_size=args.neighborhood_size+5  # or make this a command-line argument
    )

    print("\n" + "=" * 60)
    output_file = os.path.join(args.output_dir, f"{args.neighborhood_size}_cell_neighborhood.csv")
    print(f"Saving to: {output_file}")
    membership_df.to_csv(output_file, index=False)
    
    # Final check on saved file
    saved_df = pd.read_csv(output_file)
    if saved_df['cell'].duplicated().any():
        print(f"  ⚠️  WARNING: Saved file has duplicate cell IDs!")
    else:
        print(f"  ✓ Saved file verified - no duplicates")
    
    print("\n✅ Done!")
    print("=" * 60)

if __name__ == "__main__":
    main()