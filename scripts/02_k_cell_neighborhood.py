#!/usr/bin/env python3
"""
Run neighborhood LRI analysis using precomputed neighborhoods.
Updated to use adata.obs.index for cell identification.
"""

import os
import argparse
import anndata as ad
import pandas as pd
import numpy as np
from sc_neighborhood_lri_analysis import SCNeighborhoodLRIAnalyzer

def validate_membership_alignment(adata: ad.AnnData, membership: pd.DataFrame):
    """
    Validate that membership cell IDs align with adata.obs.index
    """
    print("\nValidating cell ID alignment...")
    
    # Get cell IDs from both sources
    adata_cells = set(adata.obs.index.astype(str))
    membership_cells = set(membership['cell'].astype(str))
    
    # Check for mismatches
    in_adata_not_membership = adata_cells - membership_cells
    in_membership_not_adata = membership_cells - adata_cells
    
    if in_adata_not_membership:
        print(f"  WARNING: {len(in_adata_not_membership)} cells in adata but not in membership")
        print(f"    First 5: {list(in_adata_not_membership)[:5]}")
    
    if in_membership_not_adata:
        print(f"  WARNING: {len(in_membership_not_adata)} cells in membership but not in adata")
        print(f"    First 5: {list(in_membership_not_adata)[:5]}")
    
    if not in_adata_not_membership and not in_membership_not_adata:
        print(f"  ✓ Perfect alignment: {len(adata_cells)} cells match")
    else:
        print(f"  ⚠ Partial alignment: {len(adata_cells & membership_cells)} cells match")
    
    return len(adata_cells & membership_cells) > 0

def ensure_index_alignment(adata: ad.AnnData, membership: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure membership DataFrame is properly aligned with adata.obs.index
    Returns aligned membership DataFrame.
    """
    # Convert both to string for comparison
    adata_index = adata.obs.index.astype(str).values
    membership = membership.copy()
    membership['cell'] = membership['cell'].astype(str)
    
    # Create a mapping from cell ID to neighborhood
    cell_to_nb = dict(zip(membership['cell'], membership['neighborhood_id']))
    
    # Create aligned membership DataFrame in adata order
    aligned_membership = pd.DataFrame({
        'cell': adata_index,
        'neighborhood_id': [cell_to_nb.get(cell, -1) for cell in adata_index]
    })
    
    # Check for unassigned cells
    unassigned = (aligned_membership['neighborhood_id'] == -1).sum()
    if unassigned > 0:
        print(f"  WARNING: {unassigned} cells in adata have no neighborhood assignment")
        print("  These cells will be excluded from LRI analysis")
        
        # Option 1: Remove unassigned cells
        aligned_membership = aligned_membership[aligned_membership['neighborhood_id'] != -1]
        
        # Option 2: Assign them to singleton neighborhoods (uncomment if preferred)
        # max_nb = aligned_membership['neighborhood_id'].max()
        # for i, row in aligned_membership.iterrows():
        #     if row['neighborhood_id'] == -1:
        #         max_nb += 1
        #         aligned_membership.loc[i, 'neighborhood_id'] = max_nb
    
    return aligned_membership

def main():
    p = argparse.ArgumentParser(
        description="Run neighborhood LRI analysis with precomputed neighborhoods"
    )
    p.add_argument("--data-file", required=True,
                   help="Path to processed .h5ad file")
    p.add_argument("--lr-csv", 
                   default="data/LRdatabase/CellPhoneDBv5.0.human.csv",
                   help="Path to LR database CSV")
    p.add_argument("--membership-csv", required=True,
                   help="CSV with columns: cell,neighborhood_id (cell should match adata.obs.index)")
    p.add_argument("--output-dir", default="results/random_walk",
                   help="Output directory")
    p.add_argument("--cell-type-column", default="cell_type",
                   help="Column in adata.obs for cell types")
    p.add_argument("--resource", default="cellphonedb",
                   choices=["cellphonedb", "cellchatdb"],
                   help="Which local resource parser to use")
    p.add_argument("--spliter", default="|",
                   help="Separator for LRI column naming")
    p.add_argument("--skip-validation", action="store_true",
                   help="Skip cell ID validation (not recommended)")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load data
    print("=" * 60)
    print("Loading AnnData...")
    adata = ad.read_h5ad(args.data_file)
    print(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"  First 5 cell IDs (obs.index): {list(adata.obs.index[:5])}")
    
    # Check for duplicate indices in adata
    if adata.obs.index.duplicated().any():
        n_dups = adata.obs.index.duplicated().sum()
        print(f"  ⚠ WARNING: adata has {n_dups} duplicate obs.index entries!")
        print("  This may cause issues in the analysis")

    print("\nLoading neighborhood membership...")
    membership = pd.read_csv(args.membership_csv)
    print(f"  Loaded {len(membership):,} cell-neighborhood assignments")
    print(f"  First 5 cell IDs in membership: {list(membership['cell'][:5])}")
    
    # Check for duplicates in membership
    dups = membership.duplicated(subset=['cell'], keep=False)
    if dups.sum() > 0:
        print(f"\n  ⚠ WARNING: {int(dups.sum())} duplicate cell entries in membership!")
        print("  Duplicate examples:")
        print(membership.loc[dups, ['cell', 'neighborhood_id']].head(10))
        
        # Remove duplicates (keep first occurrence)
        print("  Removing duplicates (keeping first occurrence)...")
        membership = membership.drop_duplicates(subset=['cell'], keep='first')
        print(f"  After deduplication: {len(membership):,} assignments")
    
    # Validate and align
    if not args.skip_validation:
        if not validate_membership_alignment(adata, membership):
            print("\nERROR: No matching cells between adata and membership!")
            print("Check that membership 'cell' column matches adata.obs.index")
            return
        
        print("\nAligning membership with adata.obs.index order...")
        membership = ensure_index_alignment(adata, membership)
        print(f"  Aligned membership has {len(membership):,} assignments")
    
    # Neighborhood statistics
    print("\nNeighborhood statistics:")
    nb_sizes = membership.groupby('neighborhood_id').size()
    print(f"  Number of neighborhoods: {len(nb_sizes):,}")
    print(f"  Size distribution: mean={nb_sizes.mean():.1f}, "
          f"min={nb_sizes.min()}, max={nb_sizes.max()}")
    print(f"  Singleton neighborhoods: {(nb_sizes == 1).sum():,}")
    
    print("\n" + "=" * 60)
    print("Initializing LRI analyzer...")
    
    # Initialize analyzer
    analyzer = SCNeighborhoodLRIAnalyzer(
        patch_size=50.0,  # Not used for precomputed neighborhoods
        resource_name=args.resource,
        spliter=args.spliter,
        cellchatdb_path=args.lr_csv if args.resource == "cellchatdb" else "unused",
        cellphonedb_path=args.lr_csv if args.resource == "cellphonedb" else "unused",
        cell_type_column=args.cell_type_column
    )
    
    # Feed the precomputed neighborhoods
    print("Setting precomputed neighborhoods...")
    analyzer.use_precomputed_neighborhoods(
        adata, 
        membership, 
        cell_col="cell",  # This column should now match adata.obs.index
        nb_col="neighborhood_id"
    )
    
    # Run LRI analysis
    print("\nRunning neighborhood LRI analysis...")
    results = analyzer.run_neighborhood_analysis(adata, args.output_dir)
    
    # Save additional metadata
    metadata_file = os.path.join(args.output_dir, 'analysis_metadata.txt')
    with open(metadata_file, 'w') as f:
        f.write(f"Data file: {args.data_file}\n")
        f.write(f"Membership file: {args.membership_csv}\n")
        f.write(f"LR database: {args.lr_csv}\n")
        f.write(f"Resource type: {args.resource}\n")
        f.write(f"Number of cells: {adata.n_obs}\n")
        f.write(f"Number of neighborhoods: {len(nb_sizes)}\n")
        f.write(f"Mean neighborhood size: {nb_sizes.mean():.2f}\n")
        f.write(f"Cell type column: {args.cell_type_column}\n")
    
    print("\n" + "=" * 60)
    print("✅ Analysis complete! Output files:")
    print(f"  Matrix: {os.path.join(args.output_dir, 'neighborhood_lri_matrix.npz')}")
    print(f"  Columns: {os.path.join(args.output_dir, 'neighborhood_lri_columns.csv')}")
    print(f"  Membership: {os.path.join(args.output_dir, 'cell_neighborhood_correspondence.csv')}")
    print(f"  Metadata: {metadata_file}")

if __name__ == "__main__":
    main()