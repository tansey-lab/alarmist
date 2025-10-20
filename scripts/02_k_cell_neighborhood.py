#!/usr/bin/env python3
import os
import argparse
import anndata as ad
import pandas as pd

from sc_neighborhood_lri_analysis import SCNeighborhoodLRIAnalyzer

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data-file", required=True,
                   help="Path to processed .h5ad")
    p.add_argument("--lr-csv", default="data/LRdatabase/CellPhoneDBv5.0.human.csv",
                   help="Path to CellPhoneDBv5.0.human.csv (or CellChatDB) with columns: ligand,receptor,signaling_type")
    p.add_argument("--membership-csv", required=True,
                   help="CSV with columns: cell,neighborhood_id")
    p.add_argument("--output-dir", default="results/random_walk",
                   help="Output directory")
    p.add_argument("--cell-type-column", default="cell_type",
                   help="Column in adata.obs for cell types")
    p.add_argument("--resource", default="cellphonedb",
                   choices=["cellphonedb", "cellchatdb"],
                   help="Which local resource parser to use")
    p.add_argument("--spliter", default="|",
                   help="Separator for LRI column naming")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading AnnData...")
    adata = ad.read_h5ad(args.data_file)
    print(f"  shape={adata.shape}")

    print("Loading neighborhood membership...")
    membership = pd.read_csv(args.membership_csv)

    # Instantiate with your local resource paths
    analyzer = SCNeighborhoodLRIAnalyzer(
        patch_size=50.0,                 # irrelevant here, but required by parent
        resource_name=args.resource,
        spliter=args.spliter,
        cellchatdb_path=args.lr_csv if args.resource == "cellchatdb" else "unused",
        cellphonedb_path=args.lr_csv if args.resource == "cellphonedb" else "unused",
        cell_type_column=args.cell_type_column
    )

    dups = membership.duplicated(subset=['cell'], keep=False)
    print("duplicates in membership:", int(dups.sum()))
    print(membership.loc[dups, ['cell', 'neighborhood_id']].head(10))


    # Feed the precomputed cell->neighborhood map
    analyzer.use_precomputed_neighborhoods(
        adata, membership, cell_col="cell", nb_col="neighborhood_id"
    )

    # Run the LRI counting using your logic (with minimal override for juxtacrine)
    results = analyzer.run_neighborhood_analysis(adata, args.output_dir)

    print("\nDone.")
    print(f"Matrix: {os.path.join(args.output_dir, 'neighborhood_lri_matrix.npz')}")
    print(f"Columns: {os.path.join(args.output_dir, 'neighborhood_lri_columns.csv')}")
    print(f"Membership: {os.path.join(args.output_dir, 'cell_neighborhood_correspondence.csv')}")

if __name__ == "__main__":
    main()
