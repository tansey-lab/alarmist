"""
CLI command: alarmist-project

Project single-cell data onto BPTF motifs by assigning cells to patches.
"""

import argparse
import sys

from alarmist.cli import common, log_config


def get_parser():
    """Create argument parser for project command"""
    parser = argparse.ArgumentParser(
        prog="alarmist-project",
        description="Project single-cell data onto BPTF latent motifs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  alarmist-project --adata data.h5ad --bptf-dir results/bptf \\
      --output-dir results/project

  # With custom patch-LRI directory
  alarmist-project --adata data.h5ad --bptf-dir results/bptf \\
      --output-dir results/project \\
      --patch-lri-dir results/patchify
        """,
    )

    # Add argument groups
    common.add_adata_arguments(parser)
    common.add_output_arguments(parser)

    parser.add_argument(
        "--bptf-dir",
        "-b",
        type=str,
        required=True,
        help="Directory containing BPTF results (from alarmist-bptf)",
    )
    parser.add_argument(
        "--patch-lri-dir",
        type=str,
        default=None,
        help="Directory with patch-LRI results (defaults to bptf-dir/../patchify)",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        default=True,
        help="Normalize cell loadings (default: True)",
    )
    parser.add_argument(
        "--no-normalize",
        action="store_false",
        dest="normalize",
        help="Disable normalization of cell loadings",
    )

    log_config.add_logging_args(parser)

    return parser


def main():
    """Main entry point for project command"""
    parser = get_parser()
    args = parser.parse_args()

    # Configure logging
    logger = log_config.configure_logging(args)

    # Import heavy dependencies only after argument parsing
    logger.info("Loading dependencies...")
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import scanpy as sc

    import alarmist as al

    # Load data
    logger.info(f"Loading AnnData from {args.adata}")
    adata = sc.read_h5ad(args.adata)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Load BPTF results
    logger.info(f"Loading BPTF results from {args.bptf_dir}")
    bptf_results = al.load_bptf_results(args.bptf_dir)
    patch_loadings = bptf_results["patch_loadings"]  # (n_patches, n_motifs)
    n_motifs = patch_loadings.shape[1]
    logger.info(f"Loaded {patch_loadings.shape[0]} patches, {n_motifs} motifs")

    # Determine patch-LRI directory
    if args.patch_lri_dir is None:
        bptf_path = Path(args.bptf_dir)
        patch_lri_dir = bptf_path.parent / "patchify"
        if not patch_lri_dir.exists():
            patch_lri_dir = bptf_path.parent / "patch_lri"
        logger.info(f"Using patch-LRI directory: {patch_lri_dir}")
    else:
        patch_lri_dir = Path(args.patch_lri_dir)

    # Load patch metadata
    patch_metadata_file = patch_lri_dir / "patch_metadata.parquet"
    if not patch_metadata_file.exists():
        logger.error(f"Patch metadata not found: {patch_metadata_file}")
        sys.exit(1)

    patch_metadata = pd.read_parquet(patch_metadata_file)
    logger.info(f"Loaded metadata for {len(patch_metadata)} patches")

    # Get spatial coordinates
    if "spatial" not in adata.obsm:
        logger.error("Spatial coordinates not found in adata.obsm['spatial']")
        sys.exit(1)

    cell_coords = adata.obsm["spatial"]

    # Get patch size from metadata or parameters
    params_file = patch_lri_dir / "analysis_parameters.csv"
    if params_file.exists():
        params_df = pd.read_csv(params_file)
        patch_size_row = params_df[params_df["parameter"] == "patch_size"]
        if len(patch_size_row) > 0:
            patch_size = float(patch_size_row["value"].iloc[0])
        else:
            # Infer from patch metadata
            patch_size = (
                patch_metadata["x_max"].iloc[0] - patch_metadata["x_min"].iloc[0]
            )
    else:
        patch_size = patch_metadata["x_max"].iloc[0] - patch_metadata["x_min"].iloc[0]

    logger.info(f"Patch size: {patch_size}")

    # Assign cells to patches based on spatial coordinates
    logger.info("Assigning cells to patches...")

    # Create patch boundaries array for vectorized lookup
    patch_x_min = patch_metadata["x_min"].values
    patch_x_max = patch_metadata["x_max"].values
    patch_y_min = patch_metadata["y_min"].values
    patch_y_max = patch_metadata["y_max"].values

    # For each cell, find which patch it belongs to
    cell_x = cell_coords[:, 0]
    cell_y = cell_coords[:, 1]

    # Initialize cell loadings
    cell_loadings = np.zeros((len(adata), n_motifs))
    cells_assigned = 0

    # Assign each cell to its patch
    for i in range(len(patch_metadata)):
        # Find cells in this patch
        in_patch = (
            (cell_x >= patch_x_min[i])
            & (cell_x < patch_x_max[i])
            & (cell_y >= patch_y_min[i])
            & (cell_y < patch_y_max[i])
        )
        n_cells_in_patch = in_patch.sum()
        if n_cells_in_patch > 0:
            cell_loadings[in_patch] = patch_loadings[i]
            cells_assigned += n_cells_in_patch

    logger.info(f"Assigned {cells_assigned}/{len(adata)} cells to patches")

    # Handle cells not in any patch (use nearest patch)
    unassigned = cell_loadings.sum(axis=1) == 0
    n_unassigned = unassigned.sum()
    if n_unassigned > 0:
        logger.info(f"Assigning {n_unassigned} cells to nearest patch...")
        from scipy.spatial import cKDTree

        # Build KDTree from patch centers
        patch_centers = np.column_stack(
            [(patch_x_min + patch_x_max) / 2, (patch_y_min + patch_y_max) / 2]
        )
        tree = cKDTree(patch_centers)

        # Find nearest patch for unassigned cells
        unassigned_coords = cell_coords[unassigned]
        _, nearest_patches = tree.query(unassigned_coords)
        cell_loadings[unassigned] = patch_loadings[nearest_patches]

    # Normalize if requested
    if args.normalize:
        logger.info("Normalizing cell loadings...")
        row_sums = cell_loadings.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        cell_loadings = cell_loadings / row_sums

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save results
    logger.info("Saving projection results...")

    # Save cell loadings
    cell_loadings_df = pd.DataFrame(
        cell_loadings,
        index=adata.obs_names,
        columns=[f"motif_{i}" for i in range(n_motifs)],
    )
    cell_loadings_df.to_parquet(output_dir / "cell_motif_loadings.parquet")

    # Add to adata and save
    adata.obsm["X_motif"] = cell_loadings
    adata.write_h5ad(output_dir / "projected_adata.h5ad")

    # Report results
    logger.info(f"Cell loadings shape: {cell_loadings.shape}")
    logger.info(f"Results saved to: {args.output_dir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
