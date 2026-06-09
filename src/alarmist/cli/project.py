"""
CLI command: alarmist-project

Project single-cell data onto BPTF motifs using cell-level LRI computation.
"""

import argparse
import sys

from alarmist.cli import common, log_config
from alarmist.constants import (
    COLUMN_NAME_PARAMETER,
    COLUMN_NAME_SAMPLE_ID,
    COLUMN_NAME_VALUE,
)


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
      --output-dir results/single_cell

  # With custom patch-LRI directory
  alarmist-project --adata data.h5ad --bptf-dir results/bptf \\
      --patch-lri-dir results/patchify --output-dir results/single_cell

  # Multi-sample mode
  alarmist-project --adata data.h5ad --bptf-dir results/bptf \\
      --output-dir results/single_cell --multi-sample --sample-column batch
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
        "--cellchatdb",
        type=str,
        default=None,
        help="Path to CellChatDB CSV file (uses built-in if not specified)",
    )
    parser.add_argument(
        "--resource",
        "-r",
        type=str,
        default="cellchatdb",
        choices=["cellchatdb", "cellphonedb"],
        help="LR database resource to use (default: cellchatdb). Must match the resource used in patchify.",
    )
    parser.add_argument(
        "--multi-sample",
        action="store_true",
        help="Enable multi-sample mode for concatenated AnnData",
    )
    parser.add_argument(
        "--sample-column",
        type=str,
        default=COLUMN_NAME_SAMPLE_ID,
        help="Column in adata.obs containing sample IDs (default: sample_id)",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=200,
        help="Maximum iterations for projection (default: 200)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=50000,
        help="Number of cells per chunk for memory efficiency (default: 50000)",
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
    from bptf import load_bptf

    import alarmist as al

    # Load data
    logger.info(f"Loading AnnData from {args.adata}")
    adata = sc.read_h5ad(args.adata)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Determine patch-LRI directory
    if args.patch_lri_dir is None:
        bptf_path = Path(args.bptf_dir)
        patch_lri_dir = bptf_path.parent / "patchify"
        if not patch_lri_dir.exists():
            patch_lri_dir = bptf_path.parent / "patch_lri"
        if not patch_lri_dir.exists():
            logger.error(
                f"Could not find patch-LRI directory. Tried:\n"
                f"  {bptf_path.parent / 'patchify'}\n"
                f"  {bptf_path.parent / 'patch_lri'}\n"
                f"Please specify --patch-lri-dir explicitly."
            )
            sys.exit(1)
    else:
        patch_lri_dir = Path(args.patch_lri_dir)

    logger.info(f"Using patch-LRI directory: {patch_lri_dir}")

    # Load patch-LRI results to get column names for alignment
    logger.info("Loading patch-LRI results...")
    patch_results = al.load_patch_lri_results(str(patch_lri_dir))
    column_names = patch_results["column_names"]
    logger.info(f"Loaded {len(column_names)} LRI column names for alignment")

    # Get parameters from patch analysis
    params_file = patch_lri_dir / "analysis_parameters.csv"
    if params_file.exists():
        params_df = pd.read_csv(params_file)
        patch_size_row = params_df[params_df[COLUMN_NAME_PARAMETER] == "patch_size"]
        if len(patch_size_row) > 0:
            neighborhood_size = float(patch_size_row[COLUMN_NAME_VALUE].iloc[0])
        else:
            logger.error("patch_size not found in analysis_parameters.csv")
            sys.exit(1)
    else:
        logger.error(f"Analysis parameters not found: {params_file}")
        sys.exit(1)

    logger.info(f"Using neighborhood size: {neighborhood_size}")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Compute cell-level LRI matrix
    logger.info("Computing cell-level LRI matrix...")
    cell_analyzer = al.NeighborhoodLRIAnalyzer(
        neighborhood_size=neighborhood_size,
        resource_name=args.resource,
        cell_type_column=args.cell_type_column,
        cellchatdb_path=args.cellchatdb,
    )

    if args.multi_sample:
        cell_results = cell_analyzer.run_neighborhood(
            adata,
            output_dir=str(output_dir),
            required_columns=column_names,
            multi_sample=True,
            sample_column=args.sample_column,
        )
    else:
        cell_results = cell_analyzer.run_neighborhood(
            adata,
            output_dir=str(output_dir),
            required_columns=column_names,
        )

    cell_lri_matrix = cell_results["cell_lri_matrix"]
    logger.info(f"Cell-LRI matrix shape: {cell_lri_matrix.shape}")

    # Load BPTF model
    logger.info(f"Loading BPTF model from {args.bptf_dir}")
    bptf_path = Path(args.bptf_dir)
    model_file = bptf_path / "bptf_1.npz"
    if not model_file.exists():
        # Try finding any .npz file
        npz_files = list(bptf_path.glob("*.npz"))
        if npz_files:
            model_file = npz_files[0]
        else:
            logger.error(f"No BPTF model found in {args.bptf_dir}")
            sys.exit(1)

    model = load_bptf(model_file)
    logger.info(f"Loaded BPTF model with {model.n_components} components")

    # Project cell loadings
    logger.info("Projecting cell loadings...")
    np.random.seed(42)

    cell_loadings = al.project_cell_loadings(
        model=model,
        cell_lri_matrix=cell_lri_matrix,
        model_lri_columns=column_names,
        cell_lri_columns=cell_results.get("column_names"),
        max_iter=args.max_iter,
        chunk_size=args.chunk_size,
        verbose=True,
        output_dir=str(output_dir),
    )

    logger.info(f"Cell loadings shape: {cell_loadings.shape}")

    # Save cell loadings as parquet with cell IDs
    cell_loadings_df = pd.DataFrame(
        cell_loadings,
        index=adata.obs_names,
        columns=[f"motif_{i}" for i in range(cell_loadings.shape[1])],
    )
    cell_loadings_df.to_parquet(output_dir / "cell_motif_loadings.parquet")

    # Add to adata and save
    adata.obsm["X_motif"] = cell_loadings
    adata.write_h5ad(output_dir / "projected_adata.h5ad")

    logger.info(f"Results saved to: {output_dir}")
    logger.info("  - cell_loadings.npy")
    logger.info("  - cell_motif_loadings.parquet")
    logger.info("  - projected_adata.h5ad")

    return 0


if __name__ == "__main__":
    sys.exit(main())
