"""
CLI command: alarmist-project

Project single-cell data onto BPTF motifs.
"""

import argparse
import sys

from alarmist.cli import common, log_config


def get_parser():
    """Create argument parser for project command"""
    parser = argparse.ArgumentParser(
        prog='alarmist-project',
        description='Project single-cell data onto BPTF latent motifs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  alarmist-project --adata data.h5ad --bptf-dir results/bptf --output-dir results/project

  # With custom cell type column
  alarmist-project --adata data.h5ad --bptf-dir results/bptf --output-dir results/project \\
      --cell-type-column annotation
        """
    )

    # Add argument groups
    common.add_adata_arguments(parser)
    common.add_output_arguments(parser)

    parser.add_argument(
        '--bptf-dir', '-b',
        type=str,
        required=True,
        help='Directory containing BPTF results (from alarmist-bptf)'
    )
    parser.add_argument(
        '--patch-lri-dir',
        type=str,
        default=None,
        help='Directory containing patch-LRI results (defaults to bptf-dir/../patch_lri)'
    )
    parser.add_argument(
        '--normalize',
        action='store_true',
        default=True,
        help='Normalize cell loadings (default: True)'
    )
    parser.add_argument(
        '--no-normalize',
        action='store_false',
        dest='normalize',
        help='Disable normalization of cell loadings'
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
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import alarmist as al

    # Load data
    logger.info(f"Loading AnnData from {args.adata}")
    adata = sc.read_h5ad(args.adata)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Check required columns
    if args.cell_type_column not in adata.obs.columns:
        logger.error(f"Cell type column '{args.cell_type_column}' not found in adata.obs")
        logger.error(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)

    # Load BPTF results
    logger.info(f"Loading BPTF results from {args.bptf_dir}")
    bptf_results = al.load_bptf_results(args.bptf_dir)

    # Determine patch-LRI directory
    if args.patch_lri_dir is None:
        bptf_path = Path(args.bptf_dir)
        patch_lri_dir = bptf_path.parent / 'patch_lri'
        if not patch_lri_dir.exists():
            patch_lri_dir = bptf_path.parent / 'patchify'
        logger.info(f"Using patch-LRI directory: {patch_lri_dir}")
    else:
        patch_lri_dir = Path(args.patch_lri_dir)

    # Load patch-LRI results
    patch_lri_results = al.load_patch_lri_results(str(patch_lri_dir))

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Project single cells onto motifs
    logger.info("Projecting cells onto BPTF motifs...")
    cell_loadings = al.project_cell_loadings(
        adata,
        bptf_results,
        patch_lri_results,
        cell_type_column=args.cell_type_column,
        normalize=args.normalize
    )

    # Save results
    logger.info("Saving projection results...")

    # Save cell loadings
    cell_loadings_df = pd.DataFrame(
        cell_loadings,
        index=adata.obs_names,
        columns=[f'motif_{i}' for i in range(cell_loadings.shape[1])]
    )
    cell_loadings_df.to_parquet(output_dir / 'cell_motif_loadings.parquet')

    # Add to adata and save
    adata.obsm['X_motif'] = cell_loadings
    adata.write_h5ad(output_dir / 'projected_adata.h5ad')

    # Report results
    logger.info(f"Cell loadings shape: {cell_loadings.shape}")
    logger.info(f"Results saved to: {args.output_dir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
