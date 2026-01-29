"""
CLI command: alarmist-patchify

Patchify tissue and count LRI interactions.
"""

import argparse
import sys

from alarmist.cli import common, log_config


def get_parser():
    """Create argument parser for patchify command"""
    parser = argparse.ArgumentParser(
        prog='alarmist-patchify',
        description='Patchify tissue and count ligand-receptor interactions (LRI)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  alarmist-patchify --adata data.h5ad --output-dir results/patch_lri

  # With custom patch size
  alarmist-patchify --adata data.h5ad --output-dir results/ --patch-size 80

  # Multi-sample mode
  alarmist-patchify --adata data.h5ad --output-dir results/ --multi-sample --sample-column batch
        """
    )

    # Add argument groups
    common.add_adata_arguments(parser)
    common.add_output_arguments(parser)
    common.add_lri_arguments(parser)
    common.add_multi_sample_arguments(parser)
    log_config.add_logging_args(parser)

    return parser


def main():
    """Main entry point for patchify command"""
    parser = get_parser()
    args = parser.parse_args()

    # Configure logging
    logger = log_config.configure_logging(args)

    # Validate arguments
    if args.multi_sample and args.sample_column is None:
        parser.error("--sample-column is required when --multi-sample is enabled")

    # Import heavy dependencies only after argument parsing
    logger.info("Loading dependencies...")
    import scanpy as sc
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

    if 'spatial' not in adata.obsm:
        logger.error("Spatial coordinates not found in adata.obsm['spatial']")
        sys.exit(1)

    # Create analyzer
    logger.info(f"Creating PatchLRIAnalyzer with patch_size={args.patch_size}")
    analyzer = al.PatchLRIAnalyzer(
        patch_size=args.patch_size,
        resource_name=args.resource,
        cell_type_column=args.cell_type_column,
        cellchatdb_path=args.cellchatdb_path
    )

    # Run patchify
    logger.info("Running patchify...")
    results = analyzer.run_patchify(
        adata,
        output_dir=args.output_dir,
        multi_sample=args.multi_sample,
        sample_column=args.sample_column
    )

    # Report results
    matrix = results['patch_lri_matrix']
    logger.info(f"Patch-LRI matrix shape: {matrix.shape}")
    logger.info(f"Number of patches: {matrix.shape[0]}")
    logger.info(f"Number of LRI combinations: {matrix.shape[1]}")
    logger.info(f"Results saved to: {args.output_dir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
