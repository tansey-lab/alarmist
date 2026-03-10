"""
CLI command: alarmist-glm

Run GLM differential expression analysis on BPTF results.
"""

import argparse
import sys

from alarmist.cli import common, log_config


def get_parser():
    """Create argument parser for glm command"""
    parser = argparse.ArgumentParser(
        prog='alarmist-glm',
        description='Run Poisson GLM differential expression analysis on BPTF motifs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  alarmist-glm --input-dir results/project --adata data.h5ad --output-dir results/glm

  # With condition column for differential analysis
  alarmist-glm --input-dir results/project --adata data.h5ad --output-dir results/glm \\
      --condition-column treatment
        """
    )

    # Input/output arguments
    parser.add_argument(
        '--input-dir', '-i',
        type=str,
        required=True,
        help='Directory containing projection results (from alarmist-project) or BPTF results'
    )
    parser.add_argument(
        '--adata',
        type=str,
        required=True,
        help='Path to AnnData h5ad file'
    )
    common.add_output_arguments(parser)

    # Analysis parameters
    parser.add_argument(
        '--bptf-dir',
        type=str,
        default=None,
        help='BPTF results directory (defaults to input-dir/../bptf)'
    )
    parser.add_argument(
        '--patch-lri-dir',
        type=str,
        default=None,
        help='Patch-LRI results directory (defaults to input-dir/../patch_lri)'
    )
    parser.add_argument(
        '--condition-column',
        type=str,
        default=None,
        help='Column in adata.obs for condition/group comparison'
    )
    parser.add_argument(
        '--covariates',
        type=str,
        nargs='+',
        default=None,
        help='Additional covariate columns from adata.obs'
    )
    parser.add_argument(
        '--count-layer',
        type=str,
        default='layers:counts',
        help='Which layer to use: "X", "raw", or "layers:NAME" (default: layers:counts)'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=0.05,
        help='FDR significance threshold (default: 0.05)'
    )

    common.add_seed_argument(parser)
    log_config.add_logging_args(parser)

    return parser


def main():
    """Main entry point for glm command"""
    parser = get_parser()
    args = parser.parse_args()

    # Configure logging
    logger = log_config.configure_logging(args)

    # Import heavy dependencies only after argument parsing
    logger.info("Loading dependencies...")
    from pathlib import Path
    import numpy as np
    import alarmist as al

    # Set random seed
    np.random.seed(args.seed)
    logger.info(f"Random seed: {args.seed}")

    # Determine directory structure
    input_path = Path(args.input_dir)

    # Try to find BPTF directory
    if args.bptf_dir:
        bptf_dir = args.bptf_dir
    else:
        # Try common locations
        candidates = [
            input_path / 'bptf',
            input_path.parent / 'bptf',
            input_path,
        ]
        bptf_dir = None
        for c in candidates:
            if (c / 'lri_factors.parquet').exists() or (c / 'model.pkl').exists():
                bptf_dir = str(c)
                break
        if bptf_dir is None:
            bptf_dir = str(input_path)
    logger.info(f"Using BPTF directory: {bptf_dir}")

    # Try to find patch-LRI directory
    if args.patch_lri_dir:
        patch_lri_dir = args.patch_lri_dir
    else:
        candidates = [
            input_path / 'patch_lri',
            input_path.parent / 'patch_lri',
            input_path / 'patchify',
            input_path.parent / 'patchify',
        ]
        patch_lri_dir = None
        for c in candidates:
            if c.exists():
                patch_lri_dir = str(c)
                break
        if patch_lri_dir is None:
            patch_lri_dir = str(input_path)
    logger.info(f"Using patch-LRI directory: {patch_lri_dir}")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run GLM analysis
    logger.info("Running Poisson GLM analysis...")
    try:
        results = al.run_poisson_glm_analysis(
            bptf_dir=bptf_dir,
            patch_lri_dir=patch_lri_dir,
            data_file=args.adata,
            output_dir=str(output_dir),
            count_layer=args.count_layer,
            alpha=args.alpha,
            random_state=args.seed,
        )

        logger.info(f"GLM analysis complete")
        logger.info(f"Total motif-celltype combinations analyzed: {len(results)}")
        logger.info(f"Results saved to: {args.output_dir}")

    except Exception as e:
        logger.error(f"GLM analysis failed: {e}")
        raise

    return 0


if __name__ == "__main__":
    sys.exit(main())
