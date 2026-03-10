"""
CLI command: alarmist-bptf

Run BPTF matrix factorization on patch-LRI results.
"""

import argparse
import sys

from alarmist.cli import common, log_config


def get_parser():
    """Create argument parser for bptf command"""
    parser = argparse.ArgumentParser(
        prog="alarmist-bptf",
        description="Run BPTF (Bayesian Poisson Tensor Factorization) on patch-LRI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  alarmist-bptf --input-dir results/patch_lri --output-dir results/bptf

  # With custom parameters
  alarmist-bptf --input-dir results/patch_lri --output-dir results/bptf \\
      --n-components 20 --max-iter 500 --seed 42
        """,
    )

    # Input/output arguments
    parser.add_argument(
        "--input-dir",
        "-i",
        type=str,
        required=True,
        help="Directory containing patch-LRI results (from alarmist-patchify)",
    )
    common.add_output_arguments(parser)

    # BPTF arguments
    common.add_bptf_arguments(parser)
    common.add_seed_argument(parser)

    # Logging
    log_config.add_logging_args(parser)

    return parser


def main():
    """Main entry point for bptf command"""
    parser = get_parser()
    args = parser.parse_args()

    # Configure logging
    logger = log_config.configure_logging(args)

    # Import heavy dependencies only after argument parsing
    logger.info("Loading dependencies...")
    import numpy as np

    import alarmist as al

    # Set random seed
    np.random.seed(args.seed)
    logger.info(f"Random seed: {args.seed}")

    # Load patch-LRI results
    logger.info(f"Loading patch-LRI results from {args.input_dir}")
    results = al.load_patch_lri_results(args.input_dir)

    matrix = results["patch_lri_matrix"]
    logger.info(f"Loaded matrix shape: {matrix.shape}")

    # Check BPTF availability
    if not al.BPTF_AVAILABLE:
        logger.error(
            "BPTF is not available. Install from: https://github.com/aschein/bptf"
        )
        sys.exit(1)

    # Run BPTF
    logger.info(
        f"Running BPTF with n_components={args.n_components}, max_iter={args.max_iter}"
    )
    model = al.run_bptf(
        matrix,
        n_components=args.n_components,
        max_iter=args.max_iter,
        verbose=args.verbose if hasattr(args, "verbose") else False,
        random_state=args.seed,
    )

    # Extract and save results
    logger.info("Processing BPTF results...")
    bptf_results = al.process_bptf_results(model, results, output_dir=args.output_dir)

    # Report results
    patch_loadings = bptf_results["patch_loadings"]
    lri_factors = bptf_results["lri_factors"]
    logger.info(f"Patch loadings shape: {patch_loadings.shape}")
    logger.info(f"LRI factors shape: {lri_factors.shape}")
    logger.info(f"Results saved to: {args.output_dir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
