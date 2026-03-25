#!/usr/bin/env python3
"""
Poisson GLM differential expression analysis

Based on scripts/04_poisson_glm.py
"""

import argparse
import logging
from pathlib import Path

from alarmist import log_config
from alarmist.core import run_poisson_glm_analysis

logger = logging.getLogger(__name__)


def get_parser():
    """Create argument parser"""
    parser = argparse.ArgumentParser(
        description="BPTF Differential Expression Analysis using Poisson GLM"
    )

    # Input/Output
    parser.add_argument(
        "--bptf-dir",
        "--results-dir",
        dest="bptf_dir",
        required=True,
        help="BPTF results directory",
    )
    parser.add_argument(
        "--patch-lri-dir",
        default=None,
        help="Patch-LRI results directory (if different from BPTF dir parent)",
    )
    parser.add_argument(
        "--data-file", required=True, help="AnnData file path (h5ad format)"
    )
    parser.add_argument(
        "--output-dir", required=True, help="Output directory for DE results"
    )

    # Analysis parameters
    parser.add_argument(
        "--count-layer",
        default="layers:counts",
        help='Layer to use: "X", "raw", or "layers:NAME" (default: layers:counts)',
    )
    parser.add_argument(
        "--splitter", default="|", help="Separator for LRI names (default: |)"
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="FDR significance threshold (default: 0.05)",
    )
    parser.add_argument(
        "--keep-sparse", action="store_true", help="Keep count matrix in sparse format"
    )

    # Random seed
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )

    # Logging
    log_config.add_logging_args(parser)

    return parser


def main():
    """Main function"""
    args = get_parser().parse_args()
    log_config.configure_logging(args)

    logger.info("=" * 60)
    logger.info("POISSON GLM DIFFERENTIAL EXPRESSION ANALYSIS")
    logger.info("=" * 60)
    logger.info(f"BPTF results: {args.bptf_dir}")
    logger.info(f"Data file: {args.data_file}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Random seed: {args.random_state}")
    logger.info("=" * 60)

    # Determine patch-LRI directory
    if args.patch_lri_dir is None:
        # Try to infer from BPTF directory structure
        results_path = Path(args.bptf_dir)
        if results_path.name == "bptf":
            patch_lri_dir = str(results_path.parent / "patch_lri")
        else:
            patch_lri_dir = args.bptf_dir
        logger.info(f"Inferred patch-LRI directory: {patch_lri_dir}")
    else:
        patch_lri_dir = args.patch_lri_dir

    # Run analysis
    try:
        results = run_poisson_glm_analysis(
            bptf_dir=args.bptf_dir,
            patch_lri_dir=patch_lri_dir,
            data_file=args.data_file,
            output_dir=args.output_dir,
            count_layer=args.count_layer,
            splitter=args.splitter,
            alpha=args.alpha,
            random_state=args.random_state,
            keep_sparse=args.keep_sparse,
        )

        logger.info("=" * 60)
        logger.info("ANALYSIS COMPLETED SUCCESSFULLY!")
        logger.info("=" * 60)
        logger.info(f"Results saved to: {args.output_dir}")
        logger.info(f"Total motif-celltype combinations analyzed: {len(results)}")

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()
