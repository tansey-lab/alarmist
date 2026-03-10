#!/usr/bin/env python3
"""
GLM results analysis and visualization

Based on scripts/05_glm_results.py
"""

import argparse
import logging

from alarmist import log_config
from alarmist.core import analyze_glm_results

logger = logging.getLogger(__name__)


def get_parser():
    """Create argument parser"""
    parser = argparse.ArgumentParser(
        description="GLM results analysis and visualization"
    )

    # Input/Output
    parser.add_argument(
        "--data-file", required=True, help="Input h5ad file with expression data"
    )
    parser.add_argument(
        "--results-dir", required=True, help="Directory containing GLM results CSVs"
    )
    parser.add_argument(
        "--output-dir", required=True, help="Output directory for plots and analysis"
    )

    # Analysis parameters
    parser.add_argument(
        "--n-motifs",
        type=int,
        default=20,
        help="Number of motifs to analyze (default: 20)",
    )
    parser.add_argument(
        "--min-expression-frac",
        type=float,
        default=0.02,
        help="Minimum expression fraction in cell type (default: 0.02)",
    )

    # Marker gene parameters
    parser.add_argument(
        "--marker-lfc",
        type=float,
        default=1.0,
        help="Log fold change threshold for marker genes (default: 1.0)",
    )
    parser.add_argument(
        "--marker-pvalue",
        type=float,
        default=1e-5,
        help="P-value threshold for marker genes (default: 1e-5)",
    )
    parser.add_argument(
        "--marker-subsample",
        type=int,
        default=50000,
        help="Maximum cells to use for marker gene detection (default: 50000)",
    )
    parser.add_argument(
        "--force-recompute-markers",
        action="store_true",
        help="Force recompute marker genes even if cached",
    )

    # Plotting parameters
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="FDR threshold for volcano plots (default: 0.05)",
    )
    parser.add_argument(
        "--lfc-threshold",
        type=float,
        default=0.5,
        help="Log fold change threshold for volcano plots (default: 0.5)",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=10,
        help="Number of top genes to label in volcano plots (default: 10)",
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
    logger.info("GLM RESULTS ANALYSIS AND VISUALIZATION")
    logger.info("=" * 60)
    logger.info(f"Data file: {args.data_file}")
    logger.info(f"Results directory: {args.results_dir}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Number of motifs: {args.n_motifs}")
    logger.info(f"Random seed: {args.random_state}")
    logger.info("=" * 60)

    # Run analysis
    try:
        cell_types, all_genes, exclusion_mask = analyze_glm_results(
            data_file=args.data_file,
            results_dir=args.results_dir,
            output_dir=args.output_dir,
            n_motifs=args.n_motifs,
            min_expression_frac=args.min_expression_frac,
            marker_lfc=args.marker_lfc,
            marker_pvalue=args.marker_pvalue,
            marker_subsample=args.marker_subsample,
            fdr_threshold=args.fdr_threshold,
            lfc_threshold=args.lfc_threshold,
            n_top_genes=args.n_top_genes,
            random_state=args.random_state,
            force_recompute_markers=args.force_recompute_markers,
        )

        logger.info("=" * 60)
        logger.info("ANALYSIS COMPLETED SUCCESSFULLY!")
        logger.info("=" * 60)
        logger.info(f"Cell types analyzed: {len(cell_types)}")
        logger.info(f"Total genes: {len(all_genes)}")
        logger.info(f"Results saved to: {args.output_dir}")

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()
