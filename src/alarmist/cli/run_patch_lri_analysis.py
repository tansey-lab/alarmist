#!/usr/bin/env python3
"""
01 - Patch-based Ligand-Receptor Interaction Analysis

This script runs the patch-based LRI analysis for matrix factorization approaches.
"""

import argparse
import logging
from pathlib import Path

import anndata
import numpy as np
import scipy.sparse as sp

from alarmist import log_config
from alarmist.cli.common import (
    add_data_input_arguments,
    add_lri_arguments,
    add_output_arguments,
    add_random_state_argument,
)
from alarmist.constants import COLUMN_NAME_TMA_ID
from alarmist.core import PatchLRIAnalyzer

logger = logging.getLogger(__name__)


def binomial_thinning(sparse_matrix, probs):
    """Perform binomial thinning on a sparse matrix"""
    sparse_matrix = sparse_matrix.tocoo(copy=False)
    data_thinned = np.random.binomial(n=sparse_matrix.data, p=probs[sparse_matrix.col])
    thinned_matrix = sp.coo_matrix(
        (data_thinned, (sparse_matrix.row, sparse_matrix.col)),
        shape=sparse_matrix.shape,
    )
    return thinned_matrix.tocsr()


def get_parser():
    """Create argument parser"""
    parser = argparse.ArgumentParser(
        description="Patch-based ligand-receptor interaction analysis"
    )

    # Add common arguments
    add_data_input_arguments(parser)
    add_output_arguments(parser)
    add_lri_arguments(parser)
    add_random_state_argument(parser, default=123)

    # Analysis-specific arguments
    parser.add_argument(
        "--use-batch-processing",
        action="store_true",
        help="Use batch processing for large datasets",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=500,
        help="Batch size for processing patches (default: 500)",
    )
    parser.add_argument(
        "--data-thinning",
        type=float,
        default=1.0,
        help="Thin data to which quantile (1.0 = no thinning, default: 1.0)",
    )

    log_config.add_logging_args(parser)
    return parser


def main():
    """Main patch-based LRI analysis pipeline"""
    args = get_parser().parse_args()
    log_config.configure_logging(args)

    logger.info("=" * 60)
    logger.info("01 - PATCH-BASED LIGAND-RECEPTOR INTERACTION ANALYSIS")
    logger.info("=" * 60)
    logger.info(f"Data file: {args.data_file}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Patch size: {args.patch_size} μm")
    logger.info(f"LRI database: {args.resource}")
    logger.info(f"Random state: {args.random_state}")
    logger.info("=" * 60)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info("Loading spatial transcriptomics data...")
    adata = anndata.read_h5ad(args.data_file)
    logger.info(f"Data shape: {adata.shape}")
    logger.info(
        f"Cell types: {adata.obs[args.cell_type_column].cat.categories.tolist()}"
    )

    if COLUMN_NAME_TMA_ID in adata.obs:
        logger.info(f"TMA IDs: {sorted(adata.obs[COLUMN_NAME_TMA_ID].unique())}")

    # Initialize analyzer
    analyzer = PatchLRIAnalyzer(
        patch_size=args.patch_size,
        resource_name=args.resource,
        spliter=args.spliter,
        cell_type_column=args.cell_type_column,
    )

    # Run analysis
    logger.info("Running patch-based LRI analysis...")
    results = analyzer.run_analysis(
        adata,
        str(args.output_dir),
        use_batch_processing=args.use_batch_processing,
        batch_size=args.batch_size,
    )

    # Set random seed
    np.random.seed(args.random_state)

    # Optional data thinning
    if args.data_thinning != 1.0:
        logger.info(f"Applying data thinning at {args.data_thinning} quantile...")
        column_means = np.array(results["patch_lri_matrix"].mean(axis=0)).flatten()
        t = int(args.data_thinning * 100)
        m_t = np.percentile(column_means, t)
        p_j = np.minimum(1, m_t / column_means)

        patch_lri_matrix_thinned = binomial_thinning(results["patch_lri_matrix"], p_j)
        out_path = output_dir / f"patch_lri_matrix_{t}.npz"
        sp.save_npz(out_path, patch_lri_matrix_thinned)
        logger.info(f"Saved thinned matrix to: {out_path}")

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("ANALYSIS SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total patches: {results['patch_info']['n_patches']}")
    logger.info(f"Total LRI combinations: {len(results['column_names'])}")
    logger.info(f"Matrix shape: {results['patch_lri_matrix'].shape}")
    logger.info(
        f"Matrix sparsity: {(1 - results['patch_lri_matrix'].nnz / np.prod(results['patch_lri_matrix'].shape)) * 100:.2f}%"
    )
    logger.info(f"Non-zero interactions: {results['patch_lri_matrix'].nnz:,}")

    # Cell type distribution
    logger.info("\nCell type distribution:")
    cell_type_counts = results["cell_patch_df"][args.cell_type_column].value_counts()
    for cell_type, count in cell_type_counts.items():
        logger.info(f"  {cell_type}: {count:,} cells")

    # TMA distribution
    if COLUMN_NAME_TMA_ID in adata.obs:
        logger.info("\nTMA distribution:")
        tma_counts = results["patch_tma_df"][COLUMN_NAME_TMA_ID].value_counts()
        for tma_id, count in tma_counts.items():
            logger.info(f"  TMA {tma_id}: {count} patches")

    logger.info("\n" + "=" * 60)
    logger.info("STEP 01 COMPLETED SUCCESSFULLY!")
    logger.info("=" * 60)
    logger.info(f"Results saved to: {args.output_dir}")
    logger.info("\nFiles created:")
    logger.info("  - patch_lri_matrix.npz (sparse matrix)")
    if args.data_thinning != 1.0:
        t = int(args.data_thinning * 100)
        logger.info(f"  - patch_lri_matrix_{t}.npz (thinned sparse matrix)")
    logger.info("  - patch_lri_columns.csv (column names)")
    logger.info("  - patch_tma_correspondence.csv (patch-TMA mapping)")
    logger.info("  - cell_patch_correspondence.csv (cell-patch mapping)")
    logger.info("  - analysis_parameters.csv (analysis parameters)")

    logger.info("\nNext step:")
    logger.info("  Run: alarmist-bptf --input-dir " + str(args.output_dir))


if __name__ == "__main__":
    main()
