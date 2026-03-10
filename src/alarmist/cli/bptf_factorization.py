#!/usr/bin/env python3
"""
02 - BPTF Matrix Factorization for Spatial LRI Analysis
"""

import argparse
import logging
import time
from pathlib import Path

import numpy as np

from alarmist import log_config
from alarmist.cli.common import (
    add_bptf_arguments,
    add_output_arguments,
    add_random_state_argument,
)
from alarmist.core import (
    BPTF_AVAILABLE,
    extract_factors,
    get_top_motifs,
    run_bptf,
    save_bptf_results,
)
from alarmist.data.loaders import load_patch_lri_results
from alarmist.plotting.bptf_plots import plot_bptf_diagnostics

logger = logging.getLogger(__name__)


def get_parser():
    """Create argument parser"""
    parser = argparse.ArgumentParser(
        description="BPTF matrix factorization for spatial LRI analysis"
    )

    parser.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Input directory with patch LRI results",
    )
    parser.add_argument(
        "--sparse-matrix-name",
        type=str,
        default="patch_lri_matrix.npz",
        help="Name of the sparse matrix file (default: patch_lri_matrix.npz)",
    )
    parser.add_argument(
        "--column-df-name",
        type=str,
        default="patch_lri_columns.csv",
        help="Name of column names file (default: patch_lri_columns.csv)",
    )
    parser.add_argument(
        "--meta-df-name",
        type=str,
        default="cell_patch_correspondence.csv",
        help="Name of metadata file (default: cell_patch_correspondence.csv)",
    )

    add_output_arguments(parser)
    add_bptf_arguments(parser)
    add_random_state_argument(parser, default=0)

    parser.add_argument(
        "--neighborhood",
        type=bool,
        default=False,
        help="Whether this is neighborhood-based analysis (default: False)",
    )
    parser.add_argument(
        "--single-cell",
        type=bool,
        default=False,
        help="Whether this is single-cell analysis (default: False)",
    )

    log_config.add_logging_args(parser)
    return parser


def main():
    """Main BPTF matrix factorization pipeline"""
    args = get_parser().parse_args()
    log_config.configure_logging(args)

    logger.info("=" * 60)
    logger.info("02 - BPTF MATRIX FACTORIZATION")
    logger.info("=" * 60)
    logger.info(f"Input directory: {args.input_dir}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Number of components: {args.n_components}")
    logger.info(f"Max iterations: {args.max_iter}")
    logger.info(f"Random state: {args.random_state}")
    logger.info("=" * 60)

    # Check BPTF availability
    if not BPTF_AVAILABLE:
        logger.error("BPTF not available!")
        logger.error("Install: pip install git+https://github.com/aschein/bptf.git")
        return

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "plots").mkdir(exist_ok=True)

    # Load patch LRI results
    logger.info("Loading patch-LRI results...")
    results = load_patch_lri_results(
        args.input_dir,
        sparse_matrix_name=args.sparse_matrix_name,
        column_df_name=args.column_df_name,
        meta_df_name=args.meta_df_name,
        neighborhood=args.neighborhood,
        single_cell=args.single_cell,
    )

    mat = results["patch_lri_matrix"]
    cols = results["column_names"]
    cell_meta_df = results.get("cell_patch_df", results.get("patch_tma_df"))

    logger.info(f"Matrix shape: {mat.shape}")
    logger.info(f"Matrix sparsity: {(1 - mat.nnz / np.prod(mat.shape)) * 100:.2f}%")
    logger.info(f"Number of LRIs: {len(cols)}")

    # Run BPTF
    logger.info(f"\nRunning BPTF with {args.n_components} components...")
    start_time = time.time()

    model = run_bptf(
        mat,
        n_components=args.n_components,
        max_iter=args.max_iter,
        verbose=True,
        random_state=args.random_state,
    )

    fit_time = time.time() - start_time
    logger.info(f"BPTF completed in {fit_time:.1f} seconds!")

    # Extract factors
    logger.info("Extracting factor matrices...")
    patch_loadings, lri_factors = extract_factors(model)
    logger.info(f"Patch loadings shape: {patch_loadings.shape}")
    logger.info(f"LRI factors shape: {lri_factors.shape}")

    # Save results
    logger.info("Saving BPTF results...")
    patch_metadata_df = results.get("patch_tma_df", cell_meta_df)
    save_bptf_results(
        model,
        patch_loadings,
        lri_factors,
        cols,
        patch_metadata_df,
        str(args.output_dir),
    )

    # Create diagnostic plots
    if args.make_plots:
        logger.info("Creating diagnostic plots...")
        plot_bptf_diagnostics(patch_loadings, lri_factors, str(args.output_dir))

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("BPTF SUMMARY")
    logger.info("=" * 60)

    top_motifs = get_top_motifs(patch_loadings, top_k=5)
    logger.info(f"Total motifs: {args.n_components}")
    logger.info("\nTop 5 motifs by activity:")
    for i, (motif_idx, activity, fraction) in enumerate(
        zip(
            top_motifs["motif_indices"][:5],
            top_motifs["activities"][:5],
            top_motifs["activity_fractions"][:5],
        )
    ):
        logger.info(f"  {i+1}. Motif {motif_idx}: {fraction:.1%} of total activity")

    logger.info("\n" + "=" * 60)
    logger.info("STEP 02 COMPLETED SUCCESSFULLY!")
    logger.info("=" * 60)
    logger.info(f"Results saved to: {args.output_dir}")
    logger.info("\nFiles created:")
    logger.info("  - bptf_model.npz (BPTF model)")
    logger.info("  - patch_loadings.npy, lri_factors.npy (factor matrices)")
    logger.info("  - patch_motifs.csv, lri_motifs.csv (detailed analysis)")
    logger.info("  - top_motifs.csv (motif summaries)")
    logger.info("  - factorization_parameters.csv (parameters)")
    if args.make_plots:
        logger.info("  - plots/ (diagnostic visualizations)")

    logger.info("\nNext step:")
    logger.info("  Run: alarmist-bptf-viz --bptf-dir " + str(args.output_dir))


if __name__ == "__main__":
    main()
