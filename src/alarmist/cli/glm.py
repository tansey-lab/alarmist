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
        prog="alarmist-glm",
        description="Run Poisson GLM differential expression analysis on BPTF motifs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  alarmist-glm --input-dir results/project --adata data.h5ad --output-dir results/glm
        """,
    )

    # Input/output arguments
    parser.add_argument(
        "--input-dir",
        "-i",
        type=str,
        required=True,
        help="Directory containing projection results (from alarmist-project)",
    )
    parser.add_argument(
        "--adata", type=str, required=True, help="Path to AnnData h5ad file"
    )
    parser.add_argument(
        "--cell-type-column",
        type=str,
        default="cell_type",
        help=(
            "Column in adata.obs containing cell type labels used to group the "
            "per-cell-type DE analysis (default: cell_type). The projected adata "
            "from alarmist-project preserves whatever column was passed via its "
            "own --cell-type-column flag (e.g. predicted_cell_type)."
        ),
    )
    common.add_output_arguments(parser)

    # Analysis parameters
    parser.add_argument(
        "--patch-lri-dir",
        type=str,
        default=None,
        help="Patch-LRI results directory (defaults to input-dir/../patchify)",
    )
    parser.add_argument(
        "--count-layer",
        type=str,
        default="X",
        help='Which layer to use: "X", "raw", or "layers:NAME" (default: X)',
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="FDR significance threshold (default: 0.05)",
    )
    parser.add_argument(
        "--prefilter-spearman",
        dest="prefilter_spearman",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Pre-filter genes via Spearman correlation before fitting GLMs "
            "(default: enabled). Use --no-prefilter-spearman to disable."
        ),
    )
    parser.add_argument(
        "--spearman-pval-threshold",
        type=float,
        default=0.001,
        help="Spearman p-value threshold for pre-filtering (default: 0.001)",
    )
    parser.add_argument(
        "--spearman-chunk-size",
        type=int,
        default=1000,
        help="Chunk size for Spearman correlation computation (default: 1000)",
    )
    parser.add_argument(
        "--backend",
        type=str,
        choices=["sklearn", "torch"],
        default="sklearn",
        help="GLM fitting backend. 'torch' fits all genes per motif x cell-type "
        "in one batched IRLS solve (much faster; runs on GPU if available). "
        "'sklearn' is the per-gene reference. Default: sklearn.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        help="Torch device for --backend torch: 'auto' (cuda if available else "
        "cpu), 'cpu', or 'cuda'. Default: auto.",
    )
    parser.add_argument(
        "--glm-dtype",
        type=str,
        choices=["float64", "float32"],
        default="float64",
        help="Numeric precision for the torch backend. float64 matches sklearn; "
        "float32 is ~2x faster on GPU. Default: float64.",
    )
    parser.add_argument(
        "--gene-tile",
        type=int,
        default=2048,
        help="Genes per tile for the torch backend (bounds memory: "
        "~n_cells x gene_tile x 8 bytes). Default: 2048.",
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
    import pandas as pd
    import scanpy as sc

    import alarmist as al

    # Set random seed
    np.random.seed(args.seed)
    logger.info(f"Random seed: {args.seed}")

    # Load input directory
    input_path = Path(args.input_dir)

    # Load cell loadings from project results
    cell_loadings_file = input_path / "cell_motif_loadings.parquet"
    if not cell_loadings_file.exists():
        logger.error(f"Cell loadings not found: {cell_loadings_file}")
        sys.exit(1)

    logger.info(f"Loading cell loadings from {cell_loadings_file}")
    cell_loadings_df = pd.read_parquet(cell_loadings_file)
    cell_loadings = cell_loadings_df.values
    logger.info(f"Loaded cell loadings shape: {cell_loadings.shape}")

    # Load AnnData
    logger.info(f"Loading AnnData from {args.adata}")
    adata = sc.read_h5ad(args.adata)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Find patch-LRI directory for column names
    if args.patch_lri_dir:
        patch_lri_dir = Path(args.patch_lri_dir)
    else:
        candidates = [
            input_path.parent / "patchify",
            input_path.parent / "patch_lri",
            input_path / "patchify",
            input_path / "patch_lri",
        ]
        patch_lri_dir = None
        for c in candidates:
            if c.exists() and (c / "patch_lri_columns.csv").exists():
                patch_lri_dir = c
                break
        if patch_lri_dir is None:
            logger.error("Could not find patch-LRI directory with column names")
            sys.exit(1)

    # Load LRI column names
    lri_columns_file = patch_lri_dir / "patch_lri_columns.csv"
    logger.info(f"Loading LRI column names from {lri_columns_file}")
    lri_columns_df = pd.read_csv(lri_columns_file)
    lri_column_names = lri_columns_df["column_name"]
    logger.info(f"Loaded {len(lri_column_names)} LRI columns")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run GLM analysis
    logger.info("Running Poisson GLM analysis...")
    try:
        al.run_poisson_glm_analysis(
            cell_loadings=cell_loadings,
            adata=adata,
            lri_column_names=lri_column_names,
            output_dir=str(output_dir),
            count_layer=args.count_layer,
            alpha=args.alpha,
            random_state=args.seed,
            cell_type_column=args.cell_type_column,
            prefilter_spearman=args.prefilter_spearman,
            spearman_pval_threshold=args.spearman_pval_threshold,
            spearman_chunk_size=args.spearman_chunk_size,
            backend=args.backend,
            device=args.device,
            glm_dtype=args.glm_dtype,
            gene_tile=args.gene_tile,
        )

        logger.info("GLM analysis complete")
        logger.info(f"Results saved to: {args.output_dir}")

    except Exception as e:
        logger.error(f"GLM analysis failed: {e}")
        raise

    return 0


if __name__ == "__main__":
    sys.exit(main())
