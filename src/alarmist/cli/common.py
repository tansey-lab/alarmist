"""
Common CLI argument definitions (following modulator's pattern)
"""

import argparse


def add_adata_arguments(parser: argparse.ArgumentParser):
    """Add common adata input arguments"""
    parser.add_argument(
        "--adata", type=str, required=True, help="Path to input AnnData h5ad file"
    )
    parser.add_argument(
        "--cell-type-column",
        type=str,
        default="cell_type",
        help="Column name for cell types in adata.obs (default: cell_type)",
    )


def add_output_arguments(parser: argparse.ArgumentParser):
    """Add output directory arguments"""
    parser.add_argument(
        "--output-dir",
        "-o",
        type=str,
        required=True,
        help="Output directory for results",
    )


def add_multi_sample_arguments(parser: argparse.ArgumentParser):
    """Add multi-sample support arguments"""
    parser.add_argument(
        "--multi-sample",
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Enable multi-sample mode (default: False)",
    )
    parser.add_argument(
        "--sample-column",
        type=str,
        default=None,
        help="Column name for sample IDs in adata.obs (required if --multi-sample)",
    )


def add_lri_arguments(parser: argparse.ArgumentParser):
    """Add LRI analysis arguments"""
    parser.add_argument(
        "--resource",
        type=str,
        default="cellchatdb",
        choices=["cellchatdb", "cellphonedb", "connectomedb"],
        help="LRI database resource name (default: cellchatdb)",
    )
    parser.add_argument(
        "--cellchatdb-path",
        type=str,
        default=None,
        help="Path to custom CellChatDB CSV file (optional)",
    )
    parser.add_argument(
        "--patch-size",
        type=float,
        default=50.0,
        help="Size of spatial patches in micrometers (default: 50.0)",
    )


def add_bptf_arguments(parser: argparse.ArgumentParser):
    """Add BPTF-specific arguments"""
    parser.add_argument(
        "--n-components",
        type=int,
        default=15,
        help="Number of latent factors/motifs (default: 15)",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=300,
        help="Maximum number of iterations for BPTF (default: 300)",
    )


def add_seed_argument(parser: argparse.ArgumentParser, default: int = 42):
    """Add random seed argument"""
    parser.add_argument(
        "--seed",
        type=int,
        default=default,
        help=f"Random seed for reproducibility (default: {default})",
    )
