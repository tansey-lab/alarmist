"""
Common CLI argument definitions (following modulator's pattern)
"""

import argparse
import logging

logger = logging.getLogger(__name__)


def add_data_input_arguments(parser: argparse.ArgumentParser):
    """Add common data input arguments"""
    parser.add_argument(
        '--data-file',
        type=str,
        required=True,
        help='Path to AnnData h5ad file'
    )
    parser.add_argument(
        '--cell-type-column',
        type=str,
        default='cell_type',
        help='Column name for cell types in adata.obs (default: cell_type)'
    )
    return parser


def add_output_arguments(parser: argparse.ArgumentParser):
    """Add output directory arguments"""
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--make-plots',
        default=True,
        action=argparse.BooleanOptionalAction,
        help='Generate visualization plots (default: True)'
    )
    return parser


def add_lri_arguments(parser: argparse.ArgumentParser):
    """Add LRI analysis arguments"""
    parser.add_argument(
        '--resource',
        type=str,
        default='cellchatdb',
        help='LRI database resource name (default: cellchatdb)'
    )
    parser.add_argument(
        '--patch-size',
        type=float,
        default=50.0,
        help='Size of spatial patches in micrometers (default: 50.0)'
    )
    parser.add_argument(
        '--spliter',
        type=str,
        default='|',
        help='Separator for column names (default: |)'
    )
    return parser


def add_random_state_argument(parser: argparse.ArgumentParser, default: int = 0):
    """Add random state argument"""
    parser.add_argument(
        '--random-state',
        type=int,
        default=default,
        help=f'Random seed for reproducibility (default: {default})'
    )
    return parser


def add_bptf_arguments(parser: argparse.ArgumentParser):
    """Add BPTF-specific arguments"""
    parser.add_argument(
        '--n-components',
        type=int,
        default=15,
        help='Number of latent factors/motifs (default: 15)'
    )
    parser.add_argument(
        '--max-iter',
        type=int,
        default=10000,
        help='Maximum number of iterations for BPTF (default: 10000)'
    )
    return parser
