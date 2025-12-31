#!/usr/bin/env python3
"""
01 - Create Cell-Level LRI Matrix

Based on 01_create_cell_matrix.py
"""

import argparse
import logging

from alarmist import log_config

logger = logging.getLogger(__name__)


def get_parser():
    """Create argument parser"""
    parser = argparse.ArgumentParser(
        description='Create cell-level LRI matrix'
    )

    parser.add_argument('--data-file', required=True, help='Path to AnnData file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--resource', default='cellchatdb', help='LRI database resource')
    parser.add_argument('--no-gene-expression', action='store_true', help='Do not include gene expression')

    log_config.add_logging_args(parser)
    return parser


def main():
    """Main function"""
    args = get_parser().parse_args()
    log_config.configure_logging(args)

    logger.warning("SingleCellLRIAnalyzer not yet fully implemented.")
    logger.info("Please use the original script: scripts/01_create_cell_matrix.py")


if __name__ == '__main__':
    main()
