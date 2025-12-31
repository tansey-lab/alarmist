#!/usr/bin/env python3
"""
06 - Single Cell Analysis

Based on 06_single_cell.py
"""

import argparse
import logging

from alarmist import log_config

logger = logging.getLogger(__name__)


def get_parser():
    """Create argument parser"""
    parser = argparse.ArgumentParser(description='Single cell analysis')
    log_config.add_logging_args(parser)
    return parser


def main():
    """Main function"""
    args = get_parser().parse_args()
    log_config.configure_logging(args)

    logger.warning("Single cell analysis not yet implemented.")
    logger.info("Please use the original script: scripts/06_single_cell.py")


if __name__ == '__main__':
    main()
