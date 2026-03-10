"""
Logging configuration for CLI commands
"""

import logging
import sys


def add_logging_args(parser):
    """Add logging-related arguments to parser"""
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress all output except errors"
    )


def configure_logging(args):
    """
    Configure logging based on CLI arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments containing verbose and quiet flags

    Returns
    -------
    logging.Logger
        Configured logger for alarmist
    """
    logger = logging.getLogger("alarmist")

    # Clear existing handlers
    logger.handlers = []

    # Create console handler
    handler = logging.StreamHandler(sys.stdout)

    # Set level based on flags
    if hasattr(args, "quiet") and args.quiet:
        logger.setLevel(logging.ERROR)
        handler.setLevel(logging.ERROR)
    elif hasattr(args, "verbose") and args.verbose:
        logger.setLevel(logging.DEBUG)
        handler.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        handler.setLevel(logging.INFO)

    # Create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(formatter)

    # Add handler to logger
    logger.addHandler(handler)

    return logger
