"""
Logging configuration for alarmist package
"""
import logging
import sys


def add_logging_args(parser):
    """Add logging arguments to argument parser"""
    parser.add_argument(
        '--log-level',
        type=str,
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help='Logging level'
    )
    return parser


def configure_logging(args):
    """Configure logging based on arguments"""
    log_level = getattr(logging, args.log_level.upper())

    # Configure root logger
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

    # Set level for alarmist logger
    logger = logging.getLogger('alarmist')
    logger.setLevel(log_level)
