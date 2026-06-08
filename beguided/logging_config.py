"""Logging configuration for BEguided."""

import logging
import sys


def setup_logging(level=logging.INFO):
    """Configure root logger for the package.

    Args:
        level: Logging level (default: INFO).

    Returns:
        logging.Logger: The configured logger.
    """
    logger = logging.getLogger("beguided")
    logger.setLevel(level)

    # Avoid adding duplicate handlers if already configured
    if logger.handlers:
        return logger

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%H:%M:%S",
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


def get_logger(name):
    """Get a logger under the beguided namespace."""
    return logging.getLogger(f"beguided.{name}")
