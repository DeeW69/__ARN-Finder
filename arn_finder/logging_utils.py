"""Logging helpers for ARN Finder CLI."""

from __future__ import annotations

import logging

LOGGER_NAME = "arn_finder"


def configure_logging(verbose: bool = False) -> logging.Logger:
    """Configure root logger and return the package logger."""

    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(level)
    return logger


def get_logger() -> logging.Logger:
    """Return the package logger."""

    return logging.getLogger(LOGGER_NAME)
