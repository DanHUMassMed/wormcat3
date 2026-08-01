"""
wormcat3.logger
~~~~~~~~~~~~~~~
Centralized logging configuration and utilities for wormcat3.
"""

import logging
import os
import sys
from pathlib import Path
from typing import Optional, Union

LOGGER_NAME = "wormcat3"
CONSOLE_FORMAT = "%(levelname)s [%(name)s] %(message)s"
DETAILED_FORMAT = "%(asctime)s [%(levelname)s] %(name)s standard_path=%(filename)s:%(lineno)d: %(message)s"

# Attach a NullHandler to root library logger by default
_root_logger = logging.getLogger(LOGGER_NAME)
if not _root_logger.handlers:
    _root_logger.addHandler(logging.NullHandler())


def get_logger(module_name: Optional[str] = None) -> logging.Logger:
    """
    Get a child logger under the 'wormcat3' hierarchy.

    Args:
        module_name: Name of module (e.g., __name__ or sub-component name).

    Returns:
        logging.Logger instance properly scoped under 'wormcat3'.
    """
    if module_name:
        name = module_name if module_name.startswith(LOGGER_NAME) else f"{LOGGER_NAME}.{module_name}"
        return logging.getLogger(name)
    return logging.getLogger(LOGGER_NAME)


def configure_logging(
    level: Union[int, str] = logging.INFO,
    log_file: Optional[Union[str, Path]] = None,
    disabled: bool = False,
    format_str: Optional[str] = CONSOLE_FORMAT,
) -> logging.Logger:
    """
    Configure top-level logger handlers and level.
    Respects WORMCAT_LOG_LEVEL and WORMCAT_LOG_PATH environment variables if set.

    Args:
        level: Minimum log level (e.g. "DEBUG", "INFO", "WARNING", "ERROR").
        log_file: Optional file path to record log messages. If None, respects WORMCAT_LOG_PATH.
        disabled: If True, silences all logger output.
        format_str: Custom format string for console log messages.

    Returns:
        Configured root library Logger instance.
    """
    root_logger = logging.getLogger(LOGGER_NAME)
    root_logger.handlers.clear()

    # Check environment variable override for log level
    env_level = os.getenv("WORMCAT_LOG_LEVEL")
    if env_level:
        if env_level.upper() in ["OFF", "DISABLE", "FALSE", "NONE", "0"]:
            disabled = True
        else:
            level = env_level.upper()

    # Check environment variable override for log path if log_file is not explicitly set
    if log_file is None:
        env_log_path = os.getenv("WORMCAT_LOG_PATH")
        if env_log_path:
            log_file = env_log_path

    if disabled:
        root_logger.addHandler(logging.NullHandler())
        root_logger.disabled = True
        return root_logger

    root_logger.disabled = False
    if isinstance(level, str):
        level_num = getattr(logging, level.upper(), None)
        if level_num is None:
            level_num = logging.INFO
        level = level_num

    root_logger.setLevel(level)

    # Console Stream Handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(logging.Formatter(format_str))
    root_logger.addHandler(ch)

    # Optional File Handler
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_path, encoding="utf-8")
        fh.setLevel(level)
        fh.setFormatter(logging.Formatter(DETAILED_FORMAT))
        root_logger.addHandler(fh)

    return root_logger


def set_log_level(level: Union[int, str]) -> None:
    """Easily change log level at runtime."""
    configure_logging(level=level)


def disable_logging() -> None:
    """Turn off all logging output across wormcat3."""
    configure_logging(disabled=True)


def enable_logging(level: Union[int, str] = logging.INFO) -> None:
    """Enable logging with a specified log level."""
    configure_logging(level=level, disabled=False)
