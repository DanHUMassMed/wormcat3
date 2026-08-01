import os
from pathlib import Path

import pytest

import wormcat3
from wormcat3.logger import (
    LOGGER_NAME,
    configure_logging,
    disable_logging,
    enable_logging,
    get_logger,
    set_log_level,
)


@pytest.fixture(autouse=True)
def reset_logging_state():
    """Ensure clean logging state before and after each test."""
    yield
    disable_logging()
    if "WORMCAT_LOG_LEVEL" in os.environ:
        del os.environ["WORMCAT_LOG_LEVEL"]
    if "WORMCAT_LOG_PATH" in os.environ:
        del os.environ["WORMCAT_LOG_PATH"]


def test_get_logger_naming():
    logger1 = get_logger()
    assert logger1.name == LOGGER_NAME

    logger2 = get_logger("wormcat.submodule")
    assert logger2.name == "wormcat3.wormcat.submodule"

    logger3 = get_logger("wormcat3.another")
    assert logger3.name == "wormcat3.another"


def test_default_logger_silent(capsys):
    disable_logging()
    logger = get_logger("test_module")
    logger.info("This should not appear anywhere")
    logger.error("This error should also be silent by default")

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_enable_logging_stdout(capsys):
    enable_logging(level="INFO")
    logger = get_logger("test_module")

    logger.info("Test info message")
    logger.debug("Test debug message - should be filtered out")

    captured = capsys.readouterr()
    assert "INFO [wormcat3.test_module] Test info message" in captured.out
    assert "Test debug message" not in captured.out


def test_set_log_level_debug(capsys):
    set_log_level("DEBUG")
    logger = get_logger("test_module")

    logger.debug("Debug detail message")
    captured = capsys.readouterr()
    assert "DEBUG [wormcat3.test_module] Debug detail message" in captured.out


def test_disable_logging(capsys):
    enable_logging("INFO")
    disable_logging()

    logger = get_logger("test_module")
    logger.info("Should be silenced after disable_logging")

    captured = capsys.readouterr()
    assert captured.out == ""


def test_env_var_override_disable(capsys):
    os.environ["WORMCAT_LOG_LEVEL"] = "OFF"
    configure_logging(level="DEBUG")

    logger = get_logger("test_module")
    logger.error("Error message when log level set to OFF via env")

    captured = capsys.readouterr()
    assert captured.out == ""


def test_env_var_override_debug(capsys):
    os.environ["WORMCAT_LOG_LEVEL"] = "DEBUG"
    configure_logging(level="ERROR")

    logger = get_logger("test_module")
    logger.debug("Debug message enabled by env var")

    captured = capsys.readouterr()
    assert "DEBUG [wormcat3.test_module] Debug message enabled by env var" in captured.out


def test_file_handler_logging(tmp_path: Path):
    log_file = tmp_path / "test_run.log"
    configure_logging(level="INFO", log_file=log_file)

    logger = get_logger("test_file")
    logger.info("Message saved to file")

    assert log_file.exists()
    content = log_file.read_text(encoding="utf-8")
    assert "[INFO] wormcat3.test_file" in content
    assert "Message saved to file" in content


def test_env_var_log_path(tmp_path: Path):
    env_log_file = tmp_path / "env_run.log"
    os.environ["WORMCAT_LOG_PATH"] = str(env_log_file)
    configure_logging(level="INFO")

    logger = get_logger("test_env_file")
    logger.info("Message saved via env var path")

    assert env_log_file.exists()
    content = env_log_file.read_text(encoding="utf-8")
    assert "Message saved via env var path" in content


def test_env_var_log_path_override(tmp_path: Path):
    env_log_file = tmp_path / "env_run.log"
    explicit_log_file = tmp_path / "explicit_run.log"
    os.environ["WORMCAT_LOG_PATH"] = str(env_log_file)
    configure_logging(level="INFO", log_file=explicit_log_file)

    logger = get_logger("test_override_file")
    logger.info("Message saved to explicit file")

    assert explicit_log_file.exists()
    assert not env_log_file.exists()
    content = explicit_log_file.read_text(encoding="utf-8")
    assert "Message saved to explicit file" in content


def test_exports_in_init():
    assert hasattr(wormcat3, "get_logger")
    assert hasattr(wormcat3, "set_log_level")
    assert hasattr(wormcat3, "enable_logging")
    assert hasattr(wormcat3, "disable_logging")
    assert hasattr(wormcat3, "configure_logging")
