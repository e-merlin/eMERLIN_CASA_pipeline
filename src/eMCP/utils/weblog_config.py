"""Weblog entry point selection."""


def get_weblog_function():
    """Return the supported weblog function."""
    from ..weblog.eMCP_weblog_main import start_weblog
    return start_weblog
