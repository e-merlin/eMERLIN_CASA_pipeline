"""
Simple configuration for switching between original and modern weblog
"""

# Set to True to use the modernized weblog, False for original
USE_MODERN_WEBLOG = True

def get_weblog_function():
    """Return the appropriate weblog function based on configuration"""
    if USE_MODERN_WEBLOG:
        from ..weblog.eMCP_weblog_main import start_weblog
    else:
        from ..weblog.eMCP_weblog import start_weblog
    
    return start_weblog
