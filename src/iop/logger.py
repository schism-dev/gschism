import logging
import sys
from pathlib import Path

def setup_logger(log_file: str = "gschism.log", level: int = logging.INFO) -> None:
    """
    Setup the logger with both file and console handlers.
    
    Args:
        log_file: Path to the log file
        level: Logging level (default: INFO)
    """
    # Create logger
    logger = logging.getLogger('gschism')
    logger.setLevel(level)
    
    # Create formatters
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_formatter = logging.Formatter('%(message)s')
    
    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

def log(message: str, level: int = logging.INFO) -> None:
    """
    Log a message at the specified level.
    
    Args:
        message: The message to log
        level: Logging level (default: INFO)
    """
    logger = logging.getLogger('gschism')
    logger.log(level, message)

# Convenience functions
def infoMessage(message: str) -> None:
    """Log an info message."""
    log(message, logging.INFO)

def warning(message: str) -> None:
    """Log a warning message."""
    log(message, logging.WARNING)

def error(message: str) -> None:
    """Log an error message."""
    log(message, logging.ERROR)

def debug(message: str) -> None:
    """Log a debug message."""
    log(message, logging.DEBUG) 