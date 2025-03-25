import logging
import os
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path

def setup_logging(log_level="INFO"):
    """
    Configure logging for the application.
    
    Args:
        log_level: The minimum log level to record (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    """
    # Convert string log level to actual log level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")
    
    # Create logs directory if it doesn't exist
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(numeric_level)
    
    # Clear any existing handlers to avoid duplicates
    if root_logger.handlers:
        root_logger.handlers.clear()
    
    # Configure console handler with a higher log level
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(numeric_level)
    console_format = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    console_handler.setFormatter(console_format)
    root_logger.addHandler(console_handler)
    
    # Configure file handler for all log messages
    file_handler = RotatingFileHandler(
        log_dir / "d3ai.log",
        maxBytes=10485760,  # 10MB
        backupCount=5,
        encoding="utf-8"
    )
    file_handler.setLevel(numeric_level)
    file_format = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(pathname)s:%(lineno)d - %(message)s'
    )
    file_handler.setFormatter(file_format)
    root_logger.addHandler(file_handler)
    
    # Configure error file handler for ERROR and above
    error_file_handler = RotatingFileHandler(
        log_dir / "error.log",
        maxBytes=10485760,  # 10MB
        backupCount=5,
        encoding="utf-8"
    )
    error_file_handler.setLevel(logging.ERROR)
    error_file_handler.setFormatter(file_format)
    root_logger.addHandler(error_file_handler)
    
    # Suppress excessively verbose logs from libraries
    for logger_name in ["uvicorn", "uvicorn.access", "uvicorn.error"]:
        lib_logger = logging.getLogger(logger_name)
        lib_logger.handlers.clear()
        lib_logger.addHandler(console_handler)
        lib_logger.addHandler(file_handler)
        lib_logger.propagate = False
    
    # Log that the logger was configured
    root_logger.info(f"Logging configured with level {log_level}")
    
    return root_logger 