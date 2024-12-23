import os
import logging
import sys
def setup_logging(logger_name, log_file, level=logging.INFO):

    logger = logging.getLogger(logger_name)
    
    # Avoid adding multiple handlers if the logger already exists
    if logger.hasHandlers():
        return logger
    
    # Set the logging level
    logger.setLevel(level)

    # Ensure the log file's directory exists
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Create file handler for logging to a file
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(level)

    # Create console handler for logging to the console
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)

    # Create a logging format
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)



    return logger
