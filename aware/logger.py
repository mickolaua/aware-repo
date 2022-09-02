"""
# ----------------------------------------------------------------------------- 
# Project:     AWARE 
# Name:        aware.logger 
# Purpose:     Logging facility  
# 
# Author:      npank 
# 
# Created:     2022-08-30 
# Copyright:   (c) 2004-2022 AWARE Developers 
# ----------------------------------------------------------------------------
"""
from __future__ import annotations

import atexit
import logging

import colorama

from aware.config import cfg

__all__ = ["log"]

# Config options
verbosity = logging.getLevelName(cfg["logger"]["verbosity"])
filename = cfg["logger"]["filename"]
filemode = cfg["logger"]["filemode"]


# Set logger
log = logging.getLogger("aware")
log.setLevel(verbosity)


def colored_message(msg: str, color: str) -> str:
    """
    Return colorized version of the message.
    
    Parameters
    ----------
    msg : str
        message to colorize
    color : str
        text color
    
    Returns
    -------
    str
        message with colorized text
    """
    return color + msg + colorama.Style.RESET_ALL


class ColorFormatter(logging.Formatter):
    """
    Colorful messages in console/terminal emulator via logging facility.
    """

    def format(self, record: logging.LogRecord) -> str:
        """
        Paint record message according to the logging level.
        """
        
        msg = super(ColorFormatter, self).format(record)

        lvl = record.levelno
        if lvl < logging.INFO:
            formatted_msg = colored_message(msg, colorama.Fore.BLUE)
        elif lvl < logging.WARN:
            formatted_msg = colored_message(msg, colorama.Fore.GREEN)
        elif lvl < logging.ERROR:
            formatted_msg = colored_message(msg, colorama.Fore.YELLOW)
        elif lvl < logging.CRITICAL:
            formatted_msg = colored_message(msg, colorama.Fore.LIGHTRED_EX)
        else:
            formatted_msg = colored_message(msg, colorama.Fore.RED)

        return formatted_msg


msg_format = "%(asctime)s.%(msecs)03d %(levelname)8s "
msg_format += "%(module)s.%(funcName)s(): %(message)s"

console = logging.StreamHandler()
console.setFormatter(ColorFormatter())

def _remove_handlers() -> None:
    if log.hasHandlers():
        log.handlers = []

# Remove handlers to prevent double message emmiting 
_remove_handlers()

log.addHandler(console)

logfile = logging.FileHandler(
    cfg["logger"]["filename"], mode=cfg["logger"]["filemode"]
)
logfile.setLevel(cfg["logger"]["verbosity"])
logfile.setFormatter(logging.Formatter(msg_format, "%Y-%m-%d %H:%M:%S"))
log.addHandler(logfile)

colorama.init()
atexit.register(colorama.deinit)
atexit.register(_remove_handlers)
