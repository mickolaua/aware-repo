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
from logging import LogRecord

import colorama

from .config import CfgOption

__all__ = ["log"]

# Config options
verbosity = CfgOption("verbosity", value="INFO", 
                      typ=lambda val: logging._nameToLevel[val])
filename = CfgOption("filename", value="aware.log", typ=str)
filemode = CfgOption("filemode", value="w+", typ=str)


# Set logger
log = logging.getLogger("aware")
log.setLevel(verbosity.value)


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
console.setFormatter(ColorFormatter(msg_format))

def _remove_handlers() -> None:
    if log.hasHandlers():
        log.handlers = []

# Remove handlers to prevent double message emmiting 
_remove_handlers()

log.addHandler(console)

logfile = logging.FileHandler(filename.get_value(), mode=filemode.get_value())
logfile.setLevel(verbosity.get_value())
logfile.setFormatter(logging.Formatter(msg_format, "%Y-%m-%d %H:%M:%S"))
log.addHandler(logfile)

colorama.init()
atexit.register(colorama.deinit)
atexit.register(_remove_handlers)
