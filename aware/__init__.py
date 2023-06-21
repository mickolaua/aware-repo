"""
# ----------------------------------------------------------------------------- 
# Project:     AWARE 
# Name:        aware.__init__ 
# Purpose:     AWARE package initializer 
# 
# Author:      npank 
# 
# Created:     2022-08-30 
# Copyright:   (c) 2004-2022 AWARE Developers 
# ----------------------------------------------------------------------------
"""
import collections

__modules__ = [
    "config", 
    "consumer", 
    "logger",
    "angle",
    "alert",
    "field",
    "plugin", 
    "sql", 
    "topic", 
    "planner",
    "cosmology",
    "localization",
    "visualtization",
    "voevent",
    "json",
    "cache",
    "data",
    "site",
    "credentials"
]

# Patch astroplan since Sequence was moved to collections.abc
try:
    from collections.abc import Sequence
    collections.Sequence = Sequence
except ImportError:
    pass

# Load config options
from . import config
