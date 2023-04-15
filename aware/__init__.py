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
from .plugin import register_plugins


__modules__ = [
    "config", 
    "consumer", 
    "logger", 
    "alert", 
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
    import collections
    collections.Sequence = collections.abc.Sequence
except ImportError:
    pass

