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

import asyncio
import collections
import traceback
from types import ModuleType
from ruamel.yaml import YAML

import uvloop

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
    "credentials",
    "path"
]

# Patch astroplan since Sequence was moved to collections.abc
try:
    from collections.abc import Sequence

    collections.Sequence = Sequence
except ImportError:
    pass

# Initialize all CfgOptions
import os


# Small utility to automatically load modules
# Borrowed from https://gist.github.com/dorneanu/cce1cd6711969d581873a88e0257e312
def load_module(path: str) -> ModuleType:
    from importlib import util

    name = os.path.split(path)[-1]
    spec = util.spec_from_file_location(name, path)
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_all_modules():
    path = os.path.abspath(__file__)
    dirpath = os.path.dirname(path)

    # Load all modules to intialize configuration options and plugins
    for root, dirnames, fnames in os.walk(dirpath):
        for fname in fnames:
            # Load only "real modules"
            if (
                not fname.startswith(".")
                and not fname.startswith("__")
                and fname.endswith(".py")
                and not fname.endswith("config.py") # Do not load config.py twice
            ):
                try:
                    module = load_module(os.path.join(root, fname))
                except Exception as e:
                    ...
                    # log.debug("Failed to load module: %s", e, exc_info=e)

# Load config options, doing this here, since we stuck in circular imports when
# scanning the modules from config submodule.
# ! Note: one should NOT use cfg objects directly
# from aware.config import load_config, cfg

# yaml_cfg = load_config()


def deep_update(src: dict, updates: dict):
    """Do a deep update of the src dict using updates dict.

    This will update src with values in updates, but will not delete keys in src
    not found in updates at some arbitrary depth of src. That is, updates is deeply
    merged into src.

    Parameters
    ----------
    src : dict
        original configuration dictionary
    updates : dict
        dictionary of configuration overrides

    Note: this is destructive to src, but not updates.
    Borrowed from: https://stackoverflow.com/a/52099238

    Returns
    -------
        Updates dictionary inplace
    """
    stack = [(src, updates)]
    while stack:
        src, updates = stack.pop(0)
        for k, v in updates.items():
            if not isinstance(v, dict):
                # src[k] is not a dict, nothing to merge, so just set it,
                # regardless if src[k] *was* a dict
                src[k] = v
                if src.get(k):
                    if isinstance(src[k], tuple):
                        src[k][1] = v[1]
                else:
                    src[k] = v                   
            else:
                # note: updates[k] is a dict
                if k not in src:
                    # add new key into src
                    src[k] = v
                elif not isinstance(src[k], dict):
                    # src[k] is not a dict, so just set it to updates[k],
                    # overriding whatever it was
                    if src.get(k):
                        if isinstance(src[k], tuple):
                            src[k][1] = v[1]
                    else:
                        src[k] = v   
                else:
                    # both src[k] and updates[k] are dicts, push them on the stack
                    # to merge
                    stack.append((src[k], v))


# Update configuration inplace
# deep_update(cfg, yaml_cfg)

# # It is save to remove configuration here, because it is presented it config module
# del cfg
# del yaml_cfg


asyncio.set_event_loop_policy(uvloop.EventLoopPolicy())