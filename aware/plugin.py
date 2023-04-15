"""
# ----------------------------------------------------------------------------- 
# Project:     AWARE 
# Name:        aware.plugin 
# Purpose:     Plugin manager and hooks
# 
# Author:      npank 
# 
# Created:     2022-08-30 
# Copyright:   (c) 2004-2022 AWARE Developers 
# ----------------------------------------------------------------------------
"""
from ast import Module
import importlib
from types import ModuleType
from pluggy import HookimplMarker, HookspecMarker, PluginManager
import pkgutil
from .logger import log
import sys

__all__ = ["hookimpl", "hookspec", "plugin_manager", "hook_add_specs"]

hookimpl = HookimplMarker("AWARE")
hookspec = HookspecMarker("AWARE")
plugin_manager = PluginManager("AWARE")


def _register_subclasses(mod: ModuleType) -> None:
    for name, obj in vars(mod).items():
        registered = False
        try:
            if mod.__name__ == "aware.alert":
                if issubclass(obj, mod.AlertParser) and obj is not mod.AlertParser:
                    mod.alert_parsers.add_parser(obj.topic, obj)
                    registered = True
            elif mod.__name__ ==  "aware.planner":
                if isinstance(obj, mod.Site) and obj is not mod.Site:
                    mod.sites[obj.name] = obj
                    registered = True
        except TypeError as e:
            pass
        else:
            if registered:
                log.debug("found and registered plugin '%s' in %s", name, mod.__name__)


def register_plugins(path: str) -> None:
    package = "aware"

    log.debug("\nlooking for plugins in submodules")
    for loader, submodule_name, is_pkg in pkgutil.walk_packages(path):
        # Do not lookup plugin submodule
        if submodule_name == __name__:
            continue

        module_name = package + "." + submodule_name

        try:
            mod = importlib.import_module(module_name)
        except Exception as e:
            log.exception(e)
            continue
        
        # Temporally hardcoded solution, in future Plugin class must be implemented
        _register_subclasses(mod)
            
        # Will be loaded latter
        del sys.modules[module_name]


def hook_add_specs(path: str) -> None:
    """
    Find and register hooks in modules.

    Parameters
    ----------
    path : str
        root packager
    """

    package = "aware"

    for loader, submodule_name, is_pkg in pkgutil.walk_packages(path):
        # Do not lookup plugin submodule
        if submodule_name == __name__:
            continue

        module_name = package + "." + submodule_name

        try:
            mod = importlib.import_module(module_name)
        except Exception as e:
            continue

        try:
            plugin_manager.add_hookspecs(mod)
        except Exception as e:
            continue
        else:
            log.debug("found hook specification(s) in module: '%s'", module_name)
            try:
                plugin_names = plugin_manager.register(sys.modules[module_name])
            except Exception as e:
                pass
            else:
                log.debug("hook implementation(s) from module '%s' registrated", module_name)
