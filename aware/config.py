# -*- coding: utf-8 -*-
# @Author: npank
# @Date:   2022-07-21
# @Last Modified by:   mickolaua
# @Last Modified time: 2022-10-04
from __future__ import annotations

# from importlib import util
import importlib
import inspect
import os
import sys
import traceback
from types import ModuleType
import warnings
from pathlib import Path
from typing import Any, Callable, Generic, Mapping, Sequence, TypeVar

import dotenv
import yaml

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
ABS_DIR = os.path.abspath(BASE_DIR)


def load_config(
    path: str | Path | None = None, loader: yaml.Loader = yaml.SafeLoader
) -> dict[str, Any]:
    """Load the YAML configuration from a file or enviromental variable.

    Parameters
    ----------
    path : str | Path | None, optional
        path to the configurational file, by default None
    loader : yaml.Loader, optional
        which loader is used to parse the raw config, by default
        yaml.SafeLoader

    Returns
    -------
    cfg : dict[str, Any]
        the config in the form of a dictionary
    """
    # When path is not provided, look into enviromental variable
    if not path:
        path = os.getenv("AWARE_CONFIG_FILE", "aware.yaml")

    full_path = os.path.expanduser(os.path.realpath(path))
    cfg = {}
    if not os.path.exists(full_path):
        print(
            "AWARE config file not found, using the default settings", file=sys.stderr
        )
    else:
        try:
            with open(full_path, "rb") as f:
                cfg = yaml.load(f, loader)
        except IsADirectoryError:
            print(f"Path to config file is not a file: {full_path}")
        except (OSError, IOError) as e:
            print("Error loading configuration, traceback follows")
            traceback.print_exc()

    return cfg


# Configuration is stored here
cfg = load_config()


def set_config_option_value(name: str, val: Any, typ: Callable = str):
    global cfg

    # The option is in the submodule
    names = name.split(".")
    if len(names) > 1:
        curr_d = cfg
        while names:
            n = names.pop(0)
            v = curr_d.get(n, None)
            if not v:
                curr_d[n] = {}

            if not names:
                curr_d[n] = typ(val)
            else:
                curr_d = curr_d[n]

    else:
        # The option at the root level
        cfg[name] = typ(val)


# Option return type
__T = TypeVar("__T")


class CfgOption(Generic[__T]):
    """
    YAML config option.

    Attributes
    ----------
    name : str
        a name of the option inside the YAML file
    value : Any
        a default value for the option if it is not specified in config
    typ : callable
        a function that processes and validates the value type
    """

    def __init__(self, name: str, value: Any, typ: Callable[..., __T]) -> None:
        global cfg
        self._name = name
        self.typ = typ

        # Get caller filename here or it will always return filename where is the base
        # class located
        self.__filename = inspect.currentframe().f_back.f_code.co_filename

        # ! All additional initialization MUST be done down below, since the full option
        # ! name is determined after previous line.
        self._value = self._get_val_from_nested_dict(cfg, self.name.split('.')) or value
        
        # Adding option to config dictionary, so one can dump the configuration to file
        set_config_option_value(self.name, self._value, lambda x: x)

    def __get_aware_pkg_root(self) -> str:
        """
        Get the path to the AWARE root directory
        """
        mod = __import__("aware")
        return mod.__file__

    def __get_caller_name(self) -> str:
        """
        Get the name of the file where the caller is defined
        """
        try:
            root_dir = os.path.dirname(self.__get_aware_pkg_root())
            start = len(root_dir) + 1
            return os.path.splitext(self.__filename[start:])[0]
        except AttributeError:
            return ""

    @property
    def name(self) -> str:
        """Name of the option relative to the package root"""
        try:
            name = (
                (
                    self.__get_caller_name()
                    .replace(os.path.dirname(self.__get_aware_pkg_root()), "")
                    .replace(os.sep, ".")
                )
                + "."
                + self._name
            )
        except (AttributeError, IndexError):
            name = self._name

        if name.startswith("."):
            name = name[1:]

        return name

    @property
    def section(self) -> str:
        """Get the full name of the submodule in which this option was defined"""
        try:
            return ".".join(self.name.split(".")[:-1])
        except IndexError:
            return ""

    def get_value(self) -> __T:
        warnings.warn(
            (
                "get_value() considered deprecated in next releases, use value "
                "property instead"
            ),
            DeprecationWarning,
        )
        return self.__get_value()

    def __get_value(self) -> __T:

        global cfg
        names = self.name.split(".")

        if not cfg:
            return self.typ(self._value)

        val = self._get_val_from_nested_dict(cfg, names)

        if val is not None:
            return self.typ(val)
        else:
            return self.typ(self._value)

    @property
    def value(self) -> __T:
        return self.__get_value()

    def _get_val_from_nested_dict(
        self, d: Mapping[Any, Any], names: Sequence[str]
    ) -> Any:
        """Get the value from the nested dictionary (raw YAML config is the one)"""
        names = list(names)
        while names:
            name = names.pop(0)
            d = d.get(name, None)
            if not isinstance(d, dict):
                return d

        return None


def save_cfg(path: str = ""):
    if not path:
        path = os.getenv("AWARE_CONFIG_FILE", "aware.yaml")

    with open(os.path.expanduser(os.path.realpath(path)), "w") as f:
        yaml.safe_dump(cfg, f, default_flow_style=False)


# Development switch
dev = CfgOption("dev", False, bool)

# load .env files with defined environment variables, if available
dot_env_path = CfgOption("dot_env_path", ".env", str)

if dot_env_path.value:
    if not dotenv.load_dotenv(dotenv_path=dot_env_path.get_value(), verbose=True):
        print(
            f"No environment variables found in {dot_env_path.get_value()}",
            file=sys.stderr,
        )
