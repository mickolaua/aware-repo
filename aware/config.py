# -*- coding: utf-8 -*-
# @Author: npank
# @Date:   2022-07-21
# @Last Modified by:   mickolaua
# @Last Modified time: 2022-10-04
from __future__ import annotations

import inspect
import os
import sys
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence, TypeVar

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

    if not os.path.exists(os.path.expanduser(os.path.realpath(path))):
        print(
            "AWARE config file not found, using the default settings", file=sys.stderr
        )
        return {}

    with open(os.path.expanduser(os.path.realpath(path)), "rb") as f:
        cfg = yaml.load(f, loader)

    return cfg


# Parse config
cfg = load_config()

# Option return type
__T = TypeVar("__T")


class CfgOption:
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

    def __init__(self, name: str, value: Any, typ: Callable[[Any], Any]) -> None:
        self._name = name
        self._value = value
        self.typ = typ

        # Get caller filename here or it will always return filename where is the base
        # class located
        self.__filename = inspect.currentframe().f_back.f_code.co_filename
 
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
            return (
                (
                    self.__get_caller_name()
                    .replace(os.path.dirname(self.__get_aware_pkg_root()), "")
                    .replace(os.sep, ".")
                )
                + "."
                + self._name
            )
        except (AttributeError, IndexError):
            return self._name

    @property
    def section(self) -> str:
        """Get the full name of the submodule in which this option was defined"""
        try:
            return ".".join(self.name.split(".")[:-1])
        except IndexError:
            return ""

    def get_value(self) -> Any:
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
    def value(self) -> Any:
        return self.get_value()
    

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


dev = CfgOption("dev", False, bool)
dot_env_path = CfgOption("dot_env_path", ".env", str)
if dev.get_value():
    if not dotenv.load_dotenv(dotenv_path=dot_env_path.get_value(), verbose=True):
        raise FileNotFoundError(".env file not found")
