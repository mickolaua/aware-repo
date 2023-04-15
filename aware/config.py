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
from typing import Any, Callable, Mapping, Sequence

import dotenv
import yaml


ABS_DIR = os.path.abspath(os.path.dirname(__file__))


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
        print("Config file not found!", file=sys.stderr)
        sys.exit(1)

    with open(os.path.expanduser(os.path.realpath(path)), "rb") as f:
        cfg = yaml.load(f, loader)

    return cfg


# Parse config
cfg = load_config()


class CfgOption:
    def __init__(
        self, name: str, value: Any, typ: Callable[[Any], Any]
    ) -> None:
        # Get caller's module name
        # TODO: find more beauty way of find the caller's name
        try:
            raise ValueError
        except ValueError as e:
            frame = sys.exc_info()[2].tb_frame.f_back
            caller_basename = os.path.splitext(frame.f_code.co_filename)[0]
            start_idx = caller_basename.find("aware" + os.sep)
            fullmodname = caller_basename[start_idx:].replace(
                os.sep, ".").replace("..", ".")
            self.name = fullmodname.lstrip("aware.") + "." + name

        self._value = value
        self.typ = typ

    def get_section(self) -> str:
        return ".".join(self.name.split("."))

    def get_value(self) -> Any:
        global cfg
        names = self.name.split(".")
        val = self._get_val_from_nested_dict(cfg, names)
        if val is not None:
            return self.typ(val)
        else:
            return self.typ(self._value)

    def _get_val_from_nested_dict(
            self, d: Mapping[Any, Any], names: Sequence[str]
    ) -> Any:
        names = list(names)
        while names:
            name = names.pop(0)
            d = d.get(name, None)
            if not isinstance(d, dict):
                return d

        return None


dev = CfgOption("dev", False, bool)
dot_env_path = CfgOption("dot_env_path", ".env", str)
if dev.get_value():
    if not dotenv.load_dotenv(dotenv_path=dot_env_path.get_value(), verbose=True):
        raise FileNotFoundError(".env file not found")


