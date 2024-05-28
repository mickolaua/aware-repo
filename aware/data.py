"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
data.py (c) 2023
Desc: telegram data objects and other useful stuff
Created:  2023-05-21
Modified: 2024-04-02
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
import os

from .site.main import Site
from .alert.target_info import TargetInfo
from .config import CfgOption


products_dir = CfgOption(
    "products_dir", "~/aware/products/", lambda _p: str(os.path.expanduser(_p))
)
cache_dir = CfgOption(
    "cache_dir", "~/aware/cache/", lambda _p: str(os.path.expanduser(_p))
)


def make_dirs(dirs: list[str]):
    for d in dirs:
        if not os.path.exists(d):
            # exist_ok does not required here, but for additional security i put it here
            os.makedirs(d, exist_ok=True)
        else:
            if os.path.isfile(d):
                raise NotADirectoryError(
                    1,
                    "Path must be a directory storing AWARE products, got a file",
                    d,
                )
            

make_dirs([products_dir.value, cache_dir.value])


@dataclass(frozen=True)
class AlertMessage:
    id: str
    senter: str
    body: str
    alert_type: str
    first_mention: str = ""
    is_new: bool = False
    update: bool = False
    target_info: TargetInfo = field(default_factory=lambda: TargetInfo(None))
    created: datetime = datetime.now()


@dataclass(frozen=True)
class DataPackage:
    id: str
    alert_type: str
    site: Site
    plan_filename: str
    plot_fname: str
    comment: str = ""
    target_info: TargetInfo = field(default_factory=lambda: TargetInfo(None))
    created: datetime = datetime.now()
    