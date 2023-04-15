"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
cache.py (c) 2023
Desc: description
Created:  2023-03-05
Modified: 2023-03-06
"""

from __future__ import annotations

from diskcache import Cache, JSONDisk, Disk
from .config import CfgOption


cache_dir = CfgOption("cache_dir", ".cache", str)

cache = Cache(cache_dir.get_value())
