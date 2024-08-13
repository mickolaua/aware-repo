"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
cache.py (c) 2023
Desc: Caching functions
Created:  2023-03-05
Modified: 2024-02-22
"""

from __future__ import annotations
from typing import Any

from diskcache import Cache, DEFAULT_SETTINGS, Disk
from .config import CfgOption
from .logger import log
from sqlite3 import Error


__all__ = ["cache_value"]


cache_dir = CfgOption("cache_dir", ".cache", str, comment="Cache directory")
timeout = CfgOption(
    "timeout",
    60,
    float,
    comment="Maximum number of seconds to wait for retrieval of a cached value",
)
expire = CfgOption(
    "expire", 86400, float, comment="Cached value expiration time in seconds"
)


def __access_cache():
    """
    Access the cache via calling this function, do not use global `Cache` object
    since it may be not created due to SQL errors (e.g. database is blocked)
    """
    try:
        cache = Cache(cache_dir.get_value(), timeout=timeout.value)
    except Error as e:
        log.debug("Error getting cache database: %s", e)

        class EmptyCache:
            def __init__(self, *args, **kwargs): ...

            def set(self, *args, **kwargs): ...

            def get(self, *args, **kwargs): ...

            def __enter__(self):
                return self

            def __exit__(self, *args, **kwargs): ...

        cache = EmptyCache()

    return cache


def write_cache(
    key: str,
    value: Any,
    tag: str = "",
    expire: float = expire.value,
    retry: bool = True,
):
    global __access_cache
    with __access_cache() as cache:
        # Additional arguments to cache set
        kws = {}

        # File object
        if hasattr(value, "read"):
            kws["read"] = True

        if tag is None:
            kws["tag"] = key.title()

        try:
            cache.set(key, value, expire=expire, retry=retry, **kws)
        except Exception as e:
            if isinstance(e, KeyboardInterrupt):
                raise
            else:
                log.error("Cache entry %s is not writable %s", key, e)


def read_cache(
    key: str,
) -> Any:
    global __access_cache
    value = None
    with __access_cache() as cache:
        try:
            value = cache.get(key)
        except Exception as e:
            if isinstance(e, KeyboardInterrupt):
                raise
            else:
                log.error("Cache entry %s is not readable %s", key, e)

    return value


def pop_cache(key: str):
    try:
        with __access_cache() as cache:
            cache.pop(key)
    except Exception as e:
        if isinstance(e, KeyboardInterrupt):
            raise
        else:
            log.error("Cache entry %s is not expirable %s", key, e)
