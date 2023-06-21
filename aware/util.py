from __future__ import annotations
import posixpath
from .data import cache_dir
from uuid import uuid4
from .logger import log
import requests


def download_file(
    url: str,
    timeout: float = 1.0,
    path: str = cache_dir.value,
    chunk_size: int = 2048,
    overwrite: bool = False,
    raise_file_exist: bool = False
) -> bool:
    status = False

    try:
        filename = posixpath.split(url)[1]
    except IndexError as e:
        filename = uuid4().hex
        log.error(
            "could not extract filename from URL: `%s`, assigning unique filename `%s`",
            url,
            filename,
            exc_info=e,
        )

    filename = posixpath.join(path, filename)
    if posixpath.exists(filename) and posixpath.isfile(filename) and not overwrite:
        if raise_file_exist:
            raise FileExistsError(filename)
        else:
            log.warn("file already exists: `%s`, so do not download it", filename)
            status = True
    else:
        try:
            resp = requests.get(url, timeout=timeout, stream=True)
            with open(filename, "wb") as f:
                for chunk in resp.iter_content(chunk_size=chunk_size):
                    f.write(chunk)
            status = True
        except requests.HTTPError as e:
            log.error("could not download file from URL: `%s`", url, exc_info=e)

    return status