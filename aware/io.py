from __future__ import annotations

from typing import Any
import numpy as np
from ligo.skymap.io.fits import read_sky_map
from ligo.skymap.moc import rasterize
import healpy as hp
import astropy_healpix as ah
from astropy.table import Table

from aware.topic import full_topic_name_to_short
from .data import products_dir
import os.path
from pathvalidate import sanitize_filename


def hpx_from_moc(
    moc: Table, order: int = 8, nest: bool = False
) -> tuple[np.ndarray, dict[str, Any]]:
    data = rasterize(moc, order=order)
    nside = ah.level_to_nside(order)
    area = ah.nside_to_pixel_area(nside)
    hpx = data["PROBDENSITY"] * area.to_value("sr")

    if not nest:
        nside = hp.get_nside(hpx)
        npix = hp.get_map_size(hpx)
        ring_idx = hp.ring2nest(nside, np.arange(npix))
        hpx = hpx[ring_idx]

    return hpx, moc.meta


def read_hpx_from_moc(
    path: str, header: bool = False
) -> np.ndarray | tuple[np.ndarray, dict]:
    moc = read_sky_map(path, moc=True, nest=True)
    hpx, hdr = hpx_from_moc(moc)

    if header:
        return hpx, hdr
    else:
        return hpx


def create_event_folder(
    origin: str, event: str, base_dir: str = products_dir.value
) -> str:
    folder = sanitize_filename(f"{origin}_{event}", replacement_text="_")
    local_path = os.path.join(base_dir, folder)
    os.makedirs(local_path, exist_ok=True)
    return local_path
