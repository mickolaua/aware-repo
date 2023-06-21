"""
Author: Artem Prokhorenko (prokhorenko.aiu@phystech.edu)
        Evgenii Shekotihin (z.shekotihin.2012@yandex.ru)
        Nicolai Pankov (colinsergesen@gmail.com)
mosaic.py (c) 2023
Desc: mosaic scanning of the localization skymap
Created:  2023-04-18
Modified: !date!
"""
from __future__ import annotations

from functools import partial

import healpy as hp
import numpy as np
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from healpy.projector import CartesianProj
from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.skymap.postprocess.util import find_greedy_credible_levels

from ..field import Field
from ..site import FOV
from ..angle import coord2str, coord_to_target_name
from ..localization.main import lvk_uncert_level

__all__ = ["mosaic_walker"]


def get_cur_p(
    map: np.ndarray, cur_p: tuple[float, float], bs: int
) -> tuple[float, float]:
    """Get current position of the mosaic walk in Cartesian coordinates.

    Parameters
    ----------
    map : np.ndarray
        a skymap reprojected to the Cartesian surface
    cur_p : tuple[float, float]
        a previous position (X, Y)
    bs : int
        a size of the box containing neighbors of the (X, Y)

    Returns
    -------
    tuple[float, float]
        a current position (X, Y)
    """

    box = map[cur_p[0] - bs : cur_p[0] + bs + 1, cur_p[1] - bs : cur_p[1] + bs + 1]

    if box.argmax() == 0:
        while box.argmax() == 0:
            bs += 1

            box = map[
                cur_p[0] - bs : cur_p[0] + bs + 1, cur_p[1] - bs : cur_p[1] + bs + 1
            ]
            box_max = np.array(np.unravel_index(box.argmax(), box.shape))

            cur_p = (cur_p + box_max - bs * np.ones((1, 2), dtype=int)).reshape(
                2,
            )
    else:
        box_max = np.array(np.unravel_index(box.argmax(), box.shape))
        cur_p = (cur_p + box_max - bs * np.ones((1, 2), dtype=int)).reshape(
            2,
        )

    return cur_p


def get_list(
    map: np.ndarray, box_size: int, num_epochs: int, proj: CartesianProj
) -> np.ndarray:
    """Get the list of sky field centers for mosaic scanning the localization skymap.

    Parameters
    ----------
    map : np.ndarray
        skymap reprojected to the Cartesian surface
    box_size : int
        size of the box containing the neighbors of the given pixel
    num_epochs : int
        number of sky fields to scan
    proj : CartesianProj
        Cartesian projector of the `map`

    Returns
    -------
    result: np.ndarray
        a matrix containing (ra, dec) rows of the sky field centers
    """
    # Store XY coordinates of sky fields
    res = []

    # initial pixel
    max_p = np.array(np.unravel_index(map.argmax(), map.shape))
    res.append(max_p)

    map[max_p[0]][max_p[1]] = 0
    cur_p = max_p
    bs = int((box_size - 1) / 2)

    for _ in range(num_epochs):
        try:
            cur_p = get_cur_p(map, cur_p, bs)
            res.append(cur_p)
            map[cur_p[0]][cur_p[1]] = 0
        except ValueError:
            break

    res = np.array((res))
    temp = proj.ij2xy(res[:, 0], res[:, 1])
    result = proj.xy2ang(temp, lonlat=True).T

    return result


def get_map(
    hpx: np.ndarray, fov: tuple[float, float]
) -> tuple[np.ndarray, CartesianProj]:
    """Get reprojected HEALPix skymap.

    Parameters
    ----------
    hpx : ndarray
        a probability column of the HEALPix skymap
    fov : tuple[float, float]
        a rectangular field of view in degrees

    Returns
    -------
    map: ndarray, CartesianProj
        projected skymap and Cartesian projector for this skymap
    """

    nside = hp.pixelfunc.get_nside(hpx)

    xsize = 360 / fov[0]
    ysize = 180 / fov[1]

    proj = CartesianProj(hpx, xsize=xsize, ysize=ysize)
    map = proj.projmap(map=hpx, vec2pix_func=partial(hp.vec2pix, nside))

    return map, proj


def mosaic_walker(
    hpx: np.ndarray,
    fov: FOV,
    effective_field_area: float = 0.9,
    iter_count: int | None = None,
    box_size: int = 3,
    prob: float = lvk_uncert_level.value,
) -> list[Field]:
    """
    Create the path for optimal coverage the localization region using the telescope
    with a given field of view (FOV). The localization region is split into tiles of
    size = FOV (i.e. like a mosaic). The algorithm starts from the tile with a highest
    probability to a second highest and so on. It runs several iterations `iter_count`
    until stops, i.e. it walks through `iter_count` tiles. The more area one wants to
    cover, higher `iter_count` should be. In the result, the list of the coordinates of
    the centers of the tiles is returned.

    Parameters
    ----------
    hpx : np.ndarray
        a HEALPix array of constant order, which represents the probability
    fov : FOV
        a field of view of a chosen telescope
    effective_field_area : float
        a percent of the FOV that is actually used (to account for field edge artifacts)
    iter_count : int
        a number of iterations of the algorithm

    Returns
    -------
    Field
        a list of the sky field sorted in optimal order for mosaic scanning
    """
    fovx, fovy = fov.width.to_value("deg"), fov.height.to_value("deg")
    c = find_greedy_credible_levels(hpx)
    good_hpx = hpx.copy()
    good_hpx[c > prob] = 0.0
    arr, proj = get_map(
        good_hpx, [fovx * effective_field_area, fovy * effective_field_area]
    )

    if not iter_count:
        npix = hp.get_map_size(good_hpx)
        iter_count = npix
    result = get_list(map=arr, box_size=box_size, num_epochs=iter_count, proj=proj)

    ra = result[:, 0]
    dec = result[:, 1]

    sorted_targets = [
        Field(
            SkyCoord(a * u.deg, d * u.deg),
            width=fov.width,
            height=fov.height,
            name=(
                f"Field_{fovx:.2f}dx{fovy:.2f}d_"
                f"{coord_to_target_name(SkyCoord(a*u.deg, d*u.deg))}"
            ),
        )
        for i, (a, d) in enumerate(zip(ra, dec))
    ]

    return sorted_targets
