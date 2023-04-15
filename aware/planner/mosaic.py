from __future__ import annotations

import warnings

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import SkyCoord
from healpy.visufunc import cartview

from ..site import FOV

__all__ = ["mosaic_walker"]


def _distance_matrix(
    arr: np.ndarray, y_n: int, x_n: int, scale_x: int, scale_y: int
) -> np.ndarray:
    """Helper function for constructing the distance matrix. It is not supposed for
    use directly."""
    cur_arr = np.zeros_like(arr)
    for j in range(scale_x):
        for i in range(scale_y):
            abs_x = abs(x_n-j)
            dist_x = [abs_x if (abs_x<=(scale_x/2)) else (scale_x-abs_x)][0]
            dist_y = abs(y_n-i)
            cur_arr[i][j] = (10**(5) * arr[i][j]) - (10**(-2) * (dist_y + dist_x))

    return cur_arr


def _transition_xy2ra(
    x: np.ndarray, y: np.ndarray, scale_x: int, scale_y: int
) -> tuple[np.ndarray, np.ndarray]:
    """Convert x, y coordinates to RA, Dec"""
    ra, dec = x.copy(), y.copy()
    for i, (xx, yy) in enumerate(zip(x, y)):
        cur_y = xx + 0.5
        cur_x = yy + 0.5
        d_ang_y = 180
        d_ang_x = 360
        dec[i] = (d_ang_y / scale_y * cur_y) - 90
        ra[i] = ((d_ang_x / scale_x * (scale_x - cur_x)) + 180) % 360

    return ra, dec


def mosaic_walker(
    hpx: np.ndarray, fov: FOV, effective_field_area: float = 0.9, iter_count: int = 300
) -> list[FixedTarget]:
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
    SkyCoord
        a list of the coordinates of the center of the tiles ready for telescope
        observation.
    """
    frame_x = fov.width.to_value("deg") * effective_field_area
    frame_y = fov.height.to_value("deg") * effective_field_area

    # Construct sky coordinates grid: RA: from 0 to 360 = 360 deg,
    # DEC: from -90 to +90 = 180 deg
    scale_x = 360
    scale_y = 180
    scale_x /= frame_x
    scale_y /= frame_y

    def int_r(num):
        num = int(num + (0.5 if num > 0 else -0.5))
        return num

    scale_x = int_r(scale_x)
    scale_y = int_r(scale_y)

    # Here, we only interested in returned projected map
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "get_cmap function was deprecated")
        full_map = cartview(
            hpx,
            coord=["G"],
            unit="cbar label",
            xsize=scale_x,
            ysize=scale_y,
            return_projected_map=True,
            nest=True,
        )

    # Modify copy, not original localization map
    map_copy = ma.getdata(full_map).copy()

    # Start at these indicies
    arr = map_copy.copy()
    # y_n, x_n = np.unravel_index(np.argmax(arr), arr.shape)

    # arr[y_n][x_n] = 0
    # x, y = [x_n], [y_n]
    # for _ in range(iter_count - 1):
    #     cur_arr = _distance_matrix(arr, y_n, x_n, scale_x, scale_y)
    #     y_n, x_n = np.unravel_index(np.argmax(cur_arr), cur_arr.shape)
    #     arr[y_n][x_n] = 0
    #     x.append(x_n)
    #     y.append(y_n)

    # x = np.asarray(x)
    # y = np.asarray(y)

    # print(x)

    # # TODO:Convert projected X,Y to RA, Dec
    # ra, dec = _transition_xy2ra(x, y, scale_x, scale_y)
    # print(ra)
    # print(dec)

    def cur_matrix(arr, y_n, x_n):
        cur_arr = np.zeros_like(arr)
        for j in range(scale_x):
            for i in range(scale_y):
                abs_x = abs(x_n-j)
                dist_x = [abs_x if (abs_x<=(scale_x/2)) else (scale_x-abs_x)][0]
                dist_y = abs(y_n-i)
                cur_arr[i][j] = (10**(5) * arr[i][j]) - (10**(-2) * (dist_y + dist_x))
        return cur_arr
    
    pred_xy = []
    def simple(arr, i):
        y_n, x_n = np.unravel_index(np.argmax(arr), arr.shape)
        val_n = arr[y_n][x_n]
        report = [y_n, x_n, val_n]
        arr[y_n][x_n] = 0
        pred_xy.append(report)
        for z in range(i-1):
            cur_arr = cur_matrix(arr, y_n, x_n) 
            y_n, x_n = np.unravel_index(np.argmax(cur_arr), cur_arr.shape)
            val_n = arr[y_n][x_n]
            report = [y_n, x_n, val_n]
            arr[y_n][x_n] = 0
            pred_xy.append(report)

        return pred_xy
        
    pred_xy = simple(arr, iter_count)

    def transition_xy(pred_xy):
        for i in range(len(pred_xy)):
            cur_y = pred_xy[i][0] + 0.5
            cur_x = pred_xy[i][1] + 0.5
            d_ang_y = 180
            d_ang_x = 360
            pred_xy[i][0] = round(((d_ang_y/scale_y)*cur_y) - 90, 3)
            pred_xy[i][1] = round((((d_ang_x/(scale_x))*(scale_x-cur_x)) + 180)%360, 3)

        return pred_xy
    
    pred_ang = transition_xy(pred_xy)
    coord = SkyCoord(
        [ra for (_, ra, _) in pred_ang]*u.deg, 
        [dec for (dec, _, _) in pred_ang]*u.deg
    )

    sorted_targets = [
        FixedTarget(c, name=f"field_{i+1}") for i, c in enumerate(coord)
    ]

    return sorted_targets
