"""
Author: Artem Prokhorenko (prokhorenko.aiu@phystech.edu)
        Evgenii Shekotihin (z.shekotihin.2012@yandex.ru)
        Nicolai Pankov (colinsergesen@gmail.com)
mosaic.py (c) 2023
Desc: mosaic scanning of the localization skymap
Created:  2023-04-18
Modified: 2024-03-13
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

def find_next_p(
    cur_p: tuple[int, int], 
    cur_map: np.ndarray, 
    orig_map: np.ndarray, 
    visited: np.ndarray) -> tuple[int, int]:
    """
    Given current problem pixel, current skymap 
    and original map find next pixel to exposure
    """
    # box indices
    box_ps = [[0, 0], [0, 1], [0, 2],
              [1, 0], [1, 1], [1, 2],
              [2, 0], [2, 1], [2, 2]]

    # unexposured box pixels 
    box_ps_f = [] 
    
    for p in box_ps:
        # original pixel
        op = get_p(cur_p, p, 3)
        # in visited grid
        op_visited = visited[op[0]][op[1]]
        
        # no need to exposure visited pixels 
        if op_visited == 1:
            pass
        elif p == [1, 1]:  # and current pixel 
            pass
        else:
            box_ps_f.append(p)  
    
    # from filtered pixels choose one with the biggest sum 
    # calculated for original skymap box, centered in that pixel 
    sums = []
    
    for p in box_ps_f:
        # probe pixel
        pp = get_p(p, cur_p, 3)
        # sum of a 7x7 box around it
        sums.append(get_box(orig_map, pp, 7).sum())
        
    box_max = box_ps_f[sums.index(max(sums))]
    next_p = get_p(cur_p, box_max, 3)
    
    return next_p

def get_box(
    data: np.ndarray, cur_p: tuple[int, int], box_size=3) -> np.ndarray:
    """Data box slice centered at current pixel 
    """
    bs = box_size
    return data[cur_p[0] - bs : cur_p[0] + bs + 1, cur_p[1] - bs : cur_p[1] + bs + 1]

def get_p(
    cur_p: tuple[float, float], box_p: tuple[float, float], box_size=3) -> tuple[float, float]:
    """
    Square odd box located at cur_p e.g. [121, 120],
    given box pixel index e.g. [0, 1] return pixel
    in current pixel system: 
    get_p([121, 120], [1, 2], 5) --> [120, 119]
    
    """
    if not isinstance(cur_p, np.ndarray):
        cur_p = np.array(cur_p)

    if not isinstance(box_p, np.ndarray):
        box_p = np.array(box_p)

    return (cur_p - box_size // 2  + box_p).reshape(2, )

def get_next_p(
    cur_map: np.ndarray, 
    orig_map: np.ndarray,
    visited: np.ndarray, 
    cur_p: tuple[float, float], 
    bs: int, iter: int, enlarge_box=False) -> tuple[float, float]:
    """
    Get current position of the mosaic walk in Cartesian coordinates.

    Parameters
    ----------
    cur_map : np.ndarray
        current skymap (with exposured pixels set to 0) 
    orig_map: np.ndarray
        original skymap reprojected to the Cartesian surface
    visited: np.ndarray
        array with exposured pixels set to 1
    cur_p : tuple[float, float]
        current pixel (Y, X)
    bs : int
        a size of the box containing neighbors of the (Y, X)
        
    Returns
    -------
    tuple[float, float]
        next pixel (Y, X)
    """
    box = get_box(cur_map, cur_p, bs)

    if box.max() == 0:
        # variant with enlarging box
        if enlarge_box:
            while box.max() == 0:
                print(f'Iter: {iter}, enlarging box size')
                
                bs += 1
                
                box = get_box(cur_map, cur_p, bs)
                box_max = np.array(np.unravel_index(box.argmax(), box.shape))
                next_p = get_p(cur_p, box_max, box.shape[0])
                
        # smooth variant 
        else:
            next_p = find_next_p(cur_p, cur_map, orig_map, visited)
    else:
        box_max = np.array(np.unravel_index(box.argmax(), box.shape))
        next_p = get_p(cur_p, box_max, box.shape[0])
    
    return next_p

def get_list(
    map: np.ndarray, box_size: int, num_epochs: int, proj: CartesianProj) -> np.ndarray:
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
    # grid of exposured cells
    visited = np.zeros(map.shape)
    
    # original copy
    orig_map = map.copy()
    
    # Store YX coordinates of sky fields
    res = []
    
    # initial pixel
    max_p = np.array(np.unravel_index(map.argmax(), map.shape))
    res.append(max_p)

    cur_p = max_p  # current pixel
    bs = int((box_size - 1) // 2)
    
    map[max_p[0]][max_p[1]] = 0
    visited[max_p[0]][max_p[1]] = 1
    
    for iter in range(1, num_epochs):
        try:
            new_p = get_next_p(map, orig_map, visited, cur_p, bs, iter)
        except ValueError: 
            # try with enlarged box for failures that may occur
            new_p = get_next_p(map, orig_map, visited, cur_p, bs, iter, enlarge_box=True)
            
        res.append(new_p)
        map[new_p[0]][new_p[1]] = 0
        visited[new_p[0]][new_p[1]] = 1
        cur_p = new_p
            
    res = np.array((res))
    temp = proj.ij2xy(res[:, 0], res[:, 1])
    result = proj.xy2ang(temp, lonlat=True).T

    # return result     # return for mosaic_walker
    return result, res  # return for validation notebook

def get_map(
    hpx: np.ndarray, fov: tuple[float, float]) -> tuple[np.ndarray, CartesianProj]:
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
    prob: float = lvk_uncert_level.value,) -> list[Field]:
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
