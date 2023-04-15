# -----------------------------------------------------------------------------
# Project:     AWARE
# Name:        aware/visualization.py
# Purpose:     Visualization and plotting
# Author:      N. Pankov (colinsergesen@gmail.com)
# Created:     2023-02-01
# Copyright:   (c) 2004-2023 AWARE Developers
# -----------------------------------------------------------------------------
"""
Module :mod:`aware.visualization` contains classes and functions for plotting 
objects on localization coverage map
"""
from __future__ import annotations

import gc
import os
from typing import Any, Callable, Literal, Optional, TypedDict
from contextlib import contextmanager, AbstractContextManager


try:
    from mocpy import MOC, WCS
except ImportError:
    from mocpy import MOC, World2ScreenMPL as WCS

import matplotlib as mpl
import numpy as np
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import Angle, BaseCoordinateFrame, SkyCoord
from astropy.visualization.wcsaxes import WCSAxes
from astropy.wcs import PRJ_CODES
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from ..glade import GladeGalaxy
from ..logger import log

# Select the not-interactive backend for plots
# mpl.use("agg")

# Set Serif text font, which is prettier for human reading
plt.rcParams["font.family"] = "Serif"


__all__ = ["plot_glade_galaxies", "plot_moc"]


TOTAL_SKY_AREA = 4 * np.pi * u.sr
MAX_DEPTH = 9
DELTA_DEPTH = 1
COORD_FRAME = "icrs"
PROJECTION = "MOL"

import astropy.units as u
import numpy as np
from astropy import coordinates, wcs
from matplotlib.pyplot import figure


class World2ScreenMPL(AbstractContextManager):
    """
    Create a World2ScreenMPL for vizualizing a MOC in a matplotlib axis.
    Code copied from `mocpy` and adapted for rectangular FOV.
    """

    def __init__(
        self,
        fig: Figure,
        fov: u.Unit | Angle = Angle([325, 160], "deg"),
        center=coordinates.SkyCoord(180, 0, unit="deg", frame="icrs"),
        coordsys="icrs",
        projection="MOL",
        rotation=coordinates.Angle(0, u.radian),
    ):
        self.wcs = wcs.WCS(naxis=2)

        width_px, height_px = fig.get_size_inches() * float(fig.dpi)

        cdelt_x = fov[0].to_value("deg") / float(width_px)
        cdelt_y = fov[1].to_value("deg") / float(height_px)

        self.wcs.wcs.crpix = [width_px / 2.0, height_px / 2.0]
        self.wcs.wcs.cdelt = [-cdelt_x, cdelt_y]

        if coordsys == "icrs":
            self.wcs.wcs.crval = [center.icrs.ra.deg, center.icrs.dec.deg]
            self.wcs.wcs.ctype = ["RA---" + projection, "DEC--" + projection]
        elif coordsys == "galactic":
            self.wcs.wcs.crval = [center.galactic.l.deg, center.galactic.b.deg]
            self.wcs.wcs.ctype = ["GLON-" + projection, "GLAT-" + projection]

        theta = rotation.radian
        self.wcs.wcs.pc = [
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ]

    def __enter__(self):
        return self.wcs

    def __exit__(self, __exc_type, __exc_value, __traceback):
        pass


def close_all_figures(func: Callable[[Any], Any]):
    """
    Wrapper function, which closes all opened plt.Figure figures after plotting to
    free the memory. Also calls garbage collector to avoid possible memory leaks in
    matplotlib.
    """

    def inner_func(*args: Any, **kwargs: Any):
        plt.close("all")
        gc.collect()
        res = func(*args, **kwargs)
        return res

    return inner_func


@close_all_figures
def plot_moc(
    moc: MOC,
    fov: u.Unit | Angle = Angle([325, 162], "deg"),
    center: SkyCoord = SkyCoord(180 * u.deg, 0 * u.deg, frame="icrs"),
    coordsys: Literal["icrs", "galactic"] = "icrs",
    rotation: u.Unit | Angle = Angle(0, u.degree),
    projection: Literal(PRJ_CODES) = "MOL",
    border_kws: dict[str, Any] | None = None,
    ax: plt.Axes | WCSAxes | None = None,
    fig_kws: Optional[dict[str, Any]] = None,
) -> plt.Axes | WCSAxes:
    """
    Plot the Multi-Coverage Map (MOC), e.g. the sky map of a LIGO/Virgo/KAGRO
    gravitational wave signal.

    Parameters
    ----------
    moc : MOC
        a multi-coverage map, should have at least one column with UNIQ index
    fov : u.Unit | Angle, optional
        a size WxH of the field of view , by default Angle(360, u.degree)
    center : SkyCoord, optional
        the center of the field, by default SkyCoord(180*u.deg, 0*u.deg)
    coordsys : Literal["icrs", "galactic"], optional
       A coordinate system, must be in ["icrs", "galactic"], by default "icrs"
    rotation : u.Unit | Angle, optional
        a field rotation, by default Angle(0, u.degree)
    projection : Literal(PRJ_CODES), optional
        the world base to image base projection type (see `~astropy.wcs.PRJ_CODES` for
        the list of available projections), by default "MOL"
    border_kws : dict[str, Any] | None, optional
        the keyword arguments for the MOC border formatting, by default None
    ax : Optional[plt.Axes], optional
        the keyword arguments for the axes formatting, by default None
    fig_kws : Optional[dict[str, Any]], optional
        the keyword arguments for the figure formatting, by default None

    Returns
    -------
    ax : plt.Axes
        axes with MOC plot and a possible formatting
    """

    if not ax:
        fig = plt.figure(**fig_kws if fig_kws else {})
        with World2ScreenMPL(
            fig,
            fov=fov,
            center=center,
            coordsys=coordsys,
            rotation=rotation,
            projection=projection,
        ) as wcs:
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs)
            moc.border(ax=ax, wcs=wcs, **border_kws if border_kws else {})
    else:
        moc.border(ax=ax, wcs=ax.wcs, **border_kws if border_kws else {})

    ax.set_xlabel("R.A., Hours")
    ax.set_ylabel("Dec, Degrees")

    return ax


@close_all_figures
def plot_glade_galaxies(
    galaxies: list[GladeGalaxy],
    ax: plt.Axes | WCSAxes | None = None,
    scatter_kws: Optional[dict[str, Any]] = None,
    colorbar_kws: Optional[dict[str, Any]] = None,
    fig_kws: Optional[dict[str, Any]] = None,
) -> tuple[plt.Axes | WCSAxes, plt.Axes]:
    """
    Plot

    Parameters
    ----------
    galaxies : list[GladeGalaxy]
        _description_
    ax : plt.Axes | WCSAxes | None, optional
        _description_, by default None
    scatter_kws : Optional[dict[str, Any]], optional
        _description_, by default None
    colorbar_kws : Optional[dict[str, Any]], optional
        _description_, by default None
    fig_kws : Optional[dict[str, Any]], optional
        _description_, by default None

    Returns
    -------
    tuple[plt.Axes | WCSAxes, plt.Axes]
        _description_
    """

    ra = np.asarray([o.ra.deg for o in galaxies], dtype=np.float32)
    dec = np.asarray([o.dec.deg for o in galaxies], dtype=np.float32)
    z = np.asarray([getattr(o, "z", np.inf) for o in galaxies], dtype=np.float16)

    if not ax:
        fig, ax = plt.subplots(**fig_kws, projection="MOL")

        scatter_kws = (
            {"label": "GLADE+ galaxies", "cmap": "viridis_r", "edgecolor": "k"}
            if not scatter_kws
            else scatter_kws
        )
        transform = ax.get_transform()
    else:
        fig = ax.get_figure()
        transform = ax.get_transform("world")

    # TODO: Add correct RA ticklabels in hour format
    ax.scatter(ra, dec, transform=transform, **scatter_kws)
    ax.set_xlabel("RA, Degrees")
    ax.set_ylabel("Dec, Degrees")

    # # Adapted from here https://stackoverflow.com/a/69296746
    # norm = plt.Normalize(np.min(z), np.max(z))
    # smap = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    # colorbar_kws = dict(fraction=0.1, shrink=0.8) if not colorbar_kws else colorbar_kws
    # cbar = fig.colorbar(smap, ax=ax, **colorbar_kws)

    return ax
