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
from contextlib import AbstractContextManager, contextmanager
from typing import Any, Callable, Literal, Optional, TypedDict

from matplotlib.axes import Axes

from aware.field import Field

try:
    from mocpy import MOC, WCS
except ImportError:
    from mocpy import MOC, World2ScreenMPL as WCS

import astropy.units as u
import matplotlib as mpl
import numpy as np
from astroplan.target import FixedTarget
from astropy import coordinates
from astropy import units as u
from astropy import wcs
from astropy.coordinates import Angle, BaseCoordinateFrame, SkyCoord
from astropy.visualization.wcsaxes import WCSAxes
from astropy.wcs import PRJ_CODES
from ligo.skymap.plot import cylon
from ligo.skymap.postprocess.util import find_greedy_credible_levels
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from aware.config import CfgOption

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


plot_width_inch = CfgOption("plot_width_px", 9.600, float)
plot_height_inch = CfgOption("plot_height_px", 9.600, float)
plot_dpi = CfgOption("plot_dpi", 100, float)
plot_colormap = CfgOption("plot_colormap", "cylon", str)

plt.rcParams["figure.figsize"] = (plot_width_inch.value, plot_height_inch.value)
plt.rcParams["figure.dpi"] = plot_dpi.value
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


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
    projection: Literal[PRJ_CODES] = "MOL",
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
def plot_skymap(
    hpx: np.ndarray,
    center: SkyCoord | None = None,
    radius: u.Unit | None = None,
    prob_threshold: float = 0.9,
) -> plt.Axes:
    ax_kwargs = {}
    if radius is None:
        radius = 180 * u.deg

    if center is None:
        center = SkyCoord(0.0 * u.deg, 0.0 * u.deg)

    if radius < 15 * u.deg:
        ax_kwargs["radius"] = 3 * radius
        ax_kwargs["center"] = center
        ax_kwargs["projection"] = "astro hours zoom"
    else:
        ax_kwargs["projection"] = "astro hours mollweide"

    c = find_greedy_credible_levels(hpx)
    fig = plt.figure()
    ax = plt.axes([0.05, 0.05, 0.9, 0.9], **ax_kwargs)
    fig.add_axes(ax)
    fig.tight_layout()
    cs = ax.contour_hpx(
        (100 * c, "ICRS"),
        nested=True,
        colors="k",
        linewidths=2,
        levels=[prob_threshold * 100],
    )
    fmt = r"%g\%%" if plt.rcParams["text.usetex"] else "%g%%"
    ax.clabel(cs, cs.levels, fmt=fmt, fontsize=8, inline=True)
    ax.grid(True)
    im = ax.imshow_hpx((hpx, "icrs"), nested=True, cmap="cylon")
    ax.set_xlabel("R.A.")
    ax.set_ylabel("Dec")
    ax.invert_xaxis()

    norm = plt.Normalize(np.min(hpx), np.max(hpx))
    smap = plt.cm.ScalarMappable(cmap="cylon", norm=norm)
    colorbar_kws = dict(fraction=0.1, shrink=0.5, orientation="horizontal", pad=0.05)
    cbar = fig.colorbar(smap, ax=ax, **colorbar_kws)
    cbar.set_label(r"prob. per deg$^2$")

    fig.tight_layout(pad=0.05)

    return ax


@close_all_figures
def plot_fields_on_healpix(
    hpx: np.ndarray,
    fields: list[Field],
    center: SkyCoord,
    radius: u.Unit,
    prob_threshold: float = 0.9,
) -> Axes:
    ax: WCSAxes = plot_skymap(hpx, center, radius, prob_threshold)
    fig = ax.get_figure()
    dpi = fig.get_dpi()

    for i, f in enumerate(fields):
        x0 = f.ra.deg
        y0 = f.dec.deg
        w = f.width.to_value(u.deg)
        h = f.height.to_value(u.deg)
        r = patches.Rectangle(
            xy=(x0 - w / 2, y0 - h / 2),
            width=0.9 * w,
            height=0.9 * h,
            transform=ax.get_transform("icrs"),
            facecolor="none",
            edgecolor="b",
            alpha=1.0,
            lw=0.5,
            zorder=99,
        )
        #
        # fontsize in points = 2 pix2pts * w / cdelt / dpi * inch_size
        #
        pixscale = abs(ax.wcs.pixel_scale_matrix[0, 0])
        size_inches = min(fig.get_size_inches())
        fontsize_px = 0.3 * h / pixscale
        pix2pt = 1.333
        fontsize_pts = max(3.0, fontsize_px / pix2pt)
        ax.text(
            x0,
            y0,
            str(i + 1),
            transform=ax.get_transform("icrs"),
            va="center",
            ha="center",
            size=fontsize_pts,
            color="b",
        )
        ax.add_artist(r)

    # x = [f.ra.deg for f in fields]
    # y = [f.dec.deg for f in fields]
    # U = np.diff(x)
    # V = np.diff(y)
    # pos_x = x[:-1] + U / 2
    # pos_y = y[:-1] + V / 2
    # norm = np.sqrt(U**2 + V**2)
    # ax.plot(x, y, "b-", transform=ax.get_transform("icrs"), lw=0.5, alpha=0.5)
    # ax.quiver(
    #     pos_x,
    #     pos_y,
    #     U / norm,
    #     V / norm,
    #     angles="xy",
    #     zorder=5,
    #     pivot="mid",
    #     transform=ax.get_transform("icrs"),
    # )

    return ax


@close_all_figures
def plot_targets_on_healpix(
    hpx: np.ndarray,
    targets: list[GladeGalaxy],
    center: SkyCoord,
    radius: u.Unit,
    prob_threshold: float = 0.9,
) -> Axes:
    ax = plot_skymap(hpx, center, radius, prob_threshold)
    fig = ax.get_figure()

    ra = np.asarray([t.coord.ra.deg for t in targets])
    dec = np.asarray([t.coord.dec.deg for t in targets])
    D_L = np.asarray([getattr(o, "D_L", None) for o in targets], dtype=np.float64)

    ax.scatter(
        ra,
        dec,
        ls="",
        marker="o",
        transform=ax.get_transform("icrs"),
        edgecolor="black",
        alpha=0.8,
        c=D_L,
        cmap="winter"
    )
    norm = plt.Normalize(np.min(D_L), np.max(D_L))
    smap = plt.cm.ScalarMappable(cmap="winter", norm=norm)
    colorbar_kws = dict(fraction=0.1, shrink=0.5, orientation="horizontal", pad=0.05)
    cbar = fig.colorbar(smap, ax=ax, **colorbar_kws)
    cbar.set_label(r"Luminosity distance $D_L$, Mpc")

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
