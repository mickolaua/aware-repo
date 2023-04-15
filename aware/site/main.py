"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
site.py (c) 2023
Desc: Observer sites
Created:  2023-03-03
Modified: 2023-03-13
"""
from __future__ import annotations

import warnings
from datetime import datetime, timedelta
from functools import cached_property
from typing import Any, Optional, Sequence

import numpy as np
from astroplan.observer import Observer
from astroplan.plots import plot_airmass
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import Angle

from matplotlib import pyplot as plt
from matplotlib.axes import Axes

from .. import planner
from ..config import CfgOption

__all__ = ["FOV", "Site", "Telescopes", "default_sites"]


sort_method = CfgOption("sort_method", "radec", str)
max_airmass = CfgOption("max_airmass", 3.5, float)

# Registry for all defined sites. Note, sites are theyselves added to the registry,
# a person does not need to add them manually.
Telescopes: dict[str, Site] = {}


class FOV:
    """
    Helper class representing telescope field of view.

    Attributes
    ----------
    width : float
        an angular width of the field
    height : float
        an angular height of the field
    radius : float
        an equivalent radius of the circle that contains the width x height FOV
    is_widefield : bool
        is the telescope has a wide FOV or not
    """
    def __init__(self, width: Angle, height: Angle):
        self.width = width
        self.height = height

    @cached_property
    def radius(self) -> Angle:
        return np.hypot(0.5*self.width, 0.5*self.height)
    
    @cached_property
    def is_widefield(self) -> bool:
        return (
            self.width.to_value("arcmin") >= 40 
            or self.height.to_value("arcmin") >= 40
        )


class Site(Observer):
    def __init__(
        self,
        name: str,
        latitude: u.Unit,
        longitude: u.Unit,
        elevation: u.Unit,
        timezone: str,
        fov: FOV,
        default_exposure: u.Unit,
        default_magnitude: u.Unit,
        horizon: u.Unit = 0.0 * u.degree,
        full_name: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(
            None, timezone, name, latitude, longitude, elevation, None, None, None, None
        )
        self.default_exposure = default_exposure
        self.default_magnitudee = default_magnitude
        self.horizon = horizon
        if full_name is not None:
            self.full_name = full_name
        else:
            self.full_name = self.name

        self.fov = fov
        for k, v in kwargs.items():
            try:
                setattr(self, k, v)
            except (LookupError, AttributeError, ArithmeticError):
                pass

        # Add the site to the registry
        global Telescopes
        Telescopes[self.name] = self

    def nearest_observation_window(
        self,
        time: Time | datetime | str | float,
        after_sunset: TimeDelta | timedelta | u.Unit = TimeDelta(300 * u.s),
        before_sunrise: TimeDelta | timedelta | u.Unit = TimeDelta(300 * u.s),
    ) -> tuple[Time, Time]:
        return planner.nearest_observation_window(
            self, time, self.horizon, after_sunset, before_sunrise
        )

    def observable_targets(
        self,
        targets: Sequence[FixedTarget],
        start_time: Optional[Time] = None,
        end_time: Optional[Time] = None,
        max_airmass: Optional[float] = max_airmass.get_value(),
    ) -> list[FixedTarget]:

        if not len(targets):
            return []

        if start_time is None and end_time is not None:
            time_grid = [end_time]
        elif start_time is not None and end_time is None:
            time_grid = [start_time]
        else:
            time_grid = Time(np.linspace(start_time.jd, end_time.jd, 150), format="jd")
        above_horizon, altaz = self.target_is_up(
            time_grid, targets, return_altaz=True, grid_times_targets=True
        )
        airmass = altaz.altaz.secz
        observable = np.argwhere(
            np.any((airmass < max_airmass) & (airmass > 0), axis=1)
            & np.any(above_horizon, axis=1)
        ).ravel()
        observable_targets = [targets[i] for i in observable]

        return observable_targets

    def observation_order(
        self,
        targets: Sequence[FixedTarget],
        time: Optional[Time] = None,
        start_time: Optional[Time] = None,
        end_time: Optional[Time] = None,
        method: Optional[str] = sort_method.get_value(),
    ) -> list[FixedTarget]:
        if not len(targets):
            return []

        # TODO: implement new methods and improve existing
        if method == "nearest":
            sorted_targets = planner.nn_sorter(
                self, targets, time, start_time, end_time
            )
        elif method == "radec":
            sorted_targets = planner.radec_sorter(
                self, targets, time, start_time, end_time
            )
        else:
            # Just return original list of targets
            sorted_targets = list(targets)

        return sorted_targets

    def plot_airmass(
        self,
        targets: Sequence[FixedTarget],
        time: Optional[Time] = None,
        start_time: Optional[Time] = None,
        end_time: Optional[Time] = None,
        show_moon: bool = True,
        fig_path: str = "airmass.png",
    ) -> Axes:
        if start_time is None:
            start_time = time - TimeDelta(0.5 * u.d)

        if end_time is None:
            start_time = time + TimeDelta(0.5 * u.d)

        fig, ax = plt.subplots(figsize=(8, 6), dpi=96)

        # Astroplan emits warnings on some problems with line style in graph
        # we dont want to see
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            if time is None:
                time = Time(np.linspace(start_time.jd, end_time.jd, 150), format="jd")
            ax = plot_airmass(
                targets,
                self,
                time.datetime,
                ax=ax,
                min_airmass=1.0,
                max_airmass=3.8,
            )

        if show_moon:
            airmass_time = Time(
                np.linspace(start_time.jd, end_time.jd, 150), format="jd"
            )
            moon_airmass = self.moon_altaz(airmass_time).secz
            good = moon_airmass > 0
            ax.plot(
                airmass_time[good].datetime, moon_airmass[good], "k--", label="The Moon"
            )

        secz_ticks = np.asarray(ax.get_yticks())

        def _airmass_to_alt(airmass):
            return 90 - np.rad2deg(np.arccos(1 / airmass))

        alt = _airmass_to_alt(secz_ticks)
        ax2 = ax.secondary_yaxis("right")
        ax2.set_yticks(secz_ticks)
        ax2.set_yticklabels([f"{t:.01f}" for t in alt])
        ax2.set_ylabel("Altitude, deg")

        # Add offsets of 30 minutes from each side to fit the whole
        # observational window in the graph
        llim = (start_time - TimeDelta(30 * u.min)).datetime
        rlim = (end_time + TimeDelta(30 * u.min)).datetime
        ax.set_xlim(llim, rlim)
        ax.axvline(start_time.datetime, ls="-.", label="Obs. window")
        ax.axvline(end_time.datetime, ls="-.")
        ax.grid(True, lw=1, alpha=0.3)
        ax.tick_params(axis="both", which="both", top=True, right=True, direction="in")
        ax.tick_params(axis="x", rotation=0)
        ax.set_title(
            "Site {:s}\nLongitude: {:.5f} deg, Latitude: {:.5f} deg".format(
                self.full_name,
                self.location.lon.value,
                self.location.lat.value,
            ),
            loc="left",
        )

        lgd = ax.legend(
            bbox_to_anchor=(1.15, 1),
            loc="upper left",
            borderaxespad=0,
            frameon=False,
            title="List of Objects",
        )

        # Reposition Axes to fit both Axes and legend into the canvas
        # Borrowed from https://stackoverflow.com/questions/30413789/
        fig.canvas.draw()
        inv_fig = fig.transFigure.inverted()
        lgd_pos = lgd.get_window_extent()
        lgd_coord = inv_fig.transform(lgd_pos)
        lgd_xmax = lgd_coord[1, 0]

        ax_pos = ax.get_window_extent()
        ax_coord = inv_fig.transform(ax_pos)
        ax_xmax = ax_coord[1, 0]

        # Multiplied by 0.7 here, because of the large shift between legend
        # rigiht side and right side of canvas
        shift = 1 - (lgd_xmax - ax_xmax) * 0.1
        fig.tight_layout(rect=(0, 0, shift, 1))

        if fig_path:
            fig.savefig(fig_path)

        return ax


# Define observatories here
Mondy_AZT33IK = Site(
    timezone="Asia/Irkutsk",
    name="mondy_azt33ik",
    full_name="AZT-33IK at Mondy, Irkutsk",
    latitude=51.619822 * u.deg,
    longitude=100.9214 * u.deg,
    elevation=2000 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(15 * u.arcmin, 15 * u.arcmin),
    aperture=1.5 * u.m,
    default_exposure=120 * u.s,
    default_magnitude=21.0 * u.ABmag
)
CrAO_ZTSH = Site(
    timezone="Europe/Moscow",
    name="crao_ztsh",
    full_name="ZTSh at Nauchny/Simeiz, Crimea",
    latitude=44.726667 * u.deg,
    longitude=34.015861 * u.deg,
    elevation=600 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(15 * u.arcmin, 15 * u.arcmin),
    aperture=2.65 * u.m,
    default_exposure=120 * u.s,
    default_magnitude=21.0 * u.ABmag
)
Assy_AZT20 = Site(
    timezone="Asia/Almaty",
    name="assy_azt20",
    full_name="AZT-20 at Assy-Turgen observatory, Almaty",
    latitude=43.2252778 * u.deg,
    longitude=77.8716667 * u.deg,
    elevation=2750 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(13 * u.arcmin, 13 * u.arcmin),
    aperture=1.5 * u.m,
    default_exposure=30 * u.s,
    default_magnitude=19.0 * u.ABmag
)
TSHAO_Zeiss1000 = Site(
    timezone="Asia/Almaty",
    name="tshao_zeiss1000",
    full_name="Zeiss-1000 at Tien-Shan observatory, Almaty",
    latitude=43.06 * u.deg,
    longitude=76.97 * u.deg,
    elevation=2735 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(13 * u.arcmin, 13 * u.arcmin),
    aperture=1.0 * u.m,
    default_exposure=30 * u.s,
    default_magnitude=18.0 * u.ABmag
)
Altai_Santel400 = Site(
    timezone="Asia/Barnaul",
    name="altai_santel400",
    full_name="Santel-400 at Multa observatory, Altai",
    latitude=53.43517913184 * u.deg,
    longitude=83.95733665302 * u.deg,
    elevation=600 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(82 * u.arcmin, 82 * u.arcmin),
    aperture=0.4 * u.m,
    default_exposure=120 * u.s,
    default_magnitude=16.0 * u.ABmag
)
Abastumani_AS32 = Site(
    timezone="Asia/Yerevan",
    name="abao_as32",
    full_name="AS-32 at Georgian National Astrophysical Observatory, Abastumani",
    latitude=41.7542 * u.deg,
    longitude=42.8194 * u.deg,
    elevation=1650 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(44 * u.arcmin, 44 * u.arcmin),
    aperture=0.7 * u.m,
    default_exposure=120 * u.s,
    default_magnitude=17.0 * u.ABmag
)

default_sites = CfgOption("default_sites", list(Telescopes.keys()), list)
