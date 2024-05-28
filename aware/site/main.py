"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
site.py (c) 2023
Desc: Observer sites
Created:  2023-03-03
Modified: 2024-03-31
"""

from __future__ import annotations

import os
import re
import warnings
from datetime import datetime, timedelta
from functools import cached_property
from typing import Any, Optional, Sequence

import numpy as np
import orjson
from adjustText import adjust_text
from astroplan import ObservingBlock, constraints, observability_table, scheduling
from astroplan.observer import Observer
from astroplan.plots import plot_airmass, plot_altitude, plot_schedule_airmass
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import Angle, Latitude, Longitude, SkyCoord
from astropy.time import Time, TimeDelta
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.dates import DateFormatter
from pydantic import BaseModel, ValidationError, field_validator
from glob import glob

from aware.logger import log
from aware.path import get_full_path

from .. import planning
from ..config import CfgOption

__all__ = ["FOV", "Site", "Telescopes", "default_sites"]


sort_method = CfgOption("sort_method", "radec", str)
max_airmass = CfgOption("max_airmass", 3.5, float)
min_altitude = CfgOption("min_altitude", 30, lambda x: Angle(x, unit=u.deg))
max_altitude = CfgOption("max_altitude", 90, lambda x: Angle(x, unit=u.deg))
min_altitude = CfgOption("min_altitude", 30, lambda x: Angle(x, unit=u.deg))
min_moon_separation = CfgOption(
    "min_moon_separation", 50, lambda x: Angle(x, unit=u.deg)
)
max_moon_phase = CfgOption("max_moon_illumination", 1.0, float)
site_plugin_directory = CfgOption(
    "site_plugin_directory", get_full_path("~/aware/plugins"), get_full_path
)
# Registry for all defined sites. Note, sites are theyselves added to the registry,
# a person does not need to add them manually.
Telescopes: dict[str, Site] = {}


# TODO: `is_widefield` property should be assigned to a site instead of FOV
class FOV:
    """
    Helper class representing telescope field of view.

    Attributes
    ----------
    width : Angle
        an angular width of the field
    height : Angle
        an angular height of the field
    radius :Angle
        an equivalent radius of the circle that contains the width x height FOV
    is_widefield : bool
        is the telescope has a wide FOV or not
    """

    def __init__(self, width: Angle, height: Angle):
        self.width = width
        self.height = height

    @cached_property
    def radius(self) -> Angle:
        return np.hypot(0.5 * self.width, 0.5 * self.height)

    @cached_property
    def is_widefield(self) -> bool:
        warnings.warn(
            "This method marked as deprecated and will be accessed via `Site` instance "
            "instead.",
            FutureWarning,
        )
        return (
            self.width.to_value("arcmin") >= 40 or self.height.to_value("arcmin") >= 40
        )

    def __str__(self) -> str:
        return "{:.3f} deg x {:.3f} deg".format(
            self.width.to_value(u.deg), self.height.to_value(u.deg)
        )

    __repr__ = __str__


# TODO: slew rate should be in deg per second
class Site(Observer):

    def __init__(
        self,
        name: str,
        latitude: u.Unit,
        longitude: u.Unit,
        elevation: u.Unit,
        timezone: str,
        fov: FOV,
        default_magnitude: u.Unit,
        default_exposure: u.Unit = 30 * u.s,
        default_exposure_number: int = 3,
        default_filter: str = "",
        default_target_list_fmt: str = "txt",
        default_slew_rate: u.Unit = 3 * u.min,
        horizon: u.Unit = -12.0 * u.degree,  # Nautical twilight
        full_name: Optional[str] = None,
        observer_telegram_id: str = "",
        observer_email: str = "",
        observer_socket: str = "",
        list_formatter_kws: dict | None = None,
        **kwargs: Any,
    ) -> None:
        global Telescopes
        # if name in Telescopes:
        #     raise KeyError(
        #         f"Site plugin `{name}` will not be loaded, since it overwrites "
        #         "existing site"
        #     )

        super().__init__(
            None, timezone, name, latitude, longitude, elevation, None, None, None, None
        )
        self.default_exposure = default_exposure
        self.default_exposure_number = default_exposure_number
        self.default_magnitude = default_magnitude
        self.default_filter = default_filter
        self.default_target_list_fmt = default_target_list_fmt
        self.default_slew_rate = default_slew_rate
        self.horizon = horizon
        if full_name is not None:
            self.full_name = full_name
        else:
            self.full_name = self.name

        self.observer_telegram_id = observer_telegram_id

        self.fov = fov
        for k, v in kwargs.items():
            try:
                setattr(self, k, v)
            except (LookupError, AttributeError, ArithmeticError):
                pass

        if list_formatter_kws:
            self.list_formatter_kws = list_formatter_kws
        else:
            self.list_formatter_kws = {}

        # Add the site to the registry
        Telescopes[self.name] = self

    @classmethod
    def from_json(cls, filename: str):
        class SiteModel(BaseModel):
            name: str
            latitude: float
            longitude: float
            elevation: float
            timezone: str
            fov: tuple[float, float] | dict
            default_magnitude: float
            default_exposure: float = 30
            default_exposure_number: int = 3
            default_filter: str = ""
            default_target_list_fmt: str = "txt"
            default_slew_rate: float = 3.0 * u.min
            horizon: float = 12.0  # Nautical twilight
            full_name: Optional[str] = None
            observer_telegram_id: str = ""
            observer_email: str = ""
            observer_socket: str = ""

            @field_validator("latitude")
            def check_lat(cls, v: Any):
                if isinstance(v, str):
                    if u.deg in u.Unit(v).find_equivalent_units():
                        lat = Latitude(v).to_value(u.deg)
                    else:
                        raise u.UnitsError(
                            "Latitude must have an equivalent unit to degrees"
                        )
                else:
                    lat = float(v)
                return lat

            @field_validator("longitude")
            def check_lon(cls, v: Any):
                if isinstance(v, str):
                    if u.deg in u.Unit(v).find_equivalent_units():
                        lat = Longitude(v).to_value(u.deg)
                    else:
                        raise u.UnitsError(
                            "Longitude must have an equivalent unit to degrees"
                        )
                else:
                    lon = float(v)
                return lon

            @field_validator("elevation")
            def check_elevation(cls, v: Any):
                if isinstance(v, str):
                    composite_height = u.Unit(v)
                    if u.m in composite_height.find_equivalent_units():
                        h = composite_height.to_value(u.m)
                    else:
                        raise u.UnitsError(
                            "Field `elevation` must have unit equivalent to meters"
                        )
                else:
                    h = float(v)
                return h

            @field_validator("fov")
            def check_fov(cls, v: Any):
                if isinstance(v, dict):
                    if "width" in v and "height" in v:
                        w = v["width"]
                        h = v["height"]
                    else:
                        raise KeyError(
                            "Field `fov` must provide `width` and `height` keys"
                            "if presented as a mapping"
                        )
                elif isinstance(v, tuple):
                    if len(v) == 2:
                        w, h = v
                    else:
                        raise TypeError(
                            "Field `fov` must be of length 2 if presented as an array"
                        )
                else:
                    raise TypeError(
                        "Field `fov` must be either a tuple of length 2 or a "
                        "dictionary with keys `width` and `height`"
                    )

                if isinstance(w, str):
                    if u.deg in w.find_equivalent_units():
                        w = w.to_value(u.deg)
                    else:
                        raise u.UnitsError("Width must have unit equivalent to `deg`")
                else:
                    w = float(w)

                if isinstance(h, str):
                    if u.deg in h.find_equivalent_units():
                        h = h.to_value(u.deg)
                    else:
                        raise u.UnitsError("Height must have unit equivalent to `deg`")
                else:
                    h = float(h)

                return w, h

            @field_validator("default_magnitude")
            def check_mag(cls, v):
                if isinstance(v, str):
                    mag = u.Unit(v)
                    if u.mag in mag.find_equivalent_units():
                        v = mag.to_value(u.mag)
                    else:
                        raise u.UnitsError("Mag must have unit equivalent to `mag`")
                else:
                    v = float(v)
                return v

            @field_validator("default_exposure")
            def check_exp(cls, v):
                if isinstance(v, str):
                    mag = u.Unit(v)
                    if u.s in mag.find_equivalent_units():
                        v = mag.to_value(u.s)
                    else:
                        raise u.UnitsError("Exposure must have unit equivalent to `s`")
                else:
                    v = float(v)
                return v

            # @TODO: Change unit to deg per second later, when it will be correctly
            # used in the `Site` class
            @field_validator("default_slew_rate")
            def check_slew_rate(cls, v):
                if isinstance(v, str):
                    rate = u.Unit(v)
                    if u.min in rate.find_equivalent_units():
                        v = rate.to_value(u.s)
                    else:
                        raise u.UnitsError(
                            "Exposure must have unit equivalent to `min`"
                        )
                else:
                    v = float(v)
                return v

            @field_validator("horizon")
            def check_horizon(cls, v):
                if isinstance(v, str):
                    horizon = u.Unit(v)
                    if u.deg in horizon.find_equivalent_units():
                        v = horizon.to_value(u.s)
                    else:
                        raise u.UnitsError("Horizon must have unit equivalent to `deg`")
                else:
                    v = float(v)
                return v

            @field_validator("timezone")
            def check_timezone(cls, v):
                import pytz

                if v not in pytz.common_timezones_set:
                    raise ValueError(
                        "Incorrect timezone, choose one from supported timezones:\n"
                        + "\n".join(pytz.common_timezones[:30])
                        + ("\nFull list can be retrieved from `pytz.common_timezones`")
                    )

                return v

        if os.path.exists(filename):
            with open(filename, "rb") as f:
                attrs = orjson.loads(f.read())

            if not isinstance(attrs, dict):
                raise ValueError("JSON file must contain only a single `Site` plugin")

            loaded_site = SiteModel(**attrs)

            return cls(
                name=loaded_site.name,
                latitude=loaded_site.latitude * u.deg,
                longitude=loaded_site.longitude * u.deg,
                elevation=loaded_site.elevation * u.m,
                timezone=loaded_site.timezone,
                fov=FOV(*(loaded_site.fov * u.deg)),
                default_magnitude=loaded_site.default_magnitude * u.mag,
                default_exposure=loaded_site.default_exposure * u.s,
                default_exposure_number=loaded_site.default_exposure_number,
                default_slew_rate=loaded_site.default_slew_rate * u.min,
                horizon=loaded_site.horizon * u.deg,
                full_name=loaded_site.full_name,
                observer_telegram_id=loaded_site.observer_telegram_id,
                observer_email=loaded_site.observer_email,
                observer_socket=loaded_site.observer_socket,
            )
        raise FileNotFoundError(f"Could not find Site plugin at {filename}")

    def nearest_observation_window(
        self,
        time: Time | datetime | str | float,
        after_sunset: TimeDelta | timedelta | u.Unit = TimeDelta(300 * u.s),
        before_sunrise: TimeDelta | timedelta | u.Unit = TimeDelta(300 * u.s),
    ) -> tuple[Time, Time]:
        return planning.nearest_observation_window(
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
            night_duration = end_time.jd - start_time.jd
            total_exposure = np.floor(
                self.default_exposure.to_value("day") * self.default_exposure_number
            )
            single_slew_time = total_exposure + self.default_slew_rate.to_value("day")
            Npoints = int(np.ceil(night_duration / single_slew_time))
            time_grid = Time(
                np.linspace(start_time.jd, end_time.jd, Npoints), format="jd"
            )
        above_horizon, altaz = self.target_is_up(
            time_grid,
            targets,
            return_altaz=True,
            grid_times_targets=True,
            horizon=self.horizon,
        )
        airmass = altaz.altaz.secz
        observable_targets = []
        for i, (a, above) in enumerate(zip(airmass, above_horizon)):
            good = np.where((a < max_airmass) & (a > 0))[0]
            if len(good) > 1:
                observable_targets.append(targets[i])

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
            log.debug("Sorting targets with algorithm `%s`", method)
            sorted_targets = planning.nn_sorter(
                self, targets, time, start_time, end_time
            )
        elif method == "radec":
            log.debug("Sorting targets with algorithm `%s`", method)
            sorted_targets = planning.radec_sorter(
                self, targets, time, start_time, end_time
            )
        else:
            log.debug(
                "Sorting algorithm `%s` is unknown, returning original list", method
            )
            # Just return a copy of original list of targets
            sorted_targets = list(targets)

        return sorted_targets

    def plot_airmass(
        self,
        targets: Sequence[FixedTarget],
        time: Optional[Time] = None,
        start_time: Optional[Time] = None,
        end_time: Optional[Time] = None,
        show_moon: bool = True,
        plot_kws: dict | None = None,
        fig_path: str = "",
    ) -> Axes:
        if start_time is None:
            start_time = time - TimeDelta(0.5 * u.d)

        if end_time is None:
            end_time = time + TimeDelta(0.5 * u.d)

        fig, ax = plt.subplots(figsize=(8, 6), dpi=96)

        def _airmass_to_alt(airmass):
            return 90 - np.rad2deg(np.arccos(1 / airmass))

        # Astroplan emits warnings on some problems with line style in graph
        # we dont want to see
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            if time is None:
                night_duration = end_time.jd - start_time.jd
                total_exposure = np.floor(
                    self.default_exposure.to_value("day") * self.default_exposure_number
                )
                single_slew_time = total_exposure + self.default_slew_rate.to_value(
                    "day"
                )
                Npoints = int(np.ceil(night_duration / single_slew_time))
                time = Time(
                    np.linspace(start_time.jd, end_time.jd, Npoints), format="jd"
                )

            plot_kws = plot_kws or {}
            plot_altitude(
                targets,
                self,
                time.datetime,
                ax=ax,
                min_altitude=0,
                max_altitude=90,
                airmass_yaxis=True,
                style_kwargs=plot_kws,
            )

        # ax.invert_yaxis()
        ax.xaxis.set_major_formatter(DateFormatter("%H"))

        if show_moon:
            airmass_time = Time(
                np.linspace(start_time.jd, end_time.jd, 150), format="jd"
            )
            altaz = self.moon_altaz(airmass_time)
            moon_airmass = altaz.secz
            good = moon_airmass > 0
            moon_alt = _airmass_to_alt(moon_airmass.value)
            ax.plot(
                airmass_time[good].datetime, moon_alt[good], "k--", label="The Moon"
            )

            # The Moon separation angle (in degrees) and its illumination phase
            annotate_time = Time(
                np.linspace(start_time.jd, end_time.jd, 5), format="jd"
            )
            annotations = []
            for target in targets:
                sep2d = altaz.separation(target.coord).to_value("deg")[0]
                airmass = self.altaz(annotate_time, target).secz
                alt = _airmass_to_alt(airmass.value)
                for t, h in zip(annotate_time, alt):
                    if h > 0:
                        annot = ax.text(
                            t.datetime,
                            h,
                            f"{sep2d:.1f}",
                            size=8,
                        )
                        annotations.append(annot)

            moon_brightess = self.moon_illumination(annotate_time)
            for ph, t, h in zip(
                moon_brightess,
                annotate_time,
                _airmass_to_alt(self.moon_altaz(annotate_time).secz.value),
            ):
                if h > 0:
                    annot = ax.text(
                        t.datetime,
                        h,
                        f"{ph:.1%}",
                        size=8,
                    )
                    annotations.append(annot)

            # Prevent annotations from overlapping
            adjust_text(
                annotations, arrowprops=dict(arrowstyle="->", color="r", lw=0.5)
            )

        # secz_ticks = np.asarray(ax.get_yticks())

        # alt = _airmass_to_alt(secz_ticks)
        # ax2 = ax.secondary_yaxis("right")
        # ax2.set_yticks(secz_ticks)
        # ax2.set_yticklabels([f"{t:.01f}" for t in alt])
        # ax2.set_ylabel("Altitude, deg")

        # Add offsets of 30 minutes from each side to fit the whole
        # observational window in the graph
        llim = (start_time - TimeDelta(30 * u.min)).datetime
        rlim = (end_time + TimeDelta(30 * u.min)).datetime
        ax.set_xlim(llim, rlim)
        ax.axvline(start_time.datetime, ls="-.", label="Obs. window")
        ax.axvline(end_time.datetime, ls="-.")
        ax.grid(True, lw=1, alpha=0.3)
        ax.tick_params(
            axis="both", which="both", top=False, right=False, direction="in"
        )
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


class TargetSeparation(constraints.Constraint):
    def __init__(self, min=None, max=None, boolean_constraint=True):
        self.min = min if min is not None else 0 * u.deg
        self.max = max if max is not None else 180 * u.deg
        self.boolean_constraint = boolean_constraint

    def compute_constraint(self, times, observer, targets):
        separations = np.full(len(targets), np.inf) * u.deg
        for i, t in enumerate(targets):
            # other_coord = targets[:i] + targets[i+1:]
            sep = t.separation(targets)
            separations[i] = sep[targets != t].min()

        if self.boolean_constraint:
            mask = (self.min < separations) & (separations < self.max)
            return mask

        # if we want to return a non-boolean score
        else:
            from astroplan import min_best_rescale

            rescale = min_best_rescale(
                separations.to_value(u.deg),
                self.min.to_value(u.deg),
                self.max.to_value(u.deg),
                less_than_min=1,
            )
            return rescale


# Separation from bright planets such as Mercury, Venus, Mars, Jupyter,
# Saturn, Uranus, and Neptune


# Site with better scheduler, which takes into account:
# 1) The Moon separation angle
# 2) The Moon phase
# 3) Visibility of the targets
# 4) Checks that each target must have a time slot
class BetterSite(Site):
    def observable_targets(
        self,
        targets: list[FixedTarget],
        start_time: Time,
        end_time: Time,
        min_moon_sep: Angle = min_moon_separation.value,
        max_moon_phase: float = max_moon_phase.value,
        min_altitude: Angle = min_altitude.value,
        max_altitude: Angle = max_altitude.value,
        # airmass_weight_ratio: float = 1.0,
        # curr_order_weight_ratio: float = 1000.0,
    ) -> list[FixedTarget]:

        # ! SunSeparationConstraint may be not passed, although it is a night
        cons = [
            constraints.AltitudeConstraint(min_altitude, max_altitude),
            constraints.MoonSeparationConstraint(min_moon_sep),
            constraints.MoonIlluminationConstraint(max=max_moon_phase),
            # constraints.SunSeparationConstraint(min=min_sun_separation),
            # constraints.AtNightConstraint(self.horizon)
        ]

        def calc_priorities() -> np.ndarray:
            """
            Calculate the priorities for target using the skymap probability and
            airmass. The skymap probability is implicitly used (assumed that targets
            already sorted in probability order)."""

            from astroplan import time_grid_from_range

            times = time_grid_from_range(
                (start_time, end_time), time_resolution=self.default_exposure
            )

            priorities = np.zeros(len(targets))
            for i, target in enumerate(targets):
                # altaz = self.altaz(times, target)
                # airmass = altaz.secz
                # rank =  (len(targets) - i) / len(targets)

                # combined_ratio = airmass_weight_ratio + curr_order_weight_ratio
                # airmass_weight = (2 - airmass.mean())

                # priority = (
                #     curr_order_weight_ratio / combined_ratio * rank
                #     + airmass_weight_ratio / combined_ratio * airmass_weight
                # )

                priorities[i] = (len(targets) - i) / len(targets)

            return priorities

        # Create the observation blocks
        blocks = []
        priorities = calc_priorities()[::-1]

        for i, target in enumerate(targets):
            block = ObservingBlock(
                target,
                self.default_exposure * self.default_exposure_number,
                priority=priorities[i],
            )
            blocks.append(block)

        scheduler = scheduling.PriorityScheduler(
            cons,
            self,
            gap_time=self.default_slew_rate,
            transitioner=scheduling.Transitioner(slew_rate=6 * u.deg / u.s),
            time_resolution=self.default_exposure,
        )

        schedule = scheduling.Schedule(start_time, end_time)
        scheduled: scheduling.Schedule = scheduler(blocks, schedule)

        obs_targets = []
        for block in scheduled.scheduled_blocks:
            if hasattr(block, "target"):
                obs_targets.append(block.target)

        return obs_targets


# Adapter that converts a Site to a BetterSite
class SiteAdapter:
    def __init__(self, site: Site):
        self.site = site
        self.better_site = BetterSite(
            name=site.name,
            latitude=site.latitude,
            longitude=site.longitude,
            elevation=site.elevation,
            timezone=site.timezone,
            fov=site.fov,
            default_magnitude=site.default_magnitude,
        )
        self.better_site.__dict__.update(self.site.__dict__)

    def __getattr__(self, name):
        if name != "observable_targets":
            return getattr(self.site, name)
        else:
            return getattr(self, name)

    def observable_targets(self, *args, **kwargs) -> list[FixedTarget]:
        return self.better_site.observable_targets(*args, **kwargs)


# Define observatories here
Mondy_AZT33IK = Site(
    timezone="Asia/Irkutsk",
    name="mondy_azt33ik",
    full_name="ISTP SO (Mondy) AZT-33IK",
    latitude=51.619822 * u.deg,
    longitude=100.9214 * u.deg,
    elevation=2000 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(15 * u.arcmin, 15 * u.arcmin),
    aperture=1.5 * u.m,
    default_exposure=120 * u.s,
    default_magnitude=21.0 * u.ABmag,
    default_filter="R",
    default_target_list_fmt="objlist",
    default_exposure_number=2,
)
CrAO_ZTSH = Site(
    timezone="Europe/Moscow",
    name="crao_ztsh",
    full_name="CrAO ZTSh",
    latitude=44.726667 * u.deg,
    longitude=34.015861 * u.deg,
    elevation=600 * u.m,
    fov=FOV(9.5 * u.arcmin, 9.5 * u.arcmin),
    aperture=2.65 * u.m,
    default_exposure=30 * u.s,
    default_magnitude=21.0 * u.ABmag,
    default_filter="R",
    default_target_list_fmt="objlist",
    default_exposure_number=1,
)
CrAO_Sintez = Site(
    timezone="Europe/Moscow",
    name="crao_sintez",
    full_name="CrAO Sintez-Newton",
    latitude=44.726667 * u.deg,
    longitude=34.015861 * u.deg,
    elevation=600 * u.m,
    fov=FOV(84 * u.arcmin, 56 * u.arcmin),
    aperture=2.65 * u.m,
    default_exposure=30 * u.s,
    default_magnitude=21.0 * u.ABmag,
    default_filter="R",
    default_target_list_fmt="objlist",
    default_exposure_number=1,
)
Assy_AZT20 = Site(
    timezone="Asia/Almaty",
    name="assy_azt20",
    full_name="Assy-Turgen AZT-20",
    latitude=43.2252778 * u.deg,
    longitude=77.8716667 * u.deg,
    elevation=2750 * u.m,
    fov=FOV(13 * u.arcmin, 13 * u.arcmin),
    aperture=1.5 * u.m,
    default_exposure=30 * u.s,
    default_magnitude=19.0 * u.ABmag,
    default_exposure_number=3,
    default_filter="r",
    list_formatter_kws={"header": False},
)
TSHAO_Zeiss1000 = Site(
    timezone="Asia/Almaty",
    name="tshao_zeiss1000",
    full_name="TSHAO Zeiss-1000",
    latitude=43.06 * u.deg,
    longitude=76.97 * u.deg,
    elevation=2735 * u.m,
    fov=FOV(13 * u.arcmin, 13 * u.arcmin),
    aperture=1.0 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=2,
    default_filter="R",
    list_formatter_kws={"header": False},
)
Altai_Santel400 = Site(
    timezone="Asia/Barnaul",
    name="altai_santel400",
    full_name="Multa Santel-400",
    latitude=53.43517913184 * u.deg,
    longitude=83.95733665302 * u.deg,
    elevation=600 * u.m,
    fov=FOV(324 * u.arcmin, 276 * u.arcmin),
    aperture=0.4 * u.m,
    default_exposure=120 * u.s,
    default_magnitude=16.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="objlist",
)
Abastumani_AS32 = Site(
    timezone="Asia/Yerevan",
    name="abao_as32",
    full_name="AbAO AS-32",
    latitude=41.7542 * u.deg,
    longitude=42.8194 * u.deg,
    elevation=1650 * u.m,
    fov=FOV(44 * u.arcmin, 44 * u.arcmin),
    aperture=0.7 * u.m,
    default_exposure=120 * u.s,
    default_magnitude=17.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
    default_filter="R",
)
SAO_RAS_Zeiss1000 = Site(
    timezone="Europe/Moscow",
    name="sao_ras_zeiss1000",
    full_name="SAO RAS Zeiss-1000",
    latitude=43.646825 * u.deg,
    longitude=41.440447 * u.deg,
    elevation=2100 * u.m,
    fov=FOV(7.37 * u.arcmin, 7.37 * u.arcmin),
    aperture=1.0 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
    default_filter="R",
)
SAO_RAS_BTA = Site(
    timezone="Europe/Moscow",
    name="sao_ras_bta",
    full_name="SAO RAS BTA",
    latitude=43.646825 * u.deg,
    longitude=41.440447 * u.deg,
    elevation=2100 * u.m,
    fov=FOV(6.1 * u.arcmin, 6.1 * u.arcmin),
    aperture=6.0 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
    default_filter="R",
)
UAFO_RC500 = Site(
    timezone="Asia/Vladivostok",
    name="uafo_rc500",
    full_name="UAFO RC500",
    latitude=Latitude("43:41:57 d"),
    longitude=Longitude("+132:09:56 d"),
    elevation=271 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(30 * u.arcmin, 20 * u.arcmin),
    aperture=0.5 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
)
Khureltogoot_ORI40 = Site(
    timezone="Asia/Ulaanbaatar",
    name="khureltogoot_ori40",
    full_name="Khureltogoot ORI-40",
    latitude=Latitude("47:51:55.2 d"),
    longitude=Longitude("+107:03:09 d"),
    elevation=1608 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(138 * u.arcmin, 138 * u.arcmin),
    aperture=0.4 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
)
MAO_AZT22 = Site(
    timezone="Asia/Tashkent",
    name="mao_azt22",
    full_name="MAO AZT-22",
    latitude=+38.6733 * u.deg,
    longitude=+66.8964 * u.deg,
    elevation=2593 * u.m,
    fov=FOV(18 * u.arcmin, 18 * u.arcmin),
    aperture=1.5 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
    default_filter="R",
)
Kitab_RC36 = Site(
    timezone="Asia/Tashkent",
    name="kitab_rc36",
    full_name="Kitab-ISON RC-36",
    latitude=+38.6733 * u.deg,
    longitude=+66.8964 * u.deg,
    elevation=2593 * u.m,
    fov=FOV(43.6 * u.arcmin, 43.6 * u.arcmin),
    aperture=0.36 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="objlist",
)
KGO_SAI25 = Site(
    timezone="Europe/Moscow",
    name="kgo_sai25",
    full_name="CMO SAI-25",
    latitude=+38.6733 * u.deg,
    longitude=+66.8964 * u.deg,
    elevation=2100 * u.m,
    fov=FOV(4.6 * u.arcmin, 4.6 * u.arcmin),
    aperture=2.5 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
    default_filter="r",
)
Terskol_K800 = Site(
    timezone="Europe/Moscow",
    name="terskol_k800",
    full_name="Terskol K-800",
    latitude=Latitude("43:16:28.4 d"),
    longitude=Longitude("+42:30:00 d"),
    elevation=3150 * u.m,
    fov=FOV(55.4 * u.arcmin, 55.4 * u.arcmin),
    aperture=0.8 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="objlist",
)
Terskol_Zeiss2000 = Site(
    timezone="Europe/Moscow",
    name="terskol_zeiss2000",
    full_name="Terskol Zeiss-2000",
    latitude=Latitude("43:16:28.4 d"),
    longitude=Longitude("+42:30:00 d"),
    elevation=3150 * u.m,
    fov=FOV(10.8 * u.arcmin, 10.8 * u.arcmin),
    aperture=2.0 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="objlist",
    default_filter="R",
)
Castelgrande_ORI22 = Site(
    timezone="Europe/Moscow",
    name="castelgrande_ori22",
    full_name="Castelgrande ORI-22",
    latitude=Latitude("40:49:04 d"),
    longitude=Longitude("+15:27:48 d"),
    elevation=1250 * u.m,
    fov=FOV(204 * u.arcmin, 240 * u.arcmin),
    aperture=0.22 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="objlist",
)
OPD_Zeiss_1000 = Site(
    timezone="Etc/GMT-3",
    name="opd_zeiss_1000",
    full_name="OPD Zeiss-1000",
    latitude=Latitude("-22:32:04 d"),
    longitude=Longitude("-45:34:57 d"),
    elevation=1864 * u.m,
    fov=FOV(7.94 * u.arcmin, 11.88 * u.arcmin),
    aperture=1.0 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
    default_filter="R",
)
OPD_IAG = Site(
    timezone="Etc/GMT-3",
    name="opd_aig",
    full_name="OPD IAG",
    latitude=Latitude("-22:32:04 d"),
    longitude=Longitude("-45:34:57 d"),
    elevation=1864 * u.m,
    fov=FOV(7.35 * u.arcmin, 11.00 * u.arcmin),
    aperture=0.6 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=17.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
    default_filter="R",
)
OPD_1_6m = Site(
    timezone="Etc/GMT-3",
    name="opd_1_6m",
    full_name="OPD 1.6-meter",
    latitude=Latitude("-22:32:04 d"),
    longitude=Longitude("-45:34:57 d"),
    elevation=1864 * u.m,
    fov=FOV(3.75 * u.arcmin, 5.57 * u.arcmin),
    aperture=1.6 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=20.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
    default_filter="R",
)
kubsu_05m = Site(
    timezone="Europe/Moscow",
    name="kubsu_0.5m",
    full_name="KubSU 0.5 meter",
    latitude=Latitude("+39:01:00 d"),
    longitude=Longitude("45:01:00 d"),
    elevation=71 * u.m,
    fov=FOV(15 * u.arcmin, 21 * u.arcmin),
    aperture=0.510 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=17.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
    default_filter="R",
)
assy_wfos40 = Site(
    timezone="Asia/Almaty",
    name="assy_wfos40",
    full_name="Assy-Turgen WFOS-40",
    latitude=43.2252778 * u.deg,
    longitude=77.8716667 * u.deg,
    elevation=2750 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(3.8 * u.degree, 2.5 * u.degree),
    aperture=0.4 * u.m,
    default_exposure=30 * u.s,
    default_magnitude=17.0 * u.ABmag,
    default_exposure_number=3,
    default_filter="R",
    default_target_list_fmt="txt",
    list_formatter_kws={"header": False},
)


# Hook that loads JSON plugins
def __load_site_plugins(plugin_dir: str = site_plugin_directory.value):
    global Telescopes
    for fn in glob(f"{os.path.realpath(plugin_dir)}/*.json"):
        try:
            telescope = Site.from_json(fn)
        except KeyError as e:
            log.error("Could not load JSON Site plugin: %s", e)
        else:
            Telescopes[telescope.name] = telescope


__load_site_plugins()


# TODO: in the future all sites should be already better sites
# Adapt each telescope
__TELESCOPES_ADAPTED = False


def __adapt_telescopes():
    global __TELESCOPES_ADAPTED
    global Telescopes
    if not __TELESCOPES_ADAPTED:
        for name, tel in Telescopes.items():
            Telescopes[name] = SiteAdapter(tel)
        __TELESCOPES_ADAPTED = True


__adapt_telescopes()

default_sites = CfgOption("default_sites", list(Telescopes.keys()), list)
