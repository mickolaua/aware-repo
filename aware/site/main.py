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
from astropy.coordinates import Angle, Latitude, Longitude

from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib.axes import Axes

from aware.logger import log

from .. import planning
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
        return np.hypot(0.5 * self.width, 0.5 * self.height)

    @cached_property
    def is_widefield(self) -> bool:
        return (
            self.width.to_value("arcmin") >= 40 or self.height.to_value("arcmin") >= 40
        )

    def __str__(self) -> str:
        return "{:.3f} deg x {:.3f} deg".format(
            self.width.to_value(u.deg), self.height.to_value(u.deg)
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
        default_magnitude: u.Unit,
        default_exposure: u.Unit = 1 * u.s,
        default_exposure_number: int = 1,
        default_filter: str = "",
        default_target_list_fmt: str = "txt",
        default_slew_rate: u.Unit = 3 * u.min,
        horizon: u.Unit = 0.0 * u.degree,
        full_name: Optional[str] = None,
        observer_telegram_id: str = "",
        **kwargs: Any,
    ) -> None:
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

        # Add the site to the registry
        global Telescopes
        Telescopes[self.name] = self

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
            time_grid = Time(np.linspace(start_time.jd, end_time.jd, 150), format="jd")
        above_horizon, altaz = self.target_is_up(
            time_grid,
            targets,
            return_altaz=True,
            grid_times_targets=True,
            horizon=-18 * u.deg,
        )
        airmass = altaz.altaz.secz
        # observable = np.argwhere(
        #     np.any((airmass < max_airmass) & (airmass > 0), axis=1)
        #     & np.any(above_horizon, axis=1)
        # ).ravel()
        # print(airmass)
        # print(np.any(airmass > 0, axis=0))
        # print(np.where(airmass < max_airmass))
        # print(np.where(airmass > 0))
        # print(above_horizon)
        # observable = (
        #     np.any(airmass > 0, axis=1).ravel()
        #     & np.any(airmass < max_airmass, axis=1).ravel()
        # )
        # idx = np.arange(len(airmass))
        # obs_idx = idx[observable]
        # print(observable)
        observable_targets = []
        for i, (a, above) in enumerate(zip(airmass, above_horizon)):
            good = np.where((a < max_airmass) & (a > 0))[0]
            if len(good) > 1:
                observable_targets.append(targets[i])

        # observable_targets = [targets[i] for i in obs_idx]

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

        ax.xaxis.set_major_formatter(DateFormatter("%H"))

        if show_moon:
            airmass_time = Time(
                np.linspace(start_time.jd, end_time.jd, 150), format="jd"
            )
            altaz = self.moon_altaz(airmass_time)
            moon_airmass = altaz.secz
            good = moon_airmass > 0
            ax.plot(
                airmass_time[good].datetime, moon_airmass[good], "k--", label="The Moon"
            )

            # TODO: Separation angle between the Moon and a target (in degrees)
            # for tgt in targets:

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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
    fov=FOV(13 * u.arcmin, 13 * u.arcmin),
    aperture=1.5 * u.m,
    default_exposure=30 * u.s,
    default_magnitude=19.0 * u.ABmag,
    default_exposure_number=3,
    default_filter="r",
)
TSHAO_Zeiss1000 = Site(
    timezone="Asia/Almaty",
    name="tshao_zeiss1000",
    full_name="TSHAO Zeiss-1000",
    latitude=43.06 * u.deg,
    longitude=76.97 * u.deg,
    elevation=2735 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(13 * u.arcmin, 13 * u.arcmin),
    aperture=1.0 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=18.0 * u.ABmag,
    default_exposure_number=2,
    default_filter="R",
)
Altai_Santel400 = Site(
    timezone="Asia/Barnaul",
    name="altai_santel400",
    full_name="Multa Santel-400",
    latitude=53.43517913184 * u.deg,
    longitude=83.95733665302 * u.deg,
    elevation=600 * u.m,
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    horizon=0.0 * u.degree,
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
    full_name="OPD Zeiss-1000",
    latitude=Latitude("-22:32:04 d"),
    longitude=Longitude("-45:34:57 d"),
    elevation=1864 * u.m,
    horizon=0.0 * u.degree,
    fov=FOV(3.75 * u.arcmin, 5.57 * u.arcmin),
    aperture=1.6 * u.m,
    default_exposure=60 * u.s,
    default_magnitude=20.0 * u.ABmag,
    default_exposure_number=3,
    default_target_list_fmt="txt",
    default_filter="R",
)


default_sites = CfgOption("default_sites", list(Telescopes.keys()), list)
