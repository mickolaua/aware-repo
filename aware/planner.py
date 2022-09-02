from __future__ import annotations

import warnings
from datetime import datetime, timedelta
from typing import Any, Optional, Sequence

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astroplan import FixedTarget
from astroplan.observer import Observer
from astroplan.plots import plot_airmass
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from matplotlib.axes import Axes


def _nearest_observation_window(
    site: Observer,
    time: Time | datetime | str | float,
    horizon: u.Unit = -12 * u.degree,
    after_sunset: TimeDelta | timedelta | u.Unit = TimeDelta(300*u.s),
    before_sunrise: TimeDelta | timedelta | u.Unit = TimeDelta(300*u.s)
) -> tuple[Time, Time]:
    """Get the closest time interval(s) (start and end time) when the site can 
    observe any targets. Such an interval is found inside 24-hr time window 
    centered at the specified time (see parameters above). 

    Parameters
    ----------
    site : Observer
        a site that observes targets
    time : Union[Time, datetime, str, float]
        an approximate time when targets may be observed by the site
        (used as guess to find sin rise and sun set times)
    horizon : u.Unit, optional
        an actual elevation above horizon that the site can observe
    after_sunset : Union[Time, datetime, str, float]
        time shift after sunset, to prevent observations at the sunset
    before_sunset : Union[Time, datetime, str, float]
        time shift before sunrise, to prevent observations at the sunrise

    Returns
    -------
    start_time, end_time : tuple[datetime, datetime]
        start and end time of the observational period
    """

    # Calculate the nearest observational window by finding appropriate sun 
    # rise and sun set times
    sunrise = site.sun_rise_time(time, which='nearest', horizon=horizon)
    sunset = site.sun_set_time(time, which='nearest', horizon=horizon)

    # if time >= sunrise:
    #     sunrise = site.sun_rise_time(time, which='next', horizon=horizon)

    #     if sunset >= sunrise:
    #         sunset = site.sun_set_time(time, which='previous', horizon=horizon)

    # if sunset >= sunrise:
    #     sunset = site.sun_set_time(time, which='previous', horizon=horizon)

    if sunrise > time:
        if sunset < time:
            start_time = sunset
        else:
            start_time = site.sun_set_time(
                time, which='previous', horizon=horizon
            )
        end_time = sunrise
    else:
        if sunset < time:
            start_time = site.sun_set_time(
                time, which='next', horizon=horizon
            )
        else:
            start_time = sunset
        end_time = site.sun_rise_time(
                time, which='next', horizon=horizon
            )

    # Account for time shifts after sun set and before sun rise
    start_time = start_time + after_sunset
    end_time = end_time - before_sunrise

    return start_time, end_time


def _nn_sorter(
    site: Site, 
    targets: Sequence[FixedTarget], 
    time: Optional[Time] = None,
    start_time: Optional[Time] = None,
    end_time: Optional[Time] = None
) -> list[FixedTarget]:
    """Sort targets with Nearest-Neighbor algorithm, starting from the target 
    that at the time has a maximal airmass with respect to other targets.

    Parameters
    ----------
    site: Site
        site at which targets will be observed
    targets : Sequence[FixedTarget]
        obsevational targets to sort
    time : Time
        time at which the airmass is computed

    Returns
    -------
    sorted_targets: list[FixedTarget]
        list of targets in the order they will be observed
    """
    # Time bounds provided, we can better sort targets using such approach:
    # 1. Calculate sum of airmass that a target has below airmass = 3
    # 2. Sort by this criterion
    # 3. Choose the target with the highest value
    #
    # Approach above provids better sorting, especially in cases when a target
    # from the list can be observed only in the narrow time window, but it's 
    # airmass at the certain time is not minimal compared with other targets.
    if start_time is not None and end_time is not None:
        time_grid = Time(
            np.linspace(start_time.jd, end_time.jd, 150), format="jd"
        )
        targets_by_airmass: list[FixedTarget] = sorted(
            targets,
            key=lambda _target: np.sum(site.altaz(time_grid, _target).secz < 3),
            reverse=True
        )

    else:
        if time is None:
            raise ValueError(
                "Single time should be provided when start_time and end_time "
                +"not passed"
            )

        # Start with target at highest airmass
        targets_by_airmass: list[FixedTarget] = sorted(
            targets,
            key=lambda _target: site.altaz(time, _target).secz,
            reverse=True
        )

    # Perform actual sorting with Nearest-Neighbor algorithm
    sorted_targets: list[FixedTarget]  = []
    while targets_by_airmass:
        sorted_targets.append(target := targets_by_airmass.pop(0))

        if not targets_by_airmass:
            break

        coords = SkyCoord([t.coord for t in targets_by_airmass])
        idx, _, _ = target.coord.match_to_catalog_sky(coords)
        tgt = targets_by_airmass[idx]
        targets_by_airmass.remove(tgt)
        targets_by_airmass = [tgt] + targets_by_airmass

    return sorted_targets


class Site(Observer):
    def __init__(
        self, 
        horizon: u.Unit = 0.0 * u.degree,
        full_name: Optional[str] = None,
        **kwargs: Any
    ) -> None:
        super().__init__(**kwargs)
        self.horizon = horizon
        if full_name is not None:
            self.full_name = full_name
        else:
            self.full_name = self.name

    def nearest_observation_window(
        self, 
        time: Time | datetime | str | float,
        after_sunset: TimeDelta | timedelta | u.Unit = TimeDelta(300*u.s),
        before_sunrise: TimeDelta | timedelta | u.Unit = TimeDelta(300*u.s)
    ) -> tuple[datetime, datetime]:
        return _nearest_observation_window(
            self, time, self.horizon, after_sunset, before_sunrise
        )

    def observable_targets(
        self,
        targets: Sequence[FixedTarget], 
        start_time: Optional[Time] = None,
        end_time: Optional[Time] = None,
        max_airmass: float = 3
    ) -> list[FixedTarget]:

        observable_targets = []
        if start_time is None and end_time is not None:
            time_grid = [end_time]
        elif start_time is not None and end_time is None:
            time_grid = [start_time]
        else:
            time_grid = Time(
                np.linspace(start_time.jd, end_time.jd, 150), format="jd"
            )
        for target in targets:
            above_horizon = self.target_is_up(time_grid, target)
            airmass = self.altaz(time_grid, target).secz
            if np.any(above_horizon) and np.any(airmass < max_airmass):
                observable_targets.append(target)

        return observable_targets

    def observation_order(
        self, 
        targets: Sequence[FixedTarget],
        time: Time,
        method: str = "nearest"
    ) -> list[FixedTarget]:
        
        # TODO: implement new methods and improve existing
        if method == "nearest":
            sorted_targets = _nn_sorter(self, targets, time)
        else:
            # Just return original list of targets
            sorted_targets = list(targets)
        
        return sorted_targets

    def plot_airmass(
        self, 
        targets: Sequence[FixedTarget], 
        time: Time, 
        start_time: Optional[Time] = None,
        end_time: Optional[Time] = None,
        show_moon: bool = True,
        fig_path: str = "airmass.png"
    ) -> Axes:
        if start_time is None:
            start_time = time - TimeDelta(0.5*u.d)

        if end_time is None:
            start_time = time + TimeDelta(0.5*u.d)

        fig, ax = plt.subplots(figsize=(8, 6), dpi=96)
        
        # Astroplan emits warnings on some problems with line style in graph
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",  category=UserWarning)
            plot_airmass(
                targets, 
                self, 
                time.datetime, 
                ax=ax, 
                min_airmass=0.8, 
                max_airmass=5
            )

        if show_moon:
            airmass_time = Time(
                np.linspace(start_time.jd, end_time.jd, 150), format="jd"
            )
            moon_airmass = self.moon_altaz(airmass_time.datetime).secz
            ax.plot(
                airmass_time.datetime, moon_airmass, "k--", label="The Moon"
            )
        
        # Add offsets of 30 minutes from each side to fit the whole 
        # observational window in the graph
        ax.set_xlim(
            (start_time - TimeDelta(30*u.min)).datetime, 
            (end_time + TimeDelta(30*u.min)).datetime
        )
        ax.axvline(start_time.datetime, ls="-.", label="Obs. window")
        ax.axvline(end_time.datetime, ls="-.")
        ax.grid(True, lw=1, alpha=0.3)
        ax.tick_params(
            axis="both", which="both", top=True, right=True, direction="in"
        )
        ax.tick_params(axis="x", rotation=0)

        ax.set_title(
            "Site {:s}\nLongitude: {:.5f} deg, Latitude: {:.5f} deg".format(
                self.full_name,
                self.location.lon.value,
                self.location.lat.value,
            ),
            loc="left"
        )

        lgd = ax.legend(
            bbox_to_anchor=(1.01, 1), 
            loc='upper left', 
            borderaxespad=0,
            frameon=False,
            title="List of Objects",
        )

        # Reposition Axes to fit both Axes and legend in to the canvas
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
        shift = 1 - (lgd_xmax - ax_xmax) * 0.7
        fig.tight_layout(rect=(0, 0, shift, 1))

        if fig_path:
            fig.savefig(fig_path)

        return ax


# Define observatories here
mondy = Site(
    timezone="Asia/Irkutsk", 
    name="mondy",
    full_name="Sayan Solar Observatory at Mondy, Irkutsk",
    latitude=51.619822*u.deg, 
    longitude=100.9214*u.deg,
    elevation=2000*u.m,
    horizon=0.0*u.degree
)

sites: dict[str, Site] = {}
