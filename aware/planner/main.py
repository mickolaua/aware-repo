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
import matplotlib as mpl

from aware.glade import GladeGalaxy

mpl.use("agg")

from ..config import CfgOption
from ..logger import log
from .. import site


def nearest_observation_window(
    site: Observer,
    time: Time | datetime | str | float,
    horizon: u.Unit = -12 * u.degree,
    after_sunset: TimeDelta | timedelta | u.Unit = TimeDelta(300 * u.s),
    before_sunrise: TimeDelta | timedelta | u.Unit = TimeDelta(300 * u.s),
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
    sunrise = site.sun_rise_time(time, which="nearest", horizon=horizon)
    sunset = site.sun_set_time(time, which="nearest", horizon=horizon)

    if sunrise.mask or sunset.mask:
        return None, None

    if sunrise > time:
        if sunset < time:
            start_time = sunset
        else:
            start_time = site.sun_set_time(time, which="previous", horizon=horizon)
        end_time = sunrise
    else:
        if sunset < time:
            start_time = site.sun_set_time(time, which="next", horizon=horizon)
        else:
            start_time = sunset
        end_time = site.sun_rise_time(time, which="next", horizon=horizon)

    # Account for time shifts after sun set and before sun rise
    start_time = start_time + after_sunset
    end_time = end_time - before_sunrise

    return start_time, end_time


def radec_sorter(
    site: site.Site,
    targets: Sequence[FixedTarget],
    time: Optional[Time] = None,
    start_time: Optional[Time] = None,
    end_time: Optional[Time] = None,
) -> list[FixedTarget]:
    """Sort targets firstly by R.A., then by Dec"""
    return sorted(targets, key=lambda _t: (_t.ra.deg, _t.dec.deg))


def nn_sorter(
    site: site.Site,
    targets: Sequence[FixedTarget],
    time: Optional[Time] = None,
    start_time: Optional[Time] = None,
    end_time: Optional[Time] = None,
    max_airmass: float = 3,
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
        time_grid = Time(np.linspace(start_time.jd, end_time.jd, 150), format="jd")
        at_low_airmass = np.argsort(
            np.sum(
                site.altaz(time_grid, targets, grid_times_targets=True).secz
                < max_airmass,
                axis=1,
            )
        ).ravel()
        log.info(at_low_airmass)
        targets_by_airmass = [targets[i] for i in at_low_airmass]
    else:
        if time is None:
            raise ValueError(
                "Single time should be provided when start_time and "
                "end_time not passed"
            )

        # Start with target at highest airmass
        targets_by_airmass: list[FixedTarget] = sorted(
            targets, key=lambda _target: site.altaz(time, _target).secz, reverse=True
        )

    # Perform actual sorting with Nearest-Neighbor algorithm
    sorted_targets = [targets_by_airmass[0]]
    coords = SkyCoord([t.coord for t in targets_by_airmass])

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


def mosaic_sorter(
    site: site.Site,
    targets: Sequence[FixedTarget | GladeGalaxy],
    time: Optional[Time] = None,
    start_time: Optional[Time] = None,
    end_time: Optional[Time] = None,
    max_airmass: float = 3,
) -> list[FixedTarget | GladeGalaxy]:
    


    return []
