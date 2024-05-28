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

from aware.config import CfgOption
from aware.logger import log
from aware import site


max_area_trigger = CfgOption("max_area_trigger", 50.0, float)


def nearest_observation_window(
    site: Observer,
    time: Time | datetime | str | float,
    horizon: u.Unit = -12 * u.degree,
    after_sunset: TimeDelta | timedelta | u.Unit = TimeDelta(5400 * u.s, format="sec"),
    before_sunrise: TimeDelta
    | timedelta
    | u.Unit = TimeDelta(5400 * u.s, format="sec"),
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
    log.debug("Searching for nearest night about %s", time.isot)


    sunrise = site.sun_rise_time(time, which="nearest", horizon=site.horizon)
    sunset = site.sun_set_time(time, which="nearest", horizon=site.horizon)
    log.debug("Initial night borders: since %s to %s", sunset.isot, sunrise.isot)

    if sunrise.mask or sunset.mask:
        return None, None

    def twilight_now(time, sunrise, sunset) -> bool:
        return sunset < time < sunrise

    if twilight_now(time, sunrise, sunset):
        return sunset, sunrise

    log.debug("Is twilight now?: %s", twilight_now(time, sunrise, sunset))

    if time < sunset < sunrise:
        return sunset, sunrise

    if time < sunrise < sunset:
        sunset = site.sun_set_time(sunset, which="previous", horizon=site.horizon)
        return sunset, sunrise

    if sunrise < sunset < time:
        sunrise = site.sun_rise_time(sunset, which="next", horizon=site.horizon)
        return sunset, sunrise

    if sunrise < time < sunset:
        sunrise = site.sun_rise_time(sunset, which="next", horizon=site.horizon)
        return sunset, sunrise

    return sunset, sunrise


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
        night_duration = end_time.jd - start_time.jd
        total_exposure = np.floor(
            site.default_exposure.to_value("day") * site.default_exposure_number
        )
        single_slew_time = total_exposure + site.default_slew_rate.to_value("day")
        Npoints = int(np.ceil(night_duration / single_slew_time))
        time_grid = Time(
            np.linspace(start_time.jd, end_time.jd, Npoints), format="jd"
        )
        at_low_airmass = np.argsort(
            np.sum(
                site.altaz(time_grid, targets, grid_times_targets=True).secz
                < max_airmass,
                axis=1,
            )
        ).ravel()
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

    log.debug(
        "Starting with target %s since it has highest average airmass during the night",
        targets_by_airmass[0].name,
    )
    sorted_targets: list[FixedTarget] = [targets_by_airmass[0]]

    if targets_by_airmass:
        coord = SkyCoord([t.coord for t in targets_by_airmass])
    else:
        return sorted_targets

    # Create the separation matrix
    N = len(targets_by_airmass)
    dist_matrix = np.zeros((N, N))

    # Fill the matrix
    for i, tgt in enumerate(targets_by_airmass):
        dist_matrix[i, :] = tgt.coord.separation(coord).to_value(u.arcsec)

    # Perform sorting
    order_idx = np.zeros_like(targets_by_airmass)
    mask = np.zeros_like(order_idx).astype(int)
    seq_obj_id = 0
    for i in range(N - 1):
        # Find closest neighbor id
        uniq = np.unique(order_idx).astype(int)
        mask[uniq] = 1
        dist = np.ma.masked_array(dist_matrix[seq_obj_id, :], mask=mask)
        next_obj_id = np.argmin(dist)
        
        # Assign next target id
        order_idx[i + 1] = next_obj_id
        seq_obj_id = next_obj_id

    # Combine the first target with others
    sorted_targets = [targets_by_airmass[i] for i in order_idx]

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
