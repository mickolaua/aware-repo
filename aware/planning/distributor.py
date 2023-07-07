"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
distributor.py (c) 2023
Desc: distriubte objects between telescope to observe without intersection
Created:  2023-06-04
Modified: !date!
"""
from __future__ import annotations

from datetime import datetime
from typing import Any, Sequence

import healpy as hp
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

from ..field import Field
from ..glade.main import GladeGalaxy
from ..localization.main import lvk_uncert_level
from ..site import Telescopes, default_sites
from .mosaic import mosaic_walker
from ..logger import log


def _targets_per_night(
    start: Time,
    stop: Time,
    exptime: u.Unit = 1 * u.s,
    expcount: int = 1,
    slew_time: u.Unit = 3 * u.min,
) -> int:
    night_duration = (stop.jd - start.jd) * u.day
    total_exptime = ((exptime * expcount) + slew_time).to(u.day)
    N = int(np.floor((night_duration / total_exptime)))
    return N


def _get_target_hpx_id(hpx: np.ndarray, target: GladeGalaxy) -> int:
    nside = hp.get_nside(hpx)
    ipix = hp.ang2pix(
        nside, target.coord.ra.deg, target.coord.dec.deg, lonlat=True, nest=True
    )
    return ipix


def paint_planned_healpix_from_fields(
    hpx: np.ndarray, fields: Sequence[Field]
) -> np.ndarray:
    """
    Paint the HEALPix map with zeros inside the given sky fields to mark these
    areas as already planned for observation.

    Parameters
    ----------
    hpx : np.ndarray
        an HEALPix array of per pixel probability (1D) in nest ordering
    fields : Sequence[Field]
        a sequence of sky fields returned, for example, after mosaic planning

    Returns
    -------
    hpx_copy: np.ndarray
        a copy of the original array, where pixels inside sky fields are zeroed
    """
    Npix = len(hpx)
    hpx_copy = hpx.copy()
    ipix = np.arange(Npix)
    nside = hp.get_nside(hpx)

    for field in fields:
        max_ra = field.coord.ra.deg + field.width.to_value(u.deg) * 0.5
        min_ra = field.coord.ra.deg - field.width.to_value(u.deg) * 0.5
        max_dec = field.coord.dec.deg + field.height.to_value(u.deg) * 0.5
        min_dec = field.coord.dec.deg - field.height.to_value(u.deg) * 0.5

        pix_ra, pix_dec = hp.pix2ang(nside, ipix, lonlat=True, nest=True)
        draw_mask = (
            (pix_ra <= max_ra)
            & (pix_ra >= min_ra)
            & (pix_dec >= min_dec)
            & (pix_dec <= max_dec)
        )

        hpx_copy[ipix[draw_mask]] = 0.0

    return hpx_copy


def paint_planned_healpix_from_galaxies(
    hpx: np.ndarray, galaxies: Sequence[GladeGalaxy]
) -> np.ndarray:
    """
    Paint the HEALPix map with zeros at the locations of the given Glade+ galaxies
    to mark these pixels as already planned for observation.

    Parameters
    ----------
    hpx : np.ndarray
        an HEALPix array of per pixel probability (1D) in nest ordering
    galaxies : Sequence[GladeGalaxy]
        a sequence of Glade+ galaxies returned, for example, after mosaic planning

    Returns
    -------
    hpx_copy: np.ndarray
        a copy of the original array, where pixels at the galaxies locations are zeroed
    """
    Npix = len(hpx)
    hpx_copy = hpx.copy()
    ipix = np.arange(Npix)
    nside = hp.get_nside(hpx)

    for gal in galaxies:
        draw_mask = hp.ang2pix(
            nside, gal.coord.ra.deg, gal.coord.dec.deg, nest=True, lonlat=True
        )
        hpx_copy[ipix[draw_mask]] = 0.0

    return hpx_copy


class TargetDistributor:
    def __init__(self, hpx: np.ndarray, objects: Sequence[GladeGalaxy]) -> None:
        objects = [obj for obj in objects if hpx[_get_target_hpx_id(hpx, obj)] != 0]
        self.hpx = hpx
        self._orig_hpx = hpx.copy()
        self.objects = set(objects)
        self._orig_objects = set(objects)
        self.day = 0
        self.plans: dict[int, dict[str, Sequence[GladeGalaxy]]] = {}

    def distribute(
        self,
        telescope_ids: Sequence[str] = default_sites.value,
        init_sep: u.Unit = 1 * u.deg,
        step_sep: u.Unit = 1 * u.deg,
        max_sep: u.Unit = 360 * u.deg,
        max_targets: int | None = None,
        disable_intersections: bool = True,
        initial_time: Time | datetime | None = None,
    ) -> dict[str, list[GladeGalaxy]]:
        telescopes = [Telescopes[i] for i in telescope_ids]
        groups = {}
        if list(self.objects):
            for tel in telescopes:
                it = (
                    Time(initial_time)
                    if initial_time is not None
                    else Time(datetime.now())
                )
                start, stop = tel.nearest_observation_window(it)
                if start is None or stop is None:
                    groups[tel.name] = []
                else:
                    try:
                        # How many targets per one night can a telescope observe?
                        N = _targets_per_night(
                            start,
                            stop,
                            exptime=tel.default_exposure,
                            expcount=tel.default_exposure_number,
                            slew_time=tel.default_slew_rate,
                        )
                        if max_targets is not None:
                            N = max_targets if max_targets < N else N

                        # Create a cluster of targets for further observing
                        # We are trying to keep targets in the cluster so a telescope
                        # Will not jump through targets that could be at a large
                        # separation
                        i = 0
                        sep = init_sep
                        obs_targets = tel.observable_targets(
                            list(self.objects), start_time=start, end_time=stop
                        )
                        most_probable_target = max(
                            obs_targets,
                            key=lambda t: self.hpx[_get_target_hpx_id(self.hpx, t)],
                        )
                        planned_targets = [most_probable_target]

                        # Expand separation limit while there are enough targets to
                        # include into plan
                        targets = list(set(obs_targets) - {most_probable_target})
                        coord = SkyCoord([t.coord for t in targets])
                        separations = most_probable_target.coord.separation(coord)
                        idx = np.argwhere(separations <= max_sep).ravel()
                        if idx.size:
                            sorted_neighbors = [targets[i] for i in idx]
                            sorted_neighbors.sort(
                                key=lambda t: t.coord.separation(
                                    most_probable_target.coord
                                )
                            )
                            planned_targets += sorted_neighbors[: N - 1]

                        # Sort planned targets for more efficient observation
                        ordered_targets = tel.observation_order(
                            planned_targets,
                            start_time=start,
                            end_time=stop,
                            method="nearest",
                        )
                        groups[tel.name] = ordered_targets

                        # When we want to plan only unique targets, we update master
                        # list on planning for each telescope
                        if disable_intersections:
                            self.objects.difference_update(planned_targets)

                            # Mark planned pixels
                            self.hpx = paint_planned_healpix_from_galaxies(
                                self.hpx, planned_targets
                            )
                    except (LookupError, AttributeError, ValueError) as e:
                        log.error(
                            "Observations can not be carried out with %s",
                            tel.full_name,
                            exc_info=e,
                        )

            # In the intersection mode, update master list after planning occured for
            # all the telescopes, since we allow the same galaxies in plans.
            if not disable_intersections and groups.values():
                targets_to_remove: list[GladeGalaxy] = []
                for v in groups.values():
                    targets_to_remove += v
                self.objects.difference_update(targets_to_remove)

                # Mark planned pixels
                self.hpx = paint_planned_healpix_from_galaxies(
                    self.hpx, targets_to_remove
                )

        self.day += 1
        self.plans[self.day] = groups

        return groups

    def reset(self):
        self.hpx = self._orig_hpx.copy()
        self.objects = self._orig_objects.copy()
        self.day = 0
        self.plans = {}


class FieldDistributor:
    def __init__(self, hpx: np.ndarray):
        self.hpx = hpx
        self._orig_hpx = hpx.copy()
        self.day = 0
        self.plans: dict[int, dict[str, Sequence[Field]]] = {}

    def distribute(
        self,
        telescope_ids: Sequence[str] = default_sites.value,
        prob: float = lvk_uncert_level.value,
        disable_intersections: bool = True,
        initial_time: Time | datetime | None = None,
    ) -> dict[str, list[Field]]:
        telescopes = [Telescopes[i] for i in telescope_ids]
        groups = {}
        for tel in telescopes:
            init = (
                Time(initial_time) if initial_time is not None else Time(datetime.now())
            )
            start, stop = tel.nearest_observation_window(init)

            if start is None or stop is None:
                groups[tel.name] = []
            else:
                try:
                    # First off, walk the contours and retrieve all possible sky fields
                    # of a size of the telescope FOV
                    nside = hp.get_nside(self.hpx)
                    npix = hp.get_map_size(self.hpx)
                    ring_idx = hp.ring2nest(nside, np.arange(npix))
                    fields = mosaic_walker(self.hpx[ring_idx], tel.fov, prob=prob)

                    # Retrieve the list of fields that a telescope can observe in
                    # principal
                    obs_targets = tel.observable_targets(
                        fields, start_time=start, end_time=stop
                    )

                    # How many fields per one night can a telescope observe?
                    N = _targets_per_night(
                        start,
                        stop,
                        tel.default_exposure,
                        tel.default_exposure_number,
                        tel.default_slew_rate,
                    )
                    obs_targets = list(obs_targets[:N])

                    # Add to the group
                    groups[tel.name] = obs_targets

                    # Mark planned pixels
                    if disable_intersections:
                        self.hpx = paint_planned_healpix_from_fields(
                            self.hpx, obs_targets
                        )
                except (LookupError, AttributeError, ValueError) as e:
                    log.error(
                        "Observations can not be carried out for %s",
                        tel.full_name,
                        exc_info=e,
                    )

        self.day += 1
        self.plans[self.day] = groups

        return groups

    def reset(self):
        self.hpx = self._orig_hpx.copy()
        self.day = 0
        self.plans = {}
