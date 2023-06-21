"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
distributor.py (c) 2023
Desc: distriubte objects between telescope to observe without intersection
Created:  2023-06-04
Modified: !date!
"""
from __future__ import annotations
from typing import Any, Sequence

import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from astroplan.target import FixedTarget
from astroplan.scheduling import ObservingBlock
from astroplan.constraints import MoonSeparationConstraint, MoonIlluminationConstraint
from datetime import datetime
import healpy as hp

from aware.field import Field
from aware.glade.main import GladeGalaxy

from ..glade.main import GladeGalaxy, GladeCatalog
from ..site import Telescopes, default_sites
from ..field import Field
from .mosaic import mosaic_walker
from ..localization.main import lvk_uncert_level


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
    ipix = hp.ang2pix(nside, target.coord.ra.deg, target.coord.dec.deg, lonlat=True)
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
        disable_intersections: bool = True
    ) -> dict[str, list[GladeGalaxy]]:
        telescopes = [
            Telescopes[i] for i in telescope_ids if not Telescopes[i].fov.is_widefield
        ]
        groups = {}
        for tel in telescopes:
            start, stop = tel.nearest_observation_window(Time(datetime.now()))

            if start is None or stop is None:
                groups[tel.name] = []
            else:
                # Sort master list of galaxies by probability so during observations
                # we walk around most probable pixels
                sorted_galaxies = set(
                    sorted(
                        self.objects,
                        key=lambda _t: self.hpx[_get_target_hpx_id(self.hpx, _t)],
                        reverse=True,
                    )
                )

                # Retrieve the list of galaxies that a telescope can observe
                obs_targets = tel.observable_targets(
                    list(sorted_galaxies), start_time=start, end_time=stop
                )

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
                obs_targets = list(obs_targets[:N])

                # Pack planned targets for a current telescope
                planned_targets = []
                sorted_targets = tel.observation_order(
                    obs_targets, start_time=start, end_time=stop, method="radec"
                )

                # Create a cluster of targets for further observing
                # We are trying to keep targets in the cluster so a telescope
                # Will not jump through targets that could be at a large separation
                i = 0
                sep = init_sep

                # Expand separation limit while there are enough targets to include
                # into plan
                while i < N and sep < max_sep:
                    for t in sorted_galaxies:
                        if t in sorted_targets:
                            if i and planned_targets:
                                if (
                                    planned_targets[0] is not t
                                    and planned_targets[0].coord.separation(t.coord)
                                    < sep
                                ):
                                    planned_targets.append(t)
                            else:
                                planned_targets.append(t)

                        i += 1

                    sep += step_sep

                # Sort planned targets for more efficient observation
                ordered_targets = tel.observation_order(
                    planned_targets, start_time=start, end_time=stop, method="nearest"
                )
                groups[tel.name] = ordered_targets

                if disable_intersections:
                    # Update master list, remove already planned targets
                    self.objects.difference_update(planned_targets)

                    # Mark planned pixels
                    self.hpx = paint_planned_healpix_from_galaxies(
                        self.hpx, self.objects
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
    ) -> dict[str, list[Field]]:
        telescopes = [
            Telescopes[i] for i in telescope_ids if Telescopes[i].fov.is_widefield
        ]
        groups = {}
        for tel in telescopes:
            start, stop = tel.nearest_observation_window(Time(datetime.now()))

            if start is None or stop is None:
                groups[tel.name] = []
            else:
                # First off, walk the contours and retrieve all possible sky fields
                # of a size of the telescope FOV
                nside = hp.get_nside(self.hpx)
                npix = hp.get_map_size(self.hpx)
                ring_idx = hp.ring2nest(nside, np.arange(npix))
                fields = mosaic_walker(self.hpx[ring_idx], tel.fov, prob=prob)

                # Retrieve the list of fields that a telescope can observe in principal
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
                    self.hpx = paint_planned_healpix_from_fields(self.hpx, obs_targets)

        self.day += 1
        self.plans[self.day] = groups

        return groups
