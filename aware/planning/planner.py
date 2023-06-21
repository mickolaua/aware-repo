from __future__ import annotations
from datetime import datetime

import os
from typing import Any, Sequence

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from ligo.skymap.io.fits import read_sky_map, write_sky_map
from ligo.skymap.moc import rasterize
import healpy as hp

from aware.field import Field
from aware.glade import GladeCatalog, GladeGalaxy
from aware.logger import log
from aware.planning.distributor import FieldDistributor, TargetDistributor
from aware.planning.program import create_observation_program
from aware.site.main import Telescopes

from ..io import hpx_from_moc
from ..localization.main import lvk_uncert_level


class SkymapPlanner:
    """
    Plans the optical observations of a large skymap with a telescope network. It is
    considered that skymap should be observed from most to least credible locations to
    find the astrophysical source there. There are two methods of planning implemented:
    target observations of GLADE+ galaxies with narrow-field telescopes, and mosaic
    scanning of sky fields with wide-field telescopes. Both tactics could be combined
    or used isolately. Along with it, planner provides possibility to distribute
    targets and/or sky tiles uniquely, so no pair of telescopes will observe the same
    location. All targets/sky fields are grouped in clusters of radius of several
    degrees to prevent sharp transitions between them. This approach assignes different
    sky patches to different scopes. However, this can be switched off. The bservation
    targets and sky fields are sorted in the optimal order for each telescope.

    Attributes
    ----------
    event_moc: Table
        a MOC of the event, table with four columns: UNIQ, PROBDENSITY, DISTMEAN
        and DISTSTD
    event_moc_hdr: Meta
        a meta information for the event_moc (i.e. header)
    event_hpx: ndarray
        a HEALPix probability array from event_moc (used for planning)
    event_hdr: Meta
        a meta data for event_hpx, when read from plan_fits file
    event_name: str
        the event name, default empty
    working_directory:
        a root of a directory tree, where all planning data is stored
    plan_filename_template_narrowfield: str
        a template for filename of a target list, contains four placeholders:
        telescope short name, event name, day number and file extension
    plan_filename_template_widefield: str
        the same but widefield telescope sky field lists
    plan_fits_filename: str
        a filename of the HEALPix file, that contains the event skymap painted at
        locations that were planned for observations. It is needed to achieve unique
        field assigning for each scope and next-day planning opportunity
    glade_galaxies_filename: str
        a list of GLADE+ galaxies that left unplanned, will be used in next-day planning
    prob: float
        a cummulative probability threshold for the confidence intervals, that will be
        used in planning
    day: int
        current observation day
    field_blocks: dict[str, list[Field]]
        a mapping of widefield telescope short names to their planned list of fields
    target_blocks: dict[str, list[GladeGalaxy]]
        a mapping of narrow field telescope short names to their planned list of targets
    hpx_plan: ndarray | None
        a HEALPix probability array painted at planned locations

    Methods
    -------
    plan_observations: Callable[list[str], list[str]] -> ndarray
        a main method that performs observation planning with narrow, wide-field
        telescopes or both
    save_plan_fits: Callable[list[str], list[str]]
        save the painted HEALPix array to file
    save_blocks: Callable -> dict[str, str]
        save the lists of targets/sky fields, assign their filenames to telescope short
        names
    reset: Callable
        returns the planner to original condition
    """

    def __init__(
        self,
        event_moc: Table,
        event_name: str = "",
        working_directory: str = ".",
        plan_filename_template_narrowfield: str = "galaxies_{:s}_{:s}_day{:d}.{:s}",
        plan_filename_template_widefield: str = "fields_{:s}_{:s}_day{:d}.{:s}",
        plan_fits_filename: str = "plan.fits",
        glade_galaxies_filename: str = "plan_galaxies.txt",
        prob: float = lvk_uncert_level.value,
    ) -> None:
        self.event_moc = event_moc
        self.event_moc_hdr = event_moc.meta
        self.working_directory = working_directory
        self.plan_filename_template_narrowfield = plan_filename_template_narrowfield
        self.plan_filename_template_widefield = plan_filename_template_widefield
        self.plan_fits_filename = os.path.join(working_directory, plan_fits_filename)
        if os.path.exists(self.plan_fits_filename):
            self.event_hpx, self.event_hdr = read_sky_map(plan_fits_filename, nest=True)
        else:
            self.event_hpx, self.event_hdr = hpx_from_moc(event_moc, nest=True)

        self.glade_galaxies_filename = os.path.join(
            working_directory, glade_galaxies_filename
        )
        self.prob = prob
        self.day = self.event_hdr.get("PLANDAY", 0)
        self.field_blocks: dict[str, list[Field]] = {}
        self.target_blocks: dict[str, list[GladeGalaxy]] = {}
        self.hpx_plan = None
        self.event_name = event_name

    def _load_galaxies(self) -> list[GladeGalaxy]:
        if os.path.exists(self.glade_galaxies_filename) and self.day:
            df = pd.read_csv(self.glade_galaxies_filename, delim_whitespace=True)
            galaxies = self._galaxies_from_df(df)
        else:
            dist_lo = self.event_moc_hdr["distmean"] - self.event_moc_hdr["diststd"]
            dist_hi = self.event_moc_hdr["distmean"] + 2 * self.event_moc_hdr["diststd"]
            galaxies = GladeCatalog.query_skymap_local(
                self.event_moc,
                dist_lo_bound=dist_lo,
                dist_hi_bound=dist_hi,
                prob=self.prob,
            )

        return galaxies

    def _galaxies_to_df(self, galaxies: Sequence[GladeGalaxy]) -> pd.DataFrame:
        df = pd.DataFrame(
            {
                "name": [g.name for g in galaxies],
                "ra": [g.coord.ra.deg for g in galaxies],
                "dec": [g.coord.dec.deg for g in galaxies],
                "Bmag": [g.Bmag for g in galaxies],
                "D_L": [g.D_L for g in galaxies],
            }
        )
        return df

    def _galaxies_from_df(self, df: pd.DataFrame) -> list[GladeGalaxy]:
        return [
            GladeGalaxy(
                name=row["name"],
                coord=SkyCoord(row["ra"] * u.deg, row["dec"] * u.deg),
                Bmag=row["Bmag"],
                z=None,
                D_L=row["D_L"],
            )
            for (i, row) in df.iterrows()
        ]

    def _update_gal_list(self):
        if self.target_dist.objects:
            df = self._galaxies_to_df(self.target_dist.objects)
            df.to_csv(self.glade_galaxies_filename, sep=" ", index=False)

    def _save_plan_fits(self, hpx: np.ndarray, planday: int, **kws: Any):
        hdr: dict[str, Any] = {"PLANDAY": planday}
        if kws:
            hdr.update({k: v for k, v in kws.items() if k == "HISTORY"})

        write_sky_map(self.plan_fits_filename, hpx, nest=True, **hdr)

    def plan_observations(
        self,
        wide_field_telescopes: Sequence[str] | None = None,
        narrow_field_telescopes: Sequence[str] | None = None,
        disable_intersections: bool = True,
    ):
        if wide_field_telescopes and narrow_field_telescopes:
            galaxies = self._load_galaxies()
            self.field_dist = FieldDistributor(self.event_hpx)
            self.target_dist = TargetDistributor(self.field_dist.hpx, galaxies)
            self.field_blocks = self.field_dist.distribute(
                telescope_ids=wide_field_telescopes,
                prob=self.prob,
                disable_intersections=disable_intersections,
            )
            self.target_blocks = self.target_dist.distribute(
                telescope_ids=narrow_field_telescopes,
                disable_intersections=disable_intersections,
            )
            self._update_gal_list()
            self.hpx_plan = self.target_dist.hpx
        elif wide_field_telescopes:
            self.field_dist = FieldDistributor(self.event_hpx)
            self.field_blocks = self.field_dist.distribute(
                telescope_ids=wide_field_telescopes,
                disable_intersections=disable_intersections,
            )
            self.hpx_plan = self.field_dist.hpx
        elif narrow_field_telescopes:
            galaxies = self._load_galaxies()
            self.target_dist = TargetDistributor(self.event_hpx, galaxies)
            self.target_blocks = self.target_dist.distribute(
                telescope_ids=narrow_field_telescopes,
                disable_intersections=disable_intersections,
            )
            self._update_gal_list()
            self.hpx_plan = self.target_dist.hpx
        else:
            raise ValueError("must provide altough a one list of telescopes")

        self.day += 1
        self.event_hpx = self.hpx_plan.copy()
        return self.event_hpx

    def save_plan_fits(
        self,
        wide_field_telescopes: Sequence[str] | None = None,
        narrow_field_telescopes: Sequence[str] | None = None,
    ):
        if self.hpx_plan is not None:
            telescopes = []
            if wide_field_telescopes:
                telescopes += wide_field_telescopes

            if narrow_field_telescopes:
                telescopes += narrow_field_telescopes

            self._save_plan_fits(
                self.hpx_plan, self.day, **self.event_hdr if self.event_hdr else {}
            )
            self._add_history(self.plan_fits_filename, telescopes=telescopes)

    def _add_history(self, filename, telescopes: Sequence[str]):
        with fits.open(filename, mode="update") as hdul:
            hdr = hdul[1].header
            now = datetime.now()
            hdr.append(
                ("HISTORY", f"Obs. plan for day#{self.day} created by AWARE on {now}"),
                end=True,
            )
            hdr.append(("HISTORY", "Telescopes:"), end=True)
            for t in telescopes:
                hdr.append(("HISTORY", t), end=True)

    def _save_block(self, tel_name: str, objects: list[Field | GladeGalaxy]) -> str:
        s = Telescopes[tel_name]
        program = create_observation_program(s, objects)
        log.warning("Observation program for %s is empty", s.full_name)

        if s.fov.is_widefield:
            fname = self.plan_filename_template_widefield.format(
                self.event_name, tel_name, self.day, s.default_target_list_fmt
            )
        else:
            fname = self.plan_filename_template_narrowfield.format(
                self.event_name, tel_name, self.day, s.default_target_list_fmt
            )

        output_dir = os.path.join(self.working_directory, tel_name)
        os.makedirs(output_dir, exist_ok=True)
        output_fname = os.path.join(output_dir, fname)

        with open(output_fname, "w+") as f:
            f.write(program)

        return os.path.realpath(output_fname)

    def save_blocks(self) -> dict[str, str]:
        fname_map: dict[str, str] = {}
        if not self.day:
            return fname_map

        if self.field_blocks:
            for name, fields in self.field_blocks.items():
                fname = self._save_block(name, fields)
                fname_map[name] = fname

        if self.target_blocks:
            for name, targets in self.target_blocks.items():
                fname = self._save_block(name, targets)
                fname_map[name] = fname

        return fname_map

    def reset(self):
        data = rasterize(self.event_moc, order=8)
        prob = data["PROBDENSITY"]
        nside = hp.get_nside(prob)
        area = hp.nside2pixarea(nside)
        self.event_hpx = prob * area
        self.hpx_plan = None
        self.day = 0
        self.target_blocks = {}
        self.field_blocks = {}
