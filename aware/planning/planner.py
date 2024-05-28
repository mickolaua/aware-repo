from __future__ import annotations

import asyncio
from io import StringIO
import os
import pickle
import re
import shutil
from base64 import urlsafe_b64decode, urlsafe_b64encode
from collections import UserList
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from functools import reduce
from glob import glob
from textwrap import dedent
from typing import Any, Optional, Sequence
import warnings

import aiofiles
import astropy.units as u
import healpy as hp
import numpy as np
import orjson
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from ligo.skymap.io.fits import read_sky_map, write_sky_map
from ligo.skymap.moc import rasterize
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from pathvalidate import sanitize_filename, sanitize_filepath

from aware.alert.target_info import TargetInfo
from aware.angle import coord_to_target_name
from aware.field import Field
from aware.glade import GladeCatalog, GladeGalaxy
from aware.logger import log
from aware.planning.distributor import (
    FieldDistributor,
    TargetDistributor,
    paint_planned_healpix_from_fields,
    paint_planned_healpix_from_galaxies,
)
from aware.planning.program import create_observation_program, read_observation_program
from aware.site.main import Site, Telescopes
from aware.visualization.main import plot_fields_on_healpix, plot_targets_on_healpix

from ..io import hpx_from_moc
from ..localization.main import Localization, lvk_uncert_level


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
        event_hpx: Optional[np.ndarray] = None,
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
            self.event_hpx, self.event_hdr = read_sky_map(
                self.plan_fits_filename, nest=True
            )
        else:
            # Try hpx attribute first, since it might be more accurate
            if event_hpx is not None:
                self.event_hpx = event_hpx
                self.event_hdr = event_moc.meta
            else:
                self.event_hpx, self.event_hdr = hpx_from_moc(event_moc, nest=True)

        self.glade_galaxies_filename = os.path.join(
            working_directory, glade_galaxies_filename
        )
        self.prob = prob
        self.day = self.event_hdr.get("PLANDAY", 0)
        self.field_blocks: dict[str, list[Field]] = {}
        self.target_blocks: dict[str, list[GladeGalaxy]] = {}
        self.hpx_plan = np.zeros_like(self.event_hpx)
        self.event_name = event_name
        self.field_dist = FieldDistributor(self.event_hpx)
        self.target_dist = TargetDistributor(self.field_dist.hpx, [])
        self._targets_loaded = False

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

    def _galaxies_first_load(self):
        if not self._targets_loaded:
            galaxies = self._load_galaxies()
            self.target_dist.objects.update(galaxies)
            self._targets_loaded = True

    def plan_observations(
        self,
        wide_field_telescopes: Sequence[str] | None = None,
        narrow_field_telescopes: Sequence[str] | None = None,
        disable_intersections: bool = True,
        initial_time: Optional[Time] = None,
    ):
        if wide_field_telescopes and narrow_field_telescopes:
            self._galaxies_first_load()
            self.field_blocks = self.field_dist.distribute(
                telescope_ids=wide_field_telescopes,
                prob=self.prob,
                disable_intersections=disable_intersections,
                initial_time=initial_time,
            )
            self.target_blocks = self.target_dist.distribute(
                telescope_ids=narrow_field_telescopes,
                disable_intersections=disable_intersections,
            )
            self._update_gal_list()
            self.hpx_plan = self.target_dist.hpx
        elif wide_field_telescopes:
            self.field_blocks = self.field_dist.distribute(
                telescope_ids=wide_field_telescopes,
                disable_intersections=disable_intersections,
                initial_time=initial_time,
            )
            self.hpx_plan = self.field_dist.hpx
        elif narrow_field_telescopes:
            self._galaxies_first_load()
            self.target_blocks = self.target_dist.distribute(
                telescope_ids=narrow_field_telescopes,
                disable_intersections=disable_intersections,
                initial_time=initial_time,
            )
            self._update_gal_list()
            self.hpx_plan = self.target_dist.hpx
        else:
            raise ValueError("must provide atleast a single list of telescopes")

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


@dataclass(frozen=True)
class Plan:
    site_id: str
    event_name: str
    targets: Sequence[Field | GladeGalaxy]
    plot: Axes
    epoch: Time
    comment: str = ""
    day: int = 1

    def to_program(self) -> str:
        tel = Telescopes[self.site_id]
        program = create_observation_program(tel, self.targets)
        return program

    def save(self, outdir: str) -> str:
        tel = Telescopes[self.site_id]
        if isinstance(self.targets[0], Field):
            ttype = "fields"
        else:
            ttype = "targets"
        plan_filename = (
            f"{ttype}_{self.event_name}_{tel.name}_day{self.day}"
            f".{tel.default_target_list_fmt}"
        )
        plan_filename_safe = sanitize_filename(
            plan_filename, replacement_text="_", platform="linux"
        )
        program = self.to_program()
        plan_path = os.path.join(outdir, plan_filename_safe)
        with open(plan_path, "w") as f:
            f.write(program)
        return plan_path


# Instanciate PlanCollection from a UserList to get rid of problems with overriding the
# default methods
class PlanCollection(UserList):
    def __init__(self, plans: Optional[Sequence[Plan]] = None) -> None:
        if plans is not None:
            for p in plans:
                if not isinstance(p, Plan):
                    raise TypeError(
                        "plans must contain only Plan instances, but found "
                        f"{type(p).__name__}"
                    )
        else:
            plans = []
        super().__init__(plans)
        self.data = list(plans)

    def __contains__(self, plan: Plan) -> bool:
        return plan in self.data

    def exists(self, plan: Plan) -> bool:
        return plan in self

    def append(self, __object: Plan) -> None:
        if not isinstance(__object, Plan):
            raise TypeError(
                f"object must be a Plan instance, got {type(__object).__name__}"
            )
        if not self.exists(__object):
            return self.data.append(__object)

    def __add__(self, other: Sequence[Plan]) -> PlanCollection:
        if isinstance(other, UserList):
            new_data = self.data + other.data
        else:
            new_data = self.data + list(other)
        return PlanCollection(new_data)

    @property
    def length(self):
        return len(self.data)

    def save(self, outdir: str) -> list[str]:
        filenames = []
        for plan in self:
            tel_dir = os.path.join(outdir, plan.site_id)
            os.makedirs(tel_dir, exist_ok=True)
            filename = plan.save(tel_dir)
            filenames.append(filename)
        return filenames


DEFAULT_JSON_OPTIONS = (
    # Sort keys in alphabetical order
    orjson.OPT_SORT_KEYS,
    # Use 2 space indentation
    orjson.OPT_INDENT_2,
    # Add Z to the end of UTC timestamps
    orjson.OPT_UTC_Z,
    # Omit microseconds in UTC timestamps
    orjson.OPT_OMIT_MICROSECONDS,
    # Use naive UTC timestamps, i.e. tzinfo is omitted
    orjson.OPT_NAIVE_UTC,
    # Serialize numpy arrays
    orjson.OPT_SERIALIZE_NUMPY,
)


class ObservationPlanner:
    """
    Observation planner creates plans for all kinds of skymaps.

    Attributes
    ----------
    target_info: TargetInfo
        a target information
    working_directory: str
        a root directory, where plan-related products are stored
    queue: Queue:
        a queue used to send plans to the receiver
    plans: PlanCollection:
        an object that stores plans for each telescope

    Methods
    -------
    plan_observations:
        a main method to create observation plans
    """

    def __init__(
        self,
        target_info: TargetInfo,
        working_directory: str,
        prob_threshold: float = 0.9,
    ) -> None:
        if not hasattr(target_info, "localization"):
            raise ValueError("TargetInfo must provide a localization")

        if target_info.localization is None:
            raise ValueError("Localization must be not empty")

        try:
            self.moc: Table = target_info.localization.data
        except AttributeError as e:
            raise AttributeError(
                "Localization sky map does not have a MOC table"
            ) from e

        self.target_info = target_info
        self.working_directory = working_directory
        self.queue = asyncio.Queue()
        self.plans: PlanCollection = PlanCollection()
        self.filenames: list[str] = []
        self.day = 0
        self.skymap_planner = SkymapPlanner(
            self.moc,
            event_hpx=self.localization.hpx,
            event_name=self.target_name,
            working_directory=self.event_directory,
        )
        self.prob_threshold = prob_threshold

    @property
    def target_name(self) -> str:
        return f"{self.target_info.origin} {self.target_info.event}"

    @property
    def localization(self) -> Localization:
        return self.target_info.localization

    def add_field_plot_title(self, ax: Axes, block: list[Field], tel: Site):
        ax.set_title(
            dedent(
                f"""
                Mosaic planning
                Event: {self.target_name}
                Telescope: {tel.full_name}
                Fields: {len(block):d}
                """
            )
        )

    def add_target_plot_title(self, ax: Axes, block: list[GladeGalaxy], tel: Site):
        ax.set_title(
            dedent(
                f"""
                Target observations
                Event: {self.target_name}
                Telescope: {tel.full_name}
                Targets: {len(block):d}
                """
            )
        )

    @property
    def target_name_url_safe(self) -> str:
        encoded_name = self.target_name.encode("utf-8")
        safe_name = urlsafe_b64encode(encoded_name).decode("utf-8")
        return safe_name

    def extract_target_name_from_path(self, path: str) -> str:
        if os.path.isfile(path):
            folder = os.path.dirname(path)
        else:
            folder = path

        encoded_folder = folder.encode("utf-8")
        target_name = urlsafe_b64decode(encoded_folder).decode("utf-8")
        return target_name

    @property
    def event_directory(self) -> str:
        # Here, we make path URLs safe, because event names may contain such symbols
        # as slashes `/`, pluses `+`, equal signs `=`, etc.
        outdir = os.path.join(
            self.working_directory,
            self.target_name_url_safe,
        )
        return outdir

    def telescope_directory(self, tel: Site) -> str:
        outdir = os.path.join(
            self.event_directory,
            tel.name,
        )
        return outdir

    def can_observe(self, tel: Site, epoch: Time) -> bool:
        start, stop = tel.nearest_observation_window(epoch)
        return start is not None and stop is not None and (start.jd - epoch.jd) < 1

    def is_observable_single_scan(self, tel: Site, loc: Localization) -> bool:
        return tel.fov.radius >= 0.9 * loc.error_radius()

    def select_single_scan_telescopes(self, telescopes: list[Site]) -> list[Site]:
        return [
            t
            for t in telescopes
            if self.is_observable_single_scan(t, self.localization)
        ]

    def plan_single_scan(self, tel: Site, epoch: Time, day: int, comment: str) -> Plan:
        coord = self.localization.center()
        name = self.target_name
        start_time, end_time = tel.nearest_observation_window(epoch)
        field = Field(coord, width=tel.fov.width, height=tel.fov.height, name=name)
        obs_fields = tel.observable_targets(
            [field], start_time=start_time, end_time=end_time,
        )
        ax = tel.plot_airmass(
            obs_fields,
            start_time=start_time,
            end_time=end_time,
            plot_kws={"ls": "-", "marker": "None"},
            fig_path="",
        )
        plan = Plan(
            tel.name,
            name,
            targets=obs_fields,
            plot=ax,
            epoch=epoch,
            comment=comment,
            day=day,
        )
        return plan

    def set_default_dist_constraints(self, hdr: dict[str, Any]):
        """
        Set default constraints on distances for such events as gamma-ray bursts
        """
        keys = {k.lower() for k in hdr.keys()}
        if not "distmean" in keys or not "diststd" in keys:
            hdr["distmean"] = 0
            hdr["diststd"] = 0

    def has_finite_dist_bounds(self) -> bool:
        """
        Check if there are finite distance bounds given in the skymap metadata
        """
        dist = self.skymap_planner.event_hdr.get("distmean")
        return dist > 0

    def plan_for_single_scans(
        self, telescopes: list[Site], epoch: Time, day: int, comment: str
    ) -> PlanCollection:
        plans = [self.plan_single_scan(tel, epoch, day, comment) for tel in telescopes]
        plans = [p for p in plans if p.targets]
        return PlanCollection(plans)

    def plan_for_multiple_scans(
        self,
        wide_field_telescopes: list[Site],
        narrow_field_telescopes: list[Site],
        epoch: Time,
        day: int,
        prob_threshold: float,
        disable_intersections: bool,
        narrow_field_comment: str = "",
        wide_field_comment: str = "",
    ) -> PlanCollection:
        plans = PlanCollection([])
        moc = self.moc
        skymap_planner = self.skymap_planner
        self.set_default_dist_constraints(moc.meta)
        wide_field_ids = [t.name for t in wide_field_telescopes]

        # Disable target planning if no distance constraints are specified
        if self.has_finite_dist_bounds():
            narrow_fields_ids = [t.name for t in narrow_field_telescopes]
        else:
            narrow_fields_ids = []

        planned_hpx = skymap_planner.plan_observations(
            wide_field_telescopes=wide_field_ids,
            narrow_field_telescopes=narrow_fields_ids,
            disable_intersections=disable_intersections,
            initial_time=epoch,
        )

        # Create mosaic plans
        for tel_id, block in skymap_planner.field_blocks.items():
            tel = Telescopes[tel_id]
            ax = plot_fields_on_healpix(
                skymap_planner.event_hpx,
                block,
                self.localization.center(),
                self.localization.error_radius() / 1.6,
                prob_threshold,
            )
            self.add_field_plot_title(ax, block, tel)
            plan = Plan(
                tel_id,
                self.target_name,
                targets=block,
                comment=wide_field_comment,
                plot=ax,
                epoch=epoch,
                day=day,
            )
            plans.append(plan)

        # Create target plans
        for tel_id, block in skymap_planner.target_blocks.items():
            tel = Telescopes[tel_id]
            ax = plot_targets_on_healpix(
                skymap_planner.event_hpx,
                block,
                self.localization.center(),
                self.localization.error_radius() / 1.6,
                prob_threshold,
            )
            self.add_target_plot_title(ax, block, tel)
            plan = Plan(
                tel_id,
                self.target_name,
                targets=block,
                comment=narrow_field_comment,
                plot=ax,
                epoch=epoch,
                day=day,
            )
            plans.append(plan)

        plans = PlanCollection([p for p in plans if p.targets])
        return plans

    def plan_observations(
        self,
        epoch: Time = Time(datetime(2000, 1, 1)),
        wide_field_telescopes: list[Site] | None = None,
        narrow_field_telescopes: list[Site] | None = None,
        prob_threshold: float = lvk_uncert_level.value,
        disable_intersections: bool = False,
        day: int = 0,
        narrow_field_comment: str = "",
        wide_field_comment: str = "",
        single_scan_comment: str = "",
    ) -> PlanCollection:
        if not wide_field_telescopes and not narrow_field_telescopes:
            raise ValueError("List of wide or narrow field telescopes must be provided")

        # Create event directory
        os.makedirs(self.event_directory, exist_ok=True)

        # Plan day, if not specified, plan for a next day
        if not day:
            day = self.day + 1

        self.prob_threshold = prob_threshold

        # Determine telescope that actually can observe
        telescopes = []
        if wide_field_telescopes:
            telescopes += wide_field_telescopes
        if narrow_field_telescopes:
            telescopes += narrow_field_telescopes
        good_telescopes = [t for t in telescopes if self.can_observe(t, epoch)]

        # Remove telescopes for which planning is already presented
        data = self.load_planning()
        for json_plan in data.get("plans", []):
            for i, tel in enumerate(good_telescopes):
                if tel.name == json_plan["site_id"] and day == json_plan["day"]:
                    good_telescopes.remove(tel)

        # Select those who can observe the localization in one scan
        single_scan_telescopes = self.select_single_scan_telescopes(good_telescopes)
        if single_scan_telescopes:
            single_scan_plans = self.plan_for_single_scans(
                single_scan_telescopes, epoch, day, single_scan_comment
            )
        else:
            single_scan_plans = PlanCollection([])

        # Select those telescopes, who can observe in series of scans
        multi_scan_telescopes = [
            t for t in good_telescopes if t not in single_scan_telescopes
        ]
        multi_wide_field_telescopes = [
            t for t in multi_scan_telescopes if t in wide_field_telescopes
        ]
        multi_narrow_field_telescopes = [
            t for t in multi_scan_telescopes if t in narrow_field_telescopes
        ]

        if multi_narrow_field_telescopes or multi_wide_field_telescopes:
            multi_scan_plans = self.plan_for_multiple_scans(
                multi_wide_field_telescopes,
                multi_narrow_field_telescopes,
                epoch,
                day,
                prob_threshold,
                disable_intersections,
                narrow_field_comment=narrow_field_comment,
                wide_field_comment=wide_field_comment,
            )
        else:
            multi_scan_plans = PlanCollection([])

        # Gather all plans
        plans = single_scan_plans + multi_scan_plans

        # Remove planned galaxies
        if not disable_intersections:
            # Paint the healpix already planned
            for plan in multi_scan_plans:
                if plan.targets:
                    if isinstance(plan.targets[0], Field):
                        ...
                        # self.planner.hpx_plan = paint_planned_healpix_from_fields(
                        #     self.planner.hpx_plan, plan.targets
                        # )
                        # self.planner.event_hpx = self.planner.hpx_plan
                        # self.planner.field_dist.hpx = self.planner.hpx_plan
                    else:
                        for target in plan.targets:
                            try:
                                self.skymap_planner.target_dist.objects.remove(target)
                            except LookupError:
                                ...
                        self.skymap_planner._update_gal_list()

        self.plans = plans
        self.skymap_planner.day = day
        self.day = day

        if not len(plans):
            log.warn(
                "Plan collection is empty. Probably, plans already written to disk.\n"
                "Inspect event directory at %s",
                self.event_directory,
            )

        return plans

    def save_or_update_planning(
        self,
        json_options: Sequence[int] = DEFAULT_JSON_OPTIONS,
    ) -> tuple[str, dict[str, Any]]:
        if not self.plans:
            log.warn("No plans to save. Please, run plan_observations() before saving.")
            # return {}

        outdir = self.event_directory
        os.makedirs(outdir, exist_ok=True)
        json_registry_path = os.path.join(outdir, "plan.json")
        if not os.path.exists(json_registry_path) and not self.plans:
            return json_registry_path, {}

        # Save only not empty plans
        not_empty_plans = PlanCollection([p for p in self.plans if p.targets])
        filenames = not_empty_plans.save(outdir)

        hpx_filepath = os.path.join(outdir, f"plan_{self.day}.fits")
        if not os.path.exists(hpx_filepath):
            hp.write_map(
                hpx_filepath,
                self.skymap_planner.hpx_plan,
                nest=True,
                dtype=np.float64,
                overwrite=True,
            )

        json_registry_path = os.path.join(outdir, "plan.json")
        data = {
            "event": self.target_name,
            "plans": [],
        }
        if os.path.exists(json_registry_path):
            with open(json_registry_path, "rb") as json_registry:
                try:
                    new_data = orjson.loads(json_registry.read())
                    data.update(new_data)
                except orjson.JSONEncodeError as e:
                    ...

        data["planday"] = self.day
        data["skymap"] = hpx_filepath
        data["galaxies"] = self.skymap_planner.glade_galaxies_filename
        data["prob"] = self.prob_threshold
        data["fields"] = []

        for plan, fn in zip(not_empty_plans, filenames):
            already_presented = False
            for json_plan in data["plans"]:
                if (
                    json_plan["site_id"] == plan.site_id
                    and json_plan["day"] == plan.day
                ):
                    already_presented = True
                    continue

            if not already_presented:
                tel_dir = os.path.join(outdir, plan.site_id)
                os.makedirs(tel_dir, exist_ok=True)
                fn = plan.save(tel_dir)

                plot_name = f"plot_{self.target_name}_day{self.day}.png"
                plot_name_safe = sanitize_filename(
                    plot_name,
                    replacement_text="_",
                    platform="linux",
                )
                plot_fn = os.path.join(tel_dir, plot_name_safe)
                fig = plan.plot.get_figure()
                fig.savefig(plot_fn, bbox_inches="tight", pad_inches=0.05)

                data["plans"].append(
                    {
                        "site_id": plan.site_id,
                        "day": plan.day,
                        "file": fn,
                        "epoch": plan.epoch.to_datetime(),
                        "plot_file": plot_fn,
                        "comment": plan.comment,
                    }
                )

                for site_id, fields in self.skymap_planner.field_dist.fields.items():
                    for f in fields:
                        dict_field = {
                            "ra": float(f.coord.ra.deg),
                            "dec": float(f.coord.dec.deg),
                            "site_id": site_id,
                        }
                        data["fields"].append(dict_field)

        with open(json_registry_path, "wb") as json_registry:
            json_registry.write(
                orjson.dumps(data, option=reduce(lambda a, b: a | b, json_options))
            )

        return json_registry_path, data

    def load_planning(self) -> dict[str, Any]:
        data = {}
        outdir = self.event_directory
        json_registry_path = os.path.join(outdir, "plan.json")
        if os.path.exists(json_registry_path):
            with open(json_registry_path, "rb") as json_registry:
                data.update(orjson.loads(json_registry.read()))

        try:
            hpx = hp.read_map(data["skymap"], nest=True)
            self.skymap_planner.hpx_plan = hpx
        except:
            hpx = None

        data_plans = data.get("plans", [])
        prob = data.get("prob", 0.9)
        plans = PlanCollection([])
        for data_plan in data_plans:
            site_id = data_plan["site_id"]
            fn = data_plan["file"]
            tel = Telescopes[site_id]
            epoch = data_plan["epoch"]
            bn = os.path.basename(fn)
            if bn.startswith("fields"):
                ttype = Field
                ttype_kws = {"width": tel.fov.width, "height": tel.fov.height}
                # if hpx:

                #     def plot_func(objects) -> Axes:
                #         ax = plot_fields_on_healpix(
                #             hpx,
                #             objects,
                #             self.localization.center(),
                #             self.localization.error_radius(),
                #             prob,
                #         )
                #         self.add_plot_title(ax, objects, tel)

                # else:
                #     plot_func = lambda objects: None
            else:
                ttype = GladeGalaxy
                ttype_kws = {"dL": None, "z": None, "Bmag": None}
                # start_time, end_time = tel.nearest_observation_window(Time(epoch))
                # plot_func = lambda objects: tel.plot_airmass(
                #     objects, start_time=start_time, end_time=end_time, fig_path=None
                # )

            objects = read_observation_program(fn, object_type=ttype, **ttype_kws)

            # TODO: Restore plots from JSON registry
            plot_func = lambda objects: None

            ax = plot_func(objects)
            plan = Plan(
                site_id=site_id,
                event_name=self.target_name,
                targets=objects,
                plot=ax,
                epoch=epoch,
                day=data["planday"],
            )
            plans.append(plan)

        self.plans = plans
        self.day = int(data.get("planday", 1))
        self.skymap_planner.day = self.day

        return data

    def _clear_queue(self):
        while not self.queue.empty():
            self.queue.get_nowait()
            self.queue.task_done()

    def reset(self):
        self._clear_queue()
        self.plans.clear()
        self.day = 0
        self.filenames = []

    def clear_plans(self):
        log.warn("CAUTION: This operation may harm all of your data!")
        agreement = input("Are you sure you want to continue? (y/n): ")
        if agreement.lower().strip() == "y":
            event_dir = self.event_directory
            for dirpath, dirnames, filenames in os.walk(event_dir):
                for fname in filenames:
                    # Match the filename by its basename, but remove by the full path
                    fullpath = os.path.join(dirpath, fname)

                    # Remove HEALPix skymap used for planning for a given day
                    if re.match("^plan_[0-9]+\.fits?$", fname):
                        os.remove(fullpath)

                    # Remove JSON planning registry
                    if re.match("^plan*.json$", fname):
                        os.remove(fullpath)

                # Remove telescope planning lying in separate directories
                for dirname in dirnames:
                    fullpath = os.path.join(dirpath, dirname)
                    if dirname in Telescopes:
                        shutil.rmtree(fullpath)
            log.info("Plan directory tree at %s removed", event_dir)
        else:
            log.info("No plan directory removed: you cancelled the operation.")
