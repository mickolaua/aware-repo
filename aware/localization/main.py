"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
localization.py (c) 2023
Desc: Localization maps
Created:  2023-02-05
Modified: 2023-03-15
"""
from __future__ import annotations
import asyncio
from contextlib import suppress

import functools
from io import BytesIO
import sys
from textwrap import dedent
import time
from typing import Any
import aiohttp

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import pandas as pd
from requests import Response
import requests

try:
    from mocpy import MOC, World2ScreenMPL
except ImportError:
    from mocpy import MOC, WCS as World2ScreenMPL

from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

import astropy_healpix as ah
import matplotlib as mpl
import numpy as np
import numpy.typing as npt
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import Galactic, SkyCoord
from astropy.time import Time
from astropy.table import vstack

# TODO: Implement these LIGO/Virgo functions for Windows
if sys.platform.startswith("win"):
    raise SystemExit("Sorry, but Windows is not currently supported.") 

from ligo.skymap.moc import rasterize
from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.skymap.plot import mellinger
from ligo.skymap.postprocess import find_greedy_credible_levels

from ..config import CfgOption
from ..glade import GladeCatalog, GladeGalaxy
from ..logger import log
from ..site import Site
from ..visualization import plot_moc, plot_glade_galaxies

mpl.use("agg")


__all__ = ["Localization", "CircularSkyMap", "LVCSkyMap"]


# Options
coord_frame = CfgOption("coord_frame", "icrs", str)
lvk_uncert_level = CfgOption("lvk_uncert_level", 0.9, float)
max_healpix_resolution = CfgOption("max_healpix_resolution", 9, int)
healpix_resolution_step = CfgOption("healpix_resolution_step", 2, int)
max_workers = CfgOption("max_workers", 100, int)

# Constants
TOTAL_SKY_AREA = 41252.961 * u.Unit("deg2")


class Localization:
    def __init__(self, **kwargs) -> None:
        for k, v in kwargs.items():
            setattr(self, k, v)

    def describe(self) -> str:
        raise NotImplementedError

    def area(self) -> u.Unit:
        raise NotImplementedError

    def to_dict(self) -> dict[str, Any]:
        raise NotImplementedError

    def center(self) -> SkyCoord:
        raise NotImplementedError

    def error_radius(self) -> u.Unit:
        raise NotImplementedError

    @classmethod
    def from_dict(cls, __dict: dict[str, Any]) -> Localization:
        raise NotImplementedError

    def moc(
        self,
        max_depth: int = max_healpix_resolution.get_value(),
        delta_depth: int = healpix_resolution_step.get_value(),
    ) -> MOC:
        raise NotImplementedError

    def observe(
        self, site: Site, start: Time, stop: Time
    ) -> tuple[list[FixedTarget], str]:
        """Retrieve list of targets from localization region that should be
        observed as well as commentary on observation strategy.

        Parameters
        ----------
        site : Site
            a telescope which will observe targets
        start : Time
            a start time of the observational window
        stop : Time
            an end time of the observational window

        Returns
        -------
        tuple[list[FixedTarget], str]
            a list of targets, commentary on observation strategy

        Raises
        ------
        NotImplementedError
            should be implemented in subclass
        """
        raise NotImplementedError

    def plot(
        self, site: Site, start: Time, stop: Time, targets: list[FixedTarget]
    ) -> Axes:
        raise NotImplementedError


class CircularSkyMap(Localization):
    def center(self) -> SkyCoord:
        return SkyCoord(self.ra_center, self.dec_center, unit="deg")

    def error_radius(self) -> u.Unit:
        return self.radius * u.deg

    def moc(
        self,
        max_depth: int = max_healpix_resolution.get_value(),
        delta_depth: int = healpix_resolution_step.get_value(),
    ) -> MOC:
        """Generate Multi-Order Coverage map from center and radius."""
        return MOC.from_cone(
            self.center().ra,
            self.center().dec,
            self.error_radius() * u.deg,
            max_depth,
            delta_depth=delta_depth,
        )

    def area(self) -> u.Unit:
        return np.pi * self.error_radius() ** 2

    def describe(self) -> str:
        gal_coord = self.center().transform_to(Galactic)
        ra_hour = "{:0.0f}:{:.0f}:{:.4f}".format(*self.center().ra.hms)

        # Handle lattitudes separately since they has - sign in each Hour, Min, Sec
        # value, but only signed hours needed
        dd, dm, ds = self.center().dec.dms
        dec_dms = f"{dd:+0.0f} {abs(dm):.0f} {abs(ds):.4f}"

        r_deg = self.error_radius().value
        if r_deg * 3600 < 30:
            r = f"{r_deg * 3600:.3g} arcsec"
        elif r_deg * 3600 > 30 and r_deg * 60 < 30:
            r = f"{r_deg * 60:.3g} arcmin"
        else:
            r = f"{r_deg:.3g} deg"

        return dedent(
            f"""
            Area: {self.area().value:.3g} deg^2
            RA={self.center().ra.deg:04f} d ({ra_hour}, J2000)
            Dec={self.center().dec.deg:+04f} d ({dec_dms}, J2000)
            l={gal_coord.l.deg:.04f} d
            b={gal_coord.b.deg:+.04f} d
            Error r={r}"""
        )

    async def observe(
        self, site: Site, start: Time, stop: Time
    ) -> tuple[list[FixedTarget], str]:
        sorted_targets = []
        comment = f"Can not cover localization region with {site.full_name}"
        w = 0.5 * site.fov.width.to_value("deg")
        h = 0.5 * site.fov.height.to_value("deg")

        # Note, that 10% of FOV area may be affected with artefacts due to telescope
        # optics or telescope guiding issues
        if np.hypot(w, h) >= 0.9 * self.error_radius().to_value("deg"):
            target = FixedTarget(self.center())
            obs_targets = site.observable_targets([target], start, stop)
            sorted_targets = site.observation_order(
                obs_targets, start_time=start, end_time=stop
            )
            comment = "Aim telescope at the center of the localization region"

        return sorted_targets, comment

    def plot(
        self, site: Site, start: Time, stop: Time, targets: list[FixedTarget]
    ) -> Axes:
        ax: Axes = site.plot_airmass(
            targets, start_time=start, end_time=stop, fig_path=""
        )
        return ax


class LVCSkyMap(Localization):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        level, ipix = ah.uniq_to_level_ipix(self.uniq)
        self.level = level
        self.ipix = ipix
        self.nside = ah.level_to_nside(level)

        order = 8
        data = self._rasterize_multiordered_healpix(order=order)
        nside = ah.level_to_nside(order)
        area = ah.nside_to_pixel_area(nside)
        self.hpx = data["PROBDENSITY"] * area.to_value("sr")

    @functools.lru_cache
    def lon_lat(self) -> SkyCoord:
        lon, lat = ah.healpix_to_lonlat(self.ipix, self.nside, order="nested")
        return SkyCoord(lon, lat, unit="deg")

    def center(self) -> SkyCoord:
        coord = self.lon_lat()
        ra_center = (coord.ra.deg.max() - coord.ra.deg.min()) / 2
        dec_center = (coord.dec.deg.max() - coord.dec.deg.min()) / 2
        return SkyCoord(ra_center, dec_center, unit="deg")

    def error_radius(self) -> u.Unit:
        coord = self.lon_lat()
        dlon = coord.ra.deg.max() - coord.ra.deg.min()
        dlat = coord.dec.deg.max() - coord.dec.deg.min()
        return max(dlon, dlat) / 2 * u.deg

    def area(self, uncert_level: float = lvk_uncert_level.get_value()) -> u.Unit:
        pixel_areas = ah.nside_to_pixel_area(self.nside).to("sr")
        prob = self.probdensity * pixel_areas.value
        mask = np.cumsum(prob) > uncert_level
        a = pixel_areas[mask].to_value("deg2").sum()
        return a * u.deg**2

    def describe(self) -> str:
        return dedent(
            f"""
            Area = {self.area().to_value("deg2"):.1f} deg^2
            P_BNS = {self.bns_prob:.1f}
            P_NSBH ={self.nsbh_prob:.1f}
            P_BBH = {self.bbh_prob:.1f}
            P_Terr = {self.terr_prob:.1f}
            P_hasNS = {self.has_ns_prob:.1f}
            P_hasRemnant = {self.has_remnant_prob:.1f}
            H/W injection = {self.hw_inject:d}
            GraceDB URL: {self.event_page:s}
            """
        )

    def moc(
        self,
        max_depth: int = max_healpix_resolution.get_value(),
        delta_depth: int = healpix_resolution_step.get_value(),
    ) -> MOC:
        """Generate Multi-Order Coverage map from LVK HEALPix.

        Solution was taken from mocpy examples:
        https://cds-astro.github.io/mocpy/examples/examples.html#gravitational-waves-mocs

        Returns
        -------
        moc : `class:mocpy.moc.MOC`
            multi-coverage map
        """
        area = ah.nside_to_pixel_area(self.nside).to_value(u.sr)
        prob = self.probdensity * area
        moc = MOC.from_valued_healpix_cells(
            self.uniq, prob, max_depth=max_depth, cumul_to=lvk_uncert_level.get_value()
        )
        return moc

    def _rasterize_multiordered_healpix(self, order: int | None = 8):
        hpx = rasterize(self.data, order=order)
        return hpx

    async def observe(self, site: Site, start: Time, stop: Time) -> list[FixedTarget]:
        # For wide-field scopes plan the mosaic observations and for narrow-field 
        # target observations of GLADE+ galaxies
        if site.fov.is_widefield:
            from ..planner.mosaic import mosaic_walker
            sorted_targets = mosaic_walker(self.hpx, site.fov)
            comment = "Observe these sky fields in a given order."

            return sorted_targets, comment
        else:
            data = self._rasterize_multiordered_healpix()
            coord = self.lon_lat()
            pixel_area = ah.level_to_nside(8)
            prob = self.probdensity * pixel_area#.to_value("sr")
            i = np.flipud(np.argsort(prob))
            sorted_credible_levels = np.cumsum(prob[i])
            credible_levels = np.empty_like(sorted_credible_levels)
            credible_levels[i] = sorted_credible_levels

            idx = np.argwhere(credible_levels <= lvk_uncert_level.get_value())
            pixel_size = np.sqrt(pixel_area)

            # Query Vizier for Glade+ galaxies
            galaxies: list[GladeGalaxy] = []
            log.debug(
                "Localization region is consist of %i pixels", 
                ah.nside_to_npix(ah.level_to_nside(8))
            )

            # _lock = Lock()
            # from astropy.io import votable
            # from astropy.table import Table

            max_dL = self.distmu.to_value("Mpc") + 2 * self.distsigma.to_value("Mpc")
            min_dL = self.distmu.to_value("Mpc") - self.distsigma.to_value("Mpc")

            # def _query_vizier(index: int) -> Table | None:
            #     with _lock:
            #         crd = coord[index][0]
            #         log.debug(
            #             "Querying Vizier for Glade+ galaxies within the field "
            #             "center: %s, W: %.1f deg, H: %.1f deg",
            #             crd.to_string("hmsdms"),
            #             pixel_size,
            #             pixel_size,
            #         )

            #     resp = GladeCatalog.query_field(
            #         crd,
            #         width=pixel_size * u.deg,
            #         height=pixel_size * u.deg,
            #         column_filters={"dL": f"<={max_dL} && >={min_dL}"},
            #     )
            #     try:
            #         resp_bytes_io = BytesIO(resp.content)
            #         resp_bytes_io.seek(0)
            #     except EOFError:
            #         return None

            #     try:
            #         vo_table: votable.table.tree.Table = votable.parse_single_table(
            #             resp_bytes_io
            #         )
            #         cat_table = vo_table.to_table()
            #     except IndexError as e:
            #         return None
            #     except requests.exceptions.ConnectTimeout as e:
            #         time.sleep(60)
            #         _query_vizier(index)

            #     with _lock:
            #         log.debug("%i galaxies retrieved", len(cat_table))

            #     return cat_table

            # # Very extensive process so using thread pool
            # with ThreadPoolExecutor(max_workers=10) as ex:
            #     futures = [ex.submit(_query_vizier, index) for index in idx]
            #     results = [
            #         res.to_pandas()
            #         for completed_future in as_completed(futures)
            #         if (res := completed_future.result()) is not None
            #     ]
            #     if not results:
            #         return (
            #             [], 
            #             f"Can not cover localization region with {site.full_name}"
            #         )
            #     table = pd.concat(results)

            # ra = list(table["RAJ2000"])
            # dec = list(table["DEJ2000"])
            # Bmag = list(table["Bmag"])
            # names = list(table["GLADE_"])
            # dL = list(table["dL"])
            # zz = list(table["zhelio"])

            # galaxies = [
            #     GladeGalaxy(
            #         SkyCoord(a*u.deg, d*u.deg), 
            #         Bmag=B, 
            #         z=z,
            #         D_L=DL,
            #         name=str(n)
            #     ) for a, d, B, n, DL, z in zip(ra, dec, Bmag, names, dL, zz)
            # ]

            galaxies = GladeCatalog.query_local(
                self.data, dist_lo_bound=min_dL, dist_hi_bound=max_dL
            )

            # Since the number of targets is very large, just sort them
            # firstly by R.A., then by Dec
            obs_targets = site.observable_targets(galaxies, start, stop)
            sorted_targets = site.observation_order(
                obs_targets, start_time=start, end_time=stop, method="radec"
            )
            if sorted_targets:
                frac = len(sorted_targets) / len(galaxies)
                if frac < 1.0:
                    comment = (
                        "Aim telescope at these Glade+ galaxies. Note: "
                        f"telescope covers only {frac*100}% of localization region"
                    )
                else:
                    comment = "Aim telescope at these Glade+ galaxies."
            else:
                comment = "Can not cover localization region with " f"{site.full_name}"

            return sorted_targets, comment

    def plot(
        self, site: Site, start: Time, stop: Time, targets: list[FixedTarget]
    ) -> Axes:
        c = find_greedy_credible_levels(self.hpx)

        fig = plt.figure(figsize=(10, 6), dpi=100)
        ax = plt.axes(
            [0.05, 0.05, 0.9, 0.9],
            projection='astro hours mollweide'
        )
        cs = ax.contour_hpx(
            (100 * c, 'ICRS'), 
            nested=True,
            colors='k', 
            linewidths=2,
            levels=[90],
        )
        plt.clabel(cs, fontsize=6, inline=True)

        ax.grid(True)
        ax.imshow_hpx(
            (self.hpx, "icrs"), nested=True, cmap='cylon'
        )
        ax.scatter(
            [t.coord.ra.deg for t in targets] * u.deg, 
            [t.coord.dec.deg for t in targets]* u.deg,
            transform=ax.get_transform("icrs")
        )
        # ax.set_title(f"{table.meta['OBJECT']} localization map")

        

        # ax = plot_moc(self.moc(), projection="MOL", fig_kws={"figsize": (12, 6)})

        scatter_kws = dict(label="GLADE+ galaxies", cmap="viridis_r", edgecolor="k")
        colorbar_kws = dict(fraction=0.1, shrink=0.8)
        # ax = plot_glade_galaxies(
        #     targets, ax, scatter_kws=scatter_kws, colorbar_kws=colorbar_kws
        # )

        # cax.tick_params(labelsize=11)
        # cax.set_ylabel("Redshift (z)", rotation=90, labelpad=15)
        # ax.legend()

        ax.set_xlabel("R.A.")
        ax.set_ylabel("Dec")
        ax.set_title("Localization map")

        return ax
