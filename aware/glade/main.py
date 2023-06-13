from __future__ import annotations

import asyncio
import gc
import os
import sqlite3
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from io import BytesIO
from pathlib import Path
from threading import Lock
import time
from typing import Any, Callable, Literal, Optional

import numpy as np
import pandas as pd
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.io.fits import FITS_record
from astropy.table import Table, vstack
from astropy.io import votable
from astroquery.vizier import Vizier
from ligo.skymap.postprocess.crossmatch import crossmatch
from requests import Response
import requests


from ..config import CfgOption
from ..cosmology import cosmos
from ..logger import log
from ..angle import coord2str, coord_to_target_name


vizier_url = CfgOption("vizier_url", "http://vizier.u-strasbg.fr", str)
row_limit = CfgOption("row_limit", 50, int)
catalog_path = CfgOption("catalog_path", "~/gladep.fits", str)


class CatalogNotFound(Exception):
    pass


class GladeGalaxy(FixedTarget):
    def __init__(
        self,
        coord: SkyCoord,
        Bmag: Optional[float] = None,
        z: Optional[float] = None,
        D_L: Optional[float] = None,
        name: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(coord, name, **kwargs)
        self.Bmag = Bmag
        self.z = z
        self.D_L = D_L


class GladeCatalog:
    id = "VII/291/gladep"
    field_scheme = (
        "GLADE_",
        "PGC",
        "GWGC",
        "HyperLEDA",
        "2MASS",
        "WISExSCOS",
        "SDSS-DR16Q",
        "Type",
        "RAJ2000",
        "DEJ2000",
        "Bmag",
        "e_Bmag",
        "f_Bmag",
        "BMAG",
        "Jmag",
        "e_Jmag",
        "Hmag",
        "e_Hmag",
        "Kmag",
        "e_Kmag",
        "W1mag",
        "e_W1mag",
        "W2mag",
        "e_W2mag",
        "f_W1mag",
        "Bjmag",
        "e_Bjmag",
        "zhelio",
        "zcmb",
        "f_zcmb",
        "e_z",
        "e_zhelio",
        "dL",
        "e_dL",
        "f_dL",
        "M_star",
        "e_M_star",
        "f_M_star",
        "logRate",
        "e_logRate",
    )
    table = None
    coord = None

    @staticmethod
    def lazy_load_catalog() -> tuple(Table, SkyCoord):
        if not getattr(GladeCatalog, "table", None) or not getattr(
            GladeCatalog, "coord", None
        ):
            table = Table.read(
                catalog_path.get_value(),
                format="fits",
                units=["deg", "deg", "Mpc", "mag"],
                memmap=True,
                astropy_native=True,
            )
            GladeCatalog.coord = SkyCoord(
                table["RAJ2000"], table["DEJ2000"], table["dL"]
            )
            GladeCatalog.table = table

        return GladeCatalog.table, GladeCatalog.coord

    @staticmethod
    def query_skymap_local(
        sky_map_moc: Table,
        dist_lo_bound: u.Unit = 0 * u.Mpc,
        dist_hi_bound: u.Unit = np.inf * u.Mpc,
        prob: float = 0.9,
    ) -> list[GladeGalaxy]:
        log.debug("loading local GLADE+ catalog into RAM")
        table, coord = GladeCatalog.lazy_load_catalog()

        log.debug("query GLADE+ catalog against %.1f per cent contour", prob * 100)
        match_res = crossmatch(sky_map_moc, coordinates=coord, contours=(prob,))
        galaxies = [
            GladeGalaxy(
                SkyCoord(ra * u.deg, dec * u.deg),
                Bmag=Bmag,
                z=None,
                D_L=dL,
                name=coord_to_target_name(SkyCoord(ra * u.deg, dec * u.deg)),
            )
            for (ra, dec, dL, Bmag) in table[
                (match_res.searched_prob < prob)
                & (table["dL"] > dist_lo_bound)
                & (table["dL"] < dist_hi_bound)
            ].iterrows()
        ]
        log.debug("%d galaxies loaded from the catalog", len(galaxies))

        # Free memory held by crossmatch function
        gc.collect()

        return galaxies

    @staticmethod
    def query_vizier(
        center_coord: SkyCoord,
        radius: Optional[u.Unit] = None,
        width: Optional[u.Unit] = None,
        height: Optional[u.Unit] = None,
        column_filters: Optional[dict[str, str]] = None,
    ) -> Response:
        if radius is None and (width is None or height is None):
            raise ValueError("Radius or width and height should be specified")

        # Query field in Vizier for objects
        if radius is not None:
            kws = dict(radius=radius)
        else:
            kws = dict(width=width, height=height)

        _lock = Lock()

        def _query_vizier(index: int) -> Table | None:
            with _lock:
                crd = center_coord[index][0]
                log.debug(
                    "Querying Vizier for Glade+ galaxies within the field "
                    "center: %s, W: %.1f deg, H: %.1f deg",
                    crd.to_string("hmsdms"),
                    width.to_value("deg"),
                    height.to_value("deg"),
                )

            resp: Response = Vizier(row_limit=row_limit.get_value()).query_region_async(
                center_coord,
                **kws,
                catalog=GladeCatalog.id,
                column_filters=column_filters if column_filters is not None else {},
            )

            try:
                resp_bytes_io = BytesIO(resp.content)
                resp_bytes_io.seek(0)
            except EOFError:
                log.debug()
                return None

            # resp = GladeCatalog.query_field(
            #     crd,
            #     width=pixel_size * u.deg,
            #     height=pixel_size * u.deg,
            #     column_filters={"dL": f"<={max_dL} && >={min_dL}"},
            # )

            try:
                vo_table: votable.table.tree.Table = votable.parse_single_table(
                    resp_bytes_io
                )
                cat_table = vo_table.to_table()
            except IndexError as e:
                return None
            except requests.exceptions.ConnectTimeout as e:
                time.sleep(60)
                _query_vizier(index)

            with _lock:
                log.debug("%i galaxies retrieved", len(cat_table))

            return cat_table

        # Very extensive process so using thread pool
        with ThreadPoolExecutor(max_workers=10) as ex:
            futures = [ex.submit(_query_vizier, index) for index in idx]
            results = [
                res
                for completed_future in as_completed(futures)
                if (res := completed_future.result()) is not None
            ]
            if not results:
                return None

            table: Table = vstack(results)
            ra = list(table["RAJ2000"])
            dec = list(table["DEJ2000"])
            Bmag = list(table["Bmag"])
            names = [coord_to_target_name(c) for c in SkyCoord(ra, dec)]
            dL = list(table["dL"])
            zz = list(table["zhelio"])

            galaxies = [
                GladeGalaxy(
                    SkyCoord(a * u.deg, d * u.deg), Bmag=B, z=z, D_L=DL, name=str(n)
                )
                for (a, d, B, n, DL, z) in table.iterrows()
            ]

        return galaxies

    @staticmethod
    def query_field_local_txt(
        path: str,
        center_coord: SkyCoord,
        radius: Optional[u.Unit] = None,
        width: Optional[u.Unit] = None,
        height: Optional[u.Unit] = None,
        column_filters: Optional[dict[str, Callable[[Any], Any]]] = None,
        chunk_size: int | None = None,
    ) -> list[GladeGalaxy] | None:
        LINE_SIZE_BYTES = 262
        CHUNK_SIZE = 10_000 * LINE_SIZE_BYTES
        FILE_DESCRIPTER = open(path, "r", encoding="UTF-8")
        FILE_SIZE = os.path.getsize(path)
        N_READS = (FILE_SIZE + CHUNK_SIZE) // CHUNK_SIZE

        lock = Lock()

        def find():
            with lock:
                chunk = BytesIO(
                    bytes(
                        "".join(FILE_DESCRIPTER.readlines(1_000_000)), encoding="UTF-8"
                    )
                )

            chunk.seek(0)
            try:
                df = pd.read_csv(chunk, delim_whitespace=True, header=None)
            except pd.errors.ParserError as e:
                with lock:
                    print(e.args)
                return None

            df.columns = GladeCatalog.field_scheme
            cat_coord = SkyCoord(df["RAJ2000"], df["DEJ2000"], unit=["deg"] * 2)
            idx = cat_coord.separation(center_coord) <= radius

            if np.any(idx):
                t = df[idx]
                with lock:
                    print(t)
            else:
                t = None

            return t

        with ThreadPoolExecutor(max_workers=50) as ex:
            futures = [ex.submit(find) for _ in range(N_READS)]
            tables = [
                r := completed_future.result()
                for completed_future in as_completed(futures)
                if r is not None
            ]

        FILE_DESCRIPTER.close()

        table = pd.concat(tables)

        if len(table):
            table = table[0]
        else:
            # No objects retrieved
            return []

        # Fill masked values with None to make them JSON-serializable
        table = table.filled(None)
        galaxies: list[GladeGalaxy] = []
        for (
            name,
            ra,
            dec,
            Bmag,
            z,
        ) in table.iterrows("GLADE_", "RAJ2000", "DEJ2000", "Bmag", "zhelio"):
            gal = GladeGalaxy(SkyCoord(ra * u.deg, dec * u.deg), Bmag, z, name)
            galaxies.append(gal)

        return galaxies

    @staticmethod
    def query_field_local_db(
        path: str,
        center_coord: SkyCoord,
        radius: Optional[u.Unit] = None,
        width: Optional[u.Unit] = None,
        height: Optional[u.Unit] = None,
        column_filters: Optional[dict[str, Callable[[Any], Any]]] = None,
    ) -> list[GladeGalaxy] | None:
        lock = Lock()

        connection = sqlite3.connect(path)
        cursor = connection.cursor()

        def _sql_cos_query(ang: float | str) -> str:
            # Use taylor expansion up to n=2
            return f"(1 - {ang}*{ang}/2 + {ang}*{ang}*{ang}*{ang}/24)"

        def _sql_sin_query(ang: float | str) -> str:
            # Use taylor expansion up to n=2
            return f"({ang} - {ang}*{ang}*{ang}/6 + {ang}*{ang}*{ang}*{ang}/24)"

        def _angular_distance_query(ra: float, dec: float) -> str:
            one_over_cos_theta = (
                f"({_sql_sin_query(center_coord.dec.deg)}"
                f"*{_sql_sin_query('dec')}"
                f"*{_sql_cos_query(center_coord.dec.deg)}"
                f"*{_sql_cos_query('dec')}"
                f"*{_sql_cos_query(f'({center_coord.ra.deg}-ra)')})"
            )
            theta = f"1 / ({_sql_cos_query(one_over_cos_theta)})"

            return theta

        get_objects_query = "SELECT * FROM galaxies WHERE "
        if radius is None:
            if width is not None:
                if height is None:
                    height = width
                get_objects_query += (
                    f"abs(ra-{center_coord.ra.deg}) <= {width.to_value(u.deg)} "
                    f"AND abs(ra-{center_coord.dec.deg}) <= {height.to_value(u.deg)} "
                )
        else:
            get_objects_query += (
                f"{_angular_distance_query(center_coord.ra.deg, center_coord.dec.deg)} "
                f"<= {radius.to_value(u.deg)}"
            )

        cursor.executemany(get_objects_query)
        result = cursor.fetchall()

        galaxies: list[GladeGalaxy] = [
            GladeGalaxy(SkyCoord(ra * u.deg, dec * u.deg), Bmag, z, name)
            for (name, ra, dec, Bmag, z) in result
        ]

        return galaxies
