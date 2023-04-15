from __future__ import annotations

import json
from argparse import ArgumentParser
from datetime import datetime, timedelta
from pprint import pprint
from time import time
from typing import Any, Sequence, Union

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sqlalchemy
from astroplan import FixedTarget
from astroplan.observer import Observer
from astroplan.plots import plot_airmass
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time, TimeDelta
from aware.config import cfg
from aware.logger import log
from aware.planner import sites
from aware.sql import create_session


def parse_cli_args() -> Sequence[Any]:
    descr = "Given the list of objects, sort them in the optimal way to observe"
    parser = ArgumentParser("sort_targets", description=descr)
    parser.add_argument("-i", "--input", help="Input table of objects")
    parser.add_argument("-o", "--output", help="Output file with sorted objects")
    parser.add_argument("-t", "--time", help="Desired observation time (ISOT)")
    parser.add_argument("-s", "--site", help="Observer (site) id")
    parser.add_argument("--airmass-plot", help="Where to save airmass plot", default="")

    args = parser.parse_args()

    input_filename = args.input
    output_filename = args.output
    t = Time(args.time)
    site = args.site
    airmass_plot = args.airmass_plot

    return input_filename, output_filename, t, site, airmass_plot


def try_read_table(
    filename: str, con: sqlalchemy.engine.Engine, name: str = "alert"
) -> pd.DataFrame:

    try:
        df = pd.read_sql_table(name, con)
    except Exception as e:
        pass
    else:
        return df

    try:
        df = pd.read_csv(filename)
    except Exception as e:
        pass
    else:
        return df

    try:
        df = pd.read_json(filename)
    except Exception as e:
        return df

    raise ValueError(f"Can not get table from file: {filename}")


def main():
    # CLI arguments
    input_filename, output_filename, obs_time, site, airmass_plot = parse_cli_args()
    site = sites[site]

    # Retrieving site by id
    site = sites[site]

    # Read the alert tabel
    try:
        engine, session = create_session(input_filename)
    except Exception as e:
        engine = None

    df = try_read_table(input_filename, engine)
    log.info("Successfully obtained target list from file: %s", input_filename)

    # Extract targets
    targets = [
        FixedTarget(SkyCoord(ra, dec, unit=["deg", "deg"]), event)
        for event, ra, dec in zip(
            df["event"].tolist(), df["ra_center"].tolist(), df["dec_center"].tolist()
        )
    ]

    # Find nearest time window for observations
    log.info("Estimating nearest observational window for site %s", mondy.full_name)
    start, end = site.nearest_observation_window(obs_time)
    log.info("Site can observe since %s until %s", start.datetime, end.datetime)

    # Get list of targets that can be observed by the site
    log.info("Define which targets can be observed at that period")
    obs_targets = site.observable_targets(targets, start, end)

    n_targets = len(targets)
    if obs_targets is not None:
        n_obs = len(targets)
        log.info("From %d targets, %d can be observed", n_targets, n_obs)
    else:
        raise SystemExit("No targets can be observed, quitting the program")

    # Sort the targets at the time of the observation at the center of
    # observation window
    ctr = start + (end - start) / 2
    sorted_targets = site.observation_order(obs_targets, start, end)
    # log.info("All targets sorted at epoch %s (mean time)", ctr.isot)

    if airmass_plot:
        log.info("Drawing airmass plot")
        ax = site.plot_airmass(sorted_targets, ctr, start, end, fig_path=airmass_plot)
        log.info("Saved airmass plot to file: %s", airmass_plot)

    try:
        out = {
            "observatory": {
                "name": site.name,
                "geo_lattitude": site.location.lat.value,
                "geo_longitude": site.location.lon.value,
            },
            "start_time": start.isot,
            "end_time": end.isot,
            "targets": [
                dict(name=target.name, ra=target.ra.value, dec=target.dec.value)
                for target in sorted_targets
            ],
        }
        with open(output_filename, "w+") as f:
            json.dump(out, f, indent=2)
    except Exception as e:
        log.error("output table not saved", exc_info=e)
    else:
        log.info("Saved the sorted list of targets to file: %s", output_filename)


if __name__ == "__main__":
    start_time = time()
    log.info("Started AWARE sort_targets on %s", datetime.now())
    main()
    log.info("Finished in %.3g min %.3g sec", *divmod(time() - start_time, 60))
