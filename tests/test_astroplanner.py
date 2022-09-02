from __future__ import annotations


from datetime import datetime, timedelta
from pprint import pprint
from typing import Sequence, Union
from astroplan.plots import plot_airmass
from astropy.coordinates import EarthLocation, SkyCoord
from astroplan.observer import Observer
from astroplan import FixedTarget
import astropy.units as u
from astropy.time import Time, TimeDelta
import matplotlib.pyplot as plt
import numpy as np
from aware.planner import mondy
import pandas as pd
from aware.sql import create_session


def test():
    target1 = FixedTarget(
        name="GBM Trigger 683571622", 
        coord=SkyCoord(167.3, 1.01, unit=["deg", "deg"])
    )
    target2 = FixedTarget(
        name="GBM Trigger 683644301", 
        coord=SkyCoord(220.46, 17.23, unit=["deg", "deg"])
    )
    targets = [target1, target2]

    # At 20:00 target2 should be observed first
    time1 = Time("2022-12-10T20:00:00")
    sorted_targets1 = mondy.observation_path(targets, time1)
    assert (
        sorted_targets1[0].name == target2.name 
        and sorted_targets1[1].name == target1.name 
    ), "Targets sorted incorrectly"

    # At 23:59 target1 should be observed first
    time2 = Time("2022-12-10T23:59:00")
    sorted_targets2 = mondy.observation_path(targets, time2)
    assert (
        sorted_targets2[0].name == target1.name 
        and sorted_targets2[1].name == target2.name 
    ), "Targets sorted incorrectly"


if __name__ == "__main__":
    test()

