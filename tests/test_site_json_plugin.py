"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_site_json_plugin.py (c) 2024
Desc: test loading of Site plugin presented in JSON file
Created:  2024-05-28
Modified: !date!
"""

import os
from aware.site import Site
from astropy.coordinates import Latitude, Longitude
from astropy import units as u
from pytest import approx
import pytz


def test_site_plugin(test_dir: str):
    filename = os.path.join(test_dir, "site_plugins/test_site.json")
    assert os.path.exists(filename), "File does not exist"
    assert os.path.isfile(filename), "Path is not a file"
    telescope = Site.from_json(filename)
    assert telescope is not None, "Cannot load Site plugin from JSON file"
    assert telescope.name == "test_site", "Incorrect site name"
    assert telescope.latitude.to_value(u.deg) == approx(
        Latitude(30, unit=u.deg).to_value(u.deg), rel=1e-6
    ), "Incorrect latitude"
    assert telescope.longitude.to_value(u.deg) == approx(
        Longitude(30, unit=u.deg).to_value(u.deg), rel=1e-6
    ), "Incorrect longitude"
    assert telescope.elevation.to_value(u.m) == approx(
        1600, rel=1e-6
    ), "Incorrect elevation"
    assert telescope.timezone == pytz.timezone("Europe/Moscow"), "Incorrect timezone"
    assert (
        approx(telescope.fov.width.to_value(u.deg), rel=1e-6) == 1
        and approx(telescope.fov.height.to_value(u.deg), rel=1e-6) == 1
    ), "Incorrect FOV"
    assert telescope.default_exposure_number == 3, "Incorrect exposure number"


if __name__ == "__main__":
    test_site_plugin("tests")
