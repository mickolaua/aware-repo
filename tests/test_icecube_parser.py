"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_icecube_parser.py (c) 2023
Desc: test parsers of alert messages sent by IceCube
Created:  2023-07-10
Modified: !date!
"""
from __future__ import annotations

from aware.alert.plugins.icecube import (
    ICECUBE_NU_EM_COINC_Parser,
    ICECUBE_ASTROTRACK_BRONZE_Parser,
)
import os.path


def test_bronze(test_dir: str):
    filename = os.path.join(test_dir, "alert_messages/icecube_bronze.xml")
    with open(filename, "rb") as f:
        info = ICECUBE_ASTROTRACK_BRONZE_Parser.parse_alert(f.read())

    assert info is not None, "IceCube Bronze alert message should be not empty!"
    assert info.localization.center().ra.deg == 127.1800, "Incorrect RA!"
    assert info.localization.center().dec.deg == 20.7399, "Incorrect Dec!"
    assert (
        info.localization.error_radius().to_value("deg") == 10.160
    ), "Incorrect 90\\%  error radius!"
    assert (
        info.trigger_date.isoformat(timespec="seconds") == "2023-07-07T18:56:51"
    ), "Incorrect trigger date!"
    assert info.event == "29513102", "Incorrect event id!"
    assert info.meta["FAR"] == 1.3583, "Incorrect FAR!"
    assert (
        info.meta["signalness"] == 0.4663
    ), "Incorrect probability of astrophysical origin"
    assert info.meta["energy"] == 154.2766, "Incorrect neutrino energy!"


def test_nu_em_coinc(test_dir: str):
    filename = os.path.join(test_dir, "alert_messages/icecube_nu_em_coinc.xml")
    with open(filename, "rb") as f:
        info = ICECUBE_NU_EM_COINC_Parser.parse_alert(f.read())

    assert info is not None, "Alert message must be not empty!"
    assert info.localization.center().ra.deg == 23.0100, "Incorrect RA!"
    assert info.localization.center().dec.deg == -10.1300, "Incorrect Dec!"
    assert (
        info.localization.error_radius().to_value("deg") == 0.33
    ), "Incorrect 90\\%  error radius!"
    assert (
        info.trigger_date.isoformat(timespec="seconds") == "2020-05-22T18:21:38"
    ), "Incorrect trigger date!"
    assert info.event == "48580", "Incorrect event id!"
    assert info.meta["FAR"] == 3.2255e+04, "Incorrect FAR!"
    assert (
        info.meta["deltaT"] == 17056.25
    ), "Incorrect time search window!"



if __name__ == "__main__":
    test_bronze("tests")
    test_nu_em_coinc("tests")