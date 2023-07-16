"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_maxi_parser.py (c) 2023
Desc: test MAXI alert parsers
Created:  2023-07-16
Modified: !date!
"""
from __future__ import annotations

from aware.alert.plugins.maxi import MAXIKnownAlertParser, MAXIUnknownAlertParser
import os.path


def test_known(test_dir: str):
    filename = os.path.join(test_dir, "alert_messages/maxi_known.xml")
    with open(filename, "rb") as f:
        info = MAXIKnownAlertParser.parse_alert(f.read())

    assert info is not None, "MAXI alert message on known source should be not empty!"
    assert info.localization.center().ra.deg == 284.1800, "Incorrect RA!"
    assert info.localization.center().dec.deg == 5.3099, "Incorrect Dec!"
    assert (
        info.localization.error_radius().to_value("deg") == 0.5
    ), "Incorrect 90\\%  error radius!"
    assert (
        info.trigger_date.isoformat(timespec="seconds") == "2023-07-10T11:32:56"
    ), "Incorrect trigger date!"


if __name__ == "__main__":
    test_known("tests")
