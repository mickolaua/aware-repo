"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_integral_parser.py (c) 2023
Desc: test parsers of the alert messages sent by INTEGRAL spacecraft
Created:  2023-07-10
Modified: !date!
"""
from __future__ import annotations

from aware.alert.plugins.integral import IntegralWeakParser, IntegralSPIACSParser
import os.path
import pytest
from aware.config import dev


@pytest.mark.skipif(
    not dev.value,
    reason="To test parsing of mock data events, dev switch should be turned on",
)
def test_weak(test_dir: str):
    filename = os.path.join(test_dir, "alert_messages/integral_weak.xml")
    with open(filename, "rb") as f:
        info = IntegralWeakParser.parse_alert(f.read())

    assert info is not None, "INTEGRAL weak alert message should be not empty!"
    assert info.localization.center().ra.deg == 1.2343, "Incorrect RA!"
    assert info.localization.center().dec.deg == 2.3455, "Incorrect Dec!"
    assert (
        info.localization.error_radius().to_value("deg") == 0.0169
    ), "Incorrect 90\\%  error radius!"
    assert (
        info.trigger_date.isoformat(timespec="seconds") == "2023-05-23T12:03:23"
    ), "Incorrect trigger date!"


if __name__ == "__main__":
    test_weak("tests")
