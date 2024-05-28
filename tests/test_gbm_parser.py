"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_gbm_parser.py (c) 2023
Desc: test Fermi GBM alert parsers
Created:  2023-07-16
Modified: !date!
"""
from __future__ import annotations

from aware.alert.plugins.gbm import GBMAlertAlertParser
import os.path


def test_gbm_alert(test_dir: str):
    filename = os.path.join(test_dir, "alert_messages/gbm_alert.xml")
    with open(filename, "rb") as f:
        info = GBMAlertAlertParser.parse_alert(f.read())

    assert info is not None, "Fermi GBM alert message should be not empty!"
    assert info.localization is None, "This type of alerts does not provide skymap!"
    assert (
        info.trigger_date.isoformat(timespec="seconds") == "2023-07-16T04:54:44"
    ), "Incorrect trigger date!"


if __name__ == "__main__":
    test_gbm_alert("tests")
