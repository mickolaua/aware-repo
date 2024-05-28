"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_swift_parser.py (c) 2023
Desc: description
Created:  2023-02-28
Modified: 2023-03-01
"""
from aware.alert.plugins.swift import SWIFTBATAlertParser
import os.path


def test(test_dir: str):
    filename = os.path.join(test_dir, "alert_messages/bat_alert_grb230228a.xml")
    with open(filename, "rb") as f:
        info = SWIFTBATAlertParser.parse_alert(f.read())
    
    assert info is not None, "Swift BAT alert message should be not empty!"
    assert info.localization.center().ra.deg == 18.4033, "Incorrect RA!"
    assert info.localization.center().dec.deg == 44.4885, "Incorrect Dec!"
    assert info.localization.error_radius().to_value("deg") == 0.0500, "Incorrect 90\\%  error radius!"
    assert info.trigger_date.isoformat(timespec="seconds") == "2023-02-28T05:50:50", "Incorrect trigger date!"


if __name__ == "__main__":
    test("tests")
