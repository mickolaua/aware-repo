"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_lvk_parser.py (c) 2023
Desc: description
Created:  2023-03-05
Modified: 2023-03-06
"""
import os.path

from aware.alert import AlertParsers
from aware.config import dev
import pytest


@pytest.mark.skipif(
    not dev.value,
    reason="To test parsing of mock data events, dev switch should be turned on",
)
def test(test_dir: str):
    p = AlertParsers["gcn.classic.voevent.LVC_PRELIMINARY"]
    with open(os.path.join(test_dir, "alert_messages", "gw_event.xml"), "rb") as f:
        msg = f.read()
    info = p.parse_alert(msg)

    assert info is not None, "The info on the GW event should be not empty!"
    assert info.event == "MS230223s", "Incorrect GW event name!"
    assert info.origin == "LVC", "The alert should belong to LVC!"
    assert (
        info.trigger_date.isoformat(sep="T") == "2023-02-23T18:44:06"
    ), "Incorrect trigger date!"
    assert info.localization is None, "The event skymap should be not accessed anymore!"
    assert "FAR" in info.description, "FAR should be presented in info!"
    assert "FAR" in info.description, "FAR should be presented in info!"
    assert "Instruments" in info.description, "Instruments should be presented in info!"
    assert "P_BNS" in info.description, "P_BNS should be presented in info!"
    assert "P_NSBH" in info.description, "P_NSBH should be presented in info!"
    assert "P_BBH" in info.description, "P_BBH should be presented in info!"
    assert "P_hasNS" in info.description, "P_hasNS should be presented in info!"
    assert (
        "P_hasRemnant" in info.description
    ), "P_hasRemnant should be presented in info!"


if __name__ == "__main__":
    test("tests")
