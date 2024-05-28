"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_mosaic_sorter.py (c) 2023
Desc: description
Created:  2023-05-07
Modified: !date!
"""
import asyncio
import os.path
from datetime import datetime
from matplotlib import pyplot as plt

import pytest
from astropy.time import Time

from aware.alert import AlertParsers
from aware.site import Telescopes


# TODO: Upgrade this test
@pytest.mark.skip
@pytest.mark.asyncio
async def test(requests_mock):
    requests_mock.get(
        "https://gracedb.ligo.org/api/superevents/MS230223s/files/bayestar.multiorder.fits,0",
        content=open(
            os.path.join(os.path.dirname(__file__), "bayestar.multiorder.fits,0"), "rb"
        ).read(),
    )
    p = AlertParsers["gcn.classic.voevent.LVC_PRELIMINARY"]
    msg = open(os.path.join(os.path.dirname(__file__), "gw_event.xml"), "rb").read()
    info = p.parse_alert(msg)
    loc = info.localization
    as32 = Telescopes["abao_as32"]
    start, stop = as32.nearest_observation_window(Time(datetime.now()))
    sorted_targets, _ = await loc.observe(as32, start=start, stop=stop)  # type: ignore
    assert (
        len(sorted_targets) > 0
    ), "Mosaic walker should retrieve a non-empty list of sky fields."
  

if __name__ == "__main__":
    asyncio.run(test())
