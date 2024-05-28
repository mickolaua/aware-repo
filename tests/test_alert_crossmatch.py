"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_alert_crossmatch.py (c) 2023
Desc: test the crossmatching of the alert pair
Created:  2023-05-07
Modified: 2024-03-31
"""

from __future__ import annotations

import asyncio
import os.path

import pytest

from aware.alert import AlertParsers
from aware.alert.crossmatch import crossmatch_alerts
from aware.consumer.main import dump_alert_to_db


@pytest.mark.asyncio
async def test(test_dir: str):
    parser = AlertParsers["gcn.classic.voevent.SWIFT_XRT_POSITION"]

    with open(
        os.path.join(test_dir, "alert_messages/xrt_alert_1165948_acc.xml"), "rb"
    ) as f:
        alert_msg_acc = f.read()

    with open(
        os.path.join(test_dir, "alert_messages/xrt_alert_1165948.xml"), "rb"
    ) as f:
        alert_msg = f.read()

    info_acc = parser.parse_alert(alert_msg_acc)
    info_acc.rejected = False

    info = parser.parse_alert(alert_msg)
    info.rejected = False

    await dump_alert_to_db(alert_msg_acc, info_acc)

    m = crossmatch_alerts(info)
    assert len(m) == 1, "There is only one matched event!"


if __name__ == "__main__":
    asyncio.run(test(os.path.dirname(__file__)))
